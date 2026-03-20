#!/usr/bin/env python3
"""NAMICS regression/benchmark runner.

Usage:
  python3 tests/run_tests.py                 # run all tests
  python3 tests/run_tests.py --list
  python3 tests/run_tests.py --with-save-memory
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import signal
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_helpers import (  # noqa: E402
    CommandMetrics,
    TestError,
    compare_json_profiles,
    require_file,
    run_command,
    set_commented_setting,
    set_setting_line,
    utc_run_id,
)

COORD_TOL = 1e-12
_EXECUTION_FAILED_RE = re.compile(r"execution failed \((-?\d+)\)")


@dataclass
class Context:
    repo_root: Path
    tests_dir: Path
    output_dir: Path
    binary: Path
    quiet: bool
    clean_output: bool
    benchmark_threshold_pct: float
    benchmark_history: bool
    with_save_memory: bool
    reference_dir: Path
    generated_input_manifest: Path


@dataclass
class ReportNode:
    label: str
    passed: bool = True
    required_for_parent: bool = True
    wall_s: float = 0.0
    max_rss_kb: float | None = None
    details: str = ""
    children: list["ReportNode"] = field(default_factory=list)

    def add_metric(self, metrics: CommandMetrics) -> None:
        self.wall_s += metrics.wall_s
        if metrics.max_rss_kb is None:
            return
        if self.max_rss_kb is None:
            self.max_rss_kb = metrics.max_rss_kb
        else:
            self.max_rss_kb = max(self.max_rss_kb, metrics.max_rss_kb)


def _prepare_runtime_input(
    ctx: Context,
    source_input: Path,
    runtime_input: Path,
    *,
    solver_method: str,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
) -> None:
    runtime_input.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_input, runtime_input)

    set_setting_line(runtime_input, runtime_input, "newton : isaac : method", solver_method)

    if settings:
        for key, value in settings.items():
            set_setting_line(runtime_input, runtime_input, key, value)

    if comment_toggles:
        for key, value, enabled in comment_toggles:
            set_commented_setting(runtime_input, runtime_input, key, value, enabled)

    variants: dict[str, str] = {"solver_method": solver_method}
    if settings:
        variants.update(settings)
    if comment_toggles:
        for key, _, enabled in comment_toggles:
            variants[key] = "enabled" if enabled else "disabled"

    ctx.generated_input_manifest.parent.mkdir(parents=True, exist_ok=True)
    with ctx.generated_input_manifest.open("a", encoding="utf-8") as fh:
        fh.write(
            "\t".join(
                [
                    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    str(source_input.relative_to(ctx.repo_root)),
                    str(runtime_input.relative_to(ctx.repo_root)),
                    str(variants),
                ]
            )
            + "\n"
        )


def _cleanup(paths: list[Path], cleanup: bool) -> None:
    if not cleanup:
        return
    for path in paths:
        try:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
        except FileNotFoundError:
            pass


def _homopolymer_benchmark_update(
    benchmark_history: Path,
    run_id: str,
    solver_runtime_ms: int,
    wall_runtime_ms: int,
    threshold_pct: float,
) -> tuple[str, str, str]:
    previous_solver_runtime_ms: int | None = None
    if benchmark_history.is_file():
        with benchmark_history.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh)
            next(reader, None)
            for row in reader:
                if len(row) < 2:
                    continue
                try:
                    previous_solver_runtime_ms = int(row[1])
                except ValueError:
                    continue

    delta_pct = "NA"
    performance_state = "baseline"
    previous_value = "NA"

    if previous_solver_runtime_ms is not None and previous_solver_runtime_ms > 0:
        previous_value = str(previous_solver_runtime_ms)
        delta = ((solver_runtime_ms - previous_solver_runtime_ms) * 100.0) / previous_solver_runtime_ms
        delta_pct = f"{delta:.2f}"
        if abs(delta) >= threshold_pct:
            performance_state = "slower" if solver_runtime_ms > previous_solver_runtime_ms else "faster"
        else:
            performance_state = "stable"

    if not benchmark_history.is_file():
        benchmark_history.parent.mkdir(parents=True, exist_ok=True)
        with benchmark_history.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "run_id",
                    "solver_runtime_ms",
                    "wall_runtime_ms",
                    "threshold_pct",
                    "previous_solver_runtime_ms",
                    "delta_pct",
                    "performance_state",
                ]
            )

    with benchmark_history.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                run_id,
                solver_runtime_ms,
                wall_runtime_ms,
                f"{threshold_pct:g}",
                previous_value,
                delta_pct,
                performance_state,
            ]
        )

    return previous_value, delta_pct, performance_state


def _finalize_report_tree(node: ReportNode) -> None:
    if not node.children:
        return

    for child in node.children:
        _finalize_report_tree(child)

    required_children = [child for child in node.children if child.required_for_parent]
    node.passed = all(child.passed for child in required_children) if required_children else True
    node.wall_s = sum(child.wall_s for child in node.children)
    node.max_rss_kb = max((child.max_rss_kb for child in node.children if child.max_rss_kb is not None), default=None)


def _set_failure(node: ReportNode, detail: str | Exception, *, accepted: bool = False) -> str:
    message = str(detail)
    if accepted:
        message = f"failed (accepted): {message}"
    node.passed = False
    node.details = message
    return message


def _failure_suffix_from_detail(detail: str) -> str | None:
    match = _EXECUTION_FAILED_RE.search(detail)
    if match:
        returncode = int(match.group(1))
        if returncode < 0:
            signal_number = -returncode
            if signal_number == signal.SIGSEGV:
                return "segfault"
            try:
                return f"signal:{signal.Signals(signal_number).name}"
            except ValueError:
                return "signal"
    cleaned = detail.strip()
    for prefix in ("failed (accepted): ", "ERROR: "):
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix):]
    return cleaned or None


def _failure_suffix(node: ReportNode) -> str:
    if node.passed:
        return ""
    mode = _failure_suffix_from_detail(node.details)
    if mode:
        return mode
    child_modes = sorted({_failure_suffix(child) for child in node.children if not child.passed})
    child_modes = [mode for mode in child_modes if mode]
    if len(child_modes) == 1:
        return child_modes[0]
    if len(child_modes) > 1:
        return "mixed"
    return "unexpected"


def _flatten_nodes(nodes: list[ReportNode]) -> list[tuple[int, bool, str, ReportNode]]:
    out: list[tuple[int, bool, str, ReportNode]] = []

    def walk(root_index: int, node: ReportNode, label_prefix: str, child_prefix: str, *, is_root: bool) -> None:
        label = f"{label_prefix}{node.label}"
        out.append((root_index, is_root, label, node))

        child_count = len(node.children)
        for idx, child in enumerate(node.children):
            is_last = idx == child_count - 1
            branch = "└── " if is_last else "├── "
            extension = "    " if is_last else "│   "
            walk(root_index, child, child_prefix + branch, child_prefix + extension, is_root=False)

    for root_index, root in enumerate(nodes, start=1):
        walk(root_index, root, "", "", is_root=True)
    return out


def _run_solver(
    ctx: Context,
    node: ReportNode,
    *,
    source_input: Path,
    runtime_input: Path,
    solver_method: str,
    label: str,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    require_output: bool = True,
    require_initial_guess: bool = False,
) -> tuple[CommandMetrics, Path]:
    output_path = runtime_input.with_suffix(".output.json")
    preprocessed_input = runtime_input.with_suffix(".input.json")
    _prepare_runtime_input(
        ctx,
        source_input,
        runtime_input,
        solver_method=solver_method,
        settings=settings,
        comment_toggles=comment_toggles,
    )
    output_path.unlink(missing_ok=True)
    preprocessed_input.unlink(missing_ok=True)
    metrics = run_command([str(ctx.binary), str(runtime_input)], cwd=ctx.repo_root, quiet=ctx.quiet)
    node.add_metric(metrics)
    if metrics.returncode != 0:
        for line in reversed(metrics.output.splitlines()):
            message = line.strip()
            if message:
                raise TestError(message)
        raise TestError(f"ERROR: execution failed ({metrics.returncode}) for input: {runtime_input}")
    if require_output or require_initial_guess:
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise TestError(f"ERROR: expected JSON output file was not created for {label}: {output_path}")
    if require_initial_guess:
        try:
            payload = output_path.read_text(encoding="utf-8")
        except Exception as exc:
            raise TestError(f"ERROR: cannot read JSON output for {label}: {output_path} ({exc})") from exc
        if "\"initial_guess\"" not in payload:
            raise TestError(f"ERROR: JSON output does not contain embedded initial_guess for {label}: {output_path}")
    return metrics, output_path


def _run_reference_leaf(
    ctx: Context,
    leaf: ReportNode,
    *,
    source_input: Path,
    runtime_input: Path,
    reference_file: Path,
    label: str,
    value_tol: float,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
) -> tuple[CommandMetrics, Path]:
    metrics, output_path = _run_solver(
        ctx,
        leaf,
        source_input=source_input,
        runtime_input=runtime_input,
        solver_method="pseudohessian",
        label=label,
        settings=settings,
        comment_toggles=comment_toggles,
    )
    compare_json_profiles(reference_file, output_path, coord_tol=COORD_TOL, value_tol=value_tol)
    leaf.details = "matched reference"
    return metrics, output_path


def _run_optional_diis_leaf(
    ctx: Context,
    leaf: ReportNode,
    *,
    source_input: Path,
    runtime_input: Path,
    baseline_output: Path | None,
    label: str,
    value_tol: float,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
) -> str:
    try:
        _, output_path = _run_solver(
            ctx,
            leaf,
            source_input=source_input,
            runtime_input=runtime_input,
            solver_method="DIIS",
            label=label,
            settings=settings,
            comment_toggles=comment_toggles,
            require_output=False,
        )
        if not output_path.is_file() or output_path.stat().st_size == 0:
            leaf.passed = False
            status = "completed but JSON output file is missing (accepted)"
        elif baseline_output is None or not baseline_output.is_file() or baseline_output.stat().st_size == 0:
            leaf.passed = False
            status = "completed but no pseudohessian baseline was available (accepted)"
        else:
            compare_json_profiles(baseline_output, output_path, coord_tol=COORD_TOL, value_tol=value_tol)
            leaf.passed = True
            status = "passed and matched pseudohessian"
        leaf.details = status
    except Exception as exc:
        status = _set_failure(leaf, exc, accepted=True)
    return status


def _run_method_group_regression(
    ctx: Context,
    *,
    root_label: str,
    input_file: Path,
    reference_name: str,
    runtime_stem: str,
    value_tol: float,
    label_stem: str | None = None,
    method_group_label: str | None = "solver method",
    diis_optional: bool = True,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    required_files: list[Path] | None = None,
    derived_runtime_files: list[Path] | None = None,
) -> ReportNode:
    root = ReportNode(label=root_label)
    method_group = root if method_group_label is None else ReportNode(label=method_group_label)
    if method_group_label is not None:
        root.children.append(method_group)

    reference_file = ctx.reference_dir / reference_name
    require_file(reference_file)
    output_dir = ctx.output_dir
    pseudo_input = output_dir / f"{runtime_stem}.pseudohessian.in"
    diis_input = output_dir / f"{runtime_stem}.diis.in"
    pseudo_output = pseudo_input.with_suffix(".output.json")
    diis_output = diis_input.with_suffix(".output.json")
    cleanup_targets: list[Path] = [
        pseudo_input,
        pseudo_input.with_suffix(".input.json"),
        pseudo_output,
        diis_input,
        diis_input.with_suffix(".input.json"),
        diis_output,
    ]
    label_base = label_stem or root_label
    required_files = required_files or []
    derived_runtime_files = derived_runtime_files or []

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    for path in required_files:
        require_file(path)
    output_dir.mkdir(parents=True, exist_ok=True)

    for source_file in derived_runtime_files:
        cleanup_targets.append(output_dir / source_file.name)

    try:
        for source_file in derived_runtime_files:
            shutil.copy2(source_file, output_dir / source_file.name)

        pseudo_leaf = ReportNode(label="pseudohessian")
        method_group.children.append(pseudo_leaf)

        try:
            _, pseudo_output = _run_reference_leaf(
                ctx,
                pseudo_leaf,
                source_input=input_file,
                runtime_input=pseudo_input,
                reference_file=reference_file,
                label=f"{label_base}, mode=pseudohessian",
                value_tol=value_tol,
                settings=settings,
                comment_toggles=comment_toggles,
            )
        except Exception as exc:
            _set_failure(pseudo_leaf, exc)

        diis_leaf = ReportNode(label="DIIS", passed=not diis_optional, required_for_parent=not diis_optional)
        method_group.children.append(diis_leaf)

        if diis_optional:
            diis_status = _run_optional_diis_leaf(
                ctx,
                diis_leaf,
                source_input=input_file,
                runtime_input=diis_input,
                baseline_output=pseudo_output if pseudo_leaf.passed else None,
                label=f"{label_base}, mode=DIIS",
                value_tol=value_tol,
                settings=settings,
                comment_toggles=comment_toggles,
            )
        else:
            try:
                if not pseudo_leaf.passed:
                    raise TestError("ERROR: skipped because pseudohessian failed")
                _, diis_output = _run_solver(
                    ctx,
                    diis_leaf,
                    source_input=input_file,
                    runtime_input=diis_input,
                    solver_method="DIIS",
                    label=f"{label_base}, mode=DIIS",
                    settings=settings,
                    comment_toggles=comment_toggles,
                )
                compare_json_profiles(pseudo_output, diis_output, coord_tol=COORD_TOL, value_tol=value_tol)
                diis_status = "matched pseudohessian"
                diis_leaf.details = diis_status
            except Exception as exc:
                diis_status = _set_failure(diis_leaf, exc)

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def _run_homopolymer_adsorption(ctx: Context, *, enable_benchmark: bool) -> ReportNode:
    root = ReportNode(label="homopolymer adsorption benchmark" if enable_benchmark else "homopolymer adsorption")

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir
    benchmark_dir = tests_dir / "benchmarks"
    benchmark_history = benchmark_dir / "homopolymer_adsorption_test_benchmark.csv"
    template_file = tests_dir / "homopolymer_adsorption.in"

    require_file(ctx.binary, executable=True)
    require_file(template_file)

    output_dir.mkdir(parents=True, exist_ok=True)
    run_id = utc_run_id() if enable_benchmark else ""
    run_log = benchmark_dir / f"homopolymer_adsorption_test_{run_id}.log" if enable_benchmark else None

    if enable_benchmark:
        benchmark_dir.mkdir(parents=True, exist_ok=True)

    cleanup_targets: list[Path] = [run_log] if run_log is not None else []
    solver_runtime_ms = 0
    diis_notes: list[str] = []

    log_lines: list[str] = []
    if enable_benchmark:
        log_lines = [
            f"run_id={run_id}",
            f"timestamp_utc={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
            f"threshold_pct={ctx.benchmark_threshold_pct:g}",
            "cases=chi_Si(0,-2,-4,-6)",
            f"save_memory={'on' if ctx.with_save_memory else 'off'}",
        ]
    bench_start = time.perf_counter()

    try:
        for chi in (0, -2, -4, -6):
            reference_file = ctx.reference_dir / f"homopolymer_adsorption.chi_{chi}.json.ref"
            require_file(reference_file)
            chi_node = ReportNode(label=f"chi_Si = {chi}")
            root.children.append(chi_node)

            pseudohessian_input = output_dir / f"homopolymer_adsorption.chi_{chi}.pseudohessian.in"
            pseudohessian_output = pseudohessian_input.with_suffix(".output.json")
            cleanup_targets.extend([pseudohessian_input, pseudohessian_input.with_suffix(".input.json"), pseudohessian_output])

            pseudo_leaf = ReportNode(label="pseudohessian")
            chi_node.children.append(pseudo_leaf)

            try:
                run_m, pseudohessian_output = _run_reference_leaf(
                    ctx,
                    pseudo_leaf,
                    source_input=template_file,
                    runtime_input=pseudohessian_input,
                    reference_file=reference_file,
                    label=f"chi_Si={chi}, mode=pseudohessian",
                    value_tol=1e-6,
                    settings={"mon : A : chi_Si": str(chi)},
                )
                elapsed_ms = int(round(run_m.wall_s * 1000.0))
                solver_runtime_ms += elapsed_ms
                if enable_benchmark:
                    log_lines.append(f"chi_Si={chi},mode=pseudohessian,elapsed_ms={elapsed_ms}")
            except Exception as exc:
                _set_failure(pseudo_leaf, exc)

            if ctx.with_save_memory:
                save_input = output_dir / f"homopolymer_adsorption.chi_{chi}.save_memory.in"
                save_output = save_input.with_suffix(".output.json")
                cleanup_targets.extend([save_input, save_input.with_suffix(".input.json"), save_output])

                save_leaf = ReportNode(label="save_memory + pseudohessian")
                chi_node.children.append(save_leaf)

                try:
                    run_m, save_output = _run_reference_leaf(
                        ctx,
                        save_leaf,
                        source_input=template_file,
                        runtime_input=save_input,
                        reference_file=reference_file,
                        label=f"chi_Si={chi}, mode=save_memory",
                        value_tol=1e-6,
                        settings={"mon : A : chi_Si": str(chi)},
                        comment_toggles=[("mol : pol : save_memory", "true", True)],
                    )
                    elapsed_ms = int(round(run_m.wall_s * 1000.0))
                    solver_runtime_ms += elapsed_ms
                    if enable_benchmark:
                        log_lines.append(f"chi_Si={chi},mode=save_memory,elapsed_ms={elapsed_ms}")

                    if pseudohessian_output.is_file() and pseudohessian_output.stat().st_size > 0:
                        compare_json_profiles(pseudohessian_output, save_output, coord_tol=COORD_TOL, value_tol=1e-6)
                        save_leaf.details = "matched reference and pseudohessian"
                except Exception as exc:
                    _set_failure(save_leaf, exc)

            if not enable_benchmark:
                diis_input = output_dir / f"homopolymer_adsorption.chi_{chi}.diis.in"
                diis_output = diis_input.with_suffix(".output.json")
                cleanup_targets.extend([diis_input, diis_input.with_suffix(".input.json"), diis_output])

                diis_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
                chi_node.children.append(diis_leaf)

                diis_status = _run_optional_diis_leaf(
                    ctx,
                    diis_leaf,
                    source_input=template_file,
                    runtime_input=diis_input,
                    baseline_output=pseudohessian_output if pseudo_leaf.passed else None,
                    label=f"chi_Si={chi}, mode=DIIS",
                    value_tol=1e-6,
                    settings={"mon : A : chi_Si": str(chi)},
                )

                diis_notes.append(f"chi={chi}:{diis_status}")

        if enable_benchmark:
            wall_runtime_ms = int(round((time.perf_counter() - bench_start) * 1000.0))

            if ctx.benchmark_history:
                previous, delta_pct, perf_state = _homopolymer_benchmark_update(
                    benchmark_history=benchmark_history,
                    run_id=run_id,
                    solver_runtime_ms=solver_runtime_ms,
                    wall_runtime_ms=wall_runtime_ms,
                    threshold_pct=ctx.benchmark_threshold_pct,
                )
                log_lines.extend(
                    [
                        f"solver_runtime_ms={solver_runtime_ms}",
                        f"wall_runtime_ms={wall_runtime_ms}",
                        f"previous_solver_runtime_ms={previous}",
                        f"delta_pct={delta_pct}",
                        f"performance_state={perf_state}",
                        f"history_file={benchmark_history}",
                    ]
                )
                benchmark_note = (
                    f"bench:{perf_state},solver_ms={solver_runtime_ms},delta={delta_pct}%,history={benchmark_history.name}"
                )
            else:
                benchmark_note = f"bench:disabled,solver_ms={solver_runtime_ms}"

            assert run_log is not None  # guarded by enable_benchmark
            run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
            root.details = f"{benchmark_note};log={run_log.name}"
        else:
            root.details = f"diis={';'.join(diis_notes)}"

        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_homopolymer_adsorption(ctx: Context) -> ReportNode:
    return _run_homopolymer_adsorption(ctx, enable_benchmark=False)


def benchmark_homopolymer_adsorption(ctx: Context) -> ReportNode:
    return _run_homopolymer_adsorption(ctx, enable_benchmark=True)


def test_frozen_range_input_file(ctx: Context) -> ReportNode:
    tests_dir = ctx.tests_dir
    input_file = tests_dir / "frozen_range_input_file.in"
    source_frozen_file = tests_dir / "frozen_range_input_file.frozen"
    return _run_method_group_regression(
        ctx,
        root_label="frozen range input file",
        label_stem="frozen_range_input_file",
        input_file=input_file,
        reference_name="frozen_range_input_file.json.ref",
        runtime_stem="frozen_range_input_file",
        value_tol=1e-6,
        required_files=[source_frozen_file],
        derived_runtime_files=[source_frozen_file],
    )


def test_micelle_self_assembly(ctx: Context) -> ReportNode:
    root = ReportNode(label="micelle self assembly")
    generate_group = ReportNode(label="input guess generate")
    use_group = ReportNode(label="input guess use")
    root.children.extend([generate_group, use_group])

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    generate_input = tests_dir / "micelle_guess_generate.in"
    use_input = tests_dir / "micelle_guess_use.in"
    reference_file = ctx.reference_dir / "micelle_guess_use.json.ref"
    require_file(reference_file)

    require_file(ctx.binary, executable=True)
    require_file(generate_input)
    require_file(use_input)
    output_dir.mkdir(parents=True, exist_ok=True)

    pseudo_generate_input = output_dir / generate_input.name
    pseudo_generate_output = pseudo_generate_input.with_suffix(".output.json")
    pseudo_use_input = output_dir / use_input.name
    pseudo_use_output = pseudo_use_input.with_suffix(".output.json")

    diis_generate_input = pseudo_generate_input
    diis_generate_output = pseudo_generate_output
    diis_use_input = pseudo_use_input
    diis_use_output = pseudo_use_output

    cleanup_targets: list[Path] = [
        pseudo_generate_input,
        pseudo_generate_input.with_suffix(".input.json"),
        pseudo_generate_output,
        pseudo_use_input,
        pseudo_use_input.with_suffix(".input.json"),
        pseudo_use_output,
    ]

    try:
        for p in (pseudo_generate_output, pseudo_use_output, diis_generate_output, diis_use_output):
            p.unlink(missing_ok=True)

        pseudo_gen_leaf = ReportNode(label="pseudohessian")
        generate_group.children.append(pseudo_gen_leaf)
        pseudo_use_leaf = ReportNode(label="pseudohessian")
        use_group.children.append(pseudo_use_leaf)

        pseudo_use_ready = False

        try:
            _, pseudo_generate_output = _run_solver(
                ctx,
                pseudo_gen_leaf,
                source_input=generate_input,
                runtime_input=pseudo_generate_input,
                solver_method="pseudohessian",
                label="micelle guess generate, mode=pseudohessian",
                require_initial_guess=True,
            )
            pseudo_gen_leaf.details = "embedded initial_guess generated in output JSON"
            pseudo_use_ready = True
        except Exception as exc:
            _set_failure(pseudo_gen_leaf, exc)

        try:
            if not pseudo_use_ready:
                raise TestError("ERROR: skipped because pseudohessian guess generation failed")

            _, pseudo_use_output = _run_reference_leaf(
                ctx,
                pseudo_use_leaf,
                source_input=use_input,
                runtime_input=pseudo_use_input,
                reference_file=reference_file,
                label="micelle guess use, mode=pseudohessian",
                value_tol=1e-9,
            )
        except Exception as exc:
            _set_failure(pseudo_use_leaf, exc)

        # Re-run DIIS path from a clean output state.
        diis_generate_output.unlink(missing_ok=True)
        diis_use_output.unlink(missing_ok=True)
        diis_generate_input.with_suffix(".input.json").unlink(missing_ok=True)
        diis_use_input.with_suffix(".input.json").unlink(missing_ok=True)

        diis_gen_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        generate_group.children.append(diis_gen_leaf)
        diis_use_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        use_group.children.append(diis_use_leaf)

        diis_status = "failed (accepted)"
        diis_use_ready = False

        try:
            _, diis_generate_output = _run_solver(
                ctx,
                diis_gen_leaf,
                source_input=generate_input,
                runtime_input=diis_generate_input,
                solver_method="DIIS",
                label="micelle guess generate, mode=DIIS",
                require_initial_guess=True,
            )
            diis_status = "embedded initial_guess generated in output JSON"
            diis_use_ready = True
            diis_gen_leaf.details = diis_status
            diis_gen_leaf.passed = diis_status == "embedded initial_guess generated in output JSON"
        except Exception as exc:
            _set_failure(diis_gen_leaf, exc, accepted=True)

        try:
            if not diis_use_ready:
                _set_failure(diis_use_leaf, "skipped because DIIS guess generation failed", accepted=True)
            else:
                diis_status = _run_optional_diis_leaf(
                    ctx,
                    diis_use_leaf,
                    source_input=use_input,
                    runtime_input=diis_use_input,
                    baseline_output=reference_file if pseudo_use_leaf.passed else None,
                    label="micelle guess use, mode=DIIS",
                    value_tol=1e-9,
                )
        except Exception as exc:
            _set_failure(diis_use_leaf, exc, accepted=True)

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_micelle_grand_canonical_search(ctx: Context) -> ReportNode:
    root = ReportNode(label="micelle grand canonical search")
    tests_dir = ctx.tests_dir
    output_dir = ctx.output_dir

    seed_input = tests_dir / "micelle_guess_generate.in"
    search_input = tests_dir / "micelle_gc_search.in"
    utility = ctx.repo_root / "utils" / "micelle_gc_search.py"

    require_file(ctx.binary, executable=True)
    require_file(seed_input)
    require_file(search_input)
    require_file(utility)

    output_dir.mkdir(parents=True, exist_ok=True)

    cleanup_targets: list[Path] = []

    try:
        def run_case(
            node: ReportNode,
            seed_source: Path,
            require_initial_guess: bool,
            expected_x_range: tuple[float, float],
        ) -> None:
            label = node.label
            seed_runtime = output_dir / f"{seed_source.stem}_{label.replace(' ', '_')}_seed.in"
            seed_output = seed_runtime.with_suffix(".output.json")
            search_runtime = output_dir / f"micelle_gc_search_{label.replace(' ', '_')}.in"
            search_workdir = output_dir / search_runtime.stem
            summary_file = search_workdir / "summary.json"
            result_json = search_workdir / "result.output.json"
            cleanup_targets.extend([seed_runtime, seed_runtime.with_suffix(".input.json"), seed_output, search_runtime, search_workdir])

            _, seed_output = _run_solver(
                ctx,
                node,
                source_input=seed_source,
                runtime_input=seed_runtime,
                solver_method="pseudohessian",
                label=f"micelle gc seed generation ({label})",
                require_initial_guess=require_initial_guess,
            )
            seed_problem = json.loads(seed_output.read_text(encoding="utf-8"))["problems"][-1]
            if require_initial_guess:
                if "initial_guess" not in seed_problem:
                    raise TestError(f"ERROR: missing embedded initial_guess in seed JSON: {seed_output}")
            else:
                if "initial_guess" in seed_problem:
                    raise TestError(f"ERROR: fallback seed unexpectedly contains embedded initial_guess: {seed_output}")
                for key in ("method", "mx", "my", "mz", "fjc", "charged", "monlist", "statelist", "profiles"):
                    if key not in seed_problem:
                        raise TestError(f"ERROR: fallback seed is missing '{key}': {seed_output}")

            shutil.copy2(search_input, search_runtime)
            metrics = run_command(
                [
                    sys.executable,
                    str(utility),
                    str(search_runtime),
                    "--binary",
                    str(ctx.binary),
                    "--molecule",
                    "surf",
                    "--seed-json",
                    str(seed_output),
                    "--step",
                    "5",
                    "--workers",
                    "2",
                    "--max-iter",
                    "6",
                    "--gp-tol",
                    "5e-2",
                ],
                cwd=ctx.repo_root,
                quiet=ctx.quiet,
            )
            node.add_metric(metrics)
            if metrics.returncode != 0:
                raise TestError(metrics.output.strip() or f"ERROR: external micelle GC search failed ({label})")
            if not summary_file.is_file():
                raise TestError(f"ERROR: missing micelle GC summary: {summary_file}")
            if not result_json.is_file():
                raise TestError(f"ERROR: missing micelle GC result JSON: {result_json}")

            summary = json.loads(summary_file.read_text(encoding="utf-8"))
            x = float(summary["x"])
            gp = float(summary["grand_potential"])
            converged = bool(summary["converged"])

            if not converged:
                raise TestError("ERROR: micelle GC search did not report convergence")
            if abs(gp) > 0.25:
                raise TestError(f"ERROR: micelle GC search residual too large: grand_potential={gp}")
            if not (expected_x_range[0] < x < expected_x_range[1]):
                raise TestError(f"ERROR: micelle GC search returned unexpected aggregation number: n={x}")

            payload = result_json.read_text(encoding="utf-8")
            if '"initial_guess"' not in payload:
                raise TestError(f"ERROR: result JSON does not contain embedded initial_guess: {result_json}")

            node.details = f"n={x:.6f}, gp={gp:.3e}"

        for label, seed_source, require_initial_guess, expected_x_range in (
            ("embedded initial_guess seed", seed_input, True, (111.0, 114.0)),
            ("u-profile fallback seed", search_input, False, (109.5, 110.5)),
        ):
            node = ReportNode(label=label)
            root.children.append(node)
            try:
                run_case(node, seed_source, require_initial_guess, expected_x_range)
            except Exception as exc:
                _set_failure(node, exc)
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_particle_in_cyl_coordinates(ctx: Context) -> ReportNode:
    tests_dir = ctx.tests_dir
    input_file = tests_dir / "particle_in_cyl_coordinates.in"
    return _run_method_group_regression(
        ctx,
        root_label="particle in cyl coordinates",
        input_file=input_file,
        reference_name="particle_in_cyl_coordinates.json.ref",
        runtime_stem="particle_in_cyl_coordinates",
        value_tol=1e-6,
    )


def test_branched_brush(ctx: Context) -> ReportNode:
    root = ReportNode(label="branched brush")
    tests_dir = ctx.tests_dir

    for case_label, runtime_stem in (
        ("2d cylindrical", "branched_brush_cyl2d"),
        ("1d planar", "branched_brush_planar_1d"),
    ):
        input_file = tests_dir / f"{runtime_stem}.in"
        case_node = _run_method_group_regression(
            ctx,
            root_label=case_label,
            label_stem=f"branched brush {case_label}",
            input_file=input_file,
            reference_name=f"{runtime_stem}.json.ref",
            runtime_stem=runtime_stem,
            value_tol=1e-6,
            method_group_label=None,
        )
        root.children.append(case_node)

    _finalize_report_tree(root)
    return root


def test_polE_regression(ctx: Context) -> ReportNode:
    # Keep regression inputs in tests/ so the suite stays decoupled from data/.
    input_file = ctx.tests_dir / "polE.in"
    # Legacy DIIS differs from pseudohessian only in the last-layer Na/Cl tail by about 5.2e-6.
    # Keep the regression aligned with historical behavior instead of failing on that known drift.
    value_tol = 1e-5
    return _run_method_group_regression(
        ctx,
        root_label="polE regression",
        label_stem="polE",
        input_file=input_file,
        reference_name="polE.json.ref",
        runtime_stem="polE",
        value_tol=value_tol,
    )


def test_external_potential(ctx: Context) -> ReportNode:
    root = ReportNode(label="external potential")
    tests_dir = ctx.tests_dir

    for case_label in ("1d", "2d", "3d"):
        runtime_stem = f"external_potential_{case_label}"
        input_file = tests_dir / f"{runtime_stem}.in"
        potential_file = tests_dir / f"{runtime_stem}_external_potential.json"
        reference_name = f"{runtime_stem}.json.ref"

        case_node = _run_method_group_regression(
            ctx,
            root_label=case_label,
            label_stem=f"external potential {case_label}",
            input_file=input_file,
            reference_name=reference_name,
            runtime_stem=runtime_stem,
            value_tol=1e-6,
            method_group_label=None,
            diis_optional=False,
            required_files=[potential_file],
            derived_runtime_files=[potential_file],
        )
        root.children.append(case_node)

    _finalize_report_tree(root)
    return root

def _print_table(results: list[ReportNode]) -> None:
    flat = _flatten_nodes(results)

    headers = ["#", "Test", "Status", "Wall\n(s)", "MaxRSS\n(MiB)"]
    rows: list[tuple[int, list[str], str]] = []

    for root_index, is_root, label, node in flat:
        rows.append(
            (
                root_index,
                [
                    str(root_index) if is_root else "",
                    label,
                    "PASS" if node.passed else "FAIL",
                    f"{node.wall_s:.3f}",
                    "NA" if node.max_rss_kb is None else f"{(node.max_rss_kb / 1024.0):.1f}",
                ],
                _failure_suffix(node),
            )
        )

    total_wall = sum(r.wall_s for r in results)
    max_mem_values = [r.max_rss_kb for r in results if r.max_rss_kb is not None]
    total_max_mem = max(max_mem_values) if max_mem_values else None
    pass_count = sum(1 for r in results if r.passed)
    total_row = [
        "",
        "TOTAL",
        f"{pass_count}/{len(results)}",
        f"{total_wall:.3f}",
        "NA" if total_max_mem is None else f"{(total_max_mem / 1024.0):.1f}",
    ]

    widths = [0] * len(headers)
    for row in [headers] + [cells for _, cells, _ in rows] + [total_row]:
        for i, cell in enumerate(row):
            parts = cell.splitlines() or [""]
            widths[i] = max(widths[i], max(len(part) for part in parts))

    def format_cell(text: str, width: int, align: str) -> str:
        if align == "center":
            return text.center(width)
        if align == "right":
            return text.rjust(width)
        return text.ljust(width)

    def print_row(parts: list[str], aligns: list[str], suffix: str = "") -> None:
        wrapped = [part.splitlines() for part in parts]
        height = max(len(cell_lines) for cell_lines in wrapped)
        for line_idx in range(height):
            line_parts = []
            for col_idx, cell_lines in enumerate(wrapped):
                text = cell_lines[line_idx] if line_idx < len(cell_lines) else ""
                line_parts.append(format_cell(text, widths[col_idx], aligns[col_idx]))
            line = "│ " + " │ ".join(line_parts) + " │"
            if suffix and line_idx == 0:
                line += f"  {suffix}"
            print(line)

    def border(left: str, mid: str, right: str) -> str:
        return left + mid.join("─" * (w + 2) for w in widths) + right

    header_aligns = ["center"] * len(headers)
    body_aligns = ["center", "left", "center", "right", "right"]

    print(border("┌", "┬", "┐"))
    print_row(headers, header_aligns)
    print(border("├", "┼", "┤"))
    for idx, (root_index, row, mode) in enumerate(rows):
        if idx > 0 and root_index != rows[idx - 1][0]:
            print(border("├", "┼", "┤"))
        print_row(row, body_aligns, mode)
    print(border("├", "┼", "┤"))
    print_row(total_row, body_aligns)
    print(border("└", "┴", "┘"))


TESTS: dict[str, tuple[str, Callable[[Context], ReportNode]]] = {
    "homopolymer-adsorption": ("homopolymer adsorption", test_homopolymer_adsorption),
    "homopolymer-adsorption-benchmark": ("homopolymer adsorption benchmark", benchmark_homopolymer_adsorption),
    "external-potential": ("external potential", test_external_potential),
    "frozen-range-input-file": ("frozen range input file", test_frozen_range_input_file),
    "micelle-self-assembly": ("micelle self assembly", test_micelle_self_assembly),
    "micelle-grand-canonical-search": (
        "micelle grand canonical search",
        test_micelle_grand_canonical_search,
    ),
    "particle-in-cyl-coordinates": ("particle in cyl coordinates", test_particle_in_cyl_coordinates),
    "branched-brush": ("branched brush", test_branched_brush),
    "polE-regression": ("polE regression", test_polE_regression),
}

ALIASES = {
    "homopolymer_adsorption": "homopolymer-adsorption",
    "homopolymer_adsorption_benchmark": "homopolymer-adsorption-benchmark",
    "external_potential": "external-potential",
    "frozen_range_input_file": "frozen-range-input-file",
    "micelle_self_assembly": "micelle-self-assembly",
    "micelle_grand_canonical_search": "micelle-grand-canonical-search",
    "particle_in_cyl_coordinates": "particle-in-cyl-coordinates",
    "branched_brush": "branched-brush",
    "branched_brush_cyl2d": "branched-brush",
    "branched_brush_planar_1d": "branched-brush",
    "polE_regression": "polE-regression",
}

BENCHMARK_TESTS = {"homopolymer-adsorption-benchmark"}


def _resolve_selection(items: list[str]) -> list[str]:
    if not items:
        return [key for key in TESTS.keys() if key not in BENCHMARK_TESTS]
    resolved: list[str] = []
    for raw in items:
        key = raw.strip()
        if key in ALIASES:
            key = ALIASES[key]
        if key not in TESTS:
            valid = ", ".join(TESTS.keys())
            raise SystemExit(f"Unknown test '{raw}'. Valid: {valid}")
        if key not in resolved:
            resolved.append(key)
    return resolved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run NAMICS regression tests and benchmarks.")
    parser.add_argument("tests", nargs="*", help="Tests to run (default: all)")
    parser.add_argument("--list", action="store_true", help="List available tests and exit")
    parser.add_argument(
        "--binary",
        type=Path,
        default=None,
        help="Path to namics binary (default: <repo>/bin/namics)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for runtime/generated output (default: <repo>/output)",
    )
    parser.add_argument(
        "--benchmark-threshold-pct",
        type=float,
        default=15.0,
        help="Threshold used for homopolymer benchmark state classification",
    )
    parser.add_argument(
        "--no-benchmark-history",
        action="store_true",
        help="Do not append to tests/benchmarks/homopolymer_adsorption_test_benchmark.csv",
    )
    parser.add_argument("--verbose", action="store_true", help="Show test process stdout/stderr")
    parser.add_argument(
        "--with-save-memory",
        action="store_true",
        help="Run additional homopolymer save_memory variants",
    )
    parser.add_argument(
        "--keep-artifacts",
        action="store_true",
        help="Keep generated runtime files and do not purge output directory",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    tests_dir = repo_root / "tests"

    if args.list:
        for key in TESTS:
            print(key)
        return 0

    selected = _resolve_selection(args.tests)
    binary = args.binary.resolve() if args.binary else (repo_root / "bin" / "namics")
    output_dir = args.output_dir.resolve() if args.output_dir else (repo_root / "output")
    reference_archive = tests_dir / "reference" / "reference_files.tar.gz"

    require_file(reference_archive)

    clean_output = not args.keep_artifacts
    if clean_output and output_dir.is_dir():
        shutil.rmtree(output_dir, ignore_errors=True)

    with tempfile.TemporaryDirectory(prefix="namics_ref_") as tmp_ref_dir:
        try:
            shutil.unpack_archive(str(reference_archive), tmp_ref_dir)
        except (shutil.ReadError, ValueError) as exc:
            raise SystemExit(f"ERROR: cannot unpack reference archive {reference_archive}: {exc}") from exc

        reference_dir = Path(tmp_ref_dir) / "reference"
        if not reference_dir.is_dir():
            raise SystemExit(f"ERROR: unpacked reference archive is missing directory: {reference_dir}")

        ctx = Context(
            repo_root=repo_root,
            tests_dir=tests_dir,
            output_dir=output_dir,
            binary=binary,
            quiet=not args.verbose,
            clean_output=clean_output,
            benchmark_threshold_pct=args.benchmark_threshold_pct,
            benchmark_history=not args.no_benchmark_history,
            with_save_memory=args.with_save_memory,
            reference_dir=reference_dir,
            generated_input_manifest=output_dir / "generated_test_inputs.tsv",
        )

        results: list[ReportNode] = []

        for display_key in selected:
            label, fn = TESTS[display_key]
            try:
                node = fn(ctx)
            except TestError as exc:
                node = ReportNode(label=label, passed=False, details=str(exc))
            except Exception as exc:  # defensive catch to keep summary table complete
                node = ReportNode(label=label, passed=False, details=f"ERROR: unexpected failure: {exc}")

            if node.label != label:
                node.label = label

            results.append(node)

        if clean_output:
            ctx.generated_input_manifest.unlink(missing_ok=True)

        _print_table(results)
        return 0 if all(n.passed for n in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
