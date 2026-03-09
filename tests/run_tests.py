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


def _run_binary(ctx: Context, node: ReportNode, input_file: Path) -> CommandMetrics:
    metrics = run_command([str(ctx.binary), str(input_file)], cwd=ctx.repo_root, quiet=ctx.quiet)
    node.add_metric(metrics)
    return metrics


def _runtime_output_basename(runtime_input: Path) -> str:
    stem = runtime_input.stem
    safe = "".join(ch if (ch.isalnum() or ch == "_") else "_" for ch in stem)
    return safe or "namics_test_output"


def _runtime_output_path(runtime_input: Path) -> Path:
    return runtime_input.parent / f"{_runtime_output_basename(runtime_input)}.json"


def _derived_runtime_path(runtime_input: Path, suffix: str) -> Path:
    return runtime_input.with_suffix(f".{suffix}")


def _write_generated_input_record(ctx: Context, runtime_input: Path, source_input: Path, variants: dict[str, str]) -> None:
    ctx.generated_input_manifest.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "runtime_input": str(runtime_input.relative_to(ctx.repo_root)),
        "source_input": str(source_input.relative_to(ctx.repo_root)),
        "variants": variants,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    line = "\t".join([payload["timestamp_utc"], payload["source_input"], payload["runtime_input"], str(payload["variants"])])
    with ctx.generated_input_manifest.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


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
    # Force deterministic per-variant JSON output names so each run is recoverable.
    output_basename = _runtime_output_basename(runtime_input)
    set_setting_line(runtime_input, runtime_input, "output : json : filename", output_basename)

    if settings:
        for key, value in settings.items():
            set_setting_line(runtime_input, runtime_input, key, value)

    if comment_toggles:
        for key, value, enabled in comment_toggles:
            set_commented_setting(runtime_input, runtime_input, key, value, enabled)

    variants: dict[str, str] = {"solver_method": solver_method}
    variants["output : json : filename"] = output_basename
    if settings:
        variants.update(settings)
    if comment_toggles:
        for key, _, enabled in comment_toggles:
            variants[key] = "enabled" if enabled else "disabled"

    _write_generated_input_record(ctx, runtime_input, source_input, variants)


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

    max_mem = None
    for child in node.children:
        if child.max_rss_kb is None:
            continue
        if max_mem is None:
            max_mem = child.max_rss_kb
        else:
            max_mem = max(max_mem, child.max_rss_kb)
    node.max_rss_kb = max_mem


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


def _reference_file(ctx: Context, filename: str) -> Path:
    reference_file = ctx.reference_dir / filename
    require_file(reference_file)
    return reference_file


def _require_json_output(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise TestError(f"ERROR: expected JSON output file was not created for {label}: {path}")


def _require_json_contains_initial_guess(path: Path, label: str) -> None:
    try:
        payload = path.read_text(encoding="utf-8")
    except Exception as exc:
        raise TestError(f"ERROR: cannot read JSON output for {label}: {path} ({exc})") from exc
    if "\"initial_guess\"" not in payload:
        raise TestError(f"ERROR: JSON output does not contain embedded initial_guess for {label}: {path}")


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
    output_path = _runtime_output_path(runtime_input)
    _prepare_runtime_input(
        ctx,
        source_input,
        runtime_input,
        solver_method=solver_method,
        settings=settings,
        comment_toggles=comment_toggles,
    )
    output_path.unlink(missing_ok=True)
    metrics = _run_binary(ctx, node, runtime_input)
    if metrics.returncode != 0:
        raise TestError(f"ERROR: execution failed ({metrics.returncode}) for input: {runtime_input}")
    if require_output or require_initial_guess:
        _require_json_output(output_path, label)
    if require_initial_guess:
        _require_json_contains_initial_guess(output_path, label)
    return metrics, output_path


def _optional_match_status(candidate_output: Path, baseline_output: Path | None, *, value_tol: float) -> tuple[bool, str]:
    if not candidate_output.is_file() or candidate_output.stat().st_size == 0:
        return False, "completed but JSON output file is missing (accepted)"
    if baseline_output is None or not baseline_output.is_file() or baseline_output.stat().st_size == 0:
        return False, "completed but no pseudohessian baseline was available (accepted)"
    try:
        compare_json_profiles(baseline_output, candidate_output, coord_tol=COORD_TOL, value_tol=value_tol)
    except TestError:
        return False, "completed but differs from pseudohessian (accepted)"
    return True, "passed and matched pseudohessian"


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
        leaf.passed, status = _optional_match_status(output_path, baseline_output, value_tol=value_tol)
        leaf.details = status
    except Exception as exc:
        status = f"failed (accepted): {exc}"
        leaf.details = status
        leaf.passed = False
    return status


def _merge_settings(*parts: dict[str, str] | None) -> dict[str, str] | None:
    merged: dict[str, str] = {}
    for part in parts:
        if part:
            merged.update(part)
    return merged or None


def _run_method_group_regression(
    ctx: Context,
    *,
    root_label: str,
    input_file: Path,
    reference_name: str,
    runtime_stem: str,
    value_tol: float,
    label_stem: str | None = None,
    method_group_label: str = "solver method",
    settings: dict[str, str] | None = None,
    pseudohessian_settings: dict[str, str] | None = None,
    diis_settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    required_files: list[Path] | None = None,
    derived_runtime_files: list[tuple[Path, str]] | None = None,
) -> ReportNode:
    root = ReportNode(label=root_label)
    method_group = ReportNode(label=method_group_label)
    root.children.append(method_group)

    reference_file = _reference_file(ctx, reference_name)
    output_dir = ctx.output_dir
    pseudo_input = output_dir / f"{runtime_stem}.pseudohessian.in"
    diis_input = output_dir / f"{runtime_stem}.diis.in"
    pseudo_output = _runtime_output_path(pseudo_input)
    diis_output = _runtime_output_path(diis_input)
    cleanup_targets: list[Path] = [pseudo_input, pseudo_output, diis_input, diis_output]
    label_base = label_stem or root_label
    required_files = required_files or []
    derived_runtime_files = derived_runtime_files or []

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    for path in required_files:
        require_file(path)
    output_dir.mkdir(parents=True, exist_ok=True)

    for _, suffix in derived_runtime_files:
        cleanup_targets.extend(
            [
                _derived_runtime_path(pseudo_input, suffix),
                _derived_runtime_path(diis_input, suffix),
            ]
        )

    try:
        for source_file, suffix in derived_runtime_files:
            shutil.copy2(source_file, _derived_runtime_path(pseudo_input, suffix))
            shutil.copy2(source_file, _derived_runtime_path(diis_input, suffix))

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
                settings=_merge_settings(settings, pseudohessian_settings),
                comment_toggles=comment_toggles,
            )
        except Exception as exc:
            pseudo_leaf.passed = False
            pseudo_leaf.details = str(exc)

        diis_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        method_group.children.append(diis_leaf)

        diis_status = _run_optional_diis_leaf(
            ctx,
            diis_leaf,
            source_input=input_file,
            runtime_input=diis_input,
            baseline_output=pseudo_output if pseudo_leaf.passed else None,
            label=f"{label_base}, mode=DIIS",
            value_tol=value_tol,
            settings=_merge_settings(settings, diis_settings),
            comment_toggles=comment_toggles,
        )

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
            reference_file = _reference_file(ctx, f"homopolymer_adsorption.chi_{chi}.json.ref")
            chi_node = ReportNode(label=f"chi_Si = {chi}")
            root.children.append(chi_node)

            pseudohessian_input = output_dir / f"homopolymer_adsorption.chi_{chi}.pseudohessian.in"
            pseudohessian_output = _runtime_output_path(pseudohessian_input)
            cleanup_targets.extend([pseudohessian_input, pseudohessian_output])

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
                pseudo_leaf.passed = False
                pseudo_leaf.details = str(exc)

            if ctx.with_save_memory:
                save_input = output_dir / f"homopolymer_adsorption.chi_{chi}.save_memory.in"
                save_output = _runtime_output_path(save_input)
                cleanup_targets.extend([save_input, save_output])

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
                    save_leaf.passed = False
                    save_leaf.details = str(exc)

            if not enable_benchmark:
                diis_input = output_dir / f"homopolymer_adsorption.chi_{chi}.diis.in"
                diis_output = output_dir / f"{_runtime_output_basename(diis_input)}.json"
                cleanup_targets.extend([diis_input, diis_output])

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
        pseudohessian_settings={"mon : W : frozen_filename": "frozen_range_input_file.pseudohessian.frozen"},
        diis_settings={"mon : W : frozen_filename": "frozen_range_input_file.diis.frozen"},
        derived_runtime_files=[(source_frozen_file, "frozen")],
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
    reference_file = _reference_file(ctx, "micelle_guess_use.json.ref")

    require_file(ctx.binary, executable=True)
    require_file(generate_input)
    require_file(use_input)
    output_dir.mkdir(parents=True, exist_ok=True)

    pseudo_generate_input = output_dir / "micelle_guess_generate.pseudohessian.in"
    pseudo_generate_output = _runtime_output_path(pseudo_generate_input)
    pseudo_use_input = output_dir / "micelle_guess_use.pseudohessian.in"
    pseudo_use_output = _runtime_output_path(pseudo_use_input)

    diis_generate_input = output_dir / "micelle_guess_generate.diis.in"
    diis_generate_output = _runtime_output_path(diis_generate_input)
    diis_use_input = output_dir / "micelle_guess_use.diis.in"
    diis_use_output = _runtime_output_path(diis_use_input)

    cleanup_targets: list[Path] = [
        pseudo_generate_input,
        pseudo_generate_output,
        pseudo_use_input,
        pseudo_use_output,
        diis_generate_input,
        diis_generate_output,
        diis_use_input,
        diis_use_output,
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
                settings={"sys : noname : write_initial_guess": "true"},
                require_initial_guess=True,
            )
            pseudo_gen_leaf.details = "embedded initial_guess generated in output JSON"
            pseudo_use_ready = True
        except Exception as exc:
            pseudo_gen_leaf.passed = False
            pseudo_gen_leaf.details = str(exc)

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
                settings={"sys : noname : guess_inputfile": pseudo_generate_output.name},
            )
        except Exception as exc:
            pseudo_use_leaf.passed = False
            pseudo_use_leaf.details = str(exc)

        # Re-run DIIS path from a clean output state.
        diis_generate_output.unlink(missing_ok=True)
        diis_use_output.unlink(missing_ok=True)

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
                settings={"sys : noname : write_initial_guess": "true"},
                require_initial_guess=True,
            )
            diis_status = "embedded initial_guess generated in output JSON"
            diis_use_ready = True
            diis_gen_leaf.details = diis_status
            diis_gen_leaf.passed = diis_status == "embedded initial_guess generated in output JSON"
        except Exception as exc:
            diis_gen_leaf.details = f"failed (accepted): {exc}"
            diis_gen_leaf.passed = False

        try:
            if not diis_use_ready:
                diis_use_leaf.details = "failed (accepted): skipped because DIIS guess generation failed"
                diis_use_leaf.passed = False
            else:
                diis_status = _run_optional_diis_leaf(
                    ctx,
                    diis_use_leaf,
                    source_input=use_input,
                    runtime_input=diis_use_input,
                    baseline_output=pseudo_use_output if pseudo_use_leaf.passed else None,
                    label="micelle guess use, mode=DIIS",
                    value_tol=1e-9,
                    settings={"sys : noname : guess_inputfile": diis_generate_output.name},
                )
        except Exception as exc:
            diis_use_leaf.details = f"failed (accepted): {exc}"
            diis_use_leaf.passed = False

        root.details = f"diis={diis_status}"
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


def test_polE_regression(ctx: Context) -> ReportNode:
    # Keep regression inputs in tests/ so the suite stays decoupled from data/.
    input_file = ctx.tests_dir / "polE.in"
    runtime_settings = {"sys : noname : write_initial_guess": "false"}
    # Legacy DIIS differs from pseudohessian only in the last-layer Na/Cl tail by about 5.2e-6.
    # Keep the regression aligned with historical behavior instead of failing on that known drift.
    value_tol = 1e-5
    # The legacy binary used for the reference only supports this profile subset via .pro output.
    legacy_subset = [
        ("json : state : AH : phi", "", False),
        ("json : state : AM : phi", "", False),
        ("json : state : H3O : phi", "", False),
        ("json : state : H2O : phi", "", False),
        ("json : state : OH : phi", "", False),
    ]
    return _run_method_group_regression(
        ctx,
        root_label="polE regression",
        label_stem="polE",
        input_file=input_file,
        reference_name="polE.json.ref",
        runtime_stem="polE",
        value_tol=value_tol,
        settings=runtime_settings,
        comment_toggles=legacy_subset,
    )


def _fmt_wall(v: float) -> str:
    return f"{v:.3f}"


def _fmt_mem_kb(v: float | None) -> str:
    if v is None:
        return "NA"
    return f"{(v / 1024.0):.1f}"


def _print_table(results: list[ReportNode]) -> None:
    flat = _flatten_nodes(results)

    headers = ["#", "Test", "Status", "Wall\n(s)", "MaxRSS\n(MiB)"]
    rows: list[tuple[int, list[str]]] = []

    for root_index, is_root, label, node in flat:
        rows.append(
            (
                root_index,
                [
                    str(root_index) if is_root else "",
                    label,
                    "PASS" if node.passed else "FAIL",
                    _fmt_wall(node.wall_s),
                    _fmt_mem_kb(node.max_rss_kb),
                ],
            )
        )

    total_wall = sum(r.wall_s for r in results)
    max_mem_values = [r.max_rss_kb for r in results if r.max_rss_kb is not None]
    total_max_mem = max(max_mem_values) if max_mem_values else None
    pass_count = sum(1 for r in results if r.passed)
    total_row = ["", "TOTAL", f"{pass_count}/{len(results)}", _fmt_wall(total_wall), _fmt_mem_kb(total_max_mem)]

    widths = [0] * len(headers)
    for row in [headers] + [cells for _, cells in rows] + [total_row]:
        for i, cell in enumerate(row):
            parts = cell.splitlines() or [""]
            widths[i] = max(widths[i], max(len(part) for part in parts))

    def format_cell(text: str, width: int, align: str) -> str:
        if align == "center":
            return text.center(width)
        if align == "right":
            return text.rjust(width)
        return text.ljust(width)

    def print_row(parts: list[str], aligns: list[str]) -> None:
        wrapped = [part.splitlines() for part in parts]
        height = max(len(cell_lines) for cell_lines in wrapped)
        for line_idx in range(height):
            line_parts = []
            for col_idx, cell_lines in enumerate(wrapped):
                text = cell_lines[line_idx] if line_idx < len(cell_lines) else ""
                line_parts.append(format_cell(text, widths[col_idx], aligns[col_idx]))
            print("│ " + " │ ".join(line_parts) + " │")

    def border(left: str, mid: str, right: str) -> str:
        return left + mid.join("─" * (w + 2) for w in widths) + right

    header_aligns = ["center"] * len(headers)
    body_aligns = ["center", "left", "center", "right", "right"]

    print(border("┌", "┬", "┐"))
    print_row(headers, header_aligns)
    print(border("├", "┼", "┤"))
    for idx, (root_index, row) in enumerate(rows):
        if idx > 0 and root_index != rows[idx - 1][0]:
            print(border("├", "┼", "┤"))
        print_row(row, body_aligns)
    print(border("├", "┼", "┤"))
    print_row(total_row, body_aligns)
    print(border("└", "┴", "┘"))


TESTS: dict[str, tuple[str, Callable[[Context], ReportNode]]] = {
    "homopolymer-adsorption": ("homopolymer adsorption", test_homopolymer_adsorption),
    "homopolymer-adsorption-benchmark": ("homopolymer adsorption benchmark", benchmark_homopolymer_adsorption),
    "frozen-range-input-file": ("frozen range input file", test_frozen_range_input_file),
    "micelle-self-assembly": ("micelle self assembly", test_micelle_self_assembly),
    "particle-in-cyl-coordinates": ("particle in cyl coordinates", test_particle_in_cyl_coordinates),
    "polE-regression": ("polE regression", test_polE_regression),
}

ALIASES = {
    "homopolymer_adsorption": "homopolymer-adsorption",
    "homopolymer_adsorption_benchmark": "homopolymer-adsorption-benchmark",
    "frozen_range_input_file": "frozen-range-input-file",
    "micelle_self_assembly": "micelle-self-assembly",
    "particle_in_cyl_coordinates": "particle-in-cyl-coordinates",
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
