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
    ensure_file_from_archive,
    require_file,
    run_command,
    set_commented_setting,
    set_setting_line,
    utc_run_id,
)


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
    reference_archive: Path
    reference_unpack_dir: Path
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


def _flatten_nodes(nodes: list[ReportNode]) -> list[tuple[int, ReportNode]]:
    out: list[tuple[int, ReportNode]] = []

    def walk(level: int, node: ReportNode) -> None:
        out.append((level, node))
        for child in node.children:
            walk(level + 1, child)

    for root in nodes:
        walk(0, root)
    return out


def _resolve_reference_file(ctx: Context, reference_file: Path) -> Path:
    if reference_file.is_file():
        return reference_file

    cache_target = ctx.reference_unpack_dir / reference_file.relative_to(ctx.tests_dir / "reference")
    if cache_target.is_file():
        return cache_target

    member_name = reference_file.relative_to(ctx.tests_dir).as_posix()
    recovered = ensure_file_from_archive(ctx.reference_archive, member_name, cache_target)
    if recovered and cache_target.is_file():
        return cache_target

    raise TestError(
        f"ERROR: reference file is missing and could not be recovered from archive: {reference_file} "
        f"(archive: {ctx.reference_archive})"
    )


def _check_json_output(path: Path, label: str) -> tuple[bool, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return False, f"ERROR: expected JSON output file was not created for {label}: {path}"
    return True, ""


def test_homopolymer_adsorption(ctx: Context) -> ReportNode:
    root = ReportNode(label="homopolymer adsorption")

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir
    reference_dir = tests_dir / "reference"
    benchmark_dir = tests_dir / "benchmarks"
    benchmark_history = benchmark_dir / "homopolymer_adsorption_test_benchmark.csv"
    template_file = tests_dir / "homopolymer_adsorption.in"

    require_file(ctx.binary, executable=True)
    require_file(template_file)

    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    run_id = utc_run_id()
    run_log = benchmark_dir / f"homopolymer_adsorption_test_{run_id}.log"

    cleanup_targets: list[Path] = [run_log]
    solver_runtime_ms = 0
    diis_notes: list[str] = []

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
            reference_file = _resolve_reference_file(
                ctx, reference_dir / f"homopolymer_adsorption.chi_{chi}.json.ref"
            )
            chi_node = ReportNode(label=f"chi_Si = {chi}")
            root.children.append(chi_node)

            pseudohessian_input = output_dir / f"homopolymer_adsorption.chi_{chi}.pseudohessian.in"
            pseudohessian_output = output_dir / f"{_runtime_output_basename(pseudohessian_input)}.json"
            cleanup_targets.extend([pseudohessian_input, pseudohessian_output])

            pseudo_leaf = ReportNode(label="pseudohessian")
            chi_node.children.append(pseudo_leaf)

            try:
                _prepare_runtime_input(
                    ctx,
                    template_file,
                    pseudohessian_input,
                    solver_method="pseudohessian",
                    settings={"mon : A : chi_Si": str(chi)},
                )
                pseudohessian_output.unlink(missing_ok=True)
                run_m = _run_binary(ctx, pseudo_leaf, pseudohessian_input)
                solver_runtime_ms += int(round(run_m.wall_s * 1000.0))
                log_lines.append(f"chi_Si={chi},mode=pseudohessian,elapsed_ms={int(round(run_m.wall_s * 1000.0))}")

                if run_m.returncode != 0:
                    raise TestError(
                        f"ERROR: execution failed ({run_m.returncode}) for input: {pseudohessian_input}"
                    )

                ok, msg = _check_json_output(pseudohessian_output, f"chi_Si={chi}, mode=pseudohessian")
                if not ok:
                    raise TestError(msg)

                compare_json_profiles(reference_file, pseudohessian_output, coord_tol=1e-12, value_tol=1e-6)
                pseudo_leaf.details = "matched reference"
            except Exception as exc:
                pseudo_leaf.passed = False
                pseudo_leaf.details = str(exc)

            if ctx.with_save_memory:
                save_input = output_dir / f"homopolymer_adsorption.chi_{chi}.save_memory.in"
                save_output = output_dir / f"{_runtime_output_basename(save_input)}.json"
                cleanup_targets.extend([save_input, save_output])

                save_leaf = ReportNode(label="save_memory + pseudohessian")
                chi_node.children.append(save_leaf)

                try:
                    _prepare_runtime_input(
                        ctx,
                        template_file,
                        save_input,
                        solver_method="pseudohessian",
                        settings={"mon : A : chi_Si": str(chi)},
                        comment_toggles=[("mol : pol : save_memory", "true", True)],
                    )
                    save_output.unlink(missing_ok=True)
                    run_m = _run_binary(ctx, save_leaf, save_input)
                    solver_runtime_ms += int(round(run_m.wall_s * 1000.0))
                    log_lines.append(f"chi_Si={chi},mode=save_memory,elapsed_ms={int(round(run_m.wall_s * 1000.0))}")

                    if run_m.returncode != 0:
                        raise TestError(f"ERROR: execution failed ({run_m.returncode}) for input: {save_input}")

                    ok, msg = _check_json_output(save_output, f"chi_Si={chi}, mode=save_memory")
                    if not ok:
                        raise TestError(msg)

                    compare_json_profiles(reference_file, save_output, coord_tol=1e-12, value_tol=1e-6)

                    if pseudohessian_output.is_file() and pseudohessian_output.stat().st_size > 0:
                        compare_json_profiles(pseudohessian_output, save_output, coord_tol=1e-12, value_tol=1e-6)
                    save_leaf.details = "matched reference and pseudohessian"
                except Exception as exc:
                    save_leaf.passed = False
                    save_leaf.details = str(exc)

            diis_input = output_dir / f"homopolymer_adsorption.chi_{chi}.diis.in"
            diis_output = output_dir / f"{_runtime_output_basename(diis_input)}.json"
            cleanup_targets.extend([diis_input, diis_output])

            diis_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
            chi_node.children.append(diis_leaf)

            diis_status = "failed (accepted)"
            try:
                _prepare_runtime_input(
                    ctx,
                    template_file,
                    diis_input,
                    solver_method="DIIS",
                    settings={"mon : A : chi_Si": str(chi)},
                )
                diis_output.unlink(missing_ok=True)
                run_diis = _run_binary(ctx, diis_leaf, diis_input)

                if run_diis.returncode == 0:
                    if diis_output.is_file() and diis_output.stat().st_size > 0:
                        if pseudohessian_output.is_file() and pseudohessian_output.stat().st_size > 0:
                            try:
                                compare_json_profiles(
                                    pseudohessian_output,
                                    diis_output,
                                    coord_tol=1e-12,
                                    value_tol=1e-6,
                                )
                                diis_status = "passed and matched pseudohessian"
                            except TestError:
                                diis_status = "completed but differs from pseudohessian (accepted)"
                        else:
                            diis_status = "completed but no pseudohessian baseline was available (accepted)"
                    else:
                        diis_status = "completed but JSON output file is missing (accepted)"

                diis_leaf.details = diis_status
                diis_leaf.passed = diis_status == "passed and matched pseudohessian"
                log_lines.append(f"chi_Si={chi},mode=DIIS,status={diis_status}")
            except Exception as exc:
                # DIIS failures are accepted for these cases.
                diis_leaf.details = f"failed (accepted): {exc}"
                diis_leaf.passed = False

            diis_notes.append(f"chi={chi}:{diis_status}")

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

        run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
        root.details = f"{benchmark_note};diis={';'.join(diis_notes)};log={run_log.name}"

        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_frozen_range_input_file(ctx: Context) -> ReportNode:
    root = ReportNode(label="frozen range input file")
    method_group = ReportNode(label="solver method")
    root.children.append(method_group)

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    input_file = tests_dir / "frozen_range_input_file.in"
    source_frozen_file = tests_dir / "frozen_range_input_file.frozen"
    reference_file = _resolve_reference_file(
        ctx, tests_dir / "reference" / "frozen_range_input_file.json.ref"
    )

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    require_file(source_frozen_file)
    output_dir.mkdir(parents=True, exist_ok=True)

    runtime_frozen_pseudo = output_dir / "frozen_range_input_file.pseudohessian.frozen"
    runtime_frozen_diis = output_dir / "frozen_range_input_file.diis.frozen"
    pseudo_input = output_dir / "frozen_range_input_file.pseudohessian.in"
    pseudo_output = output_dir / f"{_runtime_output_basename(pseudo_input)}.json"
    diis_input = output_dir / "frozen_range_input_file.diis.in"
    diis_output = output_dir / f"{_runtime_output_basename(diis_input)}.json"

    cleanup_targets: list[Path] = [
        runtime_frozen_pseudo,
        runtime_frozen_diis,
        pseudo_input,
        pseudo_output,
        diis_input,
        diis_output,
    ]

    try:
        shutil.copy2(source_frozen_file, runtime_frozen_pseudo)
        shutil.copy2(source_frozen_file, runtime_frozen_diis)

        pseudo_leaf = ReportNode(label="pseudohessian")
        method_group.children.append(pseudo_leaf)

        try:
            _prepare_runtime_input(ctx, input_file, pseudo_input, solver_method="pseudohessian")
            pseudo_output.unlink(missing_ok=True)
            run_m = _run_binary(ctx, pseudo_leaf, pseudo_input)

            if run_m.returncode != 0:
                raise TestError(f"ERROR: execution failed ({run_m.returncode}) for input: {pseudo_input}")

            ok, msg = _check_json_output(pseudo_output, "frozen_range_input_file, mode=pseudohessian")
            if not ok:
                raise TestError(msg)

            compare_json_profiles(reference_file, pseudo_output, coord_tol=1e-12, value_tol=1e-6)
            pseudo_leaf.details = "matched reference"
        except Exception as exc:
            pseudo_leaf.passed = False
            pseudo_leaf.details = str(exc)

        diis_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        method_group.children.append(diis_leaf)

        diis_status = "failed (accepted)"
        try:
            _prepare_runtime_input(ctx, input_file, diis_input, solver_method="DIIS")
            diis_output.unlink(missing_ok=True)
            run_diis = _run_binary(ctx, diis_leaf, diis_input)
            if run_diis.returncode == 0:
                if diis_output.is_file() and diis_output.stat().st_size > 0:
                    if pseudo_output.is_file() and pseudo_output.stat().st_size > 0:
                        try:
                            compare_json_profiles(pseudo_output, diis_output, coord_tol=1e-12, value_tol=1e-6)
                            diis_status = "passed and matched pseudohessian"
                        except TestError:
                            diis_status = "completed but differs from pseudohessian (accepted)"
                    else:
                        diis_status = "completed but no pseudohessian baseline was available (accepted)"
                else:
                    diis_status = "completed but output file is missing (accepted)"
            diis_leaf.details = diis_status
            diis_leaf.passed = diis_status == "passed and matched pseudohessian"
        except Exception as exc:
            diis_leaf.details = f"failed (accepted): {exc}"
            diis_leaf.passed = False

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_micelle_self_assembly(ctx: Context) -> ReportNode:
    root = ReportNode(label="micelle self assembly")
    generate_group = ReportNode(label="input guess generate")
    use_group = ReportNode(label="input guess use")
    root.children.extend([generate_group, use_group])

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    generate_input = tests_dir / "micelle_guess_generate.in"
    use_input = tests_dir / "micelle_guess_use.in"
    reference_file = _resolve_reference_file(ctx, tests_dir / "reference" / "micelle_guess_use.json.ref")

    guess_file = output_dir / "micelle_2.outi"

    require_file(ctx.binary, executable=True)
    require_file(generate_input)
    require_file(use_input)
    output_dir.mkdir(parents=True, exist_ok=True)

    pseudo_generate_input = output_dir / "micelle_guess_generate.pseudohessian.in"
    pseudo_generate_output = output_dir / f"{_runtime_output_basename(pseudo_generate_input)}.json"
    pseudo_use_input = output_dir / "micelle_guess_use.pseudohessian.in"
    pseudo_use_output = output_dir / f"{_runtime_output_basename(pseudo_use_input)}.json"

    diis_generate_input = output_dir / "micelle_guess_generate.diis.in"
    diis_generate_output = output_dir / f"{_runtime_output_basename(diis_generate_input)}.json"
    diis_use_input = output_dir / "micelle_guess_use.diis.in"
    diis_use_output = output_dir / f"{_runtime_output_basename(diis_use_input)}.json"

    cleanup_targets: list[Path] = [
        guess_file,
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
        for p in (guess_file, pseudo_generate_output, pseudo_use_output, diis_generate_output, diis_use_output):
            p.unlink(missing_ok=True)

        pseudo_gen_leaf = ReportNode(label="pseudohessian")
        generate_group.children.append(pseudo_gen_leaf)
        pseudo_use_leaf = ReportNode(label="pseudohessian")
        use_group.children.append(pseudo_use_leaf)

        pseudo_use_ready = False

        try:
            _prepare_runtime_input(ctx, generate_input, pseudo_generate_input, solver_method="pseudohessian")
            run_gen = _run_binary(ctx, pseudo_gen_leaf, pseudo_generate_input)
            if run_gen.returncode != 0:
                raise TestError(f"ERROR: execution failed ({run_gen.returncode}) for input: {pseudo_generate_input}")

            if not guess_file.is_file() or guess_file.stat().st_size == 0:
                raise TestError(f"ERROR: expected guess file was not created: {guess_file}")

            pseudo_gen_leaf.details = "guess file generated"
            pseudo_use_ready = True
        except Exception as exc:
            pseudo_gen_leaf.passed = False
            pseudo_gen_leaf.details = str(exc)

        try:
            if not pseudo_use_ready:
                raise TestError("ERROR: skipped because pseudohessian guess generation failed")

            _prepare_runtime_input(ctx, use_input, pseudo_use_input, solver_method="pseudohessian")
            run_use = _run_binary(ctx, pseudo_use_leaf, pseudo_use_input)
            if run_use.returncode != 0:
                raise TestError(f"ERROR: execution failed ({run_use.returncode}) for input: {pseudo_use_input}")

            ok, msg = _check_json_output(pseudo_use_output, "micelle guess use, mode=pseudohessian")
            if not ok:
                raise TestError(msg)

            compare_json_profiles(reference_file, pseudo_use_output, coord_tol=1e-12, value_tol=1e-9)
            pseudo_use_leaf.details = "matched reference"
        except Exception as exc:
            pseudo_use_leaf.passed = False
            pseudo_use_leaf.details = str(exc)

        # Re-run DIIS path from a clean guess/output state.
        guess_file.unlink(missing_ok=True)
        diis_generate_output.unlink(missing_ok=True)
        diis_use_output.unlink(missing_ok=True)

        diis_gen_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        generate_group.children.append(diis_gen_leaf)
        diis_use_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        use_group.children.append(diis_use_leaf)

        diis_status = "failed (accepted)"
        diis_use_ready = False

        try:
            _prepare_runtime_input(ctx, generate_input, diis_generate_input, solver_method="DIIS")
            run_gen_diis = _run_binary(ctx, diis_gen_leaf, diis_generate_input)
            if run_gen_diis.returncode == 0 and guess_file.is_file() and guess_file.stat().st_size > 0:
                diis_status = "guess generated"
                diis_use_ready = True
            diis_gen_leaf.details = diis_status
            diis_gen_leaf.passed = diis_status == "guess generated"
        except Exception as exc:
            diis_gen_leaf.details = f"failed (accepted): {exc}"
            diis_gen_leaf.passed = False

        try:
            if not diis_use_ready:
                diis_use_leaf.details = "failed (accepted): skipped because DIIS guess generation failed"
                diis_use_leaf.passed = False
            else:
                _prepare_runtime_input(ctx, use_input, diis_use_input, solver_method="DIIS")
                run_use_diis = _run_binary(ctx, diis_use_leaf, diis_use_input)
                if run_use_diis.returncode == 0:
                    if diis_use_output.is_file() and diis_use_output.stat().st_size > 0:
                        if pseudo_use_output.is_file() and pseudo_use_output.stat().st_size > 0:
                            try:
                                compare_json_profiles(pseudo_use_output, diis_use_output, coord_tol=1e-12, value_tol=1e-9)
                                diis_status = "passed and matched pseudohessian"
                            except TestError:
                                diis_status = "completed but differs from pseudohessian (accepted)"
                        else:
                            diis_status = "completed but no pseudohessian baseline was available (accepted)"
                    else:
                        diis_status = "completed but JSON output file is missing (accepted)"
                diis_use_leaf.details = diis_status
                diis_use_leaf.passed = diis_status == "passed and matched pseudohessian"
        except Exception as exc:
            diis_use_leaf.details = f"failed (accepted): {exc}"
            diis_use_leaf.passed = False

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def test_particle_in_cyl_coordinates(ctx: Context) -> ReportNode:
    root = ReportNode(label="particle in cyl coordinates")
    method_group = ReportNode(label="solver method")
    root.children.append(method_group)

    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    input_file = tests_dir / "particle_in_cyl_coordinates.in"
    reference_file = _resolve_reference_file(
        ctx, tests_dir / "reference" / "particle_in_cyl_coordinates.json.ref"
    )

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    output_dir.mkdir(parents=True, exist_ok=True)

    pseudo_input = output_dir / "particle_in_cyl_coordinates.pseudohessian.in"
    pseudo_output = output_dir / f"{_runtime_output_basename(pseudo_input)}.json"
    diis_input = output_dir / "particle_in_cyl_coordinates.diis.in"
    diis_output = output_dir / f"{_runtime_output_basename(diis_input)}.json"

    cleanup_targets: list[Path] = [pseudo_input, pseudo_output, diis_input, diis_output]

    try:
        pseudo_leaf = ReportNode(label="pseudohessian")
        method_group.children.append(pseudo_leaf)

        try:
            _prepare_runtime_input(ctx, input_file, pseudo_input, solver_method="pseudohessian")
            pseudo_output.unlink(missing_ok=True)
            run_pseudo = _run_binary(ctx, pseudo_leaf, pseudo_input)
            if run_pseudo.returncode != 0:
                raise TestError(f"ERROR: execution failed ({run_pseudo.returncode}) for input: {pseudo_input}")

            ok, msg = _check_json_output(pseudo_output, "particle_in_cyl_coordinates, mode=pseudohessian")
            if not ok:
                raise TestError(msg)

            compare_json_profiles(reference_file, pseudo_output, coord_tol=1e-12, value_tol=1e-6)
            pseudo_leaf.details = "matched reference"
        except Exception as exc:
            pseudo_leaf.passed = False
            pseudo_leaf.details = str(exc)

        diis_leaf = ReportNode(label="DIIS", passed=False, required_for_parent=False)
        method_group.children.append(diis_leaf)

        diis_status = "failed (accepted)"
        try:
            _prepare_runtime_input(ctx, input_file, diis_input, solver_method="DIIS")
            diis_output.unlink(missing_ok=True)
            run_diis = _run_binary(ctx, diis_leaf, diis_input)
            if run_diis.returncode == 0:
                if diis_output.is_file() and diis_output.stat().st_size > 0:
                    if pseudo_output.is_file() and pseudo_output.stat().st_size > 0:
                        try:
                            compare_json_profiles(pseudo_output, diis_output, coord_tol=1e-12, value_tol=1e-6)
                            diis_status = "passed and matched pseudohessian"
                        except TestError:
                            diis_status = "completed but differs from pseudohessian (accepted)"
                    else:
                        diis_status = "completed but no pseudohessian baseline was available (accepted)"
                else:
                    diis_status = "completed but JSON output file is missing (accepted)"
            diis_leaf.details = diis_status
            diis_leaf.passed = diis_status == "passed and matched pseudohessian"
        except Exception as exc:
            diis_leaf.details = f"failed (accepted): {exc}"
            diis_leaf.passed = False

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, cleanup=ctx.clean_output)


def _fmt_wall(v: float) -> str:
    return f"{v:.3f}"


def _fmt_mem_kb(v: float | None) -> str:
    if v is None:
        return "NA"
    return f"{(v / 1024.0):.1f}"


def _print_table(results: list[ReportNode]) -> None:
    flat = _flatten_nodes(results)

    headers = ["Test", "Status", "Wall(s)", "MaxRSS(MiB)"]
    rows = []

    for level, node in flat:
        label = f"{'  ' * level}{node.label}"
        rows.append(
            [
                label,
                "PASS" if node.passed else "FAIL",
                _fmt_wall(node.wall_s),
                _fmt_mem_kb(node.max_rss_kb),
            ]
        )

    total_wall = sum(r.wall_s for r in results)
    max_mem_values = [r.max_rss_kb for r in results if r.max_rss_kb is not None]
    total_max_mem = max(max_mem_values) if max_mem_values else None
    pass_count = sum(1 for r in results if r.passed)
    rows.append(["TOTAL", f"{pass_count}/{len(results)}", _fmt_wall(total_wall), _fmt_mem_kb(total_max_mem)])

    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def line(parts: list[str]) -> str:
        return " | ".join(parts[i].ljust(widths[i]) for i in range(len(parts)))

    sep = "-+-".join("-" * w for w in widths)
    print(line(headers))
    print(sep)
    for row in rows:
        print(line(row))


TESTS: dict[str, tuple[str, Callable[[Context], ReportNode]]] = {
    "homopolymer-adsorption": ("homopolymer adsorption", test_homopolymer_adsorption),
    "frozen-range-input-file": ("frozen range input file", test_frozen_range_input_file),
    "micelle-self-assembly": ("micelle self assembly", test_micelle_self_assembly),
    "particle-in-cyl-coordinates": ("particle in cyl coordinates", test_particle_in_cyl_coordinates),
}

ALIASES = {
    "homopolymer_adsorption": "homopolymer-adsorption",
    "frozen_range_input_file": "frozen-range-input-file",
    "micelle_self_assembly": "micelle-self-assembly",
    "particle_in_cyl_coordinates": "particle-in-cyl-coordinates",
}


def _resolve_selection(items: list[str]) -> list[str]:
    if not items:
        return list(TESTS.keys())
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
        "--clean-output",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--with-save-memory",
        action="store_true",
        help="Run additional homopolymer save_memory variants",
    )
    parser.add_argument(
        "--reference-archive",
        type=Path,
        default=None,
        help="Path to reference archive (default: tests/reference/reference_files.tar.gz)",
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

    reference_archive = args.reference_archive.resolve() if args.reference_archive else (tests_dir / "reference" / "reference_files.tar.gz")

    selected = _resolve_selection(args.tests)
    binary = args.binary.resolve() if args.binary else (repo_root / "bin" / "namics")
    output_dir = args.output_dir.resolve() if args.output_dir else (repo_root / "output")

    clean_output = not args.keep_artifacts
    if clean_output and output_dir.is_dir():
        shutil.rmtree(output_dir, ignore_errors=True)

    with tempfile.TemporaryDirectory(prefix="namics_ref_") as tmp_ref_dir:
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
            reference_archive=reference_archive,
            reference_unpack_dir=Path(tmp_ref_dir),
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
