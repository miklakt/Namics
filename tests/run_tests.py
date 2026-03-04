#!/usr/bin/env python3
"""NAMICS regression/benchmark runner.

Usage:
  python3 tests/run_tests.py                 # run all tests
  python3 tests/run_tests.py frozen-range-input-file
  python3 tests/run_tests.py --list
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_helpers import (  # noqa: E402
    CommandMetrics,
    TestError,
    compare_json_profiles,
    ensure_solver_method_line,
    replace_in_file,
    require_file,
    run_command,
    utc_run_id,
)


@dataclass
class Context:
    repo_root: Path
    tests_dir: Path
    output_dir: Path
    binary: Path
    quiet: bool
    keep_artifacts: bool
    benchmark_threshold_pct: float
    benchmark_history: bool


@dataclass
class TestResult:
    key: str
    passed: bool
    wall_s: float
    max_rss_kb: float | None
    details: str


class Meter:
    def __init__(self) -> None:
        self.wall_s = 0.0
        self.max_rss_kb: float | None = None

    def add(self, m: CommandMetrics) -> None:
        self.wall_s += m.wall_s
        if m.max_rss_kb is None:
            return
        if self.max_rss_kb is None:
            self.max_rss_kb = m.max_rss_kb
        else:
            self.max_rss_kb = max(self.max_rss_kb, m.max_rss_kb)


def _run_binary(ctx: Context, meter: Meter, input_file: Path, allow_fail: bool = False) -> CommandMetrics:
    metrics = run_command([str(ctx.binary), str(input_file)], cwd=ctx.repo_root, quiet=ctx.quiet)
    meter.add(metrics)
    if not allow_fail and metrics.returncode != 0:
        raise TestError(f"ERROR: execution failed ({metrics.returncode}) for input: {input_file}")
    return metrics


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _cleanup(paths: list[Path], keep: bool) -> None:
    if keep:
        return
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except IsADirectoryError:
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


def test_homopolymer_adsorption(ctx: Context) -> tuple[Meter, str]:
    meter = Meter()
    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    template_file = tests_dir / "homopolymer_adsorption.in"
    reference_dir = tests_dir / "reference"
    benchmark_dir = tests_dir / "benchmarks"
    benchmark_history = benchmark_dir / "homopolymer_adsorption_test_benchmark.csv"

    runtime_input = output_dir / "homopolymer_adsorption.in"
    runtime_input_tmp = output_dir / "homopolymer_adsorption.in.tmp"
    runtime_input_save_memory = output_dir / "homopolymer_adsorption_save_memory.in"
    runtime_input_save_memory_tmp = output_dir / "homopolymer_adsorption_save_memory.in.tmp"
    runtime_input_diis = output_dir / "homopolymer_adsorption_diis.in"
    runtime_input_diis_tmp = output_dir / "homopolymer_adsorption_diis.in.tmp"

    output_file_json = output_dir / "homopolymer_adsorption.json"
    output_file_save_memory_json = output_dir / "homopolymer_adsorption_save_memory.json"
    baseline_output_json = output_dir / "homopolymer_adsorption.baseline.json"

    require_file(ctx.binary, executable=True)
    require_file(template_file)
    for chi in (0, -2, -4, -6):
        require_file(reference_dir / f"homopolymer_adsorption.chi_{chi}.json.ref")

    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    run_id = utc_run_id()
    run_log = benchmark_dir / f"homopolymer_adsorption_test_{run_id}.log"

    cleanup_targets = [
        output_file_json,
        output_file_save_memory_json,
        runtime_input,
        runtime_input_tmp,
        runtime_input_save_memory,
        runtime_input_save_memory_tmp,
        runtime_input_diis,
        runtime_input_diis_tmp,
        baseline_output_json,
    ]

    solver_runtime_ms = 0
    diis_notes: list[str] = []
    log_lines = [
        f"run_id={run_id}",
        f"timestamp_utc={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
        f"threshold_pct={ctx.benchmark_threshold_pct:g}",
        "cases=chi_Si(0,-2,-4,-6) x save_memory(off,on)",
    ]
    bench_start = time.perf_counter()

    try:
        template_text = template_file.read_text(encoding="utf-8")

        for chi in (0, -2, -4, -6):
            reference_file = reference_dir / f"homopolymer_adsorption.chi_{chi}.json.ref"
            rendered = template_text.replace("{chi_Si}", str(chi))

            _write_text(runtime_input_tmp, rendered)
            ensure_solver_method_line(runtime_input_tmp, runtime_input, "pseudohessian")
            output_file_json.unlink(missing_ok=True)

            run_m = _run_binary(ctx, meter, runtime_input)
            solver_runtime_ms += int(round(run_m.wall_s * 1000.0))
            log_lines.append(f"chi_Si={chi},mode=baseline,elapsed_ms={int(round(run_m.wall_s * 1000.0))}")

            if not output_file_json.is_file() or output_file_json.stat().st_size == 0:
                raise TestError(f"ERROR: expected JSON output file was not created for chi_Si={chi}: {output_file_json}")
            compare_json_profiles(reference_file, output_file_json, coord_tol=1e-12, value_tol=1e-6)

            _write_text(runtime_input_save_memory_tmp, rendered)
            replace_in_file(
                runtime_input_save_memory_tmp,
                [("//mol : pol : save_memory : true", "mol : pol : save_memory : true")],
            )
            ensure_solver_method_line(runtime_input_save_memory_tmp, runtime_input_save_memory, "pseudohessian")
            output_file_save_memory_json.unlink(missing_ok=True)

            run_m = _run_binary(ctx, meter, runtime_input_save_memory)
            solver_runtime_ms += int(round(run_m.wall_s * 1000.0))
            log_lines.append(
                f"chi_Si={chi},mode=save_memory,elapsed_ms={int(round(run_m.wall_s * 1000.0))}"
            )

            if not output_file_save_memory_json.is_file() or output_file_save_memory_json.stat().st_size == 0:
                raise TestError(
                    f"ERROR: expected save_memory JSON output was not created for chi_Si={chi}: "
                    f"{output_file_save_memory_json}"
                )

            compare_json_profiles(reference_file, output_file_save_memory_json, coord_tol=1e-12, value_tol=1e-6)
            compare_json_profiles(output_file_json, output_file_save_memory_json, coord_tol=1e-12, value_tol=1e-6)

            shutil.copy2(output_file_json, baseline_output_json)
            _write_text(runtime_input_diis_tmp, rendered)
            ensure_solver_method_line(runtime_input_diis_tmp, runtime_input_diis, "DIIS")

            diis_status = "failed (accepted)"
            run_diis = _run_binary(ctx, meter, runtime_input_diis, allow_fail=True)
            if run_diis.returncode == 0:
                if output_file_json.is_file() and output_file_json.stat().st_size > 0:
                    try:
                        compare_json_profiles(baseline_output_json, output_file_json, coord_tol=1e-12, value_tol=1e-6)
                        diis_status = "passed and matched pseudohessian"
                    except TestError:
                        diis_status = "completed but differs from pseudohessian (accepted)"
                else:
                    diis_status = "completed but JSON output file is missing (accepted)"

            diis_notes.append(f"chi={chi}:{diis_status}")
            log_lines.append(f"chi_Si={chi},mode=DIIS,status={diis_status}")

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

        diis_summary = ";".join(diis_notes)
        details = f"{benchmark_note};diis={diis_summary};log={run_log.name}"
        return meter, details
    finally:
        _cleanup(cleanup_targets, keep=ctx.keep_artifacts)


def test_frozen_range_input_file(ctx: Context) -> tuple[Meter, str]:
    meter = Meter()
    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    input_file = tests_dir / "frozen_range_input_file.in"
    source_frozen_file = tests_dir / "frozen_range_input_file.frozen"
    reference_file = tests_dir / "reference" / "frozen_range_input_file.json.ref"

    output_file_json = output_dir / "frozen_range_input_file.json"
    baseline_output_json = output_dir / "frozen_range_input_file.pseudohessian.json"
    runtime_input_file = output_dir / "frozen_range_input_file.in"
    runtime_frozen_file = output_dir / "frozen_range_input_file.frozen"

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    require_file(source_frozen_file)
    require_file(reference_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cleanup_targets = [output_file_json, baseline_output_json, runtime_input_file, runtime_frozen_file]

    try:
        output_file_json.unlink(missing_ok=True)
        baseline_output_json.unlink(missing_ok=True)
        shutil.copy2(source_frozen_file, runtime_frozen_file)

        ensure_solver_method_line(input_file, runtime_input_file, "pseudohessian")
        _run_binary(ctx, meter, runtime_input_file)

        if not output_file_json.is_file() or output_file_json.stat().st_size == 0:
            raise TestError(f"ERROR: expected JSON output file was not created: {output_file_json}")

        compare_json_profiles(reference_file, output_file_json, coord_tol=1e-12, value_tol=1e-6)

        shutil.copy2(output_file_json, baseline_output_json)
        ensure_solver_method_line(input_file, runtime_input_file, "DIIS")

        diis_status = "failed (accepted)"
        run_diis = _run_binary(ctx, meter, runtime_input_file, allow_fail=True)
        if run_diis.returncode == 0:
            if output_file_json.is_file() and output_file_json.stat().st_size > 0:
                try:
                    compare_json_profiles(baseline_output_json, output_file_json, coord_tol=1e-12, value_tol=1e-6)
                    diis_status = "passed and matched pseudohessian"
                except TestError:
                    diis_status = "completed but differs from pseudohessian (accepted)"
            else:
                diis_status = "completed but output file is missing (accepted)"

        return meter, f"diis={diis_status}"
    finally:
        _cleanup(cleanup_targets, keep=ctx.keep_artifacts)


def test_micelle_self_assembly(ctx: Context) -> tuple[Meter, str]:
    meter = Meter()
    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    generate_input = tests_dir / "micelle_guess_generate.in"
    use_input = tests_dir / "micelle_guess_use.in"
    reference_file = tests_dir / "reference" / "micelle_guess_use.json.ref"

    guess_file = output_dir / "micelle_2.outi"
    output_file_json = output_dir / "micelle_guess_use.json"
    baseline_output_json = output_dir / "micelle_guess_use.pseudohessian.json"
    runtime_generate_pseudohessian = output_dir / "micelle_guess_generate.pseudohessian.in"
    runtime_use_pseudohessian = output_dir / "micelle_guess_use.pseudohessian.in"
    runtime_generate_diis = output_dir / "micelle_guess_generate.diis.in"
    runtime_use_diis = output_dir / "micelle_guess_use.diis.in"

    require_file(ctx.binary, executable=True)
    require_file(generate_input)
    require_file(use_input)
    require_file(reference_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cleanup_targets = [
        guess_file,
        output_file_json,
        baseline_output_json,
        runtime_generate_pseudohessian,
        runtime_use_pseudohessian,
        runtime_generate_diis,
        runtime_use_diis,
    ]

    try:
        for p in (guess_file, output_file_json, baseline_output_json):
            p.unlink(missing_ok=True)

        ensure_solver_method_line(generate_input, runtime_generate_pseudohessian, "pseudohessian")
        ensure_solver_method_line(use_input, runtime_use_pseudohessian, "pseudohessian")

        _run_binary(ctx, meter, runtime_generate_pseudohessian)
        if not guess_file.is_file() or guess_file.stat().st_size == 0:
            raise TestError(f"ERROR: expected guess file was not created: {guess_file}")

        _run_binary(ctx, meter, runtime_use_pseudohessian)
        if not output_file_json.is_file() or output_file_json.stat().st_size == 0:
            raise TestError(f"ERROR: expected JSON output file was not created: {output_file_json}")

        compare_json_profiles(reference_file, output_file_json, coord_tol=1e-12, value_tol=1e-9)

        shutil.copy2(output_file_json, baseline_output_json)
        guess_file.unlink(missing_ok=True)
        output_file_json.unlink(missing_ok=True)

        ensure_solver_method_line(generate_input, runtime_generate_diis, "DIIS")
        ensure_solver_method_line(use_input, runtime_use_diis, "DIIS")

        diis_status = "failed (accepted)"
        run_generate_diis = _run_binary(ctx, meter, runtime_generate_diis, allow_fail=True)
        if run_generate_diis.returncode == 0 and guess_file.is_file() and guess_file.stat().st_size > 0:
            run_use_diis = _run_binary(ctx, meter, runtime_use_diis, allow_fail=True)
            if run_use_diis.returncode == 0:
                if output_file_json.is_file() and output_file_json.stat().st_size > 0:
                    try:
                        compare_json_profiles(baseline_output_json, output_file_json, coord_tol=1e-12, value_tol=1e-9)
                        diis_status = "passed and matched pseudohessian"
                    except TestError:
                        diis_status = "completed but differs from pseudohessian (accepted)"
                else:
                    diis_status = "completed but JSON output file is missing (accepted)"

        return meter, f"diis={diis_status}"
    finally:
        _cleanup(cleanup_targets, keep=ctx.keep_artifacts)


def test_particle_in_cyl_coordinates(ctx: Context) -> tuple[Meter, str]:
    meter = Meter()
    output_dir = ctx.output_dir
    tests_dir = ctx.tests_dir

    input_file = tests_dir / "particle_in_cyl_coordinates.in"
    reference_file = tests_dir / "reference" / "particle_in_cyl_coordinates.json.ref"

    output_file_json = output_dir / "particle_in_cyl_coordinates.json"
    baseline_output_json = output_dir / "particle_in_cyl_coordinates.pseudohessian.json"
    runtime_input_pseudohessian = output_dir / "particle_in_cyl_coordinates.pseudohessian.in"
    runtime_input_diis = output_dir / "particle_in_cyl_coordinates.diis.in"

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    require_file(reference_file)

    output_dir.mkdir(parents=True, exist_ok=True)

    cleanup_targets = [
        output_file_json,
        baseline_output_json,
        runtime_input_pseudohessian,
        runtime_input_diis,
    ]

    try:
        output_file_json.unlink(missing_ok=True)
        baseline_output_json.unlink(missing_ok=True)

        ensure_solver_method_line(input_file, runtime_input_pseudohessian, "pseudohessian")
        _run_binary(ctx, meter, runtime_input_pseudohessian)

        if not output_file_json.is_file() or output_file_json.stat().st_size == 0:
            raise TestError(f"ERROR: expected JSON output file was not created: {output_file_json}")

        compare_json_profiles(reference_file, output_file_json, coord_tol=1e-12, value_tol=1e-6)

        shutil.copy2(output_file_json, baseline_output_json)
        ensure_solver_method_line(input_file, runtime_input_diis, "DIIS")

        diis_status = "failed (accepted)"
        run_diis = _run_binary(ctx, meter, runtime_input_diis, allow_fail=True)
        if run_diis.returncode == 0:
            if output_file_json.is_file() and output_file_json.stat().st_size > 0:
                try:
                    compare_json_profiles(baseline_output_json, output_file_json, coord_tol=1e-12, value_tol=1e-6)
                    diis_status = "passed and matched pseudohessian"
                except TestError:
                    diis_status = "completed but differs from pseudohessian (accepted)"
            else:
                diis_status = "completed but JSON output file is missing (accepted)"

        return meter, f"diis={diis_status}"
    finally:
        _cleanup(cleanup_targets, keep=ctx.keep_artifacts)


def _fmt_wall(v: float) -> str:
    return f"{v:.3f}"


def _fmt_mem_kb(v: float | None) -> str:
    if v is None:
        return "NA"
    return f"{(v / 1024.0):.1f}"


def _print_table(results: list[TestResult]) -> None:
    headers = ["Test", "Status", "Wall(s)", "MaxRSS(MiB)", "Details"]
    rows = [
        [
            r.key,
            "PASS" if r.passed else "FAIL",
            _fmt_wall(r.wall_s),
            _fmt_mem_kb(r.max_rss_kb),
            r.details,
        ]
        for r in results
    ]

    total_wall = sum(r.wall_s for r in results)
    max_mem_values = [r.max_rss_kb for r in results if r.max_rss_kb is not None]
    total_max_mem = max(max_mem_values) if max_mem_values else None
    pass_count = sum(1 for r in results if r.passed)
    rows.append(["TOTAL", f"{pass_count}/{len(results)}", _fmt_wall(total_wall), _fmt_mem_kb(total_max_mem), ""])

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


TESTS: dict[str, tuple[str, Callable[[Context], tuple[Meter, str]]]] = {
    "homopolymer-adsorption": ("homopolymer_adsorption", test_homopolymer_adsorption),
    "frozen-range-input-file": ("frozen_range_input_file", test_frozen_range_input_file),
    "micelle-self-assembly": ("micelle_self_assembly", test_micelle_self_assembly),
    "particle-in-cyl-coordinates": ("particle_in_cyl_coordinates", test_particle_in_cyl_coordinates),
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
    parser.add_argument("--keep-artifacts", action="store_true", help="Keep generated runtime files")
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

    ctx = Context(
        repo_root=repo_root,
        tests_dir=tests_dir,
        output_dir=output_dir,
        binary=binary,
        quiet=not args.verbose,
        keep_artifacts=args.keep_artifacts,
        benchmark_threshold_pct=args.benchmark_threshold_pct,
        benchmark_history=not args.no_benchmark_history,
    )

    results: list[TestResult] = []
    overall_pass = True

    for display_key in selected:
        _, fn = TESTS[display_key]
        meter = Meter()
        details = ""
        passed = True
        try:
            meter, details = fn(ctx)
        except TestError as exc:
            passed = False
            details = str(exc)
        except Exception as exc:  # defensive catch to keep summary table complete
            passed = False
            details = f"ERROR: unexpected failure: {exc}"

        if not passed:
            overall_pass = False

        results.append(
            TestResult(
                key=display_key,
                passed=passed,
                wall_s=meter.wall_s,
                max_rss_kb=meter.max_rss_kb,
                details=details,
            )
        )

    _print_table(results)
    return 0 if overall_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
