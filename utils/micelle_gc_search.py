#!/usr/bin/env python3
"""External grand-canonical micelle search driver for NAMICS."""

import argparse
import concurrent.futures
import json
import math
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
X_TOL = 1e-3
EXPAND_FACTOR = 1.5
MAX_EXPAND = 8


@dataclass
class EvalResult:
    x: float
    grand_potential: float
    residual: float
    output_json: Path
    runtime_input: Path
    wall_s: float
    converged: bool


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Search the micelle aggregation number where grand_potential = 0.")
    parser.add_argument("input", type=Path, help="Single-problem NAMICS input template.")
    parser.add_argument("--binary", type=Path, default=REPO_ROOT / "bin" / "namics", help="NAMICS binary.")
    parser.add_argument("--molecule", required=True, help="Restricted molecule to vary.")
    parser.add_argument("--quantity", default="n", choices=("n", "theta", "phibulk"), help="Molecule quantity.")
    parser.add_argument("--initial", type=float, default=None, help="Initial search value.")
    parser.add_argument("--seed-json", type=Path, default=None, help="JSON file with embedded initial_guess.")
    parser.add_argument("--seed-input", type=Path, default=None, help="Input used once to generate the first seed, copied as-is.")
    parser.add_argument("--step", type=float, default=None, help="Initial bracketing step.")
    parser.add_argument("--workers", type=int, default=1, help="Parallel workers for bracket probes.")
    parser.add_argument("--max-iter", type=int, default=10, help="Maximum polishing iterations.")
    parser.add_argument("--gp-tol", type=float, default=1e-2, help="Convergence tolerance on grand potential.")
    args = parser.parse_args()
    if args.seed_json is not None and args.seed_input is not None:
        parser.error("use either --seed-json or --seed-input, not both")
    return args


def _load_template(path: Path, molecule: str, quantity: str) -> tuple[str, str, float | None]:
    if not path.is_file():
        raise SystemExit(f"input template not found: {path}")

    lines = path.read_text(encoding="utf-8").splitlines()
    meaningful: list[tuple[int, str]] = []
    system = "noname"
    initial = None

    for index, line in enumerate(lines):
        compact = line.split("//", 1)[0].replace(" ", "").replace("\t", "").strip()
        if not compact:
            continue
        meaningful.append((index, compact))

        parts = compact.split(":")
        if len(parts) == 4 and parts[0].lower() == "sys" and system == "noname":
            system = parts[1]
        if len(parts) == 4 and parts[0].lower() == "mol" and parts[1].lower() == molecule.lower():
            if parts[2].lower() == quantity.lower():
                try:
                    initial = float(parts[3])
                except ValueError as exc:
                    raise SystemExit(
                        f"could not parse initial value for mol:{molecule}:{quantity}: {parts[3]}"
                    ) from exc

    if not meaningful:
        raise SystemExit(f"{path} is empty")

    starts = [index for index, compact in meaningful if compact.lower() == "start"]
    if len(starts) > 1:
        raise SystemExit(f"{path} contains multiple 'start' markers. Use a single-problem template.")
    if starts:
        if starts[0] != meaningful[-1][0]:
            raise SystemExit(f"{path} contains a non-trailing 'start'. Only a final single-problem 'start' is allowed.")
        del lines[starts[0]]

    return "\n".join(lines).rstrip(), system, initial


def _run_namics(binary: Path, runtime_input: Path) -> tuple[str, int, float]:
    start = time.perf_counter()
    proc = subprocess.run(
        [str(binary.resolve()), str(runtime_input.resolve())],
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return proc.stdout, proc.returncode, time.perf_counter() - start


def _read_last_problem(path: Path) -> dict:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise ValueError(f"failed to parse JSON output {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"JSON output root is not an object: {path}")
    problems = payload.get("problems")
    if not isinstance(problems, list) or not problems or not isinstance(problems[-1], dict):
        raise ValueError(f"JSON output does not contain a last problem object: {path}")
    return problems[-1]


def _nested_value(node: dict, *path: str):
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return None
        node = node[key]
    return node


def _prepare_seed_from_input(binary: Path, seed_input: Path, workdir: Path) -> Path:
    seed_input = seed_input.resolve()
    if not seed_input.is_file():
        raise SystemExit(f"seed input not found: {seed_input}")

    seed_dir = workdir / "seed"
    runtime_input = seed_dir / "seed.in"
    output_json = runtime_input.with_suffix(".output.json")
    runtime_input.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(seed_input, runtime_input)

    output, returncode, wall_s = _run_namics(binary, runtime_input)
    if returncode != 0:
        raise SystemExit(f"seed generation failed ({returncode})\n{output}")
    if not output_json.is_file():
        raise SystemExit(f"seed generation did not create {output_json}")
    if "initial_guess" not in _read_last_problem(output_json):
        raise SystemExit(f"seed generation output does not contain an embedded initial guess: {output_json}")
    print(f"seed generated in {wall_s:.2f}s -> {output_json}")
    return output_json


def _evaluate_once(
    template_text: str,
    binary: Path,
    system: str,
    molecule: str,
    quantity: str,
    worker_dir: Path,
    eval_index: int,
    x: float,
    seed_json: Path | None,
) -> EvalResult:
    output_stem = f"run_{eval_index:03d}"
    runtime_input = worker_dir / f"{output_stem}.in"
    output_json = runtime_input.with_suffix(".output.json")

    lines = [
        template_text.rstrip(),
        "output : json : append : false",
        "output : json : header_separator : _",
        f"output : json : filename : {output_stem}",
        f"json : sys : {system} : grand_potential",
        f"mol : {molecule} : {quantity} : {x:.15g}",
    ]
    if seed_json is None:
        lines.append(f"sys : {system} : initial_guess : none")
    else:
        lines.extend(
            [
                f"sys : {system} : initial_guess : file",
                f"sys : {system} : guess_inputfile : {seed_json.resolve()}",
            ]
        )

    runtime_input.parent.mkdir(parents=True, exist_ok=True)
    runtime_input.write_text("\n".join(lines) + "\n", encoding="utf-8")

    output, returncode, wall_s = _run_namics(binary, runtime_input)
    if returncode != 0:
        raise RuntimeError(f"NAMICS failed for x={x} ({returncode})\n{output}")
    if not output_json.is_file():
        raise RuntimeError(f"NAMICS did not create {output_json} for x={x}")

    problem = _read_last_problem(output_json)
    value = _nested_value(problem, "sys", system, "grand_potential")
    if value is None:
        raise ValueError(f"JSON output does not contain 'sys.{system}.grand_potential': {output_json}")

    try:
        grand_potential = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid 'sys.{system}.grand_potential' value in {output_json}") from exc

    return EvalResult(
        x=x,
        grand_potential=grand_potential,
        residual=grand_potential,
        output_json=output_json,
        runtime_input=runtime_input,
        wall_s=wall_s,
        converged="Warning: iteration not solved." not in output,
    )


def _evaluate_batch(
    template_text: str,
    binary: Path,
    system: str,
    molecule: str,
    quantity: str,
    workdir: Path,
    next_eval_index: int,
    xs: list[float],
    seed_json: Path | None,
    workers: int,
) -> tuple[list[EvalResult], int]:
    if not xs:
        return [], next_eval_index

    jobs = [(next_eval_index + offset, x, workdir / f"worker_{offset % max(1, workers):02d}") for offset, x in enumerate(xs)]
    next_eval_index += len(jobs)

    if len(jobs) == 1 or workers <= 1:
        return [
            _evaluate_once(template_text, binary, system, molecule, quantity, worker_dir, eval_index, x, seed_json)
            for eval_index, x, worker_dir in jobs
        ], next_eval_index

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _evaluate_once,
                template_text,
                binary,
                system,
                molecule,
                quantity,
                worker_dir,
                eval_index,
                x,
                seed_json,
            )
            for eval_index, x, worker_dir in jobs
        ]
        return [future.result() for future in futures], next_eval_index


def _write_summary(
    summary_file: Path,
    result_input: Path,
    result_json: Path,
    seed_json: Path | None,
    best: EvalResult,
    left: EvalResult,
    right: EvalResult,
    evaluations: int,
    converged: bool,
) -> None:
    summary_file.write_text(
        json.dumps(
            {
                "target_grand_potential": 0.0,
                "x": best.x,
                "grand_potential": best.grand_potential,
                "residual": best.residual,
                "converged": converged,
                "evaluations": evaluations,
                "search_bracket": {
                    "lower_x": left.x,
                    "lower_residual": left.residual,
                    "upper_x": right.x,
                    "upper_residual": right.residual,
                },
                "seed_json": str(seed_json) if seed_json is not None else None,
                "result_input": str(result_input.resolve()),
                "result_json": str(result_json.resolve()),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = _parse_args()
    if not args.binary.is_file():
        raise SystemExit(f"binary not found: {args.binary}")

    template_text, system, template_initial = _load_template(args.input.resolve(), args.molecule, args.quantity)
    if args.initial is None and template_initial is None:
        raise SystemExit(f"did not find 'mol:{args.molecule}:{args.quantity}:<value>' in the input template")

    workdir = (REPO_ROOT / "output" / args.input.stem).resolve()
    history_file = workdir / "history.tsv"
    result_input = workdir / "result.in"
    result_json = workdir / "result.output.json"
    summary_file = workdir / "summary.json"
    workdir.mkdir(parents=True, exist_ok=True)

    seed_json = None
    if args.seed_json is not None:
        seed_json = args.seed_json.resolve()
        print(f"seed -> {seed_json}")
    elif args.seed_input is not None:
        seed_json = _prepare_seed_from_input(args.binary, args.seed_input, workdir)

    initial = args.initial if args.initial is not None else template_initial
    step = args.step if args.step is not None else max(1.0, 0.05 * abs(initial))
    evaluations = 0
    next_eval_index = 0

    def find_bracket(results: list[EvalResult]) -> tuple[EvalResult, EvalResult] | None:
        ordered = sorted(results, key=lambda result: result.x)
        for left, right in zip(ordered, ordered[1:]):
            if left.residual == 0 or right.residual == 0:
                return left, right
            if math.copysign(1.0, left.residual) != math.copysign(1.0, right.residual):
                return left, right
        return None

    def evaluate(xs: list[float], seed_json: Path | None) -> list[EvalResult]:
        nonlocal evaluations, next_eval_index
        results, next_eval_index = _evaluate_batch(
            template_text,
            args.binary,
            system,
            args.molecule,
            args.quantity,
            workdir,
            next_eval_index,
            xs,
            seed_json,
            args.workers,
        )
        if not history_file.exists():
            history_file.write_text("x\tgrand_potential\tresidual\twall_s\tconverged\toutput_json\n", encoding="utf-8")
        with history_file.open("a", encoding="utf-8") as handle:
            for offset, result in enumerate(results):
                handle.write(
                    "\t".join(
                        [
                            f"{result.x:.15g}",
                            f"{result.grand_potential:.15g}",
                            f"{result.residual:.15g}",
                            f"{result.wall_s:.6f}",
                            "true" if result.converged else "false",
                            str(result.output_json.resolve()),
                        ]
                    )
                    + "\n"
                )
                print(
                    f"eval {evaluations + offset + 1}: "
                    f"{args.quantity}={result.x:.15g} "
                    f"gp={result.grand_potential:.15g} "
                    f"res={result.residual:.15g} "
                    f"{'ok' if result.converged else 'warn'} {result.wall_s:.2f}s"
                )
        evaluations += len(results)
        return results

    center = evaluate([initial], seed_json)[0]
    if abs(center.residual) <= args.gp_tol:
        shutil.copy2(center.runtime_input, result_input)
        shutil.copy2(center.output_json, result_json)
        _write_summary(summary_file, result_input, result_json, seed_json, center, center, center, evaluations, True)
        print(f"converged immediately at {args.quantity}={center.x:.15g}")
        return 0

    best = center
    left = center
    right = center

    for expand_round in range(MAX_EXPAND):
        delta = step * (EXPAND_FACTOR ** expand_round)
        trial_xs = []
        left_x = max(0.0, best.x - delta)
        right_x = best.x + delta
        if left_x < best.x - X_TOL:
            trial_xs.append(left_x)
        if right_x > best.x + X_TOL:
            trial_xs.append(right_x)

        trial_results = evaluate(trial_xs, seed_json)
        bracket = find_bracket([best, *trial_results])
        if bracket is not None:
            left, right = bracket
            best = min((left, right), key=lambda result: abs(result.residual))
            break

        if not trial_results:
            break
        best = min((best, *trial_results), key=lambda result: abs(result.residual))
    else:
        raise SystemExit("failed to bracket the grand-potential root")

    if left.x > right.x:
        left, right = right, left
    best = min((left, right), key=lambda result: abs(result.residual))
    converged = abs(best.residual) <= args.gp_tol

    for iteration in range(args.max_iter):
        interval = right.x - left.x
        if abs(best.residual) <= args.gp_tol or interval <= X_TOL:
            converged = True
            break

        if abs(right.residual - left.residual) < 1e-16:
            candidate_x = 0.5 * (left.x + right.x)
        else:
            candidate_x = (left.x * right.residual - right.x * left.residual) / (right.residual - left.residual)
        midpoint = 0.5 * (left.x + right.x)
        if not (left.x < candidate_x < right.x) or min(candidate_x - left.x, right.x - candidate_x) < 0.1 * interval:
            candidate_x = midpoint

        candidate = evaluate([candidate_x], seed_json)[0]
        best = min((best, candidate), key=lambda result: abs(result.residual))

        if candidate.residual == 0 or math.copysign(1.0, left.residual) != math.copysign(1.0, candidate.residual):
            right = candidate
        else:
            left = candidate
        if left.x > right.x:
            left, right = right, left

        print(
            f"iter {iteration + 1}: bracket=[{left.x:.15g}, {right.x:.15g}] "
            f"best={best.x:.15g} res={best.residual:.15g}"
        )
        converged = abs(best.residual) <= args.gp_tol

    shutil.copy2(best.runtime_input, result_input)
    shutil.copy2(best.output_json, result_json)
    _write_summary(summary_file, result_input, result_json, seed_json, best, left, right, evaluations, converged)

    if converged:
        print(f"converged: {args.quantity}={best.x:.15g} gp={best.grand_potential:.15g}")
        return 0

    print(f"did not converge: best {args.quantity}={best.x:.15g} gp={best.grand_potential:.15g}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
