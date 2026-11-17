#!/usr/bin/env python3
"""NAMICS regression suite.

Every case runs its input with each solver method and compares all numeric output against
tests/reference/<case>.npz: the same arrays, the same sizes, and max |diff| <= TOL everywhere.
See tests/README.md for where the references come from and how to add a case.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
TOL = 1e-8
SOLVER_TOLERANCE = "1e-10"


@dataclass(frozen=True)
class Case:
    name: str                                  # reference file: reference/<name>.npz
    input: str                                 # input file in tests/
    settings: tuple[str, ...] = ()             # "section : name : key : value" lines overriding the input
    files: tuple[str, ...] = ()                # extra files copied next to the input
    before: str | None = None                  # input run first in the same folder
    methods: tuple[str, ...] = ("pseudohessian", "DIIS")


CASES = [
    *[Case(f"homopolymer_adsorption.chi_{chi}", "homopolymer_adsorption.in", (f"mon : A : chi_Si : {chi}",)) for chi in (0, -2, -4, -6)],
    *[Case(f"external_potential_{d}", f"external_potential_{d}.in", files=(f"external_potential_{d}_external_potential.json",)) for d in ("1d", "2d", "3d")],
    Case("frozen_range_input_file", "frozen_range_input_file.in", files=("frozen_range_input_file.frozen",)),
    Case("micelle_guess_use", "micelle_guess_use.in", before="micelle_guess_generate.in"),
    Case("particle_in_cyl_coordinates", "particle_in_cyl_coordinates.in"),
    Case("branched_brush_planar_1d", "branched_brush_planar_1d.in"),
    Case("branched_brush_cyl2d", "branched_brush_cyl2d.in", methods=("DIIS",)),
    Case("polE", "polE.in"),
    Case("poor_solvent_brush_theta", "poor_solvent_brush_theta.in"),
    # poor solvent: pseudohessian and DIIS converge to different stable brush states (in 1563590 as well)
    Case("poor_solvent_brush_poor", "poor_solvent_brush_poor.in", methods=("pseudohessian",)),
    Case("charged_pore_salt", "charged_pore_salt.in"),
    Case("charged_pore_particle", "charged_pore_particle.in"),
    Case("angular_pore", "angular_pore.in", files=("angular_pore.frozen",)),
    Case("charged_sphere_brush", "charged_sphere_brush.in"),
    Case("sfb_empty_pore", "sfb_empty_pore.in"),
]


def _key(line: str) -> str | None:
    parts = [part.strip() for part in line.split("//")[0].split(":")]
    return ":".join(parts[:3]).lower() if len(parts) == 4 and all(parts) else None


def with_settings(text: str, settings: list[str]) -> str:
    """Drop input lines for the keys in SETTINGS and put SETTINGS before the first 'start' (or at the end)."""
    keys = {_key(setting) for setting in settings}
    lines = [line for line in text.splitlines() if _key(line) not in keys]
    start = next((i for i, line in enumerate(lines) if line.strip() == "start"), len(lines))
    return "\n".join(lines[:start] + settings + lines[start:]) + "\n"


def prepare(case: Case, method: str, work: Path) -> Path:
    """Write the runnable input for CASE and METHOD into WORK and return its path."""
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    for name in case.files:
        shutil.copy2(TESTS / name, work / name)
    settings = [*case.settings, f"newton : isaac : method : {method}", f"newton : isaac : tolerance : {SOLVER_TOLERANCE}"]
    for name in filter(None, (case.before, case.input)):
        (work / name).write_text(with_settings((TESTS / name).read_text(encoding="utf-8"), settings), encoding="utf-8")
    return work / case.input


def numeric_arrays(node: object, path: str = "") -> dict[str, np.ndarray]:
    """All numbers in a JSON tree as flat float64 arrays keyed by their ':'-joined path."""
    if isinstance(node, dict):
        return {k: v for key, child in node.items() for k, v in numeric_arrays(child, f"{path}:{key}" if path else key).items()}
    values = np.asarray(node, dtype=object).ravel()
    if values.size and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values):
        return {path: values.astype(np.float64)}
    return {}


def compare(got: dict[str, np.ndarray], reference: Path) -> str:
    with np.load(reference) as stored:
        expected = {name: stored[name] for name in stored.files}
    if got.keys() != expected.keys():
        raise AssertionError(f"quantities differ: missing {sorted(expected.keys() - got.keys())}, extra {sorted(got.keys() - expected.keys())}")
    worst = 0.0
    for name, values in expected.items():
        if got[name].shape != values.shape:
            raise AssertionError(f"{name}: {got[name].size} values, reference has {values.size}")
        diff = float(np.max(np.abs(got[name] - values))) if values.size else 0.0
        if not diff <= TOL:
            raise AssertionError(f"{name}: max |diff| {diff:.2e} > {TOL:.0e}")
        worst = max(worst, diff)
    return f"max |diff| {worst:.1e}"


def run(binary: Path, path: Path, verbose: bool) -> None:
    result = subprocess.run([str(binary), path.name], cwd=path.parent, capture_output=True, text=True)
    if verbose:
        print(result.stdout + result.stderr)
    if result.returncode:
        tail = (result.stdout + result.stderr).strip().splitlines()[-1:] or [""]
        raise AssertionError(f"{path.name} exited with {result.returncode}: {tail[0]}")


def run_case(case: Case, method: str, binary: Path, verbose: bool, work: Path) -> str:
    path = prepare(case, method, work)
    if case.before:
        run(binary, path.with_name(case.before), verbose)
    run(binary, path, verbose)
    problem = json.loads(path.with_suffix(".output.json").read_text(encoding="utf-8"))["problems"][-1]
    return compare(numeric_arrays(problem), TESTS / "reference" / f"{case.name}.npz")


def micelle_zero_grand_potential(binary: Path, work: Path) -> str:
    """Root search with utils/search_engine: the grand potential of the micelle must vanish."""
    sys.path.insert(0, str(REPO))
    from utils.search_engine import solve

    def micelle_shape(output: Path) -> str | None:
        phi = json.loads(output.read_text(encoding="utf-8"))["problems"][-1]["mol"]["surf"]["phi"]
        if not (min(phi) > 0 and phi[0] > max(0.3, phi[-1]) and max(phi) < 1.25 * phi[0]):
            return "mol.surf.phi does not describe a micelle centred at the origin"
        return None

    shutil.rmtree(work, ignore_errors=True)
    with contextlib.redirect_stdout(io.StringIO()):
        result = solve(TESTS / "micelle_zero_grand_potential_seed.in", TESTS / "micelle_zero_grand_potential_run.in",
                       "mol.surf.n", "sys.noname.grand_potential", method="secant", X0=160, x_min=40, x_max=300,
                       max_expand_step=10.0, maxiter=100, output_validation=micelle_shape, binary=binary,
                       workdir=work, return_result=True)
    if abs(result.fx) > 1e-4:
        raise AssertionError(f"grand potential stayed at {result.fx:.3e}")
    return f"n = {result.root:.6g}, grand potential {result.fx:.1e}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("filters", nargs="*", help="run only cases whose name contains one of these")
    parser.add_argument("--binary", type=Path, default=REPO / "bin" / "namics")
    parser.add_argument("--list", action="store_true", help="list the cases and exit")
    parser.add_argument("--verbose", action="store_true", help="print solver output")
    args = parser.parse_args()

    binary = args.binary.resolve()
    jobs = [(case.name, method, partial(run_case, case, method, binary, args.verbose)) for case in CASES for method in case.methods]
    jobs.append(("micelle_zero_grand_potential", "search", partial(micelle_zero_grand_potential, binary)))
    jobs = [job for job in jobs if not args.filters or any(f in job[0] for f in args.filters)]
    if args.list:
        print("\n".join(f"{name} {method}" for name, method, _ in jobs))
        return 0

    output = REPO / "output" / "tests"
    failed = 0
    for name, method, job in jobs:
        start = time.perf_counter()
        try:
            detail, status = job(output / f"{name}.{method}"), "PASS"
        except Exception as exc:  # report every failure and keep going
            detail, status = f"{type(exc).__name__}: {exc}" if not isinstance(exc, AssertionError) else str(exc), "FAIL"
            failed += 1
        print(f"{status}  {name:34s} {method:14s} {time.perf_counter() - start:7.2f} s  {detail}", flush=True)
    print(f"{len(jobs) - failed}/{len(jobs)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
