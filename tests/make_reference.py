#!/usr/bin/env python3
"""Write tests/reference/<case>.npz.

usage: make_reference.py [--json] <namics binary> [case ...]

The input is prepared exactly as run_tests.py prepares it (first solver method).
Default: the build writes .pro profiles, i.e. the fork point 1563590. The input is adapted to it: 'json' selectors
become 'pro' selectors, and the initial-guess keys added after the fork are dropped, since they only change where
the solver starts.
--json: the build writes JSON output (a current build); the input runs unchanged. Use it only where the fork point
is known to be wrong, and record why in tests/README.md.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from run_tests import CASES, TESTS, numeric_arrays, prepare

POST_FORK = re.compile(r"\s*(output\s*:|sys\s*:[^:]*:\s*(initial_guess|guess_inputfile|write_initial_guess|guess_outputfile)\s*:)")


def fork_point_input(text: str) -> str:
    lines = [re.sub(r"^\s*json\s*:", "pro :", line) for line in text.splitlines() if not POST_FORK.match(line)]
    return "\n".join(["output : pro : append : false", "output : pro : header_separator : classic", *lines]) + "\n"


def fork_point_arrays(path: Path, binary: Path, case) -> dict[str, np.ndarray]:
    path.write_text(fork_point_input(path.read_text(encoding="utf-8")), encoding="utf-8")
    for name in case.files:  # the fork point prefixes file names given in an input with the input's name
        shutil.copy2(path.parent / name, path.parent / f"{path.stem}.{name}")
    (path.parent / "output").mkdir()
    subprocess.run([str(binary), path.name], cwd=path.parent, capture_output=True)
    header, *rows = (path.parent / "output" / f"{path.stem}.pro").read_text(encoding="utf-8").splitlines()
    table = np.array([[float(cell) for cell in row.split()] for row in rows if row.strip()])
    return dict(zip(header.split(), table.T))


def json_arrays(path: Path, binary: Path, case) -> dict[str, np.ndarray]:
    if case.before:
        subprocess.run([str(binary), case.before], cwd=path.parent, capture_output=True, check=True)
    subprocess.run([str(binary), path.name], cwd=path.parent, capture_output=True, check=True)
    return numeric_arrays(json.loads(path.with_suffix(".output.json").read_text(encoding="utf-8"))["problems"][-1])


def main() -> int:
    args = sys.argv[1:]
    produce = json_arrays if args[0] == "--json" else fork_point_arrays
    args = args[1:] if args[0] == "--json" else args
    binary = Path(args[0]).resolve()
    for case in [case for case in CASES if len(args) < 2 or case.name in args[1:]]:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                arrays = produce(prepare(case, case.methods[0], Path(tmp) / "run"), binary, case)
            except (OSError, ValueError, subprocess.CalledProcessError) as exc:
                print(f"{case.name}: failed: {exc}")
                return 1
            np.savez_compressed(TESTS / "reference" / f"{case.name}.npz", **arrays)
            print(f"{case.name}: {len(arrays)} quantities")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
