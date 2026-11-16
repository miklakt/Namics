#!/usr/bin/env python3
"""Write tests/reference/<case>.npz with a build that writes .pro profiles, i.e. the fork point 1563590.

usage: make_reference.py <namics binary> [case ...]

The input is prepared exactly as run_tests.py prepares it (first solver method), then adapted to that build:
'json' selectors become 'pro' selectors, and the initial-guess keys added after the fork are dropped, since they
only change where the solver starts.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

from run_tests import CASES, TESTS, prepare

POST_FORK = re.compile(r"\s*(output\s*:|sys\s*:[^:]*:\s*(initial_guess|guess_inputfile|write_initial_guess|guess_outputfile)\s*:)")


def fork_point_input(text: str) -> str:
    lines = [re.sub(r"^\s*json\s*:", "pro :", line) for line in text.splitlines() if not POST_FORK.match(line)]
    return "\n".join(["output : pro : append : false", "output : pro : header_separator : classic", *lines]) + "\n"


def main() -> int:
    binary = Path(sys.argv[1]).resolve()
    for case in [case for case in CASES if len(sys.argv) < 3 or case.name in sys.argv[2:]]:
        with tempfile.TemporaryDirectory() as tmp:
            path = prepare(case, case.methods[0], Path(tmp) / "run")
            path.write_text(fork_point_input(path.read_text(encoding="utf-8")), encoding="utf-8")
            for name in case.files:  # the fork point prefixes file names given in an input with the input's name
                shutil.copy2(path.parent / name, path.parent / f"{path.stem}.{name}")
            (path.parent / "output").mkdir()
            result = subprocess.run([str(binary), path.name], cwd=path.parent, capture_output=True, text=True)
            profiles = path.parent / "output" / f"{path.stem}.pro"
            if result.returncode or not profiles.is_file():
                print(f"{case.name}: failed\n{result.stdout[-2000:]}{result.stderr[-2000:]}")
                return 1
            header, *rows = profiles.read_text(encoding="utf-8").splitlines()
            table = np.array([[float(cell) for cell in row.split()] for row in rows if row.strip()])
            np.savez_compressed(TESTS / "reference" / f"{case.name}.npz", **dict(zip(header.split(), table.T)))
            print(f"{case.name}: {len(header.split())} quantities, {len(table)} sites")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
