#!/usr/bin/env python3
"""Shared test helpers for NAMICS regression tests."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_METHOD_LINE_RE = re.compile(
    r"^[ \t]*newton[ \t]*:[ \t]*isaac[ \t]*:[ \t]*method[ \t]*:", re.IGNORECASE
)
_START_RE = re.compile(r"^[ \t]*start[ \t]*$", re.IGNORECASE)


class TestError(RuntimeError):
    """Raised when a test expectation fails."""


@dataclass
class CommandMetrics:
    returncode: int
    wall_s: float
    max_rss_kb: float | None


def ensure_solver_method_line(input_file: Path, output_file: Path, method: str) -> None:
    """Render INPUT to OUTPUT while forcing 'newton : isaac : method : METHOD'."""
    lines = input_file.read_text(encoding="utf-8").splitlines()
    has_method = any(_METHOD_LINE_RE.match(line) for line in lines)

    first_start: int | None = None
    for idx, line in enumerate(lines):
        if _START_RE.match(line):
            first_start = idx
            break

    out_lines: list[str] = []
    inserted = False

    for idx, line in enumerate(lines):
        if _METHOD_LINE_RE.match(line):
            out_lines.append(f"newton : isaac : method : {method}")
            continue
        if not has_method and not inserted and first_start is not None and idx == first_start:
            out_lines.append(f"newton : isaac : method : {method}")
            inserted = True
        out_lines.append(line)

    if not has_method and first_start is None:
        out_lines.append(f"newton : isaac : method : {method}")

    output_file.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise TestError(f"ERROR: cannot parse JSON {path}: {exc}") from exc


def _as_problem_object(data: Any, path: Path) -> dict[str, Any]:
    if isinstance(data, list):
        if not data:
            raise TestError(f"ERROR: JSON array is empty: {path}")
        data = data[-1]

    if not isinstance(data, dict):
        raise TestError(f"ERROR: JSON root is not an object: {path}")

    if "problems" in data:
        problems = data.get("problems")
        if not isinstance(problems, list) or not problems:
            raise TestError(f"ERROR: JSON problems section is malformed: {path}")
        data = problems[-1]
        if not isinstance(data, dict):
            raise TestError(f"ERROR: JSON problem entry is not an object: {path}")

    return data


def _extract_numeric_arrays(obj: dict[str, Any]) -> dict[str, list[float]]:
    arrays: dict[str, list[float]] = {}
    for key, value in obj.items():
        if not isinstance(value, list) or not value:
            continue
        numeric: list[float] = []
        ok = True
        for item in value:
            if isinstance(item, (dict, list, tuple)):
                ok = False
                break
            try:
                numeric.append(float(item))
            except (TypeError, ValueError):
                ok = False
                break
        if ok:
            arrays[str(key)] = numeric
    return arrays


def compare_json_profiles(left_path: Path, right_path: Path, coord_tol: float, value_tol: float) -> None:
    """Compare numeric array columns extracted from NAMICS JSON outputs."""
    left_data = _as_problem_object(_read_json(left_path), left_path)
    right_data = _as_problem_object(_read_json(right_path), right_path)
    left = _extract_numeric_arrays(left_data)
    right = _extract_numeric_arrays(right_data)

    for key, left_col in left.items():
        if key not in right:
            raise TestError(f"ERROR: output is missing column: {key}")
        right_col = right[key]
        if len(left_col) != len(right_col):
            raise TestError(
                f"ERROR: row count mismatch for column {key} "
                f"({left_path}: {len(left_col)}, {right_path}: {len(right_col)})"
            )
        tol = coord_tol if key in {"x", "y", "z"} else value_tol
        for i, (lv, rv) in enumerate(zip(left_col, right_col), start=1):
            if abs(lv - rv) > tol:
                raise TestError(
                    f"ERROR: mismatch at row {i}, column {key} "
                    f"(tol={tol}) left={lv} right={rv}"
                )

    for key in right:
        if key not in left:
            raise TestError(f"ERROR: reference is missing column: {key}")


def run_command(cmd: list[str], cwd: Path, quiet: bool = True) -> CommandMetrics:
    """Run command and return wall time + max RSS (where available)."""
    stdout = subprocess.DEVNULL if quiet else None
    stderr = subprocess.DEVNULL if quiet else None

    start = time.perf_counter()
    proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=stdout, stderr=stderr)

    if hasattr(os, "wait4"):
        _, status, rusage = os.wait4(proc.pid, 0)
        wall_s = time.perf_counter() - start
        returncode = os.waitstatus_to_exitcode(status)
        max_rss_kb = float(rusage.ru_maxrss)
        if sys.platform == "darwin":
            max_rss_kb /= 1024.0
        return CommandMetrics(returncode=returncode, wall_s=wall_s, max_rss_kb=max_rss_kb)

    returncode = proc.wait()
    wall_s = time.perf_counter() - start
    return CommandMetrics(returncode=returncode, wall_s=wall_s, max_rss_kb=None)


def require_file(path: Path, executable: bool = False) -> None:
    if executable:
        if not path.is_file() or not os.access(path, os.X_OK):
            raise TestError(f"ERROR: required executable not found: {path}")
        return
    if not path.is_file():
        raise TestError(f"ERROR: required file not found: {path}")


def replace_in_file(path: Path, replacements: list[tuple[str, str]]) -> None:
    text = path.read_text(encoding="utf-8")
    for old, new in replacements:
        text = text.replace(old, new)
    path.write_text(text, encoding="utf-8")


def utc_run_id() -> str:
    ns = time.time_ns()
    sec = ns // 1_000_000_000
    ms = (ns // 1_000_000) % 1000
    return time.strftime("%Y%m%dT%H%M%S", time.gmtime(sec)) + f"{ms:03d}Z"
