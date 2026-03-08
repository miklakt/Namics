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

_START_RE = re.compile(r"^[ \t]*start[ \t]*$", re.IGNORECASE)
_LINE_COMMENT_RE = re.compile(r"^[ \t]*//[ \t]*")


class TestError(RuntimeError):
    """Raised when a test expectation fails."""


@dataclass
class CommandMetrics:
    returncode: int
    wall_s: float
    max_rss_kb: float | None


def _normalize_key(key: str) -> str:
    return ":".join(part.strip().lower() for part in key.split(":"))


def _split_key_value(line: str) -> tuple[str, str] | None:
    raw = line.strip()
    if not raw:
        return None
    uncommented = _LINE_COMMENT_RE.sub("", raw, count=1)
    if ":" not in uncommented:
        return None
    parts = [part.strip() for part in uncommented.split(":")]
    if len(parts) < 2:
        return None
    key = ":".join(part.lower() for part in parts[:-1])
    value = parts[-1]
    return key, value


def _find_first_start(lines: list[str]) -> int | None:
    for idx, line in enumerate(lines):
        if _START_RE.match(line):
            return idx
    return None


def set_setting_line(input_file: Path, output_file: Path, setting_key: str, value: str) -> None:
    """Render INPUT to OUTPUT while forcing 'SETTING_KEY : VALUE'."""
    lines = input_file.read_text(encoding="utf-8").splitlines()
    target_key = _normalize_key(setting_key)
    rendered = f"{setting_key} : {value}"

    out_lines: list[str] = []
    replaced = False
    first_start = _find_first_start(lines)

    for idx, line in enumerate(lines):
        split = _split_key_value(line)
        if split is not None and split[0] == target_key:
            if not replaced:
                out_lines.append(rendered)
                replaced = True
            continue

        if not replaced and first_start is not None and idx == first_start:
            out_lines.append(rendered)
            replaced = True

        out_lines.append(line)

    if not replaced:
        out_lines.append(rendered)

    output_file.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def set_commented_setting(
    input_file: Path,
    output_file: Path,
    setting_key: str,
    value: str,
    enabled: bool,
) -> None:
    """Toggle a commented setting line while preserving runnable output files."""
    lines = input_file.read_text(encoding="utf-8").splitlines()
    target_key = _normalize_key(setting_key)
    active_line = f"{setting_key} : {value}"
    rendered = active_line if enabled else f"//{active_line}"

    out_lines: list[str] = []
    replaced = False
    first_start = _find_first_start(lines)

    for idx, line in enumerate(lines):
        raw = line.strip()
        uncommented = _LINE_COMMENT_RE.sub("", raw, count=1).strip()
        split = _split_key_value(line)
        if (split is not None and split[0] == target_key) or (
            uncommented and _normalize_key(uncommented) == target_key
        ):
            if not replaced:
                out_lines.append(rendered)
                replaced = True
            continue

        if not replaced and first_start is not None and idx == first_start:
            out_lines.append(rendered)
            replaced = True

        out_lines.append(line)

    if not replaced:
        out_lines.append(rendered)

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


def utc_run_id() -> str:
    ns = time.time_ns()
    sec = ns // 1_000_000_000
    ms = (ns // 1_000_000) % 1000
    return time.strftime("%Y%m%dT%H%M%S", time.gmtime(sec)) + f"{ms:03d}Z"
