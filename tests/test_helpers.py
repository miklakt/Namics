#!/usr/bin/env python3
"""Shared test helpers for NAMICS regression tests."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
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
    output: str = ""


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


def set_setting_line(input_file: Path, output_file: Path, setting_key: str, value: str) -> None:
    """Render INPUT to OUTPUT while forcing 'SETTING_KEY : VALUE'."""
    lines = input_file.read_text(encoding="utf-8").splitlines()
    target_key = _normalize_key(setting_key)
    rendered = f"{setting_key} : {value}"

    out_lines: list[str] = []
    replaced = False
    first_start = next((idx for idx, line in enumerate(lines) if _START_RE.match(line)), None)

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
    first_start = next((idx for idx, line in enumerate(lines) if _START_RE.match(line)), None)

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


def compare_json_profiles(left_path: Path, right_path: Path, coord_tol: float, value_tol: float) -> None:
    """Compare profile series while ignoring wrapper structure."""
    left_data = _as_problem_object(_read_json(left_path), left_path)
    right_data = _as_problem_object(_read_json(right_path), right_path)

    def numeric_list(values: list[Any]) -> list[float] | None:
        out: list[float] = []
        for item in values:
            if isinstance(item, dict):
                return None
            if isinstance(item, (list, tuple)):
                nested = numeric_list(list(item))
                if nested is None:
                    return None
                out.extend(nested)
                continue
            try:
                out.append(float(item))
            except (TypeError, ValueError):
                return None
        return out

    def collect_series(value: Any, path: list[str], out: dict[str, tuple[str, list[float]]]) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                if key == "initial_guess":
                    continue
                collect_series(child, path + [key], out)
            return
        if isinstance(value, list):
            numbers = numeric_list(value)
            if numbers is not None:
                if not path:
                    raise TestError("ERROR: encountered numeric JSON list at root")
                key = "_".join(path)
                out[key] = (path[-1], numbers)
                return
            for idx, child in enumerate(value):
                if isinstance(child, (dict, list)):
                    collect_series(child, path + [str(idx)], out)

    left_series: dict[str, tuple[str, list[float]]] = {}
    right_series: dict[str, tuple[str, list[float]]] = {}
    collect_series(left_data, [], left_series)
    collect_series(right_data, [], right_series)

    for key in left_series.keys() | right_series.keys():
        left_entry = left_series.get(key)
        right_entry = right_series.get(key)
        if left_entry is None or right_entry is None:
            raise TestError(f"ERROR: structure mismatch at {key}")
        label, left_numbers = left_entry
        _, right_numbers = right_entry
        if len(left_numbers) != len(right_numbers):
            raise TestError(
                f"ERROR: row count mismatch for column {key} "
                f"({left_path}: {len(left_numbers)}, {right_path}: {len(right_numbers)})"
            )
        tol = coord_tol if label in {"x", "y", "z"} else value_tol
        max_diff = 0.0
        for lv, rv in zip(left_numbers, right_numbers):
            diff = abs(lv - rv)
            if diff > max_diff:
                max_diff = diff
            if diff > tol:
                raise TestError(
                    f"ERROR: numerical drift at {key} "
                    f"(max diff {max_diff:.3e} > tol {tol:.3e})"
                )


def run_command(cmd: list[str], cwd: Path, quiet: bool = True) -> CommandMetrics:
    """Run command and return wall time + max RSS (where available)."""
    with tempfile.TemporaryFile() as stream:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, cwd=str(cwd), stdout=stream, stderr=stream)

        if hasattr(os, "wait4"):
            _, status, rusage = os.wait4(proc.pid, 0)
            wall_s = time.perf_counter() - start
            returncode = os.waitstatus_to_exitcode(status)
            max_rss_kb = float(rusage.ru_maxrss)
            if sys.platform == "darwin":
                max_rss_kb /= 1024.0
        else:
            returncode = proc.wait()
            wall_s = time.perf_counter() - start
            max_rss_kb = None

        stream.seek(0)
        output = stream.read().decode("utf-8", errors="replace")
        if not quiet and output:
            print(output, end="")

    return CommandMetrics(returncode=returncode, wall_s=wall_s, max_rss_kb=max_rss_kb, output=output)


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
