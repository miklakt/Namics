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
    """Compare nested JSON profiles without flattening the object tree."""
    left_data = _as_problem_object(_read_json(left_path), left_path)
    right_data = _as_problem_object(_read_json(right_path), right_path)
    left_data = {key: value for key, value in left_data.items() if key != "initial_guess"}
    right_data = {key: value for key, value in right_data.items() if key != "initial_guess"}
    path_key = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")

    def render_path(path: list[str]) -> str:
        out = ""
        for part in path:
            if part.isdigit():
                out += f"[{part}]"
            elif path_key.match(part):
                out += f".{part}" if out else part
            else:
                out += f"[{json.dumps(part)}]"
        return out or "<root>"

    def numeric_list(values: list[Any]) -> list[float] | None:
        out: list[float] = []
        for item in values:
            if isinstance(item, (dict, list, tuple)):
                return None
            try:
                out.append(float(item))
            except (TypeError, ValueError):
                return None
        return out

    def meaningful(value: Any) -> bool:
        if isinstance(value, dict):
            return any(meaningful(child) for child in value.values())
        if isinstance(value, list):
            return numeric_list(value) is not None or any(
                meaningful(child) for child in value if isinstance(child, (dict, list))
            )
        return False

    def walk(left: Any, right: Any, path: list[str], label: str | None) -> None:
        if isinstance(left, dict) and isinstance(right, dict):
            for key in left.keys() | right.keys():
                left_value = left.get(key)
                right_value = right.get(key)
                left_meaningful = meaningful(left_value)
                right_meaningful = meaningful(right_value)
                if not left_meaningful and not right_meaningful:
                    continue
                if not left_meaningful or not right_meaningful:
                    raise TestError(f"ERROR: structure mismatch at {render_path(path + [key])}")
                walk(left_value, right_value, path + [key], key)
            return
        if isinstance(left, list) and isinstance(right, list):
            left_numbers = numeric_list(left)
            right_numbers = numeric_list(right)
            if left_numbers is not None and right_numbers is not None:
                if len(left_numbers) != len(right_numbers):
                    raise TestError(
                        f"ERROR: row count mismatch for column {render_path(path)} "
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
                            f"ERROR: numerical drift at {render_path(path)} "
                            f"(max diff {max_diff:.3e} > tol {tol:.3e})"
                        )
                return
            left_meaningful = meaningful(left)
            right_meaningful = meaningful(right)
            if not left_meaningful and not right_meaningful:
                return
            if not left_meaningful or not right_meaningful:
                raise TestError(f"ERROR: structure mismatch at {render_path(path)}")
            if any(isinstance(item, (dict, list)) for item in left) or any(
                isinstance(item, (dict, list)) for item in right
            ):
                if len(left) != len(right):
                    raise TestError(
                        f"ERROR: row count mismatch for column {render_path(path)} "
                        f"({left_path}: {len(left)}, {right_path}: {len(right)})"
                    )
                for idx, (left_item, right_item) in enumerate(zip(left, right)):
                    if not meaningful(left_item) and not meaningful(right_item):
                        continue
                    if not meaningful(left_item) or not meaningful(right_item):
                        raise TestError(f"ERROR: structure mismatch at {render_path(path + [str(idx)])}")
                    walk(left_item, right_item, path + [str(idx)], label)

    walk(left_data, right_data, [], None)


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
