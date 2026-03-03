#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def load_table(path: Path):
    text = path.read_text(encoding="utf-8")
    stripped = text.lstrip()
    is_json = path.suffix.lower() == ".json" or stripped.startswith("{") or stripped.startswith("[")
    if is_json:
        data = json.loads(text)
        if isinstance(data, list):
            if len(data) == 0:
                raise ValueError(f"JSON array is empty: {path}")
            data = data[-1]
        if not isinstance(data, dict):
            raise ValueError(f"JSON root is not an object: {path}")
        if "problems" in data:
            problems = data.get("problems", [])
            if not isinstance(problems, list) or len(problems) == 0:
                raise ValueError(f"JSON problems section is malformed: {path}")
            data = problems[-1]
            if not isinstance(data, dict):
                raise ValueError(f"JSON problem entry is not an object: {path}")

            header = []
            columns = []
            expected_len = None
            for key, value in data.items():
                if not isinstance(value, list):
                    continue
                try:
                    numeric = [float(v) for v in value]
                except (TypeError, ValueError):
                    continue
                if expected_len is None:
                    expected_len = len(numeric)
                elif len(numeric) != expected_len:
                    raise ValueError(f"JSON column array lengths are inconsistent: {path}")
                header.append(str(key))
                columns.append(numeric)

            if expected_len is None:
                raise ValueError(f"JSON problems entry does not contain numeric arrays: {path}")

            rows = []
            for i in range(expected_len):
                rows.append([columns[c][i] for c in range(len(columns))])
            return header, rows

        profiles = data.get("profiles", {})
        header = profiles.get("header", [])
        rows = profiles.get("rows", [])
        if not isinstance(header, list) or not isinstance(rows, list):
            raise ValueError(f"JSON profiles section is malformed: {path}")
        numeric_rows = []
        for row in rows:
            if not isinstance(row, list):
                raise ValueError(f"JSON profile row is not a list: {path}")
            numeric_rows.append([float(v) for v in row])
        return [str(h) for h in header], numeric_rows

    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"Tabulated file has too few lines: {path}")
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        rows.append([float(v) for v in parts])
    return header, rows


def main():
    parser = argparse.ArgumentParser(description="Compare profile content between tabulated reference and JSON output.")
    parser.add_argument("--left", required=True, help="Left file path (tabulated reference)")
    parser.add_argument("--right", required=True, help="Right file path (.json or tabulated)")
    parser.add_argument("--coord-tol", type=float, default=1e-12, help="Tolerance for x/y/z coordinates")
    parser.add_argument("--value-tol", type=float, default=1e-6, help="Tolerance for profile values")
    args = parser.parse_args()

    left_path = Path(args.left)
    right_path = Path(args.right)
    if not left_path.is_file():
        raise SystemExit(f"ERROR: left file not found: {left_path}")
    if not right_path.is_file():
        raise SystemExit(f"ERROR: right file not found: {right_path}")

    left_header, left_rows = load_table(left_path)
    right_header, right_rows = load_table(right_path)

    if len(left_rows) != len(right_rows):
        raise SystemExit(
            f"ERROR: row count mismatch ({left_path}: {len(left_rows)}, {right_path}: {len(right_rows)})"
        )

    right_indices = {name: idx for idx, name in enumerate(right_header)}
    missing = [name for name in left_header if name not in right_indices]
    if missing:
        raise SystemExit(f"ERROR: output is missing columns: {', '.join(missing)}")

    for r, (left_row, right_row) in enumerate(zip(left_rows, right_rows), start=2):
        for c, name in enumerate(left_header):
            if c >= len(left_row):
                raise SystemExit(f"ERROR: malformed left row at line {r} in {left_path}")
            right_idx = right_indices[name]
            if right_idx >= len(right_row):
                raise SystemExit(f"ERROR: malformed right row at line {r} in {right_path}")
            tol = args.coord_tol if name in {"x", "y", "z"} else args.value_tol
            if abs(left_row[c] - right_row[right_idx]) > tol:
                raise SystemExit(
                    f"ERROR: mismatch at row {r}, column '{name}' (tol={tol}) "
                    f"left={left_row[c]} right={right_row[right_idx]}"
                )


if __name__ == "__main__":
    main()
