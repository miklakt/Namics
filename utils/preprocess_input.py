#!/usr/bin/env python3
"""Mirror the NAMICS input preprocessing flow without running the solver."""

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


MASK_FILE_KEYS = ("file", "filename", "path")
MASK_COORDINATES_KEY = "coordinates"
LEGACY_RANGE_KEYWORDS = {
    "firstlayer",
    "lastlayer",
    "var_pos",
    "firstlayer_x",
    "firstlayer_y",
    "firstlayer_z",
    "lastlayer_x",
    "lastlayer_y",
    "lastlayer_z",
    "lowerbound",
    "upperbound",
    "lowerbound_x",
    "lowerbound_y",
    "lowerbound_z",
    "upperbound_x",
    "upperbound_y",
    "upperbound_z",
}
INT_RE = re.compile(r"[+-]?\d+")
FLOAT_RE = re.compile(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")
MASK_KEYS = (("pinned_range", "pinned_filename"), ("frozen_range", "frozen_filename"))


class PreprocessError(RuntimeError):
    """Raised when preprocessing fails."""


@dataclass
class ProblemShape:
    gradients: int = 1
    layers: list[int] = field(default_factory=lambda: [0, 1, 1])
    bc: list[str] = field(default_factory=lambda: ["mirror"] * 6)


@dataclass
class NormalizeOptions:
    emitted_path: Path
    link_data: bool = False


def normalize_path(path: Path) -> Path:
    return Path(os.path.normpath(str(path)))


def json_path(base_dir: Path, stem: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("._") or "data"
    return base_dir / f"{safe}.json"


def sidecar_dir(options: NormalizeOptions) -> Path:
    return options.emitted_path.parent / f"{options.emitted_path.stem}.data"


def sidecar_link(options: NormalizeOptions, stem: str, payload: Any) -> str:
    target_dir = sidecar_dir(options)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = json_path(target_dir, stem)
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return os.path.relpath(target, options.emitted_path.parent)


def inline_or_link(options: NormalizeOptions, stem: str, payload: Any) -> Any:
    return sidecar_link(options, stem, payload) if options.link_data else payload


def read_text(path: Path) -> str:
    if not path.is_file():
        raise PreprocessError(f"Inputfile {path} is not found.")
    return path.read_text(encoding="utf-8")


def compact_legacy_line(line: str) -> str:
    return line.replace(" ", "").replace("\t", "")


def split_legacy(line: str, delim: str) -> list[str]:
    parts: list[str] = []
    for item in line.split(delim):
        compact = item.replace(" ", "")
        comment = compact.find("//")
        parts.append(compact[:comment] if comment >= 0 else compact)
    return parts


def parse_legacy_value(value: str) -> Any:
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if INT_RE.fullmatch(value):
        return int(value)
    if FLOAT_RE.fullmatch(value):
        return float(value)
    return value


def clone_json_like(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clone_json_like(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clone_json_like(item) for item in value]
    return value


def resolve_relative_to(base_path: Path, other_path: str) -> Path:
    candidate = Path(other_path)
    if candidate.is_absolute():
        return normalize_path(candidate)
    parent = base_path.parent if base_path.parent != Path() else Path(".")
    return normalize_path(parent / candidate)


def strip_json_comments(text: str) -> str:
    out: list[str] = []
    in_string = False
    escaped = False
    i = 0
    while i < len(text):
        ch = text[i]
        if in_string:
            out.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == '"':
                in_string = False
            i += 1
            continue
        if ch == '"':
            in_string = True
            out.append(ch)
            i += 1
            continue
        if ch == "/" and i + 1 < len(text) and text[i + 1] == "/":
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def read_json_file(path: Path) -> Any:
    try:
        return json.loads(strip_json_comments(read_text(path)))
    except Exception as error:  # pragma: no cover - exercised via CLI
        raise PreprocessError(f"Failed to parse JSON input file {path}: {error}") from error


def looks_like_json_input(path: Path) -> bool:
    if not path.is_file():
        return False
    text = read_text(path)
    i = 0
    while i < len(text):
        ch = text[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "/" and i + 1 < len(text) and text[i + 1] == "/":
            i += 2
            while i < len(text) and text[i] != "\n":
                i += 1
            continue
        return ch in "{["
    return False


def read_legacy_entries(path: Path) -> list[str]:
    entries: list[str] = []
    for raw_line in read_text(path).splitlines():
        line = compact_legacy_line(raw_line)
        if not line or line.startswith("//"):
            continue
        if line.startswith("include:"):
            include_spec = line[8:]
            if include_spec.endswith("::"):
                include_spec = include_spec[:-2]
            include_path = resolve_relative_to(path, include_spec)
            try:
                include_text = read_text(include_path)
            except PreprocessError as error:
                raise PreprocessError(f"'include : {include_path}' not working because file is not found") from error
            for include_raw in include_text.splitlines():
                include_line = compact_legacy_line(include_raw)
                if include_line and not include_line.startswith("//"):
                    entries.append(include_line)
            continue
        entries.append(line)
    if entries and entries[-1] != "start":
        entries.append("start")
    return entries


def append_json_selector(problem: dict[str, Any], section: str, entry: str, prop: str) -> None:
    existing = problem.setdefault("json", {}).setdefault(section, {}).get(entry)
    if existing is None:
        problem["json"][section][entry] = prop
        return
    if isinstance(existing, str):
        problem["json"][section][entry] = [existing, prop]
        return
    if not isinstance(existing, list):
        existing = problem["json"][section][entry] = []
    existing.append(prop)


def convert_legacy_entries(entries: list[str]) -> list[dict[str, Any]]:
    problems: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    for entry in entries:
        if entry == "start":
            problems.append(clone_json_like(current))
            continue
        parts = split_legacy(entry, ":")
        if not parts:
            continue
        if parts[0] == "json":
            if len(parts) != 4:
                raise PreprocessError(f"Malformed legacy json selector: {entry}")
            append_json_selector(current, parts[1], parts[2], parts[3])
            continue
        if len(parts) != 4:
            raise PreprocessError(f"Malformed legacy entry: {entry}")
        current.setdefault(parts[0], {}).setdefault(parts[1], {})[parts[2]] = parse_legacy_value(parts[3])
    return problems


def emitted_json_path(source_path: Path) -> Path:
    if source_path.name.endswith(".input.json"):
        return source_path
    return source_path.with_name(f"{source_path.stem or source_path.name}.input.json")


def read_sanitized_text_file(path: Path) -> list[str]:
    tokens: list[str] = []
    for raw_line in read_text(path).splitlines():
        line = compact_legacy_line(raw_line)
        if not line or line.startswith("//"):
            continue
        for token in line.split("#"):
            if token:
                tokens.append(token)
    return tokens


def find_member(obj: Any, keys: tuple[str, ...] | list[str]) -> Any | None:
    return next((obj[key] for key in keys if isinstance(obj, dict) and key in obj), None)


def is_int_value(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def is_number_value(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def parse_coordinate_list(value: Any, dimensions: int) -> list[list[int]] | None:
    if not isinstance(value, list):
        return None
    if any(not isinstance(point, list) or len(point) != dimensions for point in value):
        return None
    if any(not is_int_value(coord) for point in value for coord in point):
        return None
    return [list(point) for point in value]


def load_problem_shape(problem: dict[str, Any]) -> ProblemShape:
    lat_section = problem.get("lat")
    if not isinstance(lat_section, dict) or not lat_section:
        raise PreprocessError("Unable to infer lattice shape for mask normalization.")
    parameters = next(iter(lat_section.values()))
    if not isinstance(parameters, dict):
        raise PreprocessError("Invalid lattice parameters for mask normalization.")

    try:
        gradients = parameters.get("gradients", 1)
        if not is_int_value(gradients) or gradients < 1 or gradients > 3:
            raise PreprocessError("Invalid gradients for mask normalization.")

        layers = [
            parameters.get("n_layers", -1) if gradients == 1 else parameters.get("n_layers_x", -1),
            parameters.get("n_layers_y", -1) if gradients >= 2 else 1,
            parameters.get("n_layers_z", -1) if gradients >= 3 else 1,
        ]
        if any(not is_int_value(layer) for layer in layers[:gradients]):
            raise PreprocessError("Invalid lattice dimensions for mask normalization.")
        if layers[0] < 1 or (gradients >= 2 and layers[1] < 1) or (gradients >= 3 and layers[2] < 1):
            raise PreprocessError("Invalid lattice dimensions for mask normalization.")

        bc = ["mirror"] * 6
        if gradients == 1:
            bc[0] = str(parameters.get("lowerbound", "mirror"))
            bc[3] = str(parameters.get("upperbound", "mirror"))
        else:
            bc[0] = str(parameters.get("lowerbound_x", "mirror"))
            bc[3] = str(parameters.get("upperbound_x", "mirror"))
            bc[1] = str(parameters.get("lowerbound_y", "mirror"))
            bc[4] = str(parameters.get("upperbound_y", "mirror"))
            if gradients == 3:
                bc[2] = str(parameters.get("lowerbound_z", "mirror"))
                bc[5] = str(parameters.get("upperbound_z", "mirror"))
        return ProblemShape(gradients=gradients, layers=layers, bc=bc)
    except PreprocessError:
        raise
    except Exception as error:  # pragma: no cover - defensive
        raise PreprocessError(f"Invalid lattice configuration for mask normalization: {error}") from error


def parse_legacy_mask_file(path: Path, shape: ProblemShape) -> list[list[int]]:
    coordinates: list[list[int]] = []
    for token in read_sanitized_text_file(path):
        parts = split_legacy(token, ",")
        if len(parts) != shape.gradients:
            raise PreprocessError(f"Invalid coordinate count in mask file {path}.")
        point: list[int] = []
        for part in parts:
            if not INT_RE.fullmatch(part):
                raise PreprocessError(f"Invalid coordinate '{part}' in mask file {path}.")
            point.append(int(part))
        coordinates.append(point)
    return coordinates


def is_legacy_range_text(value: str) -> bool:
    compact = value.strip()
    return ";" in compact or "," in compact or compact.lower() in LEGACY_RANGE_KEYWORDS


def expand_legacy_range_alias(alias: str, kind: str, shape: ProblemShape) -> tuple[list[int], list[int]] | None:
    first = [1] * shape.gradients
    last = [shape.layers[axis] for axis in range(shape.gradients)]
    lowered = alias.lower()
    if kind == "pinned":
        if lowered in {"firstlayer", "firstlayer_x"}:
            first[0] = last[0] = 1
        elif lowered in {"lastlayer", "lastlayer_x"}:
            first[0] = last[0] = shape.layers[0]
        elif lowered in {"firstlayer_y", "lastlayer_y"} and shape.gradients >= 2:
            first[1] = last[1] = 1 if lowered == "firstlayer_y" else shape.layers[1]
        elif lowered in {"firstlayer_z", "lastlayer_z"} and shape.gradients >= 3:
            first[2] = last[2] = 1 if lowered == "firstlayer_z" else shape.layers[2]
        else:
            return None
        return first, last

    def set_surface(axis: int, upper: bool) -> tuple[list[int], list[int]] | None:
        first[axis] = last[axis] = shape.layers[axis] + 1 if upper else 0
        boundary = shape.bc[axis + (3 if upper else 0)]
        return (first, last) if boundary == "surface" else None

    if lowered in {"lowerbound", "lowerbound_x"}:
        return set_surface(0, False)
    if lowered in {"upperbound", "upperbound_x"}:
        return set_surface(0, True)
    if lowered in {"lowerbound_y", "upperbound_y"} and shape.gradients >= 2:
        return set_surface(1, lowered == "upperbound_y")
    if lowered in {"lowerbound_z", "upperbound_z"} and shape.gradients >= 3:
        return set_surface(2, lowered == "upperbound_z")
    return None


def parse_legacy_range(value: str, kind: str, shape: ProblemShape, var_pos: int) -> list[list[int]]:
    compact = value.strip().replace(" ", "").rstrip(";")
    expanded = expand_legacy_range_alias(compact, kind, shape)
    if expanded is None:
        endpoints = split_legacy(compact, ";")
        if len(endpoints) != 2:
            raise PreprocessError(f"Invalid legacy range: {value}")

        def parse_endpoint(text: str) -> list[int]:
            parts = split_legacy(text, ",")
            if len(parts) != shape.gradients:
                raise PreprocessError(f"Invalid legacy range endpoint: {text}")
            point: list[int] = []
            for axis, raw in enumerate(parts):
                token = raw.lower()
                if token == "var_pos":
                    point.append(var_pos)
                elif token == "firstlayer":
                    point.append(1)
                elif token == "lastlayer":
                    point.append(shape.layers[axis])
                elif kind == "frozen" and token == "lowerbound":
                    point.append(0)
                elif kind == "frozen" and token == "upperbound":
                    point.append(shape.layers[axis] + 1)
                elif INT_RE.fullmatch(raw):
                    point.append(int(raw))
                else:
                    raise PreprocessError(f"Invalid legacy range endpoint token: {raw}")
            return point

        first = parse_endpoint(endpoints[0])
        last = parse_endpoint(endpoints[1])
    else:
        first, last = expanded

    if any(a > b for a, b in zip(first, last)):
        raise PreprocessError(f"Invalid legacy range order: {value}")

    coordinates: list[list[int]] = []
    for x in range(first[0], last[0] + 1):
        if shape.gradients == 1:
            coordinates.append([x])
            continue
        for y in range(first[1], last[1] + 1):
            if shape.gradients == 2:
                coordinates.append([x, y])
                continue
            for z in range(first[2], last[2] + 1):
                coordinates.append([x, y, z])
    return coordinates


def extract_mask_coordinates_from_json_file(
    path: Path,
    key: str,
    shape: ProblemShape,
    var_pos: int,
) -> list[list[int]]:
    if not looks_like_json_input(path):
        return parse_legacy_mask_file(path, shape)
    document = read_json_file(path)
    if isinstance(document, dict) and key in document:
        return extract_mask_coordinates(document[key], key, shape, path, var_pos)
    return extract_mask_coordinates(document, key, shape, path, var_pos)


def extract_mask_coordinates(
    source: Any,
    key: str,
    shape: ProblemShape,
    base_path: Path,
    var_pos: int,
) -> list[list[int]]:
    if isinstance(source, str):
        if is_legacy_range_text(source):
            kind = "pinned" if key == "pinned_range" else "frozen"
            return parse_legacy_range(source, kind, shape, var_pos)
        return extract_mask_coordinates_from_json_file(resolve_relative_to(base_path, source), key, shape, var_pos)

    coordinates = parse_coordinate_list(source, shape.gradients)
    if coordinates is not None:
        return coordinates

    if not isinstance(source, dict):
        raise PreprocessError(f"Unable to normalize mask '{key}'.")

    path_value = find_member(source, MASK_FILE_KEYS)
    if isinstance(path_value, str):
        return extract_mask_coordinates_from_json_file(resolve_relative_to(base_path, path_value), key, shape, var_pos)
    if MASK_COORDINATES_KEY in source:
        coordinates = parse_coordinate_list(source[MASK_COORDINATES_KEY], shape.gradients)
        if coordinates is None:
            raise PreprocessError(f"Invalid coordinates for mask '{key}'.")
        return coordinates
    if key in source:
        return extract_mask_coordinates(source[key], key, shape, base_path, var_pos)
    raise PreprocessError(f"Unable to normalize mask '{key}'.")


def canonical_mask(coordinates: list[list[int]]) -> dict[str, Any]:
    return {MASK_COORDINATES_KEY: coordinates}


def parse_number_list(value: Any) -> list[float | int] | None:
    if not isinstance(value, list) or any(not is_number_value(item) for item in value):
        return None
    return list(value)


def find_last_problem_object(document: Any) -> dict[str, Any] | None:
    if not isinstance(document, dict):
        return None
    problems = document.get("problems")
    if isinstance(problems, dict):
        return problems
    if not isinstance(problems, list):
        return None
    return next((problem for problem in reversed(problems) if isinstance(problem, dict)), None)


def find_external_potential_node(document: Any) -> list[float | int] | None:
    values = parse_number_list(document)
    if values is not None:
        return values
    if not isinstance(document, dict):
        return None
    values = parse_number_list(document.get("external_potential"))
    if values is not None:
        return values
    profiles = document.get("profiles")
    if isinstance(profiles, dict):
        values = parse_number_list(profiles.get("external_potential"))
        if values is not None:
            return values
    problem = find_last_problem_object(document)
    return find_external_potential_node(problem) if problem is not None else None


def load_external_potential_source(source: Any, base_path: Path) -> list[float | int]:
    if isinstance(source, str):
        path = resolve_relative_to(base_path, source)
        values = find_external_potential_node(read_json_file(path))
        if values is None:
            raise PreprocessError(f"No embedded external_potential found in {source}.")
        return values
    if isinstance(source, dict):
        file_source = find_member(source, MASK_FILE_KEYS)
        if isinstance(file_source, str):
            return load_external_potential_source(file_source, base_path)
    values = find_external_potential_node(source)
    if values is None:
        raise PreprocessError("Invalid external_potential payload.")
    return values


def normalize_external_potentials(
    problem: dict[str, Any],
    base_path: Path,
    options: NormalizeOptions,
    problem_index: int,
) -> None:
    mon_section = problem.get("mon")
    if not isinstance(mon_section, dict):
        return
    for mon_name, mon_value in mon_section.items():
        if not isinstance(mon_value, dict):
            continue
        source = mon_value.get("external_potential", mon_value.get("external_potential_filename"))
        if source is None:
            continue
        values = load_external_potential_source(source, base_path)
        if options.link_data:
            mon_value["external_potential_filename"] = sidecar_link(
                options,
                f"problem_{problem_index}.{mon_name}.external_potential",
                {"external_potential": values},
            )
            mon_value.pop("external_potential", None)
            continue
        mon_value["external_potential"] = values
        mon_value.pop("external_potential_filename", None)


def normalize_masks(
    problem: dict[str, Any],
    shape: ProblemShape,
    base_path: Path,
    options: NormalizeOptions,
    problem_index: int,
) -> None:
    mon_section = problem.get("mon")
    if not isinstance(mon_section, dict):
        return
    for mon_name, mon_value in mon_section.items():
        if not isinstance(mon_value, dict):
            continue
        freedom = mon_value.get("freedom")
        if freedom == "free":
            mon_value.pop("pinned_range", None)
            mon_value.pop("pinned_filename", None)
            mon_value.pop("frozen_range", None)
            mon_value.pop("frozen_filename", None)
        elif freedom == "pinned":
            mon_value.pop("frozen_range", None)
            mon_value.pop("frozen_filename", None)
        elif freedom == "frozen":
            mon_value.pop("pinned_range", None)
            mon_value.pop("pinned_filename", None)
        var_pos = mon_value.get("var_pos", 0)
        if not is_int_value(var_pos):
            raise PreprocessError("Invalid var_pos for mask normalization.")
        for key, filename_key in MASK_KEYS:
            if filename_key in mon_value:
                mon_value[key] = mon_value[filename_key]
                del mon_value[filename_key]
            if key not in mon_value:
                continue
            coordinates = extract_mask_coordinates(mon_value[key], key, shape, base_path, var_pos)
            mon_value[key] = inline_or_link(
                options, f"problem_{problem_index}.{mon_name}.{key}", canonical_mask(coordinates)
            )


def find_initial_guess_object(document: Any) -> dict[str, Any] | None:
    if isinstance(document, dict) and isinstance(document.get("profiles"), dict):
        return document
    if not isinstance(document, dict):
        return None
    initial_guess = document.get("initial_guess")
    if isinstance(initial_guess, dict):
        return initial_guess
    problem = find_last_problem_object(document)
    return find_initial_guess_object(problem) if problem is not None else None


def load_initial_guess_source(source: Any, base_path: Path) -> dict[str, Any]:
    if isinstance(source, str):
        document = read_json_file(resolve_relative_to(base_path, source))
        embedded = find_initial_guess_object(document)
        if embedded is None:
            raise PreprocessError(f"No embedded initial_guess found in {source}.")
        return embedded
    if not isinstance(source, dict):
        raise PreprocessError("Invalid initial_guess source.")
    file_source = find_member(source, MASK_FILE_KEYS)
    if isinstance(file_source, str):
        return load_initial_guess_source(file_source, base_path)
    embedded = find_initial_guess_object(source)
    if embedded is None:
        raise PreprocessError("Invalid initial_guess payload.")
    return embedded


def normalize_initial_guess(
    problem: dict[str, Any],
    base_path: Path,
    options: NormalizeOptions,
    problem_index: int,
) -> None:
    if "initial_guess" in problem:
        problem["initial_guess"] = inline_or_link(
            options,
            f"problem_{problem_index}.initial_guess",
            load_initial_guess_source(problem["initial_guess"], base_path),
        )

    sys_section = problem.get("sys")
    if not isinstance(sys_section, dict):
        return

    sys_names = [name for name, sys_value in sys_section.items() if isinstance(sys_value, dict)]
    for sys_name in sys_names:
        sys_value = problem["sys"][sys_name]
        if "initial_guess" not in sys_value:
            continue
        if isinstance(sys_value["initial_guess"], str) and sys_value["initial_guess"] in {"previous_result", "file", "none"}:
            continue

        problem["initial_guess"] = inline_or_link(
            options,
            f"problem_{problem_index}.initial_guess",
            load_initial_guess_source(sys_value["initial_guess"], base_path),
        )
        problem["sys"][sys_name]["initial_guess"] = "file"


def iter_problems(document: Any) -> list[tuple[int, Any]]:
    if isinstance(document, list):
        return list(enumerate(document, start=1))
    if not isinstance(document, dict):
        return []
    problems = document.get("problems", document)
    return list(enumerate(problems if isinstance(problems, list) else [problems], start=1))


def prepare_input_document(requested_path: Path, options: NormalizeOptions) -> Any:
    if looks_like_json_input(requested_path):
        document = read_json_file(requested_path)
    else:
        document = {"problems": convert_legacy_entries(read_legacy_entries(requested_path))}
    for index, problem in iter_problems(document):
        if not isinstance(problem, dict):
            continue
        normalize_external_potentials(problem, requested_path, options, index)
        mon_section = problem.get("mon")
        if isinstance(mon_section, dict) and any(
            isinstance(mon_value, dict) and any(key in mon_value for key, filename_key in MASK_KEYS for key in (key, filename_key))
            for mon_value in mon_section.values()
        ):
            normalize_masks(problem, load_problem_shape(problem), requested_path, options, index)
        normalize_initial_guess(problem, requested_path, options, index)
    return document


def write_preprocessed_input(requested_path: Path, output_path: Path | None = None, link_data: bool = False) -> Path:
    emitted = output_path if output_path is not None else emitted_json_path(requested_path)
    emitted.parent.mkdir(parents=True, exist_ok=True)
    document = prepare_input_document(requested_path, NormalizeOptions(emitted_path=emitted, link_data=link_data))
    emitted.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return emitted


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Mirror the NAMICS input preprocessing flow.")
    parser.add_argument("inputs", nargs="+", type=Path, help="Input file(s) to preprocess.")
    parser.add_argument("-o", "--output", type=Path, help="Explicit output path for a single input.")
    parser.add_argument(
        "--link-data",
        action="store_true",
        help="Write normalized masks, external potentials, and initial_guess payloads to sidecar JSON files and keep links in the emitted input.",
    )
    parser.add_argument("--stdout", action="store_true", help="Write the preprocessed JSON for a single input to stdout.")
    args = parser.parse_args()
    if args.output is not None and len(args.inputs) != 1:
        parser.error("--output only supports a single input")
    if args.stdout and len(args.inputs) != 1:
        parser.error("--stdout only supports a single input")
    if args.output is not None and args.stdout:
        parser.error("--output and --stdout are mutually exclusive")
    if args.stdout and args.link_data:
        parser.error("--stdout and --link-data are mutually exclusive")
    return args


def main() -> int:
    args = parse_args()
    try:
        if args.stdout:
            document = prepare_input_document(args.inputs[0], NormalizeOptions(emitted_path=emitted_json_path(args.inputs[0])))
            sys.stdout.write(json.dumps(document, indent=2) + "\n")
            return 0

        for input_path in args.inputs:
            output_path = args.output if args.output is not None else None
            emitted = write_preprocessed_input(input_path, output_path, args.link_data)
            print(emitted)
        return 0
    except PreprocessError as error:
        print(error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
