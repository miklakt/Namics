"""Grand-canonical micelle root search helper."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import root_scalar


@dataclass
class SearchResult:
    root: float
    fx: float
    evaluations: int
    seed_json: Path
    output_json: Path
    method: str


def solve(
    seed_script,
    run_scrtipt,
    what_is_X_provided_as_nested_key,
    what_is_Fx_provided_as_nested_key,
    method="brentq",
    X0=None,
    X1=None,
    xtol=1e-3,
    rtol=0.0,
    maxiter=50,
    binary=Path("bin/namics"),
    workdir=None,
    return_result=False,
):
    here = Path(__file__).resolve().parent.parent / "data" / "micelle_grand_canonical"

    def resolve_source(source: Any) -> Path:
        path = source if isinstance(source, Path) else Path(str(source))
        if path.is_file():
            return path.resolve()
        if not path.is_absolute():
            candidate = here / path
            if candidate.is_file():
                return candidate.resolve()
        raise FileNotFoundError(path)

    def split_path(value: Any) -> tuple[str, ...]:
        if callable(value):
            return ()
        if isinstance(value, (list, tuple)):
            return tuple(str(part) for part in value)
        text = str(value)
        if ":" in text:
            return tuple(part for part in text.split(":") if part)
        if "." in text:
            return tuple(part for part in text.split(".") if part)
        return (text,)

    def read_json(path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))

    def write_json(path: Path, payload: dict) -> None:
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    def set_nested(node: dict, path: tuple[str, ...], value: Any) -> None:
        for key in path[:-1]:
            child = node.get(key)
            if not isinstance(child, dict):
                child = {}
                node[key] = child
            node = child
        node[path[-1]] = value

    def get_nested(node: Any, path: tuple[str, ...]) -> Any:
        if callable(what_is_Fx_provided_as_nested_key):
            return what_is_Fx_provided_as_nested_key(node)
        for key in path:
            if not isinstance(node, dict):
                return None
            node = node.get(key)
        return node

    def get_x(problem: dict) -> Any:
        if not x_path:
            return None
        node: Any = problem
        for key in x_path:
            if not isinstance(node, dict):
                return None
            node = node.get(key)
        return node

    seed_source = resolve_source(seed_script)
    run_source = resolve_source(run_scrtipt)
    x_path = split_path(what_is_X_provided_as_nested_key)
    fx_path = split_path(what_is_Fx_provided_as_nested_key)
    rtol = max(float(rtol), 4.0 * np.finfo(float).eps)
    if not x_path and not callable(what_is_X_provided_as_nested_key):
        raise ValueError("X path/delegate is required")
    if not fx_path and not callable(what_is_Fx_provided_as_nested_key):
        raise ValueError("F(X) path/delegate is required")

    base_workdir = Path(workdir) if workdir is not None else here / ".tmp"
    base_workdir.mkdir(parents=True, exist_ok=True)
    workdir = Path(tempfile.mkdtemp(prefix="micelle_gc_", dir=base_workdir))
    binary = Path(binary)
    (workdir / "output").mkdir(exist_ok=True)
    print(f"[micelle] workspace: {workdir}", flush=True)

    seed_input = workdir / seed_source.name
    run_input = workdir / run_source.name
    shutil.copy2(seed_source, seed_input)
    shutil.copy2(run_source, run_input)

    seed_payload = read_json(seed_input)
    run_payload = read_json(run_input)
    seed_problem = seed_payload["problems"][-1]

    seed_x = X0 if X0 is not None else X1
    if seed_x is None:
        seed_x = get_x(seed_problem)
    if seed_x is None:
        seed_x = 1.0

    def set_x(problem: dict, x: float) -> None:
        if callable(what_is_X_provided_as_nested_key):
            what_is_X_provided_as_nested_key(problem, x)
        else:
            set_nested(problem, x_path, x)

    def run_once(input_path: Path, payload: dict, x: float) -> tuple[float, Path]:
        nonlocal last_output_json
        print(f"[micelle] evaluate x={x:.15g}", flush=True)
        set_x(payload["problems"][-1], x)
        write_json(input_path, payload)
        output_name = (
            payload["problems"][-1]
            .get("output", {})
            .get("json", {})
            .get("filename", input_path.stem)
        )
        output_json = input_path.parent / f"{output_name}.output.json"
        proc = subprocess.run(
            [str(binary.resolve()), str(input_path.resolve())],
            cwd=str(workdir),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if proc.returncode != 0:
            if proc.stdout:
                print(proc.stdout, end="" if proc.stdout.endswith("\n") else "\n", flush=True)
            raise subprocess.CalledProcessError(proc.returncode, proc.args, output=proc.stdout)
        value = get_nested(read_json(output_json)["problems"][-1], fx_path)
        if value is None:
            raise ValueError(f"missing F(X) at {fx_path}")
        last_output_json = output_json
        return float(value), output_json

    evaluations = 0
    last_output_json = seed_input.with_name(f"{seed_payload['problems'][-1].get('output', {}).get('json', {}).get('filename', seed_input.stem)}.output.json")
    seed_fx, seed_output = run_once(seed_input, seed_payload, float(seed_x))
    evaluations += 1
    print(f"[micelle] seed f(x)={seed_fx:.15g}", flush=True)

    if method == "secant":
        if X0 is None or X1 is None:
            raise ValueError("secant method requires X0 and X1")

        def f(x: float) -> float:
            nonlocal evaluations, run_payload
            fx, _ = run_once(run_input, run_payload, float(x))
            evaluations += 1
            return fx

        result = root_scalar(f, method="secant", x0=float(X0), x1=float(X1), xtol=xtol, rtol=rtol, maxiter=maxiter)
        root = float(result.root)
        fun = getattr(result, "fun", None)
        fx = float(fun) if fun is not None else f(root)
        return SearchResult(root, fx, evaluations, seed_output, last_output_json, "secant") if return_result else root

    x0 = float(seed_x if X0 is None else X0)
    x1 = float(x0 + 1.0 if X1 is None else X1)

    f0 = seed_fx if x0 == float(seed_x) else run_once(run_input, run_payload, x0)[0]
    evaluations += 0 if x0 == float(seed_x) else 1
    if f0 == 0.0:
        print(f"[micelle] converged at x={x0:.15g}", flush=True)
        return SearchResult(x0, 0.0, evaluations, seed_output, last_output_json, "brentq") if return_result else x0

    f1, _ = run_once(run_input, run_payload, x1)
    evaluations += 1
    print(f"[micelle] bracket start x0={x0:.15g} f0={f0:.15g} x1={x1:.15g} f1={f1:.15g}", flush=True)
    for _ in range(maxiter):
        if f0 * f1 < 0.0:
            break
        step = (x1 - x0) or (1.0 if x1 == 0.0 else abs(x1) * 0.5)
        x0, f0 = x1, f1
        x1 = x1 + step * 2.0
        f1, _ = run_once(run_input, run_payload, x1)
        evaluations += 1
        print(f"[micelle] expand bracket x1={x1:.15g} f1={f1:.15g}", flush=True)
    else:
        result = root_scalar(
            lambda x: run_once(run_input, run_payload, float(x))[0],
            method="secant",
            x0=x0,
            x1=x1,
            xtol=xtol,
            rtol=rtol,
            maxiter=maxiter,
        )
        root = float(result.root)
        fun = getattr(result, "fun", None)
        fx = float(fun) if fun is not None else float(run_once(run_input, run_payload, root)[0])
        return SearchResult(root, fx, evaluations, seed_output, last_output_json, "secant") if return_result else root

    result = root_scalar(
        lambda x: run_once(run_input, run_payload, float(x))[0],
        method="brentq",
        bracket=(x0, x1),
        xtol=xtol,
        rtol=rtol,
        maxiter=maxiter,
    )
    root = float(result.root)
    fun = getattr(result, "fun", None)
    fx = float(fun) if fun is not None else float(run_once(run_input, run_payload, root)[0])
    print(f"[micelle] done root={root:.15g} fx={fx:.15g}", flush=True)
    return SearchResult(root, fx, evaluations, seed_output, last_output_json, "brentq") if return_result else root
