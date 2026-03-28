"""Root search helper."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

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


class SearchAbort(RuntimeError):
    """Raised when the search must stop early."""


def solve(
    seed_script,
    run_scrtipt,
    x_key,
    fx_key,
    method="brentq",
    X0=None,
    X1=None,
    xtol=1e-3,
    rtol=0.0,
    x_min=None,
    x_max=None,
    max_expand_step=None,
    output_validation: Callable[[Path], str | None] | None = None,
    maxiter=50,
    binary=Path("bin/namics"),
    workdir=None,
    return_result=False,
):
    def resolve_source(source: Any) -> Path:
        path = source if isinstance(source, Path) else Path(str(source))
        if path.is_file():
            return path.resolve()
        raise FileNotFoundError(path)

    def preprocess_source(source: Path) -> Path:
        if source.suffix != ".in":
            return source
        output = workdir / f"{source.stem}.input.json"
        created_paths.add(output)
        subprocess.run(
            [sys.executable, str(Path(__file__).resolve().parent / "preprocess_input.py"), str(source), "--output", str(output)],
            check=True,
        )
        return output

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

    def set_nested(node: dict, path: tuple[str, ...], value: Any) -> None:
        for key in path[:-1]:
            child = node.get(key)
            if not isinstance(child, dict):
                child = {}
                node[key] = child
            node = child
        node[path[-1]] = value

    x_path = split_path(x_key)
    fx_path = split_path(fx_key)
    rtol = max(float(rtol), 4.0 * np.finfo(float).eps)
    x_min = None if x_min is None else float(x_min)
    x_max = None if x_max is None else float(x_max)
    max_expand_step = None if max_expand_step is None else float(max_expand_step)
    if x_min is not None and x_max is not None and x_min >= x_max:
        raise ValueError("x_min must be smaller than x_max")
    if max_expand_step is not None and max_expand_step <= 0.0:
        raise ValueError("max_expand_step must be positive")
    if not x_path and not callable(x_key):
        raise ValueError("X path/delegate is required")
    if not fx_path and not callable(fx_key):
        raise ValueError("F(X) path/delegate is required")

    workdir = Path(workdir) if workdir is not None else Path.cwd() / ".tmp"
    workdir.mkdir(parents=True, exist_ok=True)
    binary = Path(binary)
    print(f"[search] workspace: {workdir}", flush=True)
    created_paths: set[Path] = set()

    seed_source_raw, run_source_raw = resolve_source(seed_script), resolve_source(run_scrtipt)
    seed_source, run_source = preprocess_source(seed_source_raw), preprocess_source(run_source_raw)
    created_paths.discard(run_source)
    seed_payload, run_payload = json.loads(seed_source.read_text(encoding="utf-8")), json.loads(run_source.read_text(encoding="utf-8"))
    for payload, source in ((seed_payload, seed_source_raw), (run_payload, run_source_raw)):
        for problem in payload.get("problems", []):
            if not isinstance(problem, dict):
                continue
            json_section = problem.get("output", {}).get("json")
            if isinstance(json_section, dict) and json_section.pop("filename", None) is not None:
                print(f"[search] warning: ignoring output:json:filename in {source}", file=sys.stderr, flush=True)
    seed_problem = seed_payload["problems"][-1]

    if X0 is not None:
        seed_x = X0
    elif X1 is not None:
        seed_x = X1
    elif x_path:
        node: Any = seed_problem
        for key in x_path:
            if not isinstance(node, dict):
                node = None
                break
            node = node.get(key)
        seed_x = node
    else:
        seed_x = None
    if seed_x is None:
        seed_x = 1.0

    def ensure_x_in_bounds(x: float, phase: str) -> None:
        if x_min is not None and x <= x_min:
            raise SearchAbort(f"{phase} x={x:.15g} is not greater than lower bound {x_min:.15g}")
        if x_max is not None and x >= x_max:
            raise SearchAbort(f"{phase} x={x:.15g} is not less than upper bound {x_max:.15g}")

    def clamp_x(x: float, anchor: float) -> float:
        if x_min is not None and x <= x_min:
            x = 0.5 * (float(np.nextafter(x_min, np.inf)) + anchor)
        if x_max is not None and x >= x_max:
            x = 0.5 * (float(np.nextafter(x_max, -np.inf)) + anchor)
        return x

    def line_root(points: list[tuple[float, float]]) -> float | None:
        sample = np.asarray(points[-4:], dtype=float)
        if sample.shape[0] >= 3:
            xs = sample[:, 0]
            fs = sample[:, 1]
            slope, intercept = np.linalg.lstsq(np.column_stack((xs, np.ones_like(xs))), fs, rcond=None)[0]
            if slope != 0.0:
                return float(-intercept / slope)
        x0, f0 = sample[-2]
        x1, f1 = sample[-1]
        return None if f1 == f0 else float(x1 - f1 * (x1 - x0) / (f1 - f0))

    def run_once(template_source: Path, payload: dict, x: float, *, phase: str) -> tuple[float, Path]:
        nonlocal last_output_json
        nonlocal evaluations

        ensure_x_in_bounds(x, phase)
        run_dir = Path(tempfile.mkdtemp(prefix=f"{phase}_{evaluations + 1:04d}_", dir=workdir))
        created_paths.add(run_dir)
        (run_dir / "output").mkdir(exist_ok=True)
        input_path = run_dir / template_source.name
        created_paths.add(input_path)
        output_section = payload["problems"][-1].get("json", {})
        if isinstance(output_section, dict):
            sys_outputs = output_section.get("sys")
            if isinstance(sys_outputs, dict):
                sys_section = payload["problems"][-1].get("sys")
                if not isinstance(sys_section, dict):
                    sys_section = {}
                    payload["problems"][-1]["sys"] = sys_section
                for sys_name in sys_outputs:
                    sys_section.setdefault(sys_name, {})
        if callable(x_key):
            x_key(payload["problems"][-1], x)
        else:
            set_nested(payload["problems"][-1], x_path, x)
        if phase != "seed" and last_output_json.is_file():
            sys_section = payload["problems"][-1].get("sys", {})
            for sys_value in sys_section.values():
                if not isinstance(sys_value, dict):
                    continue
                if sys_value.get("initial_guess") not in (None, "file"):
                    continue
                sys_value["initial_guess"] = "file"
                sys_value["guess_inputfile"] = str(last_output_json.resolve())
        base_name = input_path.name[: -len(".input.json")] if input_path.name.endswith(".input.json") else input_path.stem
        output_json = run_dir / "output" / f"{base_name}.output.json"
        input_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        proc = subprocess.run(
            [str(binary.resolve()), str(input_path.resolve())],
            cwd=str(run_dir),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if proc.returncode != 0:
            stdout = proc.stdout or ""
            if stdout:
                print(stdout, end="" if stdout.endswith("\n") else "\n", flush=True)
            raise subprocess.CalledProcessError(proc.returncode, proc.args, output=proc.stdout)

        if output_validation is not None and phase != "seed":
            decision = output_validation(output_json)
            if decision is not None:
                raise SearchAbort(f"{phase} x={x:.15g}: {decision}")

        problem = json.loads(output_json.read_text(encoding="utf-8"))["problems"][-1]
        if callable(fx_key):
            fx = float(fx_key(problem))
        else:
            node: Any = problem
            for key in fx_path:
                if not isinstance(node, dict):
                    node = None
                    break
                node = node.get(key)
            fx = float(node)
        evaluations += 1
        print(f"[search] iteration {evaluations}: x={x:.15g} fx={fx:.15g}", flush=True)
        last_output_json = output_json
        return fx, output_json

    evaluations = 0
    last_output_json = workdir / "result.output.json"
    try:
        seed_fx, seed_output = run_once(seed_source, seed_payload, float(seed_x), phase="seed")
        print(f"[search] seed f(x)={seed_fx:.15g}", flush=True)

        if method == "secant":
            x0 = float(seed_x if X0 is None else X0)
            if X1 is None:
                x1 = x0 + 1.0
                if x_max is not None and x1 >= x_max:
                    x1 = x0 - 1.0
            else:
                x1 = float(X1)
            points = [
                (x0, seed_fx if x0 == float(seed_x) else run_once(run_source, run_payload, x0, phase="run")[0]),
                (x1, run_once(run_source, run_payload, x1, phase="run")[0]),
            ]
            for _ in range(maxiter):
                if points[-1][1] == 0.0:
                    break
                root = line_root(points)
                if root is None:
                    break
                if max_expand_step is not None:
                    delta = root - points[-1][0]
                    if abs(delta) > max_expand_step:
                        root = points[-1][0] + np.copysign(max_expand_step, delta)
                root = clamp_x(root, points[-1][0])
                if abs(root - points[-1][0]) <= xtol + rtol * abs(root):
                    break
                fx, _ = run_once(run_source, run_payload, root, phase="run")
                points.append((root, fx))
            root, fx = points[-1]
            final_output_json = workdir / last_output_json.name
            shutil.copy2(last_output_json, final_output_json)
            return SearchResult(root, fx, evaluations, seed_output, final_output_json, "secant") if return_result else root

        x0 = float(seed_x if X0 is None else X0)
        if X1 is None:
            x1 = x0 + 1.0
            if x_max is not None and x1 >= x_max:
                x1 = x0 - 1.0
        else:
            x1 = float(X1)

        f0 = seed_fx if x0 == float(seed_x) else run_once(run_source, run_payload, x0, phase="run")[0]
        if f0 == 0.0:
            print(f"[search] converged at x={x0:.15g}", flush=True)
            final_output_json = workdir / last_output_json.name
            shutil.copy2(last_output_json, final_output_json)
            return SearchResult(x0, 0.0, evaluations, seed_output, final_output_json, "brentq") if return_result else x0

        f1, _ = run_once(run_source, run_payload, x1, phase="run")
        if f1 == 0.0:
            print(f"[search] converged at x={x1:.15g}", flush=True)
            final_output_json = workdir / last_output_json.name
            shutil.copy2(last_output_json, final_output_json)
            return SearchResult(x1, 0.0, evaluations, seed_output, final_output_json, "brentq") if return_result else x1
        print(f"[search] bracket start x0={x0:.15g} f0={f0:.15g} x1={x1:.15g} f1={f1:.15g}", flush=True)

        if x0 <= x1:
            xl, fl = x0, f0
            xr, fr = x1, f1
        else:
            xl, fl = x1, f1
            xr, fr = x0, f0

        if fl * fr >= 0.0:
            if fl == fr:
                expand_side = "right" if abs(fr) <= abs(fl) else "left"
            elif fl < fr:
                expand_side = "right" if fr < 0.0 else "left"
            else:
                expand_side = "right" if fr > 0.0 else "left"

            for _ in range(maxiter):
                if fl * fr < 0.0:
                    break

                width = xr - xl
                step = width if width > 0.0 else max(1.0, abs(xl), abs(xr))
                if max_expand_step is not None:
                    step = min(step, max_expand_step)
                jump = 1.2 * step

                if expand_side == "right":
                    x_new = xr + jump
                    if x_max is not None and x_new >= x_max:
                        expand_side = "left"
                        continue
                    f_new, _ = run_once(run_source, run_payload, x_new, phase="run")
                    print(f"[search] expand bracket right xr={x_new:.15g} fr={f_new:.15g}", flush=True)
                    if f_new == 0.0:
                        print(f"[search] converged at x={x_new:.15g}", flush=True)
                        final_output_json = workdir / last_output_json.name
                        shutil.copy2(last_output_json, final_output_json)
                        return SearchResult(x_new, 0.0, evaluations, seed_output, final_output_json, "brentq") if return_result else x_new
                    if fr * f_new < 0.0:
                        xl, fl = xr, fr
                        xr, fr = x_new, f_new
                        break
                    xr, fr = x_new, f_new
                else:
                    x_new = xl - jump
                    if x_min is not None and x_new <= x_min:
                        expand_side = "right"
                        continue
                    f_new, _ = run_once(run_source, run_payload, x_new, phase="run")
                    xl, fl = x_new, f_new
                    print(f"[search] expand bracket left xl={xl:.15g} fl={fl:.15g}", flush=True)
                    if fl == 0.0:
                        print(f"[search] converged at x={xl:.15g}", flush=True)
                        final_output_json = workdir / last_output_json.name
                        shutil.copy2(last_output_json, final_output_json)
                        return SearchResult(xl, 0.0, evaluations, seed_output, final_output_json, "brentq") if return_result else xl
            else:
                points = [(xl, fl), (xr, fr)]
                for _ in range(maxiter):
                    root = line_root(points)
                    if root is None:
                        break
                    root = clamp_x(root, points[-1][0])
                    if abs(root - points[-1][0]) <= xtol + rtol * abs(root):
                        break
                    fx, _ = run_once(run_source, run_payload, root, phase="run")
                    points.append((root, fx))
                    if fx == 0.0:
                        break
                root, fx = points[-1]
                final_output_json = workdir / last_output_json.name
                shutil.copy2(last_output_json, final_output_json)
                return SearchResult(root, fx, evaluations, seed_output, final_output_json, "secant") if return_result else root

        result = root_scalar(
            lambda x: run_once(run_source, run_payload, float(x), phase="run")[0],
            method="brentq",
            bracket=(xl, xr),
            xtol=xtol,
            rtol=rtol,
            maxiter=maxiter,
        )
        root = float(result.root)
        fun = getattr(result, "fun", None)
        fx = float(fun) if fun is not None else float(run_once(run_source, run_payload, root, phase="run")[0])
        print(f"[search] done root={root:.15g} fx={fx:.15g}", flush=True)
        final_output_json = workdir / last_output_json.name
        shutil.copy2(last_output_json, final_output_json)
        return SearchResult(root, fx, evaluations, seed_output, final_output_json, "brentq") if return_result else root
    finally:
        for path in sorted(created_paths, key=lambda p: len(p.parts), reverse=True):
            shutil.rmtree(path, ignore_errors=True) if path.is_dir() else path.unlink(missing_ok=True)
