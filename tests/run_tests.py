#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import shutil
import signal
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_helpers import (  # noqa: E402
    CommandMetrics,
    TestError,
    compare_json_profiles,
    require_file,
    run_command,
    set_commented_setting,
    set_setting_line,
    utc_run_id,
)

COORD_TOL = 1e-12
EXECUTION_FAILED_RE = re.compile(r"execution failed \((-?\d+)\)")


@dataclass
class Context:
    repo_root: Path
    tests_dir: Path
    output_dir: Path
    binary: Path
    quiet: bool
    clean_output: bool
    with_save_memory: bool
    reference_dir: Path
    generated_input_manifest: Path


@dataclass
class ReportNode:
    label: str
    passed: bool = True
    required_for_parent: bool = True
    wall_s: float = 0.0
    max_rss_kb: float | None = None
    details: str = ""
    children: list["ReportNode"] = field(default_factory=list)

    def add_metric(self, metrics: CommandMetrics) -> None:
        self.wall_s += metrics.wall_s
        if metrics.max_rss_kb is None:
            return
        self.max_rss_kb = metrics.max_rss_kb if self.max_rss_kb is None else max(self.max_rss_kb, metrics.max_rss_kb)


def add_child(parent: ReportNode, label: str, **kwargs: object) -> ReportNode:
    node = ReportNode(label=label, **kwargs)
    parent.children.append(node)
    return node

def last_problem(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))["problems"][-1]

def solver_artifacts(path: Path) -> list[Path]:
    return [path, path.with_suffix(".input.json"), path.with_suffix(".output.json")]

def elapsed_ms(metrics: CommandMetrics) -> int:
    return int(round(metrics.wall_s * 1000.0))

def expect_profiles(left: Path, right: Path, tol: float, detail: str = "matched reference") -> str:
    compare_json_profiles(left, right, coord_tol=COORD_TOL, value_tol=tol)
    return detail


def _prepare_runtime_input(
    ctx: Context,
    source_input: Path,
    runtime_input: Path,
    solver_method: str,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
) -> None:
    runtime_input.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_input, runtime_input)
    set_setting_line(runtime_input, runtime_input, "newton : isaac : method", solver_method)
    for key, value in (settings or {}).items():
        set_setting_line(runtime_input, runtime_input, key, value)
    for key, value, enabled in comment_toggles or ():
        set_commented_setting(runtime_input, runtime_input, key, value, enabled)
    variants = {"solver_method": solver_method, **(settings or {})}
    variants.update({key: "enabled" if enabled else "disabled" for key, _, enabled in comment_toggles or ()})
    ctx.generated_input_manifest.parent.mkdir(parents=True, exist_ok=True)
    with ctx.generated_input_manifest.open("a", encoding="utf-8") as fh:
        fh.write(
            "\t".join(
                [
                    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    str(source_input.relative_to(ctx.repo_root)),
                    str(runtime_input.relative_to(ctx.repo_root)),
                    str(variants),
                ]
            )
            + "\n"
        )


def _cleanup(paths: list[Path], enabled: bool) -> None:
    if not enabled:
        return
    for path in paths:
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        else:
            path.unlink(missing_ok=True)

def _finalize_report_tree(node: ReportNode) -> None:
    if not node.children:
        return
    for child in node.children:
        _finalize_report_tree(child)
    required = [child for child in node.children if child.required_for_parent]
    node.passed = all(child.passed for child in required) if required else True
    node.wall_s = sum(child.wall_s for child in node.children)
    node.max_rss_kb = max((child.max_rss_kb for child in node.children if child.max_rss_kb is not None), default=None)


def _set_failure(node: ReportNode, detail: str | Exception, *, accepted: bool = False) -> str:
    node.passed = False
    node.details = f"failed (accepted): {detail}" if accepted else str(detail)
    return node.details

def _failure_suffix_from_detail(detail: str) -> str | None:
    detail = detail.strip()
    match = EXECUTION_FAILED_RE.search(detail)
    if match:
        returncode = int(match.group(1))
        if returncode < 0:
            signal_number = -returncode
            if signal_number == signal.SIGSEGV:
                return "segfault"
            try:
                return f"signal:{signal.Signals(signal_number).name}"
            except ValueError:
                return "signal"
    for prefix in ("failed (accepted): ", "ERROR: ", "RuntimeError: "):
        detail = detail.removeprefix(prefix)
    lines = [line.strip() for line in detail.splitlines() if line.strip()]
    filtered = [
        line.removeprefix("RuntimeError: ").removeprefix("ERROR: ")
        for line in lines
        if not line.startswith(
            (
                "Traceback (most recent call last):",
                "File \"",
                "sys.exit(",
                "raise ",
                "seed -> ",
                "Problem nr ",
            )
        )
        and set(line) != {"^"}
    ]
    for line in reversed(filtered or lines):
        if line.startswith("No initial_guess found in ") and "Read guess for initial guess failed" in line:
            return "Read guess for initial guess failed"
        if line:
            return line
    return None


def _failure_suffix(node: ReportNode) -> str:
    if node.passed:
        return ""
    mode = _failure_suffix_from_detail(node.details)
    if mode:
        return mode
    child_modes = sorted({mode for child in node.children if not child.passed for mode in [_failure_suffix(child)] if mode})
    return child_modes[0] if len(child_modes) == 1 else ("mixed" if child_modes else "unexpected")

def _flatten_nodes(nodes: list[ReportNode]) -> list[tuple[int, bool, str, ReportNode]]:
    out: list[tuple[int, bool, str, ReportNode]] = []

    def walk(root_index: int, node: ReportNode, label_prefix: str, child_prefix: str, is_root: bool) -> None:
        out.append((root_index, is_root, f"{label_prefix}{node.label}", node))
        for idx, child in enumerate(node.children):
            is_last = idx == len(node.children) - 1
            walk(
                root_index,
                child,
                child_prefix + ("└── " if is_last else "├── "),
                child_prefix + ("    " if is_last else "│   "),
                False,
            )

    for root_index, root in enumerate(nodes, start=1):
        walk(root_index, root, "", "", True)
    return out


def _run_solver(
    ctx: Context,
    node: ReportNode,
    source_input: Path,
    runtime_input: Path,
    solver_method: str,
    label: str,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    require_output: bool = True,
    require_initial_guess: bool = False,
) -> tuple[CommandMetrics, Path]:
    output_path = runtime_input.with_suffix(".output.json")
    _prepare_runtime_input(ctx, source_input, runtime_input, solver_method, settings, comment_toggles)
    output_path.unlink(missing_ok=True)
    runtime_input.with_suffix(".input.json").unlink(missing_ok=True)
    metrics = run_command([str(ctx.binary), str(runtime_input)], cwd=ctx.repo_root, quiet=ctx.quiet)
    node.add_metric(metrics)
    if metrics.returncode:
        message = next((line.strip() for line in reversed(metrics.output.splitlines()) if line.strip()), "")
        raise TestError(message or f"ERROR: execution failed ({metrics.returncode}) for input: {runtime_input}")
    if (require_output or require_initial_guess) and (not output_path.is_file() or not output_path.stat().st_size):
        raise TestError(f"ERROR: expected JSON output file was not created for {label}: {output_path}")
    if require_initial_guess:
        initial_guess = last_problem(output_path)["initial_guess"]
        if not isinstance(initial_guess, dict) or set(initial_guess) != {"profiles"} or not initial_guess["profiles"]:
            raise TestError(f"ERROR: initial_guess is missing profile data: {output_path}")
    return metrics, output_path


def _run_leaf(
    ctx: Context,
    node: ReportNode,
    *,
    source_input: Path,
    runtime_input: Path,
    solver_method: str,
    label: str,
    detail: str = "",
    compare: Callable[[Path], str | None] | None = None,
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    accepted: bool = False,
    require_output: bool = True,
    require_initial_guess: bool = False,
) -> tuple[Path | None, CommandMetrics | None]:
    try:
        metrics, output = _run_solver(
            ctx,
            node,
            source_input,
            runtime_input,
            solver_method,
            label,
            settings=settings,
            comment_toggles=comment_toggles,
            require_output=require_output,
            require_initial_guess=require_initial_guess,
        )
        node.details = compare(output) if compare else detail
        return output, metrics
    except Exception as exc:
        _set_failure(node, exc, accepted=accepted)
        return None, None


def _leaf_stem(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def _run_method_group(
    ctx: Context,
    *,
    root_label: str,
    input_file: Path,
    reference_name: str,
    runtime_stem: str,
    value_tol: float,
    label_stem: str | None = None,
    method_group_label: str | None = "solver method",
    settings: dict[str, str] | None = None,
    comment_toggles: list[tuple[str, str, bool]] | None = None,
    copied_inputs: list[Path] | None = None,
    leaves: list[tuple[str, str, str | None, bool, bool]] | None = None,
) -> ReportNode:
    root = ReportNode(root_label)
    group = root if method_group_label is None else add_child(root, method_group_label)
    reference_file = ctx.reference_dir / reference_name
    leaf_specs = leaves or []
    cleanup_targets = [path for label, *_ in leaf_specs for path in solver_artifacts(ctx.output_dir / f"{runtime_stem}.{_leaf_stem(label)}.in")]
    label_base = label_stem or root_label

    require_file(ctx.binary, executable=True)
    require_file(input_file)
    require_file(reference_file)
    for path in copied_inputs or ():
        require_file(path)
    ctx.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        for source in copied_inputs or ():
            target = ctx.output_dir / source.name
            shutil.copy2(source, target)
            cleanup_targets.append(target)

        outputs: dict[str, Path | None] = {}
        for leaf_label, solver_method, compare_to, required_for_parent, accepted in leaf_specs:
            leaf = add_child(group, leaf_label, required_for_parent=required_for_parent)
            runtime_input = ctx.output_dir / f"{runtime_stem}.{_leaf_stem(leaf_label)}.in"
            try:
                _, output = _run_solver(
                    ctx,
                    leaf,
                    input_file,
                    runtime_input,
                    solver_method,
                    f"{label_base}, mode={leaf_label}",
                    settings=settings,
                    comment_toggles=comment_toggles,
                    require_output=not accepted,
                )
            except Exception as exc:
                outputs[leaf_label] = None
                _set_failure(leaf, exc, accepted=accepted)
                continue

            outputs[leaf_label] = output
            try:
                if accepted and (not output.is_file() or not output.stat().st_size):
                    leaf.passed = False
                    leaf.details = "completed but JSON output file is missing (accepted)"
                elif compare_to == "reference":
                    leaf.details = expect_profiles(reference_file, output, value_tol)
                elif compare_to:
                    baseline = outputs.get(compare_to)
                    if baseline is None or not baseline.is_file() or not baseline.stat().st_size:
                        leaf.details = (
                            f"completed but no {compare_to} baseline was available (accepted)"
                            if accepted
                            else _set_failure(leaf, f"ERROR: skipped because {compare_to} failed")
                        )
                        leaf.passed = False if accepted else leaf.passed
                    else:
                        leaf.details = expect_profiles(
                            baseline,
                            output,
                            value_tol,
                            ("passed and matched " if accepted else "matched ") + compare_to,
                        )
            except Exception as exc:
                _set_failure(leaf, exc, accepted=accepted)
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, ctx.clean_output)


def _run_homopolymer_adsorption(ctx: Context, benchmark: bool) -> ReportNode:
    root = ReportNode("homopolymer adsorption benchmark" if benchmark else "homopolymer adsorption")
    template_file = ctx.tests_dir / "homopolymer_adsorption.in"
    benchmark_dir = ctx.tests_dir / "benchmarks"
    run_id = utc_run_id() if benchmark else ""
    run_log = benchmark_dir / f"homopolymer_adsorption_test_{run_id}.log" if benchmark else None
    cleanup_targets = [run_log] if run_log else []
    solver_runtime_ms = 0
    diis_notes: list[str] = []
    log_lines = (
        [
            f"run_id={run_id}",
            f"timestamp_utc={time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}",
            "cases=chi_Si(0,-2,-4,-6)",
            f"save_memory={'on' if ctx.with_save_memory else 'off'}",
        ]
        if benchmark
        else []
    )

    require_file(ctx.binary, executable=True)
    require_file(template_file)
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    if benchmark:
        benchmark_dir.mkdir(parents=True, exist_ok=True)

    start = time.perf_counter()
    try:
        for chi in (0, -2, -4, -6):
            reference_file = ctx.reference_dir / f"homopolymer_adsorption.chi_{chi}.json.ref"
            require_file(reference_file)
            case = add_child(root, f"chi_Si = {chi}")
            settings = {"mon : A : chi_Si": str(chi)}

            pseudo_input = ctx.output_dir / f"homopolymer_adsorption.chi_{chi}.pseudohessian.in"
            cleanup_targets.extend(solver_artifacts(pseudo_input))
            pseudo_leaf = add_child(case, "pseudohessian")
            pseudo_output, pseudo_metrics = _run_leaf(
                ctx,
                pseudo_leaf,
                source_input=template_file,
                runtime_input=pseudo_input,
                solver_method="pseudohessian",
                label=f"chi_Si={chi}, mode=pseudohessian",
                settings=settings,
                compare=lambda out, ref=reference_file: expect_profiles(ref, out, 1e-6),
            )
            if pseudo_metrics:
                solver_runtime_ms += elapsed_ms(pseudo_metrics)
                if benchmark:
                    log_lines.append(f"chi_Si={chi},mode=pseudohessian,elapsed_ms={elapsed_ms(pseudo_metrics)}")

            if ctx.with_save_memory:
                save_input = ctx.output_dir / f"homopolymer_adsorption.chi_{chi}.save_memory.in"
                cleanup_targets.extend(solver_artifacts(save_input))
                save_leaf = add_child(case, "save_memory + pseudohessian")

                def compare_save(out: Path, ref: Path = reference_file, baseline: Path | None = pseudo_output) -> str:
                    expect_profiles(ref, out, 1e-6)
                    if baseline and baseline.is_file() and baseline.stat().st_size:
                        return expect_profiles(baseline, out, 1e-6, "matched reference and pseudohessian")
                    return "matched reference"

                _, save_metrics = _run_leaf(
                    ctx,
                    save_leaf,
                    source_input=template_file,
                    runtime_input=save_input,
                    solver_method="pseudohessian",
                    label=f"chi_Si={chi}, mode=save_memory",
                    settings=settings,
                    comment_toggles=[("mol : pol : save_memory", "true", True)],
                    compare=compare_save,
                )
                if save_metrics:
                    solver_runtime_ms += elapsed_ms(save_metrics)
                    if benchmark:
                        log_lines.append(f"chi_Si={chi},mode=save_memory,elapsed_ms={elapsed_ms(save_metrics)}")

            if not benchmark:
                diis_input = ctx.output_dir / f"homopolymer_adsorption.chi_{chi}.diis.in"
                cleanup_targets.extend(solver_artifacts(diis_input))
                diis_leaf = add_child(case, "DIIS", required_for_parent=False)
                diis_output, _ = _run_leaf(
                    ctx,
                    diis_leaf,
                    source_input=template_file,
                    runtime_input=diis_input,
                    solver_method="DIIS",
                    label=f"chi_Si={chi}, mode=DIIS",
                    settings=settings,
                    accepted=True,
                    require_output=False,
                )
                diis_status = diis_leaf.details
                if diis_output:
                    if not diis_output.is_file() or not diis_output.stat().st_size:
                        diis_leaf.passed = False
                        diis_status = diis_leaf.details = "completed but JSON output file is missing (accepted)"
                    elif not pseudo_leaf.passed or pseudo_output is None or not pseudo_output.is_file() or not pseudo_output.stat().st_size:
                        diis_leaf.passed = False
                        diis_status = diis_leaf.details = "completed but no pseudohessian baseline was available (accepted)"
                    else:
                        try:
                            diis_status = diis_leaf.details = expect_profiles(
                                pseudo_output,
                                diis_output,
                                1e-6,
                                "passed and matched pseudohessian",
                            )
                        except Exception as exc:
                            diis_status = _set_failure(diis_leaf, exc, accepted=True)
                diis_notes.append(
                    f"chi={chi}:{diis_status}"
                )

        if benchmark:
            wall_runtime_ms = int(round((time.perf_counter() - start) * 1000.0))
            log_lines.extend([f"solver_runtime_ms={solver_runtime_ms}", f"wall_runtime_ms={wall_runtime_ms}"])
            root.details = f"bench:solver_ms={solver_runtime_ms}"
            run_log.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
            root.details += f";log={run_log.name}"
        else:
            root.details = f"diis={';'.join(diis_notes)}"

        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, ctx.clean_output)


def test_homopolymer_adsorption(ctx: Context) -> ReportNode:
    return _run_homopolymer_adsorption(ctx, False)

def benchmark_homopolymer_adsorption(ctx: Context) -> ReportNode:
    return _run_homopolymer_adsorption(ctx, True)

def test_frozen_range_input_file(ctx: Context) -> ReportNode:
    return _run_method_group(
        ctx,
        root_label="frozen range input file",
        label_stem="frozen_range_input_file",
        input_file=ctx.tests_dir / "frozen_range_input_file.in",
        reference_name="frozen_range_input_file.json.ref",
        runtime_stem="frozen_range_input_file",
        value_tol=1e-6,
        copied_inputs=[ctx.tests_dir / "frozen_range_input_file.frozen"],
        leaves=[("pseudohessian", "pseudohessian", "reference", True, False), ("DIIS", "DIIS", "pseudohessian", False, True)],
    )


def test_micelle_self_assembly(ctx: Context) -> ReportNode:
    root = ReportNode("micelle self assembly")
    generate_group = add_child(root, "input guess generate")
    use_group = add_child(root, "input guess use")
    generate_input = ctx.tests_dir / "micelle_guess_generate.in"
    use_input = ctx.tests_dir / "micelle_guess_use.in"
    reference_file = ctx.reference_dir / "micelle_guess_use.json.ref"
    runtime_generate = ctx.output_dir / generate_input.name
    runtime_use = ctx.output_dir / use_input.name
    cleanup_targets = [*solver_artifacts(runtime_generate), *solver_artifacts(runtime_use)]

    require_file(ctx.binary, executable=True)
    require_file(generate_input)
    require_file(use_input)
    require_file(reference_file)
    ctx.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        pseudo_gen = add_child(generate_group, "pseudohessian")
        pseudo_output, _ = _run_leaf(
            ctx,
            pseudo_gen,
            source_input=generate_input,
            runtime_input=runtime_generate,
            solver_method="pseudohessian",
            label="micelle guess generate, mode=pseudohessian",
            detail="embedded initial_guess generated in output JSON",
            require_initial_guess=True,
        )

        pseudo_use = add_child(use_group, "pseudohessian")
        if not pseudo_output:
            _set_failure(pseudo_use, "ERROR: skipped because pseudohessian guess generation failed")
        else:
            _run_leaf(
                ctx,
                pseudo_use,
                source_input=use_input,
                runtime_input=runtime_use,
                solver_method="pseudohessian",
                label="micelle guess use, mode=pseudohessian",
                compare=lambda out: expect_profiles(reference_file, out, 1e-8),
            )

        diis_gen = add_child(generate_group, "DIIS", required_for_parent=False)
        diis_generate_output, _ = _run_leaf(
            ctx,
            diis_gen,
            source_input=generate_input,
            runtime_input=runtime_generate,
            solver_method="DIIS",
            label="micelle guess generate, mode=DIIS",
            detail="embedded initial_guess generated in output JSON",
            accepted=True,
            require_initial_guess=True,
        )

        diis_use = add_child(use_group, "DIIS", required_for_parent=False)
        if not diis_generate_output:
            diis_status = _set_failure(diis_use, "skipped because DIIS guess generation failed", accepted=True)
        else:
            diis_use_output, _ = _run_leaf(
                ctx,
                diis_use,
                source_input=use_input,
                runtime_input=runtime_use,
                solver_method="DIIS",
                label="micelle guess use, mode=DIIS",
                accepted=True,
                require_output=False,
            )
            diis_status = diis_use.details
            if diis_use_output is None:
                pass
            elif not diis_use_output.is_file() or not diis_use_output.stat().st_size:
                diis_use.passed = False
                diis_status = diis_use.details = "completed but JSON output file is missing (accepted)"
            elif not pseudo_use.passed or not reference_file.is_file() or not reference_file.stat().st_size:
                diis_use.passed = False
                diis_status = diis_use.details = "completed but no pseudohessian baseline was available (accepted)"
            else:
                try:
                    diis_status = diis_use.details = expect_profiles(reference_file, diis_use_output, 1e-9, "passed and matched pseudohessian")
                except Exception as exc:
                    diis_status = _set_failure(diis_use, exc, accepted=True)

        root.details = f"diis={diis_status}"
        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, ctx.clean_output)


def test_micelle_grand_canonical_search(ctx: Context) -> ReportNode:
    root = ReportNode("micelle grand canonical search")
    search_input = ctx.tests_dir / "micelle_gc_search.in"
    utility = ctx.repo_root / "utils" / "micelle_gc_search.py"
    cleanup_targets: list[Path] = []

    require_file(ctx.binary, executable=True)
    require_file(search_input)
    require_file(utility)
    ctx.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        def run_case(node: ReportNode, require_initial_guess: bool, settings: dict[str, str] | None = None) -> None:
            case_name = node.label.replace(" ", "_")
            seed_runtime = ctx.output_dir / f"{search_input.stem}_{case_name}_seed.in"
            search_runtime = ctx.output_dir / f"micelle_gc_search_{case_name}.in"
            search_workdir = ctx.output_dir / search_runtime.stem
            summary_file = search_workdir / "summary.json"
            result_json = search_workdir / "result.output.json"
            cleanup_targets.extend([*solver_artifacts(seed_runtime), search_runtime, search_workdir])

            _, seed_output = _run_solver(
                ctx,
                node,
                search_input,
                seed_runtime,
                "pseudohessian",
                f"micelle gc seed generation ({node.label})",
                settings=settings,
                require_initial_guess=require_initial_guess,
            )
            if not require_initial_guess:
                seed_problem = last_problem(seed_output)
                if "initial_guess" in seed_problem:
                    raise TestError(f"ERROR: seed JSON unexpectedly contains embedded initial_guess: {seed_output}")
                for key in ("method", "mx", "my", "mz", "fjc", "charged", "monlist", "statelist", "profiles"):
                    if key in seed_problem:
                        raise TestError(f"ERROR: unexpected fallback seed field '{key}' in output JSON: {seed_output}")

            shutil.copy2(search_input, search_runtime)
            metrics = run_command(
                [
                    sys.executable,
                    str(utility),
                    str(search_runtime),
                    "--binary",
                    str(ctx.binary),
                    "--molecule",
                    "surf",
                    "--seed-json",
                    str(seed_output),
                    "--step",
                    "5",
                    "--workers",
                    "2",
                    "--max-iter",
                    "6",
                    "--gp-tol",
                    "5e-2",
                ],
                cwd=ctx.repo_root,
                quiet=ctx.quiet,
            )
            node.add_metric(metrics)
            if metrics.returncode:
                raise TestError(metrics.output.strip() or f"ERROR: external micelle GC search failed ({node.label})")
            if not summary_file.is_file():
                raise TestError(f"ERROR: missing micelle GC summary: {summary_file}")
            if not result_json.is_file():
                raise TestError(f"ERROR: missing micelle GC result JSON: {result_json}")

            summary = json.loads(summary_file.read_text(encoding="utf-8"))
            x = float(summary["x"])
            gp = float(summary["grand_potential"])
            if not bool(summary["converged"]):
                raise TestError("ERROR: micelle GC search did not report convergence")
            if abs(gp) > 0.25:
                raise TestError(f"ERROR: micelle GC search residual too large: grand_potential={gp}")
            if not 109.5 < x < 110.5:
                raise TestError(f"ERROR: micelle GC search returned unexpected aggregation number: n={x}")
            node.details = f"n={x:.6f}, gp={gp:.3e}"

        for label, require_initial_guess, settings in (
            ("embedded initial_guess seed", True, {"sys : noname : write_initial_guess": "true"}),
            ("u-profile seed", False, None),
        ):
            node = add_child(root, label)
            try:
                run_case(node, require_initial_guess, settings)
            except Exception as exc:
                _set_failure(node, exc)

        _finalize_report_tree(root)
        return root
    finally:
        _cleanup(cleanup_targets, ctx.clean_output)


def test_particle_in_cyl_coordinates(ctx: Context) -> ReportNode:
    return _run_method_group(
        ctx,
        root_label="particle in cyl coordinates",
        input_file=ctx.tests_dir / "particle_in_cyl_coordinates.in",
        reference_name="particle_in_cyl_coordinates.json.ref",
        runtime_stem="particle_in_cyl_coordinates",
        value_tol=1e-6,
        leaves=[("pseudohessian", "pseudohessian", "reference", True, False), ("DIIS", "DIIS", "pseudohessian", False, True)],
    )


def test_branched_brush(ctx: Context) -> ReportNode:
    root = ReportNode("branched brush")
    leaves = [("pseudohessian", "pseudohessian", "reference", True, False), ("DIIS", "DIIS", "pseudohessian", False, True)]
    for case_label, runtime_stem in (("2d cylindrical", "branched_brush_cyl2d"), ("1d planar", "branched_brush_planar_1d")):
        root.children.append(
            _run_method_group(
                ctx,
                root_label=case_label,
                label_stem=f"branched brush {case_label}",
                input_file=ctx.tests_dir / f"{runtime_stem}.in",
                reference_name=f"{runtime_stem}.json.ref",
                runtime_stem=runtime_stem,
                value_tol=2e-6,
                method_group_label=None,
                leaves=leaves,
            )
        )
    _finalize_report_tree(root)
    return root


def test_polE_regression(ctx: Context) -> ReportNode:
    return _run_method_group(
        ctx,
        root_label="polE regression",
        label_stem="polE",
        input_file=ctx.tests_dir / "polE.in",
        reference_name="polE.json.ref",
        runtime_stem="polE",
        value_tol=1e-5,
        leaves=[("pseudohessian", "pseudohessian", "reference", True, False), ("DIIS", "DIIS", "pseudohessian", False, True)],
    )


def test_external_potential(ctx: Context) -> ReportNode:
    root = ReportNode("external potential")
    leaves = [("pseudohessian", "pseudohessian", "reference", True, False), ("DIIS", "DIIS", "pseudohessian", True, False)]
    for axis in ("1d", "2d", "3d"):
        stem = f"external_potential_{axis}"
        root.children.append(
            _run_method_group(
                ctx,
                root_label=axis,
                label_stem=f"external potential {axis}",
                input_file=ctx.tests_dir / f"{stem}.in",
                reference_name=f"{stem}.json.ref",
                runtime_stem=stem,
                value_tol=1e-6,
                method_group_label=None,
                copied_inputs=[ctx.tests_dir / f"{stem}_external_potential.json"],
                leaves=leaves,
            )
        )
    _finalize_report_tree(root)
    return root


def _print_table(results: list[ReportNode]) -> None:
    rows = [
        (
            root_index,
            [
                str(root_index) if is_root else "",
                label,
                "PASS" if node.passed else "FAIL",
                f"{node.wall_s:.3f}",
                "NA" if node.max_rss_kb is None else f"{node.max_rss_kb / 1024.0:.1f}",
            ],
            _failure_suffix(node),
        )
        for root_index, is_root, label, node in _flatten_nodes(results)
    ]
    headers = ["#", "Test", "Status", "Wall\n(s)", "MaxRSS\n(MiB)"]
    total_row = [
        "",
        "TOTAL",
        f"{sum(node.passed for node in results)}/{len(results)}",
        f"{sum(node.wall_s for node in results):.3f}",
        "NA"
        if not [node.max_rss_kb for node in results if node.max_rss_kb is not None]
        else f"{max(node.max_rss_kb for node in results if node.max_rss_kb is not None) / 1024.0:.1f}",
    ]
    widths = [max(len(part) for cell in column for part in cell.splitlines()) for column in zip(headers, *[row for _, row, _ in rows], total_row)]

    def print_row(parts: list[str], aligns: list[str], suffix: str = "") -> None:
        wrapped = [part.splitlines() for part in parts]
        for line_idx in range(max(len(cell) for cell in wrapped)):
            cells = []
            for col_idx, cell in enumerate(wrapped):
                text = cell[line_idx] if line_idx < len(cell) else ""
                cells.append(
                    text.center(widths[col_idx])
                    if aligns[col_idx] == "center"
                    else text.rjust(widths[col_idx])
                    if aligns[col_idx] == "right"
                    else text.ljust(widths[col_idx])
                )
            line = "│ " + " │ ".join(cells) + " │"
            print(line + (f"  {suffix}" if suffix and line_idx == 0 else ""))

    def border(left: str, mid: str, right: str) -> str:
        return left + mid.join("─" * (width + 2) for width in widths) + right

    print(border("┌", "┬", "┐"))
    print_row(headers, ["center"] * len(headers))
    print(border("├", "┼", "┤"))
    for idx, (root_index, row, mode) in enumerate(rows):
        if idx and root_index != rows[idx - 1][0]:
            print(border("├", "┼", "┤"))
        print_row(row, ["center", "left", "center", "right", "right"], mode)
    print(border("├", "┼", "┤"))
    print_row(total_row, ["center", "left", "center", "right", "right"])
    print(border("└", "┴", "┘"))


TEST_SPECS = [
    ("homopolymer-adsorption", "homopolymer adsorption", test_homopolymer_adsorption, True, ("homopolymer_adsorption",)),
    ("homopolymer-adsorption-benchmark", "homopolymer adsorption benchmark", benchmark_homopolymer_adsorption, False, ("homopolymer_adsorption_benchmark",)),
    ("external-potential", "external potential", test_external_potential, True, ("external_potential",)),
    ("frozen-range-input-file", "frozen range input file", test_frozen_range_input_file, True, ("frozen_range_input_file",)),
    ("micelle-self-assembly", "micelle self assembly", test_micelle_self_assembly, True, ("micelle_self_assembly",)),
    ("micelle-grand-canonical-search", "micelle grand canonical search", test_micelle_grand_canonical_search, True, ("micelle_grand_canonical_search",)),
    ("particle-in-cyl-coordinates", "particle in cyl coordinates", test_particle_in_cyl_coordinates, True, ("particle_in_cyl_coordinates",)),
    ("branched-brush", "branched brush", test_branched_brush, True, ("branched_brush", "branched_brush_cyl2d", "branched_brush_planar_1d")),
    ("polE-regression", "polE regression", test_polE_regression, True, ("polE_regression",)),
]

TESTS = {key: (label, fn) for key, label, fn, _, _ in TEST_SPECS}
ALIASES = {alias: key for key, _, _, _, aliases in TEST_SPECS for alias in aliases}
ENABLED_TESTS = [key for key, _, _, enabled, _ in TEST_SPECS if enabled]


def _resolve_selection(items: list[str]) -> list[str]:
    if not items:
        return ENABLED_TESTS
    selected: list[str] = []
    for raw in items:
        key = ALIASES.get(raw.strip(), raw.strip())
        if key not in TESTS:
            raise SystemExit(f"Unknown test '{raw}'. Valid: {', '.join(TESTS)}")
        if key not in selected:
            selected.append(key)
    return selected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run NAMICS regression tests and benchmarks.")
    parser.add_argument("tests", nargs="*", help="Tests to run (default: enabled list)")
    parser.add_argument("--list", action="store_true", help="List available tests and exit")
    parser.add_argument("--binary", type=Path, default=None, help="Path to namics binary")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for runtime/generated output")
    parser.add_argument("--verbose", action="store_true", help="Show solver stdout/stderr")
    parser.add_argument("--with-save-memory", action="store_true", help="Run homopolymer save_memory variants")
    parser.add_argument("--keep-artifacts", action="store_true", help="Keep generated runtime files")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    tests_dir = repo_root / "tests"

    if args.list:
        for key in TESTS:
            print(key)
        return 0

    output_dir = args.output_dir.resolve() if args.output_dir else repo_root / "output"
    if not args.keep_artifacts and output_dir.is_dir():
        shutil.rmtree(output_dir, ignore_errors=True)

    reference_archive = tests_dir / "reference" / "reference_files.tar.gz"
    require_file(reference_archive)

    with tempfile.TemporaryDirectory(prefix="namics_ref_") as tmp_ref_dir:
        shutil.unpack_archive(str(reference_archive), tmp_ref_dir)
        ctx = Context(
            repo_root=repo_root,
            tests_dir=tests_dir,
            output_dir=output_dir,
            binary=args.binary.resolve() if args.binary else repo_root / "bin" / "namics",
            quiet=not args.verbose,
            clean_output=not args.keep_artifacts,
            with_save_memory=args.with_save_memory,
            reference_dir=Path(tmp_ref_dir) / "reference",
            generated_input_manifest=output_dir / "generated_test_inputs.tsv",
        )

        results: list[ReportNode] = []
        for key in _resolve_selection(args.tests):
            label, fn = TESTS[key]
            try:
                node = fn(ctx)
            except TestError as exc:
                node = ReportNode(label, passed=False, details=str(exc))
            except Exception as exc:
                node = ReportNode(label, passed=False, details=f"ERROR: unexpected failure: {exc}")
            if node.label != label:
                node.label = label
            results.append(node)

        if ctx.clean_output:
            ctx.generated_input_manifest.unlink(missing_ok=True)
        _print_table(results)
        return 0 if all(node.passed for node in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
