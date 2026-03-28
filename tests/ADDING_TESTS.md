# Adding A Test

This file is the preferred pattern for extending [run_tests.py](/home/ml/Namics/tests/run_tests.py).

## Rules

- Keep runnable inputs and support assets in `tests/`.
- Keep the test honest: compare against stored references or explicit scalar checks.
- Do not add hidden behavior to the input. Only use `settings` or `comment_toggles` when the test is explicitly about that variant.
- Prefer less code. Reuse the existing runner shape before adding new shared helpers.
- Keep custom logic local to the test that needs it.

## If You Have An Input File And A Reference Output

1. Add the input file under `tests/`, for example `tests/my_case.in`.
2. Add any extra runtime asset under `tests/`, for example `tests/my_case_table.json`.
3. Add the reference JSON as `reference/my_case.json.ref` inside `tests/reference/reference_files.tar.gz`.
4. Add a small `test_*` function in `tests/run_tests.py`.
5. Register it in `TEST_SPECS`.

## Default Pattern

If the test is just "run NAMICS on one input and compare the output JSON profiles to a reference", use `_run_method_group`.

Important: the core does not inject `pseudohessian` or `DIIS` leaves for you anymore. The leaf list belongs to the test unit.

```python
def test_my_case(ctx: Context) -> ReportNode:
    return _run_method_group(
        ctx,
        root_label="my case",
        input_file=ctx.tests_dir / "my_case.in",
        reference_name="my_case.json.ref",
        runtime_stem="my_case",
        value_tol=1e-6,
        copied_inputs=[ctx.tests_dir / "my_case_table.json"],  # only if needed
        leaves=[
            ("pseudohessian", "pseudohessian", "reference", True, False),
            ("DIIS", "DIIS", "pseudohessian", False, True),
        ],
    )
```

Leaf tuple fields:

- printed leaf label
- solver method written into the runtime input
- what to compare against: `"reference"`, another earlier leaf label, or `None`
- whether the leaf is required for the parent result
- whether failures are accepted at the leaf level

What `_run_method_group` still gives you:

- cleanup of generated runtime files
- execution and comparison for the leaves you asked for
- the same report-tree style as the existing simple tests

Common patterns:

```python
leaves=[("pseudohessian", "pseudohessian", "reference", True, False)]
```

Only one required leaf.

```python
leaves=[
    ("pseudohessian", "pseudohessian", "reference", True, False),
    ("DIIS", "DIIS", "pseudohessian", False, True),
]
```

Required pseudohessian plus optional accepted DIIS.

```python
leaves=[
    ("pseudohessian", "pseudohessian", "reference", True, False),
    ("DIIS", "DIIS", "pseudohessian", True, False),
]
```

Required pseudohessian plus required DIIS.

## When To Use A Local Custom Test

Do not force everything into `_run_method_group`.

Use a local custom test when the workflow has:

- multiple phases
- generated seeds or chained runs
- scalar checks instead of profile-reference checks
- benchmark logging
- a custom subtest hierarchy

Current examples:

- `test_micelle_self_assembly`
- `test_micelle_zero_grand_potential`
- `_run_homopolymer_adsorption`

For those, keep helper code inside the test or next to it. Do not widen shared helpers unless at least two tests need the same thing.

The intended split is:

- `_run_method_group` for "one input, local leaf list, standard compare flow"
- local test code for anything with phases, seeds, custom extraction, or custom acceptance rules

## If The Reference Is A Few Scalars

If you only need to check a few values, keep the comparison local.

```python
out, _ = _run_leaf(...)
if out:
    value = float(last_problem(out)["sys"]["noname"]["grand_potential"])
    if abs(value) > 1e-3:
        _set_failure(leaf, f"ERROR: unexpected grand_potential={value}")
    else:
        leaf.details = f"grand_potential={value:.3e}"
```

Do not add a new generic framework just for one scalar-style test.

## Registering The Test

Add one tuple to `TEST_SPECS`:

```python
("my-case", "my case", test_my_case, True, ("my_case",))
```

Field order:

- CLI key
- printed root label
- function
- enabled by default or not
- aliases

## Naming And Style

- `root_label` should match the printed table label.
- `runtime_stem` should be short and filesystem-safe.
- `reference_name` should normally match `runtime_stem + ".json.ref"`.
- leaf labels should be short and readable because they become report-tree labels
- If a test is only one case, return one `ReportNode`.
- If a test is a family of cases, make one root node and append child cases with `add_child` or `_run_method_group(..., method_group_label=None)`.

## Support Files

If the input depends on another file at runtime, pass it through `copied_inputs`.

Example:

```python
copied_inputs=[ctx.tests_dir / "external_table.json"]
```

That means:

- require the file before the run
- copy it into `output/`
- clean it up with the rest of the generated runtime files

## Reference Archive

The runner unpacks `tests/reference/reference_files.tar.gz` and expects the archive to contain a top-level `reference/` directory.

When adding a new reference, keep that layout:

- `reference/my_case.json.ref`

## What Not To Do

- Do not add new makefile hooks.
- Do not store active test assets outside `tests/`.
- Do not leave dead compatibility paths behind.
- Do not add a new shared helper if one local `if` block is enough.
