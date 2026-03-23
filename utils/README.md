# Utilities

## `preprocess_input.py`

`utils/preprocess_input.py` mirrors the NAMICS input preprocessing flow and writes compact `.input.json` files without running the solver.

## Usage

Write the default sibling `.input.json` file:

```bash
python3 utils/preprocess_input.py tests/homopolymer_adsorption.in
```

Write to an explicit path:

```bash
python3 utils/preprocess_input.py tests/frozen_range_input_file.in \
  --output /tmp/frozen_range_input_file.input.json
```

Print the preprocessed JSON to stdout:

```bash
python3 utils/preprocess_input.py tests/external_potential_1d.in --stdout
```

Keep normalized masks, external potentials, and `initial_guess` payloads as links to sidecar JSON files:

```bash
python3 utils/preprocess_input.py tests/micelle_guess_use.in --link-data
```

The utility mirrors the current C++ preprocessor behavior for:

- legacy `.in` to compact JSON conversion
- legacy pinned/frozen range unwrapping and mask file normalization
- embedded `initial_guess` normalization

In addition, it can normalize monomer `external_potential_filename` inputs into embedded `external_potential` arrays.

By default the output embeds normalized masks, external potentials, and `initial_guess` payloads directly in the emitted JSON. With `--link-data`, those normalized payloads are written under `<output_stem>.data/` and referenced from the emitted input instead.

## `micelle_gc_search.py`

`utils/micelle_gc_search.py` runs NAMICS repeatedly, varies one restricted molecule quantity such as `mol : surf : n`, and searches for the value where the grand potential becomes zero.

## Usage

The script needs:

- a single-problem search template
- `--molecule`
- either `--seed-json` or `--seed-input`
- `--seed-json` must point to a JSON file with embedded `initial_guess`

The system name is taken from the input if present, otherwise `noname` is used. The work directory is always `output/<input_stem>/`.
The `--seed-input` file is copied into a temporary seed run without extra edits.

Generate the first seed from an input:

```bash
python3 utils/micelle_gc_search.py tests/micelle_gc_search.in \
  --binary bin/namics \
  --molecule surf \
  --seed-input tests/micelle_guess_generate.in \
  --step 5 \
  --workers 2
```

Reuse an existing JSON seed:

```bash
python3 utils/micelle_gc_search.py tests/micelle_gc_search.in \
  --binary bin/namics \
  --molecule surf \
  --seed-json output/micelle_gc_seed.pseudohessian.json \
  --step 5 \
  --workers 2
```

Start from an explicit value:

```bash
python3 utils/micelle_gc_search.py tests/micelle_gc_search.in \
  --binary bin/namics \
  --molecule surf \
  --seed-input tests/micelle_guess_generate.in \
  --initial 120 \
  --step 5
```

The search template may end with a single trailing `start`. Any other `start` layout is rejected.

## Search Procedure

1. The script gets the initial guess from `--seed-json` or by running `--seed-input` once as provided.
2. It takes the first search value from `--initial`, otherwise from `mol : <molecule> : <quantity>` in the search template.
3. Each evaluation writes a small runtime input, asks NAMICS to emit `json : sys : <system> : grand_potential`, and reads `problems[-1]["sys_<system>_grand_potential"]` from the JSON output.
4. It expands left and right until the residual changes sign, then refines the bracket with a safeguarded false-position step and midpoint fallback.
5. It stops when `abs(grand_potential) <= gp_tol` or the bracket becomes very small.

## Files Written

The work directory is `output/<input_stem>/`. The script keeps:

- `history.tsv`
- `result.in`
- `result.output.json`
- `summary.json`

Each evaluation also leaves `run_XXX.in` and `run_XXX.output.json` inside worker subdirectories.
