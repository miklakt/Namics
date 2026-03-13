# `micelle_gc_search.py`

`utils/micelle_gc_search.py` runs NAMICS repeatedly, varies one restricted molecule quantity such as `mol : surf : n`, and searches for the value where the grand potential becomes zero.

## Usage

The script needs:

- a single-problem search template
- `--molecule`
- either `--seed-json` or `--seed-input`

The system name is taken from the input if present, otherwise `noname` is used. The work directory is always `output/<input_stem>/`.

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

1. The script gets the first `initial_guess` by copying `--seed-json` or by running `--seed-input` once.
2. It takes the first search value from `--initial`, otherwise from `mol : <molecule> : <quantity>` in the search template.
3. Each evaluation writes a small runtime input, asks NAMICS to emit `json : sys : <system> : grand_potential`, and reads `problems[-1]["sys_<system>_grand_potential"]` from the JSON output.
4. After every evaluation it promotes the best JSON file to `current_seed.json`, so the next outer iteration starts from the best available guess.
5. It expands left and right until the residual changes sign, then refines the bracket with a safeguarded false-position step and midpoint fallback.
6. It stops when `abs(grand_potential) <= gp_tol` or the bracket becomes very small.

## Files Written

The work directory is `output/<input_stem>/`. The script keeps:

- `history.tsv`
- `current_seed.json`
- `result.in`
- `result.json`
- `summary.json`

Each evaluation also leaves `run_XXX.in` and `run_XXX.json` inside worker subdirectories.
