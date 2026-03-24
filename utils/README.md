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
