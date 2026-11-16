## Namics

Self-consistent field simulation code.

This repository is a minimalist experimental version of Namics. It keeps only the core functionality needed for refactoring, cleanup, and experimentation, and does not aim to preserve the full original features.

For the original project, see: <https://github.com/leermakers/Namics>

## Build

```bash
git submodule update --init --recursive
make
```

## Run

```bash
./bin/namics <input-file>
```

A JSON input has one shape; paths inside it resolve relative to the input file:

```json
{
  "output": { "file": "results/run1", "append": false, "write_bounds": false,
              "select": { "mol": { "pol": "phi" }, "mon": { "*": ["phi", "u"] }, "sys": { "*": "grand_potential" } } },
  "problems": [ { "lat": {...}, "mon": {...}, "mol": {...}, "newton": {...} } ]
}
```

- Results go to `<input dir>/<stem>.output.json` (`foo.in`, `foo.json` and `foo.input.json` all give `foo.output.json`), or to `output.file` with `.output.json` appended if missing. Missing directories are created.
- The output is `{"input": <source path>, "problems": [...]}`; `append: true` extends `problems` of an existing file.
- JSON inputs are never rewritten. A legacy `.in` file is converted (`output : json : key : value` → `output.key`, `json : sec : name : prop` → `output.select`) and the result is saved as `<stem>.input.json`.

## Tests

```bash
make test-all
make test-homopolymer-adsorption
make benchmark-homopolymer-adsorption
make test-frozen-range-input-file
make test-micelle-self-assembly
make test-particle-in-cyl-coordinates

# optional variants / maintenance
python3 tests/run_tests.py homopolymer-adsorption-benchmark
python3 tests/run_tests.py homopolymer-adsorption --with-save-memory
python3 tests/run_tests.py --keep-artifacts
```

Reference files are stored packed in `tests/reference/reference_files.tar.gz` and are temporarily unpacked during test runs.
