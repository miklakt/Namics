## Namics

Self-consistent field simulation code.

## Build

```bash
git submodule update --init --recursive
make
```

## Run

```bash
./bin/namics <input-file>
```

## Tests

```bash
make test-all
make test-homopolymer-adsorption
make test-frozen-range-input-file
make test-micelle-self-assembly
make test-particle-in-cyl-coordinates

# optional variants / maintenance
python3 tests/run_tests.py homopolymer-adsorption --with-save-memory
python3 tests/run_tests.py --keep-artifacts
```

Reference files are stored packed in `tests/reference/reference_files.tar.gz` and are temporarily unpacked during test runs.
