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
