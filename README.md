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
make test-homopolymer-adsorption
make test-frozen-range-input-file
```
