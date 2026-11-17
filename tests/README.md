# Regression tests

```bash
python3 tests/run_tests.py              # everything
python3 tests/run_tests.py charged polE # cases whose name contains a filter
```

Each entry of `CASES` in `run_tests.py` is one input in this folder. It runs once per solver method with
`newton : isaac : method` and `tolerance : 1e-10` set by the runner, and every number in its output (all arrays
under `problems[-1]`, coordinates included) must match `reference/<case>.npz`: the same quantities, the same sizes,
and max |diff| <= `TOL = 1e-8`. One tolerance for the whole suite; nothing is fitted, resampled or rescaled.
`micelle_zero_grand_potential` additionally checks the root search in `utils/search_engine.py`.

## References

The `.npz` files hold flat float64 arrays named `key:name:prop` (`x`, `y`, `z` for coordinates), exactly as computed.

- Default: written by commit `1563590`, the fork point from upstream `leermakers/Namics`, built unchanged:
  ```bash
  git worktree add /tmp/namics-1563590 1563590 && make -C /tmp/namics-1563590 -j8
  python3 tests/make_reference.py /tmp/namics-1563590/bin/namics <case> ...
  ```
- `external_potential_*` (feature added after the fork) and `micelle_guess_use` (needs the initial guess, which the
  fork point cannot read) keep the references of their original tests, from the builds of March 2026.
- `polE` comes from a current build (`make_reference.py --json`), because the fork point is wrong there: its planar
  lattice gave the boundary ghost sites nonzero neighbour weights, which leaks into the solution once surface
  states are involved (4.6e-6 in the ion profiles). Check: the same system as a cylinder of radius R, whose code is
  unchanged since the fork, converges to the current planar result as 1/R (2.7e-7 at R = 1e7, 2.7e-8 at R = 1e8)
  and stays 4.6e-6 away from the fork point's.

Case origins: `poor_solvent_brush_*` are identical to records in `pore_poor_solvent/data` and also agree with
them to 1e-7; `charged_pore_*`, `angular_pore` and `sfb_empty_pore` are scaled-down inputs built with the
generators of `charged_pore`, `scf_pore_angular` and `sfb_empty_pore` (see the header of each input);
`charged_sphere_brush` covers spherical electrostatics.

## Adding a case

1. Put the input (with its `json :` selectors) and any files it reads in this folder.
2. Add a `Case` to `CASES`.
3. Write the reference with `make_reference.py` and the fork-point build; if that build cannot run the input,
   say in this file where the reference comes from.
