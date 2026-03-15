# Pruning Review 2026-03-11

Timestamp:
- UTC: `2026-03-11T14:29:12Z`
- Europe/Paris: `2026-03-11T15:29:12+0100`

Baseline verification:
- `make -j8` passed.
- `python3 tests/run_tests.py` passed.
- The default test selection excludes the homopolymer benchmark (`tests/run_tests.py:985-1010`).

Goal of this review:
- Identify the smallest code and build surface that still passes the current default regression suite.
- Record what is actually exercised so later pruning work does not mistake compiled legacy code for tested behavior.

## Findings

### 1. `main()` treats many setup and input failures as success, and it still relies on manual raw-pointer cleanup.

Evidence:
- [src/namics.cpp](../src/namics.cpp) returns `0` for many failures during input parsing and object setup (`lines 95-221`, especially `96-99`, `111`, `117-120`, `147`, `169-183`, `191-221`, `271-273`).
- The same function owns `Segment`, `State`, `Reaction`, and `Molecule` instances via raw-pointer vectors and only deletes them in the happy-path tail (`lines 150-211`, `351-371`).

Why this matters for pruning:
- If a pruning change breaks setup logic, the binary can still exit with code `0`, which weakens failure detection outside the JSON-producing tests.
- Early returns skip the manual delete loops, so cleanup behavior depends on staying on the success path. That makes aggressive feature removal harder to reason about safely.

Recommendation:
- Before large-scale code deletion, convert `main()` to RAII-owned containers and return non-zero on setup/input failure.

### 2. The regression suite intentionally accepts several DIIS failures, so “tests pass” is weaker than “DIIS still works broadly”.

Evidence:
- Optional DIIS leaves are explicitly marked non-required in the generic runner (`tests/run_tests.py:355-392`, `395-470`).
- Homopolymer DIIS is non-gating (`tests/run_tests.py:608-627`).
- Micelle DIIS generation/use is non-gating and accepted on failure (`tests/run_tests.py:783-824`).
- Current baseline run still passes overall while reporting:
  - homopolymer `chi_Si = -2/-4/-6`: DIIS failed with `Detected invalid n in computed solver state for molecule pol.`
  - micelle guess use: DIIS failed with `numerical drift, failed at tol < 2.562921097570836e-09`

Why this matters for pruning:
- A pruning pass that keeps only the DIIS behavior needed by `external-potential` can still leave the default suite green.
- If the long-term product goal still includes robust DIIS for polymer and micelle cases, the current default suite is not sufficient as the sole gate.

Recommendation:
- Treat pseudohessian as mandatory.
- Treat DIIS outside the external-potential cases as “currently tolerated but semantically live” unless the user explicitly agrees to prune it.

### 3. The Makefile compiles every `src/*.cpp`, so untested legacy code remains part of the minimal build until it is removed or the source list is made explicit.

Evidence:
- [Makefile](../Makefile) discovers sources with `find` and compiles all of them (`lines 100-101`).
- Branched-molecule support is still always built via [src/mol_branched.cpp](../src/mol_branched.cpp) and selected from [src/namics.cpp](../src/namics.cpp) `203-208`.
- Branched parsing still lives in [src/molecule.cpp](../src/molecule.cpp) `531-590`, where `[` / `]` switch the molecule type to `branched`.
- No current test input uses square-bracket composition syntax; the test inputs exercise monomer and linear molecules only.

Why this matters for pruning:
- “Unused by tests” does not reduce build size or compile time unless the code is actually deleted or omitted from the build list.
- Branched support is one of the clearest current candidates for removal in a strict “default-tests-only” build.

Recommendation:
- If the target is truly the smallest build that passes the default suite, branched molecules are a high-value candidate for complete removal.
- Make that change coherently: parser, factory selection, and source list together.

### 4. `Input::MakeLists()` still contains copy-paste scaffolding that is a low-risk cleanup target.

Evidence:
- [src/input.cpp](../src/input.cpp) calls the same `newton` cardinality check twice at `473` and `475`.

Why this matters for pruning:
- This is not a functional blocker, but it shows there is still low-risk legacy control flow in core setup code.
- Similar duplicated checks and compatibility leftovers are good follow-up cleanup targets after correctness/ownership fixes.

Recommendation:
- Remove obviously duplicated validation first when simplifying parser/setup code.

## Exercised Surface

The current default `python3 tests/run_tests.py` suite exercises:

- Lattice variants
  - 1-gradient planar simple cubic
  - 1-gradient spherical hexagonal with `FJC_choices`
  - 2-gradient planar simple cubic
  - 2-gradient cylindrical hexagonal
  - 3-gradient simple cubic through `LGrad3`
- Solver/output behavior
  - pseudohessian
  - DIIS, but only strictly required by `external-potential`
  - JSON output writing and JSON append handling
  - initial guess write/read through JSON
- Segment/system features
  - `frozen_range`
  - `frozen_filename`
  - `external_potential_filename`
  - `pinned_range`
  - charged states / `valence`
  - reactions / equilibria
  - `neutralizer`
  - `lowerbound : surface`
  - `bondlength`
- Molecule shapes
  - monomer molecules
  - linear molecules

## Not Exercised By The Default Suite

These surfaces appear untested by the default regression command and are therefore top pruning candidates or at least poor places to rely on tests alone:

- branched molecules (`mol_branched`, square-bracket composition parsing)
- ring molecule behavior
- `save_memory` paths
  - only covered by `python3 tests/run_tests.py --with-save-memory`
- `range_restricted`
- `fill_range`
- `restricted_range`
- `delta_range`
- `delta_inputfile`
- `delta_molecules`
- `find_local_solution`
- `compute_Gibbs_excess`
- `compute_kJ0`
- non-JSON output formats
  - input parsing already rejects non-JSON outputs in [src/input.cpp](../src/input.cpp) `419-449`

## Suggested Pruning Order

1. Fix `main()` ownership and exit-status behavior first.
2. Decide whether the target is:
   - the smallest build that passes the current default suite, or
   - the smallest build that also preserves broader DIIS behavior.
3. If the target is the default suite only, remove branched-molecule support as one coherent change.
4. Remove default-suite-unused system/molecule options next (`delta_*`, `find_local_solution`, ring-specific paths, etc.).
5. Replace the Makefile’s `find`-based source discovery with an explicit source list once the feature set is trimmed.
6. Only prune DIIS further if the user explicitly accepts losing non-external-potential DIIS behavior, because the current suite will not fully protect it.

## Practical Notes For Future Agents

- Read this file before deleting any lattice, molecule, or solver code for “minimal build” work.
- Do not infer that a feature is required just because it compiles; check whether it is exercised above.
- Do not infer that DIIS is fully protected by the default suite; it is not.
- After any structural prune, re-run:
  - `make -j8`
  - `python3 tests/run_tests.py`
