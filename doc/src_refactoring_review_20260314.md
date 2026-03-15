# `src/` Refactoring Review

Date: 2026-03-14

Scope:
- Static review of `src/`
- Focus on duplication, dead code, stale code, refactoring candidates, module responsibilities, and responsibility chains
- No production code changes were made for this review

Inputs reviewed:
- `src/*.h`
- `src/*.cpp`
- `doc/pruning_review_20260311.md`

Clarified goal for the next phase:
- smallest possible build that still passes the default regression suite
- reduce total line count aggressively
- avoid architecture bloat
- do de-duplication and storage cleanup before large-scale interface design
- prefer deletion and consolidation over adding new abstraction layers

## 1. Refactoring Plan

The immediate target is not a cleaner architecture in the abstract. The immediate target is:

1. remove everything not needed for the default tested build
2. de-duplicate obvious repeated code
3. normalize storage and ownership with the smallest possible code
4. use standard algorithms and views where they reduce code
5. postpone new interfaces/base classes unless they clearly shrink the codebase

The code already has recognizable domains, but the next phase should treat that as a pruning and consolidation problem first, not as an invitation to add framework layers.

### Minimal-First Rules

- Prefer deletion over encapsulation.
- Prefer one helper or one local utility over a new class hierarchy.
- Prefer `std::vector`, `std::array`, `std::span`, and standard algorithms over custom storage wrappers.
- Do not introduce interfaces “for future backends” until there is a real second backend to support.
- Do not split a class unless the split reduces code or deletes unsupported feature branches.
- Keep hot numerical kernels straightforward; use standard algorithms where they make code shorter and clearer, not where they hide stencil logic.

### Phase 1: Shrink To The Default Tested Surface

This should use [`doc/pruning_review_20260311.md`](../doc/pruning_review_20260311.md) as the baseline.

Highest-value targets from that review:
- branched-molecule support if the default suite does not require it
- ring-specific paths if not covered
- `save_memory` paths if not part of default tests
- `delta_*`, `find_local_solution`, `compute_Gibbs_excess`, `compute_kJ0`, and related optional branches if they are not in the default tested surface
- any non-JSON output remnants

Rule for this phase:
- if a feature is not required by `make -j8` plus default `python3 tests/run_tests.py`, deletion is preferable to abstracting it

### Phase 2: Remove Stale APIs And Dead Scaffolding

Do this before any structural redesign.

Strong candidates already identified:
- remove stale declarations with no definitions or call sites
- remove empty lifecycle hooks that are never used
- remove duplicate aliases like `Lat` and `lat`
- remove unused state fields that have writes but no reads
- collapse duplicated control flow in input and output lookup paths

This phase should reduce line count immediately and make the remaining design easier to reason about.

### Phase 3: De-Duplicate Repeated Patterns

Only after pruning the feature surface should shared helpers be introduced.

Best early targets:

#### A. Parameter/output registry de-duplication

The repeated `push` / `PushOutput` / `GetValue` / `GetPointer` pattern is a strong consolidation target, but the consolidation should stay small.

Preferred approach:
- shared helper functions or a small utility type
- no big reflection framework
- no new inheritance tree just for output

#### B. Range parsing de-duplication

`pinned_range`, `frozen_range`, `restricted_range`, and `delta_range` all perform related parsing and mask construction.

Preferred approach:
- one compact parser for coordinate/range syntax
- one mask builder path
- small policy differences layered on top

#### C. Output traversal de-duplication

`Output::GetValue`, `GetPointer`, and `GetPointerInt` repeat the same object-graph traversal pattern.

Preferred approach:
- one shared lookup path
- keep string-driven lookup if it is smaller than introducing typed descriptors at this stage

### Phase 4: Storage Model Cleanup

Storage should be cleaned up with standard library types first.

Recommended storage policy:

- `std::vector<Real>` for owned dynamic numeric buffers
- `std::vector<int>` for owned dynamic integer buffers
- `std::array<int, 6>` for fixed-size coordinate/range tuples
- `std::span<Real>` / `std::span<const Real>` for views into buffers
- `std::unique_ptr<T[]>` only when a vector would meaningfully complicate layout or semantics

Avoid in this phase:
- a general `Field` class
- allocator hierarchies
- backend abstraction for GPU/shared memory

Those can come later if the code still needs them after pruning.

### Phase 5: Standard Algorithm Cleanup

After storage normalization, replace obvious hand-written loops only where the result is shorter or clearer.

Good candidates:
- repeated clear/fill/copy patterns
- reductions like sums and norms
- elementwise transforms
- duplicate append/filter/find logic in parser/output code

Do not force standard algorithms into:
- complex stencil kernels
- propagation code where explicit indexing is the clearest form

### Phase 6: Larger Interface Design, Only If Still Needed

Only after phases 1 through 5 should the project revisit bigger abstractions such as:
- more general lattice composition
- solver/problem interfaces
- reader/writer interfaces
- propagator objects

At that point the design problem will be smaller, because unsupported and duplicate code will already be gone.

### Recommended Execution Order

1. Prune to the default tested surface.
2. Remove stale declarations, no-op hooks, and dead fields.
3. Collapse duplicated parser/output/property plumbing.
4. Normalize storage with standard containers and spans.
5. Replace obvious utility loops with standard algorithms.
6. Re-evaluate whether larger interfaces still buy a net line-count reduction.

This order is the lowest-risk path to a smaller build and a smaller codebase.

## 2. Problems And Findings

### 2.1 `main()` is still a manual factory, manual lifecycle manager, and orchestration hub

Evidence:
- [`src/namics.cpp`](../src/namics.cpp#L95-L221)
- [`src/namics.cpp`](../src/namics.cpp#L260-L375)

Problems:
- constructs a temporary lattice just to discover `gradients` and `geometry`
- manually dispatches molecule subclasses
- manually loops through `Segment`, `State`, `Reaction`, `Molecule`, `System`, `Solve_scf`, `Output`
- owns many failure paths
- still mixes `unique_ptr` with raw-pointer vectors and manual `delete`

Refactoring implication:
- `main()` knows too much about every module’s construction order and subclass selection
- this blocks interface-driven design because all variation is hard-coded at the top level

### 2.2 `System` is a god object

Evidence:
- header size and breadth: [`src/system.h`](../src/system.h#L10-L138)
- residual path: [`src/system.cpp`](../src/system.cpp#L1516-L1614)
- density assembly path: [`src/system.cpp`](../src/system.cpp#L1616-L1903)
- thermodynamics path: [`src/system.cpp`](../src/system.cpp#L1905-L2387)

Current responsibilities include:
- owning many work buffers
- building iteration lists
- applying constraints
- electrostatics
- calling molecule propagation
- normalizing free/restricted/solvent/neutralizer molecules
- block heuristics for `find_local_solution`
- computing free energy, grand potential, curvature-related quantities
- exposing output reflection data

Refactoring implication:
- this class is the main reason interfaces are currently too large
- it should be split into at least:
  - `ScfState`
  - `ResidualEvaluator`
  - `ConstraintHandler`
  - `ElectrostaticsModel`
  - `ThermodynamicsCalculator`

### 2.3 The lattice abstraction is oversized and encodes too many variation axes in inheritance

Evidence:
- base interface: [`src/lattice.h`](../src/lattice.h#L104-L137)
- concrete families: [`src/LGrad1.h`](../src/LGrad1.h), [`src/LGrad2.h`](../src/LGrad2.h), [`src/LGrad3.h`](../src/LGrad3.h), [`src/LG1Planar.h`](../src/LG1Planar.h), [`src/LG2Planar.h`](../src/LG2Planar.h)

Problems:
- one interface mixes shape, metric, stencil, electrostatics, range parsing, and boundary mutation
- concrete classes are split by both dimension and geometry, which is why so much code repeats
- several combinations are explicitly “not implemented” or “not certified”, which shows the current inheritance tree does not scale cleanly

Refactoring implication:
- use composition of shape + metric + boundary + propagation policies
- treat “planar vs cylindrical vs spherical” and “1D vs 2D vs 3D” as independent axes where possible

### 2.4 `Molecule` mixes parsing, topology, workspace allocation, propagation, and reporting

Evidence:
- broad interface: [`src/molecule.h`](../src/molecule.h#L94-L121)
- topology parser: [`src/molecule.cpp`](../src/molecule.cpp#L368-L616)
- reporting/output: [`src/molecule.cpp`](../src/molecule.cpp#L726-L893)
- propagation helpers: [`src/molecule.cpp`](../src/molecule.cpp#L1047-L1317)

Problems:
- the base class already knows about linear and branched parsing
- topology parsing is tightly coupled to the execution format used by propagators
- propagation functions expose storage details and traversal counters (`s`, `generation`, `M`) directly
- output statistics are mixed with computational state

Refactoring implication:
- parse topology once into an explicit graph/tree representation
- make propagators consume that representation
- keep statistics in a separate analyzer/report object

### 2.5 Strong cross-cutting duplication in parameter and output plumbing

Evidence:
- same registry pattern appears in most core classes:
  - [`src/state.cpp`](../src/state.cpp#L118-L200)
  - [`src/reaction.cpp`](../src/reaction.cpp#L143-L221)
  - [`src/segment.cpp`](../src/segment.cpp#L777-L1016)
  - [`src/molecule.cpp`](../src/molecule.cpp#L703-L946)
  - [`src/lattice.cpp`](../src/lattice.cpp#L659-L762)
  - [`src/system.cpp`](../src/system.cpp#L933-L1168)
  - [`src/solve_scf.cpp`](../src/solve_scf.cpp#L205-L305)

Problems:
- same `push` overloads
- same clear-and-rebuild `PushOutput()`
- same linear lookup in `GetValue()`
- mostly null or ad hoc `GetPointer()` and `GetPointerInt()`

Refactoring implication:
- the project is paying the maintenance cost of a reflection system without actually centralizing it
- this should become one shared typed inspection API

### 2.6 Output is stringly-typed and tightly coupled to every module’s internals

Evidence:
- [`src/output.cpp`](../src/output.cpp#L230-L438)
- [`src/output.cpp`](../src/output.cpp#L442-L659)

Problems:
- output selection works by string keys, string names, and string-encoded profile descriptors
- output must know whether a value lives in `System`, `Molecule`, `Segment`, `State`, `Lattice`, or `Solve_scf`
- profile addressing uses encoded strings such as `"profile;0"`
- wildcard molecule expansion is custom logic inside `Output::Load()`

Refactoring implication:
- output currently acts as a second reflective traversal engine layered on top of the first
- typed descriptors or visitor-style inspection would drastically reduce coupling and runtime-only errors

### 2.7 There is substantial parser duplication and duplicated control flow

Evidence:
- duplicate `newton` list check: [`src/input.cpp`](../src/input.cpp#L470-L476)
- `pinned_range` parsing and keyword expansion: [`src/segment.cpp`](../src/segment.cpp#L95-L217)
- `frozen_range` parsing and keyword expansion: [`src/segment.cpp`](../src/segment.cpp#L226-L458)
- many repeated key/name/property switches in output lookup: [`src/output.cpp`](../src/output.cpp#L230-L438)

Problems:
- `pinned_range` and `frozen_range` are structurally similar but maintained separately
- output lookup repeats nearly the same traversal for ints, reals, bools, strings, and profiles
- input parsing still uses many manual token loops where a typed intermediate representation would be simpler

Refactoring implication:
- range parsing should become one reusable parser plus policy differences
- output lookup should become data-driven, not switch-driven

### 2.8 Memory ownership is inconsistent and harder than necessary

Evidence:
- `malloc/free` in `System`, `Segment`, `Molecule`, `SFNewton`
- `new[]/delete[]` in `Lattice` and `Solve_scf`
- alias pairs like `H_phi`/`phi`, `H_psi`/`psi`, `H_BETA`/`BETA`, `H_beta`/`beta`

Problems:
- ownership is implicit
- lifetime is manual
- aliasing hides which pointer actually owns storage
- this blocks backend abstraction for shared memory or GPU

Refactoring implication:
- storage should be explicit and value-like
- views should be non-owning and clearly named

### 2.9 Dead code and stale API candidates

These are the strongest static candidates for removal or consolidation.

#### Clear stale declarations

- `Solve_scf` declares `ComputePhis(bool)`, `PutU()`, and `Put_U()` in the header, but there are no matching implementations or call sites for those zero-argument variants:
  - [`src/solve_scf.h`](../src/solve_scf.h#L98-L103)
- `Molecule` declares `Theta(void)` but there is no implementation or call site:
  - [`src/molecule.h`](../src/molecule.h#L101-L104)

#### No-op lifecycle hooks

- `State::AllocateMemory()` and `State::PrepareForCalculations()` are empty and unused:
  - [`src/state.cpp`](../src/state.cpp#L15-L23)
- `Reaction::AllocateMemory()` and `Reaction::PrepareForCalculations()` are empty and unused:
  - [`src/reaction.cpp`](../src/reaction.cpp#L16-L24)

These look like legacy interface remnants rather than active behavior.

#### Parsed but not used configuration

- `StoreFileGuess` and `ReadFileGuess` are parsed but not otherwise used:
  - declaration: [`src/solve_scf.h`](../src/solve_scf.h#L49-L55)
  - parsing only: [`src/solve_scf.cpp`](../src/solve_scf.cpp#L178-L179)

#### Likely dead state fields

- `Segment::prepared` is written in the constructor and reset in `main()`, but no reads were found:
  - declaration: [`src/segment.h`](../src/segment.h#L26-L30)
  - writes: [`src/segment.cpp`](../src/segment.cpp#L5), [`src/namics.cpp`](../src/namics.cpp#L322-L323)
- `Segment::used_in_mol_nr` is initialized but no read/write after construction was found:
  - declaration: [`src/segment.h`](../src/segment.h#L55-L60)
  - initialization: [`src/segment.cpp`](../src/segment.cpp#L26)

These should be verified before deletion, but they are strong dead-code candidates.

### 2.10 Duplicate aliases and naming drift increase cognitive load

Evidence:
- duplicate lattice pointers in several classes:
  - [`src/system.h`](../src/system.h#L20-L25)
  - [`src/solve_scf.h`](../src/solve_scf.h#L40-L46)
  - [`src/segment.h`](../src/segment.h#L13-L16)
  - [`src/molecule.h`](../src/molecule.h#L15-L18)

Problems:
- both `Lat` and `lat` exist with the same meaning
- storage aliases use inconsistent prefixes (`H_`, none, uppercase duplicate)
- naming drift makes it hard to tell ownership from view semantics

Refactoring implication:
- one canonical member name per dependency
- one naming convention for owned storage vs views

### 2.11 Several features remain partially implemented inside the mainline architecture

Evidence:
- unsupported combinations are kept in normal code paths:
  - [`src/molecule.cpp`](../src/molecule.cpp#L320-L348)
  - [`src/LGrad2.cpp`](../src/LGrad2.cpp)
  - [`src/LGrad3.cpp`](../src/LGrad3.cpp)

Problems:
- “not implemented”, “untested territory”, and “not certified” branches are still part of the core abstraction
- this makes interfaces larger than the actually trusted behavior

Refactoring implication:
- move experimental combinations behind explicit capability checks
- reduce the default abstraction surface to tested, supported behavior

## 3. Module Responsibility Reference

### Entry and shared types

- `src/namics.h`
  - global numeric types/constants and enums
  - should eventually hold only small shared type definitions
- `src/namics.cpp`
  - current composition root, factory selection, iteration over `start` problems
  - should become only “read config -> build problem -> solve -> write result”

### Parsing and I/O support

- `src/input.h`, `src/input.cpp`
  - custom input parsing, includes, validation, list extraction, path resolution, output-item parsing
  - should become a reader that emits typed problem config
- `src/output.h`, `src/output.cpp`
  - resolves requested quantities and writes JSON payloads
  - should become a serializer over typed result objects
- `src/json_writer.h`, `src/json_writer.cpp`
  - file append/merge mechanics for JSON output
  - already reasonably focused
- `src/io_utils.h`
  - low-level JSON parsing helpers for initial guesses and external profiles
  - should remain low-level utility code, not domain logic
- `src/debug_log.h`
  - debug logging macro and formatting

### Solver layer

- `src/sfnewton.h`, `src/sfnewton.cpp`
  - generic nonlinear solver algorithms, Hessian updates, DIIS, line search
  - should remain solver-infrastructure code
- `src/solve_scf.h`, `src/solve_scf.cpp`
  - adapter from SCF problem to solver engine plus solver-specific config
  - should become a thin controller over `IScfProblem` and `INonlinearSolver`

### Domain model

- `src/segment.h`, `src/segment.cpp`
  - monomer species, site masks, per-site fields, state metadata, constraints
  - should focus on segment properties and fields, not range parsing or output reflection
- `src/state.h`, `src/state.cpp`
  - internal state metadata for segments, `alphabulk`, `valence`, `chi`
  - currently more metadata than active computational object
- `src/reaction.h`, `src/reaction.cpp`
  - reaction equation parsing and bulk-state closure logic
  - should be a chemistry constraint object, not a mini-reflection object
- `src/molecule.h`, `src/molecule.cpp`
  - molecule composition parsing, topology, workspace allocation, propagation scaffolding, output stats
  - should split into topology + propagator + report data
- `src/mol_linear.h`, `src/mol_linear.cpp`
  - linear-chain propagation algorithm
- `src/mol_branched.h`, `src/mol_branched.cpp`
  - branched-chain propagation algorithm
- `src/system.h`, `src/system.cpp`
  - SCF problem assembly, residual evaluation, normalization, electrostatics, thermodynamics
  - should be decomposed into smaller services

### Lattice and low-level kernels

- `src/lattice.h`, `src/lattice.cpp`
  - base lattice config, indexing, memory, some generic helpers, large abstract interface
  - should become a composition root for shape/stencil/metric/boundary collaborators
- `src/LGrad1.h`, `src/LGrad1.cpp`
  - 1-gradient non-planar lattice implementation
- `src/LGrad2.h`, `src/LGrad2.cpp`
  - 2-gradient lattice implementation
- `src/LGrad3.h`, `src/LGrad3.cpp`
  - 3-gradient lattice implementation
- `src/LG1Planar.h`, `src/LG1Planar.cpp`
  - 1-gradient planar specialization
- `src/LG2Planar.h`, `src/LG2Planar.cpp`
  - 2-gradient planar specialization
- `src/tools_host.h`
  - host-side low-level array kernels: boundary copying, DIIS helpers, block distribution

## 4. Responsibility Chains

### 4.1 Build / startup chain

1. `main()` creates `Input`
2. `Input` parses the file, resolves `include`, validates keywords, and builds per-start name lists
3. `main()` constructs a provisional lattice to discover `gradients` and `geometry`
4. `main()` selects the final concrete lattice type
5. `main()` constructs `Segment`, then `State`, then `Reaction`, then `Molecule`
6. `State::CheckInput()` feeds state definitions back into `Segment::AddState()`
7. `main()` constructs `System`
8. `main()` constructs `Solve_scf`
9. `main()` optionally constructs `Output`

Why this matters:
- startup order is not derived from interfaces; it is encoded procedurally in `main()`
- several objects depend on side effects from earlier objects rather than explicit build products

### 4.2 Iteration / SCF residual chain

1. `Solve_scf::Solve()` chooses a nonlinear method
2. solver engine calls `Solve_scf::residuals()`
3. `Solve_scf::residuals()` delegates classical SCF work to `System::Classical_residual()`
4. `System::Classical_residual()` copies iteration variables, calls `ComputePhis(x, ...)`, and then assembles the residual
5. `System::ComputePhis(x, ...)` unpacks iteration variables through `PutU()`
6. `System::PrepareForCalculations()` prepares lattice, segments, and molecules
7. `System::ComputePhis(residual)` triggers molecule propagation and normalization
8. each `Molecule::ComputePhi()` calls propagation helpers that eventually delegate to lattice propagation kernels
9. densities are accumulated back into segments and the system total
10. electrostatics and constraints are applied
11. the solver receives the residual vector and continues iterating

Why this matters:
- the actual “SCF problem” is spread across `Solve_scf`, `System`, `Molecule`, `Segment`, and `Lattice`
- this is the exact chain that should be hidden behind `IScfProblem`

### 4.3 State / reaction chain

1. `State::CheckInput()` validates a state and registers it into its owning segment
2. `Reaction::CheckInput()` parses equations in terms of state names and maps them back to segments and state slots
3. before SCF solving, `Solve_scf::Solve()` may run a weak nonlinear solve over reactions
4. `Reaction::GuessAlpha()` and `Reaction::PutAlpha()` mutate segment state fractions
5. segment state fractions then affect `Segment::PrepareForCalculations()` and downstream propagation

Why this matters:
- chemistry closure is currently implemented as mutations on segment internals rather than a separate chemistry model
- the ownership direction is backwards: reactions know too much about segment storage layout

### 4.4 Output chain

1. after solving, `Solve_scf::PushOutput()` forces every major object to populate its output registry
2. `Output::Load()` parses user-selected output quantities
3. `Output::GetValue()` / `GetPointer()` / `GetPointerInt()` walk through the object graph using string keys
4. `Output::WriteOutput()` converts scalars and profiles into JSON arrays
5. `JsonWriter::WriteProblem()` appends or rewrites the target file

Why this matters:
- there is no typed result model
- output works only because each runtime object exposes reflective string registries

### 4.5 Data ownership chain

Current practical ownership looks like this:

1. `main()` owns object lifetimes
2. `System` allocates global SCF work buffers
3. `Segment` allocates per-segment fields
4. `Molecule` allocates per-molecule propagation buffers
5. `Lattice` allocates geometry/stencil arrays
6. `Solve_scf` allocates iteration vectors

Why this matters:
- ownership is split by implementation convenience, not by stable abstraction boundaries
- that is exactly why raw pointers leak into every public interface

## 5. Best-Practice Minimalist Refactoring Candidates

These are the highest-value changes under the clarified “smallest build, smallest code” goal.

### 5.1 Delete unsupported or untested features instead of abstracting them

Purpose:
- reduce build surface immediately

Benefit:
- maximum line-count reduction
- fewer branches to preserve during later cleanup

### 5.2 Replace repeated registries with one very small shared utility

Purpose:
- remove duplicated `push` / `GetValue` / `PushOutput` logic

Benefit:
- large de-duplication win without adding a heavyweight reflection system

### 5.3 Unify range parsing and mask construction

Purpose:
- stop maintaining separate but similar parsers for pinned/frozen/restricted/delta syntax

Benefit:
- reduces parser complexity and duplicated validation code

### 5.4 Normalize storage with standard containers, not custom wrapper classes

Purpose:
- make ownership explicit with less code

Benefit:
- safer code and lower line count than the current raw-pointer + alias model

### 5.5 Use standard algorithms in utility/control code, not as a blanket rule

Purpose:
- shorten repeated loops in setup, parser, output, and reductions

Benefit:
- cleaner code without obscuring the numerical kernels

### 5.6 Introduce bigger abstractions only if they delete more code than they add

Purpose:
- prevent “refactoring into frameworks”

Benefit:
- keeps the project aligned with the minimal-build goal

## 6. Immediate Priorities

Given the clarified direction, the next practical sequence should be:

1. use the pruning review to remove default-suite-unused features
2. remove stale declarations, empty hooks, and dead fields
3. de-duplicate output/property plumbing
4. de-duplicate range parsing and mask creation
5. convert owned raw buffers to standard containers where this reduces code
6. only then revisit whether any new interface or base class is still justified

The key rule is simple:
- if a proposed abstraction does not produce a net reduction in code size or tested complexity, defer it
