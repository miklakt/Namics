# `src/` Deletion Log

Date: 2026-03-14

Purpose:
- explain every deletion made in the first `src/` shrink pass
- record why each removed item was considered safe
- keep the reasoning separate from the broader refactoring review

Verification after the deletion pass:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

Net change for this pass:
- 101 deleted lines
- 24 added lines

## `src/input.cpp`

- Deleted the second `test_count(NewtonList, "newton", ...)` call in `Input::MakeLists`.
- Reason: it was an exact duplicate of the immediately preceding validation call and had no distinct effect.

## `src/molecule.h`

- Deleted `Lattice* Lat;`.
- Reason: `Molecule` only used `lat`; `Lat` was a duplicate alias with no independent role.

- Deleted `int n_mol;`.
- Reason: no reads or writes were found outside the declaration; it was dead state.

- Deleted `int mol_nr;`.
- Reason: no reads or writes were found outside the declaration; it was dead state.

- Deleted `string composition;`.
- Reason: the class never stored composition in this member; it always read it back through `GetValue("composition")`.

- Deleted `void PutTheta(Real);`.
- Reason: no call sites were found.

- Deleted `Real Theta(void);`.
- Reason: it had no implementation and no call sites.

## `src/molecule.cpp`

- Deleted the constructor assignment `Lat = Lat_;`.
- Reason: paired with deletion of the dead `Lat` member.

- Deleted the constructor assignment `lat = Lat;`.
- Reason: after removing the duplicate alias, the constructor can initialize `lat` directly from `Lat_`.

- Deleted the heap allocation `int* r = (int*) malloc(6*sizeof(int));` in the `restricted_range` path.
- Reason: replaced by `std::array<int, 6>` because the size is fixed and known at compile time.

- Deleted the matching `free(r);` in the `restricted_range` path.
- Reason: no longer needed after replacing the raw buffer with `std::array`.

- Deleted `Molecule::PutTheta`.
- Reason: it was unused, and its behavior was trivial enough that keeping a dead method only added API surface.

## `src/namics.cpp`

- Deleted the loop that reset `all_segments->prepared = false;`.
- Reason: `Segment::prepared` was removed because it had no readers and therefore no runtime effect.

## `src/reaction.h`

- Deleted `void DeAllocateMemory();`.
- Reason: the method body was empty and had no call sites except the destructor.

- Deleted `void AllocateMemory(int,int);`.
- Reason: the method body was empty and had no call sites.

- Deleted `void PrepareForCalculations();`.
- Reason: the method body was empty and had no call sites.

## `src/reaction.cpp`

- Deleted the destructor body that only called `DeAllocateMemory()`.
- Reason: after removing the empty helper, the destructor could be defaulted.

- Deleted `Reaction::DeAllocateMemory`.
- Reason: empty function, no behavior.

- Deleted `Reaction::AllocateMemory`.
- Reason: empty function, no behavior.

- Deleted `Reaction::PrepareForCalculations`.
- Reason: empty function, no behavior.

## `src/segment.h`

- Deleted `Lattice* Lat;`.
- Reason: `Segment` only used `lat`; `Lat` was a duplicate alias.

- Deleted `bool prepared;`.
- Reason: it was written but never read.

- Deleted `int* r;`.
- Reason: replaced by `std::array<int, 6> r` because the buffer size is fixed.

- Deleted `int used_in_mol_nr;`.
- Reason: it was initialized once and never used afterward.

## `src/segment.cpp`

- Deleted the constructor assignment `Lat = Lat_;`.
- Reason: paired with deletion of the dead `Lat` member.

- Deleted the constructor assignment `prepared = 0;`.
- Reason: paired with deletion of the dead `prepared` member.

- Deleted the constructor assignment `used_in_mol_nr = -2;`.
- Reason: paired with deletion of the dead `used_in_mol_nr` member.

- Deleted `free(r);` from `Segment::DeAllocateMemory`.
- Reason: `r` is no longer heap-allocated.

- Deleted the allocation `r = (int*) malloc(6*sizeof(int));` from `Segment::AllocateMemory`.
- Reason: replaced by `r.fill(0)` on a fixed-size `std::array`.

- Deleted the allocation `r = (int*) malloc(6*sizeof(int));` from `Segment::CheckInput`.
- Reason: replaced by `r.fill(0)` on the member `std::array`.

- Deleted the matching `free(r);` from `Segment::CheckInput`.
- Reason: no longer needed after removing heap allocation.

- Deleted the raw-pointer calls `lat->ReadRange(r, ...)` and `lat->CreateMASK(..., r, ...)`.
- Reason: these deletions are the pointer-form call sites that became unnecessary after switching to `r.data()`.

## `src/solve_scf.h`

- Deleted `Lattice* Lat;`.
- Reason: `Solve_scf` already had `lat`; the capitalized alias only duplicated that dependency.

- Deleted `string StoreFileGuess;`.
- Reason: the field was only assigned during parsing and never used afterward.

- Deleted `string ReadFileGuess;`.
- Reason: the field was only assigned during parsing and never used afterward.

- Deleted `void ComputePhis(bool);`.
- Reason: it had no implementation and no call sites.

- Deleted `bool PutU();`.
- Reason: it had no implementation and no call sites.

- Deleted `bool Put_U();`.
- Reason: it had no implementation and no call sites.

## `src/solve_scf.cpp`

- Deleted `Lat{Lat_}` from the constructor initializer list.
- Reason: paired with deletion of the duplicate `Lat` member.

- Deleted the constructor statement `lat = Lat;`.
- Reason: after removing the duplicate alias, `lat` is initialized directly from `Lat_`.

- Deleted parsing of `store_guess`.
- Reason: the parsed value was never used, and the key was not part of the accepted `KEYS` list anyway.

- Deleted parsing of `read_guess`.
- Reason: same as `store_guess`: parsed but never used, and not actually accepted by `CheckParameters`.

- Deleted the use of `Lat` when constructing `SCF_LBFGS`.
- Reason: after removing the duplicate member, the code uses `lat` directly.

## `src/state.h`

- Deleted `void DeAllocateMemory();`.
- Reason: the method body was empty and had no useful call sites.

- Deleted `void AllocateMemory(int,int);`.
- Reason: the method body was empty and had no call sites.

- Deleted `void PrepareForCalculations();`.
- Reason: the method body was empty and had no call sites.

## `src/state.cpp`

- Deleted the destructor body that only called `DeAllocateMemory()`.
- Reason: after removing the empty helper, the destructor could be defaulted.

- Deleted `State::DeAllocateMemory`.
- Reason: empty function, no behavior.

- Deleted `State::AllocateMemory`.
- Reason: empty function, no behavior.

- Deleted `State::PrepareForCalculations`.
- Reason: empty function, no behavior.

## `src/system.h`

- Deleted `Lattice* Lat;`.
- Reason: `System` only used `lat`; `Lat` was a duplicate alias.

- Deleted `bool prepared;`.
- Reason: it was only assigned and never read.

## `src/system.cpp`

- Deleted the constructor assignment `Lat = Lat_;`.
- Reason: paired with deletion of the dead `Lat` alias.

- Deleted the constructor assignment `lat = Lat;`.
- Reason: after removing the alias, `lat` is initialized directly from `Lat_`.

- Deleted the constructor assignment `prepared = false;`.
- Reason: paired with deletion of the dead `prepared` field.

- Deleted the assignment `prepared = true;` in `PrepareForCalculations`.
- Reason: paired with deletion of the dead `prepared` field.

- Deleted the heap allocation `int* bc = (int*) malloc(6*sizeof(int));`.
- Reason: replaced by `std::array<int, 6>` because the size is fixed.

- Deleted the matching `free(bc);`.
- Reason: no longer needed after replacing the raw buffer with `std::array`.

## `src/namics.h`

- No deletion happened in this file.
- The only change was adding `<array>` to support deletion of small fixed-size heap allocations elsewhere.

## Notes

- This log explains the first conservative shrink pass only.
- It does not cover future pruning work such as removing untested features or larger de-duplication across the output/property plumbing.

## Second Sweep

Verification after the second sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This second pass stayed focused on declaration-only methods and fields that had become obviously disconnected from the live code graph.

### `src/molecule.h`

- Deleted `Real n_range;`.
- Reason: it was assigned once in the `range_restricted` path and never read afterward.

- Deleted `Real GN1, GN2;`.
- Reason: no reads or writes were found outside the declaration; they were dead state.

- Deleted `int GetChainlength(void);`.
- Reason: no call sites were found.

### `src/molecule.cpp`

- Deleted the assignment `n_range = theta_range/chainlength;`.
- Reason: paired with deletion of the dead `n_range` member.

- Deleted `Molecule::GetChainlength`.
- Reason: it only returned `chainlength` and had no call sites.

### `src/segment.h`

- Deleted `void PutContraintBC();`.
- Reason: no call sites were found.

- Deleted `bool IsFree();`.
- Reason: no call sites were found.

- Deleted `bool IsFrozen();`.
- Reason: no call sites were found.

### `src/segment.cpp`

- Deleted `Segment::PutContraintBC`.
- Reason: it was uncalled boundary-writing scaffolding with no effect on the tested execution paths.

- Deleted `Segment::IsFree`.
- Reason: no call sites were found.

- Deleted `Segment::IsFrozen`.
- Reason: no call sites were found.

### `src/system.h`

- Deleted `string GetMonName(int);`.
- Reason: no call sites were found.

- Deleted `bool Put_U(Real*);`.
- Reason: no call sites were found.

- Deleted `Real GetSpontaneousCurvature();`.
- Reason: no call sites were found.

- Deleted `Real GetKBar();`.
- Reason: no call sites were found.

### `src/system.cpp`

- Deleted `System::GetMonName`.
- Reason: it only forwarded `Seg[index]->name` and had no callers.

- Deleted `System::Put_U`.
- Reason: it was an older inverse of `PutU` with no remaining call sites; removing it also eliminated duplicated field-adjustment logic.

- Deleted `System::GetSpontaneousCurvature`.
- Reason: it had no callers; the current output path computes `kJ0` and `kbar` directly where needed.

- Deleted `System::GetKBar`.
- Reason: it had no callers.

## Third Sweep

Verification after the third sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass stayed in the same lane as the second one: remove dead accessors, dead fallback code, and dead temporary state without changing live output or solver behavior.

### `src/segment.h`

- Deleted `Real* GetMASK();`.
- Reason: no call sites were found.

- Deleted `Real* GetPhi();`.
- Reason: no call sites were found.

### `src/segment.cpp`

- Deleted `Segment::GetMASK`.
- Reason: it was an unused accessor for the segment mask buffer.

- Deleted `Segment::GetPhi`.
- Reason: it was an unused accessor for the segment phi buffer.

### `src/solve_scf.h`

- Deleted the `rescue` enum.
- Reason: it only existed to support an unused DIIS rescue path.

- Deleted `bool attempt_DIIS_rescue();`.
- Reason: no call sites were found.

- Deleted `rescue rescue_status;`.
- Reason: paired with deletion of the unused rescue code.

### `src/solve_scf.cpp`

- Deleted the constructor assignment `rescue_status = NONE;`.
- Reason: paired with deletion of the dead rescue state.

- Deleted `Solve_scf::attempt_DIIS_rescue`.
- Reason: the method had no callers, so its DIIS rescue ladder was dead code.

### `src/system.h`

- Deleted `Real volume;`.
- Reason: it was only used as temporary scratch state inside `generate_mask()` and was never read as stored object state.

### `src/system.cpp`

- Deleted writes to the member `volume` inside `System::generate_mask()`.
- Reason: after removing the dead field, the code computes the same quantity in a local `accessible_volume` variable and stores it directly into `lat->Accesible_volume`.

## Fourth Sweep

Verification after the fourth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass moved one step beyond dead declarations and removed one isolated feature branch that the default regression suite does not exercise: `compute_Gibbs_excess`.

### `src/molecule.h`

- Deleted `Real R_Gibbs;`.
- Reason: paired with removal of the `compute_Gibbs_excess` feature path.

- Deleted `Real theta_Gibbs;`.
- Reason: paired with removal of the `compute_Gibbs_excess` feature path.

- Deleted `Real ComputeGibbs(Real);`.
- Reason: it was only used by the removed `compute_Gibbs_excess` branch.

### `src/molecule.cpp`

- Deleted output of `theta_Gibbs`.
- Reason: the value only existed for the removed `compute_Gibbs_excess` feature.

- Deleted conditional output of `R_Gibbs`.
- Reason: same as `theta_Gibbs`: it only existed for the removed feature.

- Deleted `Molecule::ComputeGibbs`.
- Reason: it was the molecule-side implementation of the removed feature.

### `src/system.cpp`

- Deleted `KEYS.push_back("compute_Gibbs_excess");`.
- Reason: the feature is no longer supported in the minimal default-tested build.

- Deleted the `GetValue("compute_Gibbs_excess")` branch in `System::CheckResults`.
- Reason: this was the only caller of `Molecule::ComputeGibbs`, and it belonged entirely to the removed feature path.

## Fifth Sweep

Verification after the fifth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed another untested output branch and another dead output path, while also shrinking some persistent object state into locals.

### `src/molecule.h`

- Deleted `Real Delta_MU;`.
- Reason: no assignments were found beyond constructor initialization, so it only fed a dead always-zero output path.

### `src/molecule.cpp`

- Deleted the constructor assignment `Delta_MU=0;`.
- Reason: paired with removal of the dead `Delta_MU` member.

- Deleted output of `DeltaMu`.
- Reason: it always reflected the dead `Delta_MU` field.

### `src/segment.h`

- Deleted `Real M1, M2, Fl;`.
- Reason: these values were only temporary scratch for `PushOutput()` and did not need to persist as object state.

### `src/segment.cpp`

- Deleted writes to the members `M1`, `M2`, and `Fl` in `PushOutput()`.
- Reason: they were replaced by local variables `first_moment`, `second_moment`, and `fluctuations`, producing the same output with less stored state.

### `src/system.cpp`

- Deleted `KEYS.push_back("compute_kJ0");`.
- Reason: `compute_kJ0` was identified in the pruning review as outside the default-tested surface, and its implementation lived only in output handling.

- Deleted the `GetValue("compute_kJ0")` output branch.
- Reason: this branch only conditionally produced `kJ0` / `kbar` output and was not required by the default test suite.

- Deleted `Sprod` accumulation and output.
- Reason: it depended only on `Mol[i]->Delta_MU`, which was never assigned anywhere beyond zero initialization, so `Sprod` was dead output.

## Sixth Sweep

Verification after the sixth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed another cluster of stale output-only state that no longer had any producer, plus one now-dead system field.

### `src/segment.h`

- Deleted `Real phi_LB_X;`.
- Reason: it was only zero-initialized and emitted as output; no writes beyond initialization remained.

- Deleted `Real phi_UB_X;`.
- Reason: same as `phi_LB_X`: no writes beyond initialization remained.

- Deleted `Real phi_LB_Y;`.
- Reason: same as `phi_LB_X`: no writes beyond initialization remained.

- Deleted `Real phi_UB_Y;`.
- Reason: same as `phi_LB_X`: no writes beyond initialization remained.

- Deleted `Real B;`.
- Reason: it was only initialized to `1` and emitted as output; no meaningful updates were present.

- Deleted `Real J;`.
- Reason: it was only initialized to `0` and emitted as output; no meaningful updates were present.

### `src/segment.cpp`

- Deleted constructor initialization of `phi_LB_X`, `phi_UB_X`, `phi_LB_Y`, and `phi_UB_Y`.
- Reason: paired with deletion of the dead fields.

- Deleted constructor initialization of `B` and `J`.
- Reason: paired with deletion of the dead fields.

- Deleted output of `B`, `J`, `phi_LB_x`, `phi_UB_x`, `phi_LB_y`, and `phi_UB_y`.
- Reason: these values no longer had any producer and only reflected stale constant state.

### `src/molecule.h`

- Deleted `Real J;`.
- Reason: it only accumulated the removed dead segment `J` values and was otherwise output-only state.

### `src/molecule.cpp`

- Deleted constructor initialization of `J`.
- Reason: paired with deletion of the dead `J` field.

- Deleted accumulation of molecule `J` from segment `J`.
- Reason: after removing the dead segment-side `J` field, this aggregation path had no live input.

- Deleted output of molecule `J`.
- Reason: it only reflected the removed dead aggregation path.

### `src/system.h`

- Deleted `Real pos_interface;`.
- Reason: it became dead object state after the earlier removal of the `compute_kJ0` branch.

### `src/system.cpp`

- Deleted constructor initialization of `pos_interface`.
- Reason: paired with deletion of the dead field.

## Seventh Sweep

Verification after the seventh sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `find_local_solution` / block-normalization branch and the scratch state that existed only to support it.

### `src/system.h`

- Deleted `bool local_solution;`.
- Reason: it only gated the removed `find_local_solution` feature path.

- Deleted `bool do_blocks;`.
- Reason: it only tracked activation of the removed block-normalization fallback.

- Deleted `int split;`.
- Reason: it was only parsed for the removed block-normalization path.

- Deleted `Real old_residual;`.
- Reason: it only tracked solver progress for the removed `find_local_solution` heuristic.

- Deleted `int progress;`.
- Reason: same as `old_residual`: it only existed for the removed heuristic.

### `src/system.cpp`

- Deleted `KEYS.push_back("find_local_solution");` and `KEYS.push_back("split");`.
- Reason: both parameters belonged only to the removed feature.

- Deleted constructor and allocation-time initialization of `local_solution`, `do_blocks`, `split`-related state, `old_residual`, and `progress`.
- Reason: paired with removal of the dead feature state.

- Deleted `find_local_solution` / `split` parsing and validation from `System::CheckInput`.
- Reason: the default suite does not exercise this feature, and the branch was explicitly called out as a pruning candidate in the review.

- Deleted the residual-progress tracking block at the start of `System::ComputePhis`.
- Reason: it only existed to trigger the removed fallback mode.

- Deleted the `do_blocks` normalization block in `System::ComputePhis`.
- Reason: it only invoked the removed molecule-side block helpers.

- Deleted the `prepare_for_blocks` activation block at the end of `System::ComputePhis`.
- Reason: it was the final hook that turned on the removed feature.

### `src/molecule.h`

- Deleted `vector<Real> block;`.
- Reason: it only stored per-block totals for the removed `find_local_solution` path.

- Deleted declarations of `SetThetaBlocks(int)` and `NormPerBlock(int)`.
- Reason: both helpers were only called from the removed system-side block fallback.

### `src/molecule.cpp`

- Deleted `Molecule::NormPerBlock`.
- Reason: it only renormalized restricted molecules for the removed block fallback.

- Deleted `Molecule::SetThetaBlocks`.
- Reason: it only precomputed per-block totals for the removed block fallback.

## Eighth Sweep

Verification after the eighth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `compute_width_interface` feature and the width-related output state that existed only for that optional calculation.

### `src/molecule.h`

- Deleted `Real phi1, phiM, width, Dphi, pos_interface, phi_av;`.
- Reason: these fields only held temporary state for `compute_width_interface` and its output keys.

- Deleted declaration of `ComputeWidth()`.
- Reason: the helper belonged entirely to the removed feature.

### `src/molecule.cpp`

- Deleted `KEYS.push_back("compute_width_interface");`.
- Reason: the minimal default-tested build no longer accepts that optional molecule flag.

- Deleted constructor initialization of `phi1`, `phiM`, `width`, `Dphi`, `pos_interface`, and `phi_av`.
- Reason: paired with removal of the dead feature state.

- Deleted the conditional `ComputeWidth()` call in `PushOutput()`.
- Reason: it was the only execution site for the removed feature.

- Deleted output of `width`, `phi1`, `phiM`, `Dphi`, `pos_interface`, and `phi_average`.
- Reason: those keys existed only to expose the removed optional calculation.

- Deleted `Molecule::ComputeWidth`.
- Reason: it was the implementation of the removed `compute_width_interface` feature.

## Ninth Sweep

Verification after the ninth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed untested ring-molecule behavior and the special propagation code that existed only to support closed-chain handling.

### `src/molecule.h`

- Deleted `bool ring;`.
- Reason: ring behavior is outside the default-tested surface and the flag only enabled the removed closed-chain path.

### `src/molecule.cpp`

- Deleted `KEYS.push_back("ring");`.
- Reason: the minimal default-tested build no longer accepts the optional ring flag.

- Deleted constructor initialization of `ring`.
- Reason: paired with removal of the dead feature flag.

- Deleted ring parsing and validation in `Molecule::CheckInput`.
- Reason: that whole block existed only to validate and activate the removed ring mode.

- Simplified `Molecule::fraction`.
- Reason: it no longer needs to skip the first backbone segment to compensate for ring closure.

### `src/mol_linear.cpp`

- Deleted the ring-specific branch in `mol_linear::ComputePhi`.
- Reason: closed-chain propagation for linear molecules was only reachable through the removed ring option.

- Deleted the ring-only scratch objects and pinned-distance bookkeeping from `mol_linear::ComputePhi`.
- Reason: those allocations and contour-distance calculations were only needed by the removed branch.

### `src/mol_branched.cpp`

- Deleted ring-specific branches from `BackwardBra2ndO`, `ForwardBra2ndO`, `BackwardBra`, and `ForwardBra`.
- Reason: those special cases only changed propagation at the backbone ends for the removed ring mode.

- Deleted the ring-specific top-level branch in `mol_branched::ComputePhi`.
- Reason: the normal open-chain branched propagation is now the only supported path.

- Deleted the ring-only scratch objects, contour-distance bookkeeping, and pinned-reachability scanning from `mol_branched::ComputePhi`.
- Reason: they were only needed to enumerate valid closure points for ring molecules.

- Deleted one now-unused local `N`.
- Reason: it became dead after the ring-specific backward-propagation block was removed.

## Tenth Sweep

Verification after the tenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `save_memory` mode and the compressed propagation-storage machinery that existed only to support it.

### `src/molecule.h`

- Deleted `bool save_memory;`.
- Reason: the flag only enabled the removed optional storage mode.

- Deleted `vector<int> memory;` and `vector<int> last_stored;`.
- Reason: both vectors only tracked compressed forward-storage positions for the removed mode.

- Deleted `Real* Gs;`.
- Reason: this scratch buffer only existed for compressed forward/backward reconstruction in the removed mode.

### `src/molecule.cpp`

- Deleted `KEYS.push_back("save_memory");`.
- Reason: the minimal default-tested build no longer accepts the optional `save_memory` flag.

- Deleted constructor initialization and parsing of `save_memory`.
- Reason: paired with removal of the feature flag.

- Deleted conditional `free(Gs)` in `DeAllocateMemory`.
- Reason: `Gs` no longer exists after removing the optional storage mode.

- Deleted compressed-storage setup in `AllocateMemory`.
- Reason: the `memory` / `last_stored` bookkeeping and reduced-`N` allocation path belonged only to the removed mode.

- Simplified `AllocateMemory` to always use full forward storage length.
- Reason: the default-tested build now supports only the normal full-storage propagation path.

- Deleted allocation of `Gs`.
- Reason: the scratch buffer was only needed for the removed compressed-storage propagation.

- Deleted the `save_memory` branches in both `propagate_forward` overloads.
- Reason: they were the implementations of compressed forward propagation and sparse snapshot storage.

- Deleted the `save_memory` branches in both `propagate_backward` overloads.
- Reason: they reconstructed intermediate forward states only for the removed compressed-storage mode.

### `src/mol_linear.cpp`

- Deleted the `save_memory` initialization branch in `mol_linear::ComputePhi`.
- Reason: the normal backward propagator already handles initialization for the remaining full-storage path.

### `src/mol_branched.cpp`

- Deleted `save_memory` endpoint selection in `BackwardBra2ndO` and `BackwardBra`.
- Reason: branch recursion now always reads the normal full-storage endpoints from `last_s`.

- Deleted `save_memory` copy-back branches in `ForwardBra2ndO` and `ForwardBra`.
- Reason: branch joins now always write directly into the normal full-storage `Gg_f` sequence.

- Deleted the top-level `save_memory` initialization branch in `mol_branched::ComputePhi`.
- Reason: the remaining full-storage backward propagation path does not need that extra setup.

## Eleventh Sweep

Verification after the eleventh sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `range_restricted` / `restricted_range` molecule option and the system-side normalization logic that existed only to support it.

### `src/molecule.h`

- Deleted `Real theta_range;`.
- Reason: it only stored the target amount for the removed `range_restricted` mode.

- Deleted `Real* R_mask;`.
- Reason: this mask pointer was only allocated and used by the removed `restricted_range` path.

### `src/molecule.cpp`

- Deleted `KEYS.push_back("restricted_range");`.
- Reason: the minimal default-tested build no longer accepts the optional `restricted_range` parameter.

- Deleted the warning branch for `restricted_range` without `freedom = range_restricted`.
- Reason: both the parameter and the freedom mode were removed together.

- Removed `range_restricted` from accepted molecule freedom lists.
- Reason: the mode is no longer supported in the minimal build.

- Updated multiple validation/error messages to refer only to `restricted`.
- Reason: those messages no longer need to mention the removed `range_restricted` option.

- Deleted the `freedom == "range_restricted"` parsing branch.
- Reason: it was the only code that read `restricted_range`, allocated `R_mask`, and stored `theta_range`.

### `src/system.cpp`

- Deleted the `Mol[i]->freedom == "range_restricted"` block in `System::ComputePhis`.
- Reason: this was the only system-side consumer of `R_mask` and `theta_range`, and it existed only to renormalize the removed molecule mode.

## Twelfth Sweep

Verification after the twelfth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `fill_range` molecule option and the system-side fill-mask plumbing that existed only to support it.

### `src/molecule.h`

- Deleted `vector<int> FillRangesList;`.
- Reason: it only stored the pinned segment matches used by the removed `fill_range` preparation path.

- Deleted `bool Filling;`.
- Reason: it only flagged molecules participating in the removed `fill_range` feature.

### `src/molecule.cpp`

- Deleted constructor initialization of `FillRangesList` and `Filling`.
- Reason: paired with removal of the feature state.

- Removed `fill_range` from the accepted pinned-molecule freedom list.
- Reason: the minimal default-tested build no longer supports that optional mode.

- Deleted the `freedom == "fill_range"` parser branch.
- Reason: it only rewrote the freedom to `restricted` and activated the removed fill-range preparation path.

### `src/system.h`

- Deleted `vector<int> FillList;`.
- Reason: it only tracked the segment list used by the removed fill-mask setup.

- Deleted `Real* FILL;`.
- Reason: this mask buffer only existed for the removed feature.

- Deleted `bool Filling;`.
- Reason: it only tracked whether the removed fill-range mode was active for any molecule.

### `src/system.cpp`

- Deleted `free(FILL);` from `DeAllocateMemory`.
- Reason: the fill mask buffer no longer exists.

- Deleted allocation and initialization of `FILL` in `AllocateMemory`.
- Reason: that buffer only served the removed feature.

- Deleted the `Filling` / `FillList` / `FILL` preparation block in `System::PrepareForCalculations`.
- Reason: this was the entire system-side implementation of `fill_range`, including theta reconstruction from pinned masks.

- Deleted the `FILL` masking applied to free segments in `System::PrepareForCalculations`.
- Reason: it only enforced the removed fill-range exclusion mask.

## Thirteenth Sweep

Verification after the thirteenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the untested `delta_*` / system constraint-field path while leaving the separate segment-point-constraint machinery intact.

### `src/system.h`

- Deleted `vector<int> DeltaMolList;`, `px`, `py`, and `pz`.
- Reason: those vectors only stored delta-constraint molecules and delta-range coordinates.

- Deleted `bool constraintfields;`.
- Reason: it only enabled the removed delta-constraint feature path.

- Deleted `Real* H_BETA;`, `Real* BETA;`, `Real* H_beta;`, and `Real* beta;`.
- Reason: these arrays only held the removed delta-constraint field and its transformed mask.

- Deleted `Real phi_ratio;`.
- Reason: it only parameterized the removed delta two-molecule constraint.

- Deleted `string delta_inputfile;`.
- Reason: it only named the file-based delta mask input for the removed feature.

### `src/system.cpp`

- Deleted `constraint`, `delta_range`, `delta_range_units`, `delta_inputfile`, `delta_molecules`, and `phi_ratio` from `System::KEYS`.
- Reason: those parameters belonged entirely to the removed delta-constraint path.

- Deleted constructor initialization of `constraintfields`.
- Reason: paired with removal of the dead feature flag.

- Deleted conditional allocation, initialization, aliasing, and freeing of `H_beta` / `H_BETA`.
- Reason: those buffers only served the removed delta field.

- Deleted the `exp(-beta)` preparation step in `PrepareForCalculations`.
- Reason: it transformed the removed delta field before SCF iterations.

- Deleted the `constraint` / `delta_*` parser branch in `System::CheckInput`.
- Reason: this was the entire input path for the removed feature, including delta range parsing, optional file loading, delta molecule lookup, and `phi_ratio` validation.

- Deleted output of `delta_range` and `phi_ratio`.
- Reason: those keys only exposed the removed optional feature state.

- Deleted loading of `BETA` from the solver vector in `ComputePhis(Real*, bool, Real)`.
- Reason: the solver no longer carries a delta-field block.

- Deleted the `constraintfields` residual contribution in `Classical_residual`.
- Reason: it was the only residual equation contributed by the removed delta feature.

- Simplified `System::ComputePhis` to call `Mol[i]->ComputePhi()` directly for every molecule.
- Reason: the removed feature was the only reason to route molecules through the `ComputePhi(BETA, id)` overload.

- Deleted the beta-dependent terms in `GetFreeEnergy` and `GetGrandPotential`.
- Reason: those thermodynamic corrections existed only for the removed delta constraint.

### `src/molecule.h`

- Deleted declaration of `ComputePhi(Real*, int)`.
- Reason: the overload existed only for the removed delta-field pathway.

### `src/molecule.cpp`

- Deleted `Molecule::ComputePhi(Real*, int)`.
- Reason: it was the molecule-side implementation that applied and removed the delta field around `ComputePhi()`.

### `src/solve_scf.cpp`

- Deleted the extra `iv += M` branch for `Sys->constraintfields`.
- Reason: the SCF vector no longer includes a delta-field block.

## Fourteenth Sweep

Verification after the fourteenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the now-dead lattice `FillMask` API and its `LGrad1` / `LGrad2` / `LGrad3` implementations after the delta-constraint feature was removed.

### `src/lattice.h`

- Deleted the pure virtual `FillMask(Real*, vector<int>, vector<int>, vector<int>, string)`.
- Reason: the delta-constraint sweep removed its last caller, so the virtual API became dead.

### `src/LGrad1.h`

- Deleted declaration of `LGrad1::FillMask`.
- Reason: paired with removal of the dead lattice API.

### `src/LGrad2.h`

- Deleted declaration of `LGrad2::FillMask`.
- Reason: paired with removal of the dead lattice API.

### `src/LGrad3.h`

- Deleted declaration of `LGrad3::FillMask`.
- Reason: paired with removal of the dead lattice API.

### `src/LGrad1.cpp`

- Deleted `LGrad1::FillMask`.
- Reason: it only built the removed delta mask for 1-gradient systems.

### `src/LGrad2.cpp`

- Deleted `LGrad2::FillMask`.
- Reason: it only built the removed delta mask for 2-gradient systems.

### `src/LGrad3.cpp`

- Deleted `LGrad3::FillMask`.
- Reason: it only built the removed delta mask for 3-gradient systems.

## Fifteenth Sweep

Verification after the fifteenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass removed the remaining segment-point-constraint feature path, including its segment-side parser state and the extra SCF residual slots that existed only to support it.

### `src/segment.h`

- Deleted `vector<int> constraint_z;`, `vector<Real> constraint_phi;`, and `vector<Real> constraint_beta;`.
- Reason: these vectors only stored the removed per-site segment constraint data.

- Deleted `bool constraints;`.
- Reason: it only flagged whether the removed segment constraint feature was active.

- Deleted declarations of `Get_g(int)` and `Put_beta(int, Real)`.
- Reason: both helpers existed only for the removed SCF constraint slots.

### `src/segment.cpp`

- Deleted `KEYS.push_back("phi");`.
- Reason: the minimal default-tested build no longer accepts the segment `phi` input used for the removed constraint feature.

- Deleted constructor initialization of `constraints`.
- Reason: paired with removal of the feature flag.

- Deleted the `GetValue("phi")` parser/validation block in `Segment::CheckInput`.
- Reason: this was the entire input path for segment point constraints, including bracket parsing and `(z, phi)` validation.

- Deleted the diagnostic block that printed parsed segment constraints.
- Reason: it only reported state for the removed feature.

- Deleted `Segment::Get_g`.
- Reason: it computed the removed extra residual equations for segment constraints.

- Deleted `Segment::Put_beta`.
- Reason: it applied the removed extra SCF variables back onto constrained segment sites.

### `src/system.h`

- Deleted `int extra_constraints;`.
- Reason: it only tracked the number of extra SCF variables contributed by the removed segment constraint feature.

### `src/system.cpp`

- Deleted constructor initialization and `generate_mask()` accumulation of `extra_constraints`.
- Reason: paired with removal of the feature counter.

- Deleted loading of extra constraint variables in `ComputePhis(Real*, bool, Real)`.
- Reason: those SCF vector slots no longer exist.

- Deleted the `extra_constraints` residual block in `Classical_residual`.
- Reason: those residual equations existed only for the removed feature.

### `src/solve_scf.cpp`

- Deleted the `iv += Seg[i]->constraint_z.size()` accumulation in `AllocateMemory`.
- Reason: the SCF vector no longer reserves space for segment-point-constraint variables.

## Sixteenth Sweep

Verification after the sixteenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass stayed in the minimal-build lane but shifted from feature pruning to de-duplication and small API cleanup.

### `src/molecule.h`

- Deleted `int GetPinnedSeg(void);`.
- Reason: after removing `fill_range`, no call sites remained; it had become dead query scaffolding.

### `src/molecule.cpp`

- Deleted `Molecule::GetPinnedSeg`.
- Reason: paired with removal of the dead declaration.

- Deleted duplicate restricted-input validation branches in `Molecule::CheckInput`.
- Reason: pinned and non-pinned molecules still had separate copies of the same `theta`/`n` validation logic after the earlier pruning passes. The function now parses `freedom` once and shares the remaining restricted-path validation.

- Deleted impossible pinned checks from the non-pinned `freedom` branch.
- Reason: once the parser branch is entered with `pinned == false`, checks rejecting pinned `solvent` or pinned `neutralizer` states are unreachable.

- Deleted `GetValue("composition")` re-reads in `Molecule::CheckInput`.
- Reason: the value is now read once into a local and reused, which shortens the function and makes the flow easier to follow.

- Simplified `Molecule::IsPinned`.
- Reason: it now returns immediately on the first pinned segment match instead of carrying a temporary flag to the end of the loop.

### `src/segment.h`

- Deleted `string GetFreedom();`.
- Reason: `Segment::freedom` is already a public field, so the wrapper had become an unnecessary duplicate API.

### `src/segment.cpp`

- Deleted `Segment::GetFreedom`.
- Reason: paired with removal of the redundant wrapper declaration.

- Deleted the remaining `GetFreedom()` call sites in molecule code.
- Reason: they were replaced by direct reads of the already-public `Segment::freedom` field.

## Seventeenth Sweep

Verification after the seventeenth sweep:
- `make -j8` passed
- `python3 tests/run_tests.py` passed

This pass continued the same cleanup style: remove tiny single-use helpers and trivial alias wrappers when direct code was shorter and clearer.

### `src/molecule.h`

- Deleted `int GetMonNr(string);`.
- Reason: it had a single call site in `Interpret`, so keeping a dedicated helper added more API surface than value.

### `src/molecule.cpp`

- Deleted `Molecule::GetMonNr`.
- Reason: paired with removal of the single-use declaration.

- Deleted the `GetMonNr()` call inside `Interpret`.
- Reason: the short segment-name lookup loop is now spelled out where it is used.

### `src/system.h`

- Deleted `int GetMonNr(string);`.
- Reason: it had a single remaining use in `CheckChi_values`.

- Deleted `bool IsCharged();`.
- Reason: it had a single remaining use in `CheckInput` and only wrapped a short loop over molecules.

### `src/system.cpp`

- Deleted `System::GetMonNr`.
- Reason: paired with removal of the single-use declaration.

- Deleted `System::IsCharged`.
- Reason: paired with removal of the single-use declaration.

- Deleted the `GetMonNr()` call in `CheckChi_values`.
- Reason: the direct scan over segment names is shorter than keeping a dedicated helper.

- Deleted the `IsCharged()` call in `CheckInput`.
- Reason: the charged-molecule scan now lives directly next to the neutralizer logic that depends on it.

### `src/segment.h`

- Deleted `bool IsPinned();`.
- Reason: no call sites remained after earlier pruning and wrapper cleanup.

### `src/segment.cpp`

- Deleted `Segment::IsPinned`.
- Reason: paired with removal of the dead declaration.
