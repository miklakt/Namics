# Removed Features List

Date: 2026-03-14

Purpose:
- list features intentionally removed during the pruning passes
- exclude pure dead-code cleanup, unused accessors, and no-op scaffolding
- record the user-visible or behavior-level features that existed before pruning

Interpretation note:
- "worked before these passes" means the feature had a real parser/runtime/output path before removal
- it does not mean the feature was covered by the default `python3 tests/run_tests.py` suite

## Output And Analysis Features

### `compute_Gibbs_excess`

- Removed in: Fourth Sweep
- Type: optional output/analysis feature
- User-facing surface:
  - system result key `compute_Gibbs_excess`
  - molecule outputs `theta_Gibbs`
  - conditional molecule output `R_Gibbs`
- What it did:
  - ran after convergence in `System::CheckResults`; it did not change the SCF solve itself
  - only worked for one-gradient systems
  - used the solvent profile to compute a Gibbs dividing-surface position `R_Gibbs` from the solvent excess amount
  - used that same `R_Gibbs` to compute `theta_Gibbs` for each non-solvent molecule, meaning the molecule excess relative to the Gibbs dividing surface rather than relative to the raw box average
  - warned if the total non-solvent Gibbs excess was not close to zero

### `compute_kJ0`

- Removed in: Fifth Sweep
- Type: optional output/analysis feature
- User-facing surface:
  - system result key `compute_kJ0`
  - conditional outputs `kJ0` and `kbar`
- What it did:
  - ran only in output generation; it did not affect the computed fields
  - only worked for planar one-gradient systems
  - computed `kJ0` as the first moment of the grand-potential-density profile around the interface position, or the box midpoint if no interface position had been set
  - computed `kbar` as the second moment of the same profile
  - was essentially a curvature post-processing branch on top of an already converged system

### `compute_width_interface`

- Removed in: Eighth Sweep
- Type: optional molecule analysis feature
- User-facing surface:
  - molecule key `compute_width_interface`
  - molecule outputs `width`, `phi1`, `phiM`, `Dphi`, `pos_interface`, `phi_average`
- What it did:
  - only worked for one-gradient systems and only when the key was explicitly set to `true`
  - looked at the molecule total profile near the two sides of the box to define `phi1`, `phiM`, and `Dphi`
  - scanned the profile for the steepest local step and treated that as the interface location
  - converted the inverse of that maximum normalized slope into an interface width
  - reported both the width and the midpoint concentration at the detected interface

## Solver And Search Features

### `Markov == 2` semiflexible-chain mode

- Removed in: current pruning pass
- Type: legacy molecule and lattice propagation mode
- User-facing surface:
  - molecule input keys `Markov` and `k_stiff`
  - Markov-dependent output fields `Markov`, `k_stiff`, and `P[...]`
- What it did:
  - turned the chain model from a single propagation state into a multi-state bond-direction model
  - used a stiffness weight table to bias transitions between contour directions
  - carried extra forward and backward propagator arrays for the different direction states
  - specialized the 1d, 2d, and 3d gradient kernels so they could combine those states differently depending on geometry and lattice type
  - had special handling for branched molecules, where the branch bookkeeping had to propagate the stiff-state information through both backbone and side branches

- What was deleted:
  - all input parsing for `Markov` and `k_stiff`
  - molecule-side stiffness bookkeeping, including the old `P` table
  - buffer allocation for the multi-state forward/backward propagators
  - the Markov-2 overloads of `propagate_forward`, `propagate_backward`, `ComputeGN`, `AddPhiS`, `Initiate`, and `Terminate`
  - lattice-side state buffers and coefficients such as `l1`, `l11`, `l_1`, `l_11`, `LABDA`, and `LABDA_1`
  - the `propagateF` and `propagateB` implementations in `LGrad1`, `LGrad2`, and `LGrad3`
  - the Markov-2 branches inside `ComputeLambdas` that filled the stiffness coefficient tables
  - all `Markov == 2` output writing and the old warning branches that tried to keep unsupported combinations running

- Why this matters for restoration:
  - restoring the mode means rebuilding the full state-space model across the solver, not just accepting the old input keys again
  - the old code was tightly coupled to geometry-specific and lattice-specific branch logic, so the reimplementation should be treated as a fresh feature with dedicated tests
  - any future restoration should first define the intended physics for each geometry and lattice combination, then reintroduce the supporting buffers and propagation kernels in a consistent way

### `find_local_solution`

- Removed in: Seventh Sweep
- Type: optional solver/search heuristic
- User-facing surface:
  - system keys `find_local_solution` and `split`
- What it did:
  - only worked for 3-gradient lattices with power-of-two box sizes
  - monitored SCF progress through the residual history
  - if progress stalled for long enough, it switched on a block-based fallback controlled by `split`
  - that fallback stored the amount of each restricted molecule in each spatial block and then renormalized each block separately instead of only normalizing the molecule globally
  - the goal was to keep a spatially structured local solution from collapsing into a different solution during SCF

### `save_memory`

- Removed in: Tenth Sweep
- Type: optional propagation/storage mode
- User-facing surface:
  - molecule key `save_memory`
- What it did:
  - replaced the normal full forward-propagator storage with a compressed checkpoint scheme
  - stored only selected forward contour states for each block instead of every state
  - reconstructed missing forward states during backward propagation when they were needed for density accumulation
  - traded simpler code and more memory for more control flow and less storage
  - aimed to keep the same physics while reducing the propagation memory footprint

## Molecule And Chain Features

### Ring Molecules

- Removed in: Ninth Sweep
- Type: optional molecule topology
- User-facing surface:
  - molecule key `ring`
- What it did:
  - changed molecule propagation from open-chain behavior to closed-chain behavior
  - enumerated candidate closure sites in the lattice, started the chain at a chosen site, propagated around the contour, and only counted configurations that closed back onto the same site
  - for pinned systems, filtered candidate closure sites using contour-distance reachability checks so that pinned monomers could still be reached
  - added special first/last-segment handling in both linear and branched propagation code so the ring closed consistently

### `range_restricted` With `restricted_range`

- Removed in: Eleventh Sweep
- Type: optional molecule freedom mode
- User-facing surface:
  - molecule freedom value `range_restricted`
  - molecule parameter `restricted_range`
- What it did:
  - was a variant of `freedom = restricted` that targeted a selected region rather than the whole box
  - built a mask from `restricted_range`
  - computed the molecule amount inside that mask only
  - renormalized the whole molecule so that the amount inside the masked region matched the requested `theta` or `n`
  - derived `phibulk`, `n`, and `theta` from that masked normalization, so the user constrained occupancy inside a region rather than total occupancy everywhere

### `fill_range`

- Removed in: Twelfth Sweep
- Type: optional pinned-molecule freedom mode
- User-facing surface:
  - molecule freedom value `fill_range`
- What it did:
  - was a special pinned-molecule mode that auto-computed the amount of a molecule needed to fill a pinned region
  - identified pinned segments whose masks matched the reference pinned segment of the filling molecule
  - built a combined fill mask from those pinned ranges
  - estimated how much of that volume was already occupied by other molecules and set the filling molecule `theta` to the remainder
  - then masked that region out of the normal free-segment preparation so the filled region behaved as reserved

## Constraint Features

### `delta_*` Constraint Field Path

- Removed in: Thirteenth Sweep
- Type: optional system-level constraint feature
- User-facing surface:
  - system keys `constraint`, `delta_range`, `delta_range_units`, `delta_inputfile`, `delta_molecules`, `phi_ratio`
  - output keys `delta_range` and `phi_ratio`
- What it did:
  - defined a constrained spatial region from `delta_range` coordinates or from a mask file
  - selected exactly two molecules through `delta_molecules`
  - added an extra field `BETA` on the constrained region as additional SCF unknowns
  - biased one selected molecule by multiplying its segment weights by `BETA` and biased the other by dividing by `BETA`
  - enforced a local composition relation between the two selected molecules through an extra residual block parameterized by `phi_ratio` or `critical_ratio`
  - added matching extra terms to the free-energy and grand-potential calculations

### Segment Point Constraints

- Removed in: Fifteenth Sweep
- Type: optional per-segment constraint feature
- User-facing surface:
  - segment key `phi`
- What it did:
  - only worked for one-gradient systems
  - accepted explicit local constraints of the form `(z, phi)` on a segment type
  - added one extra SCF unknown per constrained point
  - applied that extra unknown as an additive local field on the constrained lattice site
  - enforced the target local segment density through an extra residual equation at each constrained point
  - effectively let the user pin the value of a specific segment profile at specific positions

## Not Counted As Removed Features

The following were intentionally not counted here because they were implementation fallout or dead/stale code rather than distinct supported features:
- unused accessors, duplicate aliases, empty lifecycle hooks, and dead fields
- stale constant outputs like `J`, `B`, `phi_LB_*`, `phi_UB_*`, `DeltaMu`, and `Sprod`
- the lattice `FillMask` API, which was internal support code for the removed `delta_*` feature rather than a standalone feature

## Summary

Feature-level removals from these pruning passes were:
- `compute_Gibbs_excess`
- `compute_kJ0`
- `find_local_solution`
- `compute_width_interface`
- ring molecules
- `save_memory`
- `range_restricted` with `restricted_range`
- `fill_range`
- `delta_*` constraint fields
- segment point constraints via segment `phi`
