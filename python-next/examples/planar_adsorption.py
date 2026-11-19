"""Homopolymer adsorption at a user-defined frozen surface."""
import matplotlib.pyplot as plt

from namics import FrozenSegment, Molecule, PlanarLattice, Segment, Species, System, solve

a = Segment('A', chi={'surface': -2})
w, surface = Segment('W'), Segment('surface')
n = 10000
polymer = Molecule('polymer', segments=(a,) * n, connections=tuple((i, i + 1) for i in range(n - 1)))
system = System(
    PlanarLattice(200),
    species=(
        Species(polymer, bulk_fraction=1e-4),
        Species(Molecule('water', segments=(w,)), bulk_fraction=1 - 1e-4),
    ),
    frozen=(FrozenSegment(surface, sites=(0,)),),
)
result = solve(system)
if not result.converged:
    raise RuntimeError(f'SCF did not converge: {result.residual}')
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(result.molecules['polymer'], label='Polymer')
axes[0].set_ylabel('Polymer volume fraction')
axes[0].legend()
axes[1].plot(result.molecules['water'], label='Water', color='tab:blue')
axes[1].set_xlabel('Lattice site (surface at 0)')
axes[1].set_ylabel('Water volume fraction')
axes[1].legend()
fig.suptitle('Homopolymer adsorption')
fig.tight_layout()
plt.show()
