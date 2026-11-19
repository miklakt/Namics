"""Linear AB copolymer near a user-defined frozen surface."""
import matplotlib.pyplot as plt

from namics import FrozenSegment, Molecule, PlanarLattice, Segment, Species, System, solve

a, b = Segment('A', chi={'surface': -.6}), Segment('B')
water, surface = Segment('W'), Segment('surface')
polymer = Molecule('AB', segments=(a, b), connections=((0, 1),))
system = System(
    PlanarLattice(16),
    species=(
        Species(polymer, bulk_fraction=.1),
        Species(Molecule('water', segments=(water,)), bulk_fraction=.9),
    ),
    frozen=(FrozenSegment(surface, sites=(0,)),),
)
result = solve(system)
if not result.converged:
    raise RuntimeError(f'SCF did not converge: {result.residual}')
plt.plot(result.molecules['AB'], label='AB molecule')
plt.plot(result.segments['A'], label='A segments')
plt.plot(result.segments['B'], label='B segments')
plt.xlabel('Lattice site (surface at 0)')
plt.ylabel('Volume fraction')
plt.title('Copolymer adsorption')
plt.legend()
plt.tight_layout()
plt.show()
