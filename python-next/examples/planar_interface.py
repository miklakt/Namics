"""Periodic A/B interface with a small concentration of AB copolymer."""
import matplotlib.pyplot as plt
import numpy as np
from namics import Molecule, PlanarLattice, Segment, Species, System, solve
from _benchmark import coexistence, reference

a, b = Segment('A', chi={'B': 3}), Segment('B')
system = System(
    PlanarLattice(40, boundary='periodic'),
    species=(
        Species(Molecule('A_mol', segments=(a,)), bulk_fraction=.49),
        Species(Molecule('B_mol', segments=(b,)), bulk_fraction=.49),
        Species(Molecule('AB', segments=(a,) * 5 + (b,) * 5,
                         connections=tuple((i, i + 1) for i in range(9))), bulk_fraction=.02),
    ),
)
seed = np.full(40, .9)
seed[20:] = .1
result = solve(system, initial_fields={'A': -np.log(2 * seed), 'B': -np.log(2 * (1 - seed))},
               max_iterations=15000)
if not result.converged:
    raise RuntimeError(f'SCF did not converge: {result.residual}')
cpp = reference('ternary')
low, high, bulk_ab = coexistence(3, monomer_bulk=.49, copolymer_bulk=.02, block_length=5)
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(result.segments['A'], label='A, Python/C++')
axes[0].plot(result.segments['B'], label='B, Python/C++')
axes[0].plot(cpp['mon']['A']['phi'], 'o', markevery=3, fillstyle='none', label='A, Namics C++')
axes[0].plot(cpp['mon']['B']['phi'], 's', markevery=3, fillstyle='none', label='B, Namics C++')
axes[0].axhline(low, color='black', linestyle='--', label='Analytical bulk coexistence')
axes[0].axhline(high, color='black', linestyle='--')
axes[0].set_ylabel('Segment volume fraction')
axes[0].legend()
axes[1].plot(result.molecules['AB'], label='AB, Python/C++', color='tab:green')
axes[1].plot(cpp['mol']['AB']['phi'], 'o', markevery=3, fillstyle='none', label='AB, Namics C++')
axes[1].axhline(bulk_ab, color='black', linestyle='--', label='Analytical bulk AB')
axes[1].set_xlabel('Lattice site')
axes[1].set_ylabel('AB volume fraction')
axes[1].legend()
fig.suptitle('Copolymer at an A/B interface')
fig.tight_layout()
plt.show()
