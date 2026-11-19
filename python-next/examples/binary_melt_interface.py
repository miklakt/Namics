"""Periodic interface between two monomer species A and B."""
import matplotlib.pyplot as plt
import numpy as np
from namics import Molecule, PlanarLattice, Segment, Species, System, solve
from _benchmark import coexistence, reference

a, b = Segment('A', chi={'B': 3}), Segment('B')
system = System(
    PlanarLattice(40, boundary='periodic'),
    species=(
        Species(Molecule('A_mol', segments=(a,)), bulk_fraction=.5),
        Species(Molecule('B_mol', segments=(b,)), bulk_fraction=.5),
    ),
)
seed = np.full(40, .9)
seed[20:] = .1
result = solve(system, initial_fields={'A': -np.log(2 * seed), 'B': -np.log(2 * (1 - seed))},
               max_iterations=15000)
if not result.converged:
    raise RuntimeError(f'SCF did not converge: {result.residual}')
cpp = reference('binary')
low, high, _ = coexistence(3, monomer_bulk=.5)
fig, ax = plt.subplots()
ax.plot(result.segments['A'], label='A, Python/C++')
ax.plot(result.segments['B'], label='B, Python/C++')
ax.plot(cpp['mon']['A']['phi'], 'o', markevery=3, fillstyle='none', label='A, Namics C++')
ax.plot(cpp['mon']['B']['phi'], 's', markevery=3, fillstyle='none', label='B, Namics C++')
ax.axhline(low, color='black', linestyle='--', label='Analytical bulk coexistence')
ax.axhline(high, color='black', linestyle='--')
ax.set(xlabel='Lattice site', ylabel='Segment volume fraction', title='Binary monomer interface')
ax.legend()
fig.tight_layout()
plt.show()
