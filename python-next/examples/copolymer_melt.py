"""Periodic, incompatible symmetric AB diblock melt."""
import matplotlib.pyplot as plt
import numpy as np
from namics import Molecule, PlanarLattice, Segment, Species, System, solve
from _benchmark import reference, strong_segregation_profiles

a, b = Segment('A', chi={'B': 1.2}), Segment('B')
n = 10
segments = (a,) * n + (b,) * n
chain = Molecule('AB', segments=segments, connections=tuple((i, i + 1) for i in range(2*n - 1)))
system = System(PlanarLattice(20, boundary='periodic'), species=(Species(chain, bulk_fraction=1),))
x = np.arange(system.lattice.layers)
# A reproducible random initial field, with no prescribed spatial wavelength.
perturbation = np.random.default_rng(2).normal(0, .5, system.lattice.layers)
perturbation -= perturbation.mean()
result = solve(system, initial_fields={'A': -perturbation, 'B': perturbation}, max_step=.05)
if not result.converged:
    raise RuntimeError(f'SCF did not converge: {result.residual}')
cpp = reference('copolymer_melt')
print(f"Converged: {result.converged}; iterations: {result.iteration}; "
      f"maximum residual: {result.residual:.3g}")
print("Maximum A-profile difference from Namics C++:",
      f"{np.max(np.abs(result.segments['A'] - cpp['mon']['A']['phi'])):.3g}")

# Approximation for an isolated interface between long homopolymers. The
# numerical crossing locations set its position; it is not an exact diblock
# solution, and this short-chain, one-layer-wide interface is a severe test.
analytical_x, analytical_a, analytical_b = strong_segregation_profiles(result.segments['A'], 1.2)
fig, ax = plt.subplots()
ax.plot(x, result.segments['A'], label='A, Python/C++')
ax.plot(x, result.segments['B'], label='B, Python/C++')
ax.plot(x, cpp['mon']['A']['phi'], 'o', markevery=2, fillstyle='none', label='A, Namics C++')
ax.plot(x, cpp['mon']['B']['phi'], 's', markevery=2, fillstyle='none', label='B, Namics C++')
ax.plot(analytical_x, analytical_a, '--', color='tab:blue', label='A, strong-segregation approximation')
ax.plot(analytical_x, analytical_b, '--', color='tab:orange', label='B, strong-segregation approximation')
ax.set(xlabel='Lattice site', ylabel='Segment volume fraction',
       title='Incompatible A10B10 melt, χAB = 1.2')
ax.legend()
fig.tight_layout()
plt.show()
