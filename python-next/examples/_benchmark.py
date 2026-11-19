"""Reference data and bulk mean-field predictions for the example scripts."""
import json
import math
from pathlib import Path


def reference(case):
    path = Path(__file__).parent / 'reference' / f'{case}.output.json'
    with path.open(encoding='utf-8') as stream:
        return json.load(stream)['problems'][0]


def coexistence(chi, *, monomer_bulk, copolymer_bulk=0.0, block_length=1):
    """Bulk coexistence for symmetric A/B monomers and an optional A_n B_n chain.

    Write m = phi_A - phi_B and c = phi_AB. In a uniform bulk phase,
    incompressibility and ideal chain statistics give
        c = K * (((1-c)**2 - m**2)/4)**n,
        m = (1-c) * tanh(chi*m/2),
    where K = copolymer_bulk / monomer_bulk**(2*n). The equations are
    solved by bisection. They predict bulk plateaus, not an interface shape.
    """
    k = copolymer_bulk / monomer_bulk**(2 * block_length)

    def copolymer_fraction(m):
        if k == 0:
            return 0.0
        low, high = 0.0, 1.0 - abs(m)
        for _ in range(60):
            c = (low + high) / 2
            product = max(0.0, ((1-c)**2 - m*m) / 4)
            if c < k * product**block_length:
                low = c
            else:
                high = c
        return (low + high) / 2

    def equation(m):
        return (1 - copolymer_fraction(m)) * math.tanh(chi*m/2) - m

    low, high = 1e-8, 1.0
    if equation(low) <= 0:
        raise ValueError('Parameters have no nonzero coexistence root in this model')
    for _ in range(80):
        mid = (low + high) / 2
        if equation(mid) > 0:
            low = mid
        else:
            high = mid
    m = (low + high) / 2
    return (1 - m) / 2, (1 + m) / 2, copolymer_fraction(m)


def strong_segregation_profiles(profile, chi, *, samples_per_site=40):
    """Helfand–Tagami single-interface approximation, aligned to profile crossings.

    The width is 1/sqrt(6*chi), taking the statistical segment length as one
    lattice unit. This is a long-homopolymer, isolated-interface approximation,
    not an exact solution for a finite diblock melt. Only the locations and
    orientations of the interfaces are taken from the numerical profile.
    """
    import numpy as np

    profile = np.asarray(profile, dtype=float)
    crossings = []
    for i in range(len(profile)):
        left, right = profile[i], profile[(i + 1) % len(profile)]
        if (left < .5 <= right) or (right < .5 <= left):
            center = (i + (.5 - left) / (right - left)) % len(profile)
            crossings.append((center, np.sign(right - left)))
    if not crossings:
        raise ValueError('A segregated profile must cross 0.5')
    x = np.linspace(0, len(profile), len(profile) * samples_per_site + 1)
    centers = np.array([center for center, _ in crossings])
    orientations = np.array([orientation for _, orientation in crossings])
    distances = (x[:, None] - centers + len(profile) / 2) % len(profile) - len(profile) / 2
    nearest = np.argmin(np.abs(distances), axis=1)
    distance = distances[np.arange(len(x)), nearest]
    width = 1 / np.sqrt(6 * chi)
    a = .5 * (1 + np.tanh(orientations[nearest] * distance / width))
    return x, a, 1 - a
