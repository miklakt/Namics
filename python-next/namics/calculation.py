"""SCF calculation state and solver for the supported one-dimensional lattice."""
from __future__ import annotations

from typing import Mapping

import numpy as np

from ._kernel import load_kernel
from .model import Molecule, System


def _linear_order(molecule: Molecule) -> list[int]:
    if not molecule.is_linear:
        raise NotImplementedError("The current solver supports only linear molecules")
    neighbors = [[] for _ in molecule.segments]
    for a, b in molecule.connections:
        neighbors[a].append(b)
        neighbors[b].append(a)
    current = next(i for i, row in enumerate(neighbors) if len(row) <= 1)
    order, previous = [], -1
    while True:
        order.append(current)
        following = [i for i in neighbors[current] if i != previous]
        if not following:
            return order
        previous, current = current, following[0]


def _readonly(array):
    view = array.view()
    view.flags.writeable = False
    return view


class Calculation:
    """Own a single SCF calculation and its live working profiles.

    Profiles are exposed as read-only views. They change in place after step().
    At every public boundary, fields, densities, and residual describe the same
    evaluated state. Copy a profile explicitly to retain an earlier iteration.
    """
    def __init__(self, system: System, *, initial_fields: Mapping[str, np.ndarray] | None = None,
                 tolerance: float = 1e-7, max_step: float = 0.0001):
        if not isinstance(system, System):
            raise TypeError("Calculation expects a System")
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("tolerance must be positive and finite")
        if not np.isfinite(max_step) or max_step <= 0:
            raise ValueError("max_step must be positive and finite")
        self.system = system
        self.tolerance = tolerance
        self.max_step = max_step
        self.iteration = 0
        self.evaluations = 0
        self._segments = system.segment_types()
        if any(segment.states for segment in self._segments):
            raise NotImplementedError("Segment states are not supported by the current solver")
        self._indices = {segment.name: i for i, segment in enumerate(self._segments)}
        mobile_names = {segment.name for item in system.species for segment in item.molecule.segments}
        self._mobile = [i for i, segment in enumerate(self._segments) if segment.name in mobile_names]
        self._m = system.lattice.layers
        self._types = len(self._segments)
        self._chi = system.interactions()
        self._fixed = np.zeros((self._types, self._m))
        for item in system.frozen:
            self._fixed[self._indices[item.segment.name], list(item.sites)] = 1
        self._accessible = self._fixed.sum(axis=0) == 0
        self._count = int(self._accessible.sum())
        self._bulk = np.zeros(self._types)
        self._chains = []
        for item in system.species:
            mol = item.molecule
            order = _linear_order(mol)
            indices = np.array([self._indices[mol.segments[i].name] for i in order], dtype=np.int32)
            for index in indices:
                self._bulk[index] += item.bulk_fraction / len(indices)
            self._chains.append((item, indices, np.empty(self._m * len(indices))))
        self._kernel = load_kernel()
        self._scratch = np.empty(2 * self._m)
        self._weights = np.zeros((self._types, self._m))
        self._density = self._fixed.copy()
        self._molecule_density = {item.molecule.name: np.empty(self._m) for item in system.species}
        self._field_storage = np.zeros((len(self._mobile), self._m))
        self.segments = {segment.name: _readonly(self._density[i])
                         for i, segment in enumerate(self._segments)}
        self.molecules = {name: _readonly(values) for name, values in self._molecule_density.items()}
        self.fields = {self._segments[index].name: _readonly(self._field_storage[row])
                       for row, index in enumerate(self._mobile)}
        self._x = np.zeros(len(self._mobile) * self._count)
        if initial_fields is not None:
            unknown = set(initial_fields) - set(self.fields)
            if unknown:
                raise ValueError(f"Unknown initial field segment names: {sorted(unknown)}")
            for row, index in enumerate(self._mobile):
                name = self._segments[index].name
                if name in initial_fields:
                    values = np.asarray(initial_fields[name], dtype=np.float64)
                    if values.shape != (self._m,) or not np.isfinite(values).all():
                        raise ValueError(f"Initial field {name!r} must be a finite profile of length {self._m}")
                    self._x.reshape(len(self._mobile), self._count)[row] = values[self._accessible]
        self._history_x = []
        self._history_r = []
        self._residual_vector = self._evaluate()

    def _evaluate(self):
        self.evaluations += 1
        u = self._x.reshape(len(self._mobile), self._count)
        self._field_storage[:, self._accessible] = u
        self._weights.fill(0)
        with np.errstate(over="ignore", invalid="ignore"):
            self._weights[np.ix_(self._mobile, self._accessible)] = np.exp(-u)
            self._density[:] = self._fixed
            for item, indices, forward in self._chains:
                molecule_density = self._molecule_density[item.molecule.name]
                molecule_density.fill(0)
                if self._kernel(self._weights.ravel(), indices, self.system.lattice.neighbors.ravel(),
                                self.system.lattice.transitions.ravel(), self.system.lattice.neighbors.shape[1],
                                self._types, self._m, len(indices), item.bulk_fraction, forward,
                                self._density.ravel(), molecule_density, self._scratch):
                    raise RuntimeError("Propagation kernel rejected its inputs")
        total = self._density[:, self._accessible].sum(axis=0)
        if not np.all(np.isfinite(total)) or np.any(total <= 0):
            raise FloatingPointError("Non-finite or zero density during solve")
        sides = np.sum(self._density[:, self.system.lattice.neighbors] * self.system.lattice.transitions, axis=2)
        effective = u - self._chi[self._mobile] @ (sides[:, self._accessible] / total - self._bulk[:, None])
        residual = (effective - effective.mean(axis=0) + (1 / total - 1)).ravel()
        if not np.all(np.isfinite(residual)):
            raise FloatingPointError("Non-finite residual during solve")
        return residual

    @property
    def residual(self):
        return float(np.max(np.abs(self._residual_vector)))

    @property
    def converged(self):
        return bool(self.residual <= self.tolerance)

    def step(self, n: int = 1) -> Calculation:
        """Perform up to n field updates; stop early if already converged."""
        if type(n) is not int or n < 0:
            raise ValueError("n must be a non-negative integer")
        for _ in range(n):
            if self.converged:
                break
            # Preserve the last accepted state if the candidate evaluation fails.
            previous_x = self._x.copy()
            previous_history_x = self._history_x.copy()
            previous_history_r = self._history_r.copy()
            self._history_x.append(previous_x)
            self._history_r.append(self._residual_vector.copy())
            if len(self._history_x) > 8:
                self._history_x.pop(0)
                self._history_r.pop(0)
            count = len(self._history_x)
            matrix = np.empty((count + 1, count + 1))
            matrix[:count, :count] = np.array(self._history_r) @ np.array(self._history_r).T
            matrix[:count, count] = matrix[count, :count] = 1
            matrix[count, count] = 0
            rhs = np.zeros(count + 1)
            rhs[count] = 1
            try:
                coeff = np.linalg.solve(matrix, rhs)[:count]
                candidate = sum(c * (v - 0.1 * g) for c, v, g in zip(coeff, self._history_x, self._history_r))
            except np.linalg.LinAlgError:
                candidate = self._x - 0.1 * self._residual_vector
            delta = candidate - self._x
            max_delta = np.max(np.abs(delta))
            if max_delta > self.max_step:
                delta *= self.max_step / max_delta
            self._x += delta
            try:
                residual = self._evaluate()
            except Exception:
                self._x[:] = previous_x
                self._history_x = previous_history_x
                self._history_r = previous_history_r
                self._residual_vector = self._evaluate()
                raise
            self._residual_vector = residual
            self.iteration += 1
        return self

    def run(self, *, max_iterations: int = 10000) -> Calculation:
        """Run up to max_iterations additional updates, or until converged."""
        if type(max_iterations) is not int or max_iterations < 1:
            raise ValueError("max_iterations must be a positive integer")
        return self.step(max_iterations)


def solve(system: System, *, initial_fields: Mapping[str, np.ndarray] | None = None,
          tolerance: float = 1e-7, max_step: float = 0.0001, max_iterations: int = 10000) -> Calculation:
    """Create and run a calculation; the returned object owns its profiles."""
    return Calculation(system, initial_fields=initial_fields, tolerance=tolerance,
                       max_step=max_step).run(max_iterations=max_iterations)
