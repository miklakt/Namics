"""Immutable physical model definitions and their user-facing views."""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np


@dataclass(frozen=True)
class State:
    """A named internal state of a segment type (thermodynamics not yet modeled)."""
    name: str

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("A state must have a non-empty name")


@dataclass(frozen=True)
class Segment:
    """A segment type and its symmetric chi interactions with named types."""
    name: str
    chi: Mapping[str, float] = field(default_factory=dict, kw_only=True)
    states: tuple[State, ...] = field(default=(), kw_only=True)

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("A segment must have a non-empty name")
        states = tuple(self.states)
        if not all(isinstance(state, State) for state in states):
            raise ValueError("states must contain State objects")
        if len({state.name for state in states}) != len(states):
            raise ValueError("State names must be unique within a segment")
        object.__setattr__(self, "states", states)
        values = dict(self.chi)
        for other, value in values.items():
            if not isinstance(other, str) or not other:
                raise ValueError("Chi targets must be non-empty segment names")
            if not isinstance(value, (int, float)) or not np.isfinite(value):
                raise ValueError("Chi values must be finite numbers")
            if other == self.name and value != 0:
                raise ValueError("A segment's self-interaction chi must be zero")
        object.__setattr__(self, "chi", MappingProxyType(values))


@dataclass(frozen=True)
class Molecule:
    """Segment occurrences and undirected connections indexed from zero.

    Repeated entries in segments are distinct occurrences of the same type.
    Connected rings and branches may be represented; the current solver only
    supports linear molecules.
    """
    name: str
    segments: tuple[Segment, ...] = field(kw_only=True)
    connections: tuple[tuple[int, int], ...] = field(default=(), kw_only=True)

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("A molecule must have a non-empty name")
        segments = tuple(self.segments)
        if not segments or not all(isinstance(segment, Segment) for segment in segments):
            raise ValueError("A molecule must contain Segment occurrences")
        types = {}
        for segment in segments:
            if segment.name in types and types[segment.name] != segment:
                raise ValueError("Segment occurrences with the same name must have identical definitions")
            types[segment.name] = segment
        neighbors = [set() for _ in segments]
        connections = []
        for edge in self.connections:
            if not isinstance(edge, (list, tuple)) or len(edge) != 2:
                raise ValueError("Each connection must contain two segment indices")
            a, b = edge
            if type(a) is not int or type(b) is not int or not (0 <= a < len(segments) and 0 <= b < len(segments)) or a == b:
                raise ValueError("Connection indices must identify two different segment occurrences")
            if b in neighbors[a]:
                raise ValueError("Duplicate connection")
            neighbors[a].add(b)
            neighbors[b].add(a)
            connections.append((a, b))
        visited, pending = set(), [0]
        while pending:
            index = pending.pop()
            if index not in visited:
                visited.add(index)
                pending.extend(neighbors[index] - visited)
        if len(visited) != len(segments):
            raise ValueError("Molecule connections must form a connected graph")
        object.__setattr__(self, "segments", segments)
        object.__setattr__(self, "connections", tuple(connections))

    @property
    def is_linear(self) -> bool:
        if len(self.connections) != len(self.segments) - 1:
            return False
        degrees = [0] * len(self.segments)
        for a, b in self.connections:
            degrees[a] += 1
            degrees[b] += 1
        return max(degrees) <= 2


@dataclass(frozen=True)
class PlanarLattice:
    """Regular 1D simple-cubic lattice with mirror or periodic boundaries."""
    layers: int
    boundary: str = field(default="mirror", kw_only=True)
    neighbors: np.ndarray = field(init=False, repr=False, compare=False)
    transitions: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        if type(self.layers) is not int or self.layers < 2:
            raise ValueError("A planar lattice needs at least two layers")
        if self.boundary not in ("mirror", "periodic"):
            raise ValueError("boundary must be mirror or periodic")
        i = np.arange(self.layers, dtype=np.int32)
        left = (i - 1) % self.layers if self.boundary == "periodic" else np.maximum(i - 1, 0)
        right = (i + 1) % self.layers if self.boundary == "periodic" else np.minimum(i + 1, self.layers - 1)
        neighbors = np.ascontiguousarray(np.column_stack((left, i, right)), dtype=np.int32)
        transitions = np.tile(np.array([1/6, 4/6, 1/6]), (self.layers, 1))
        neighbors.flags.writeable = transitions.flags.writeable = False
        object.__setattr__(self, "neighbors", neighbors)
        object.__setattr__(self, "transitions", transitions)


@dataclass(frozen=True)
class Species:
    """A molecule at a prescribed bulk volume fraction."""
    molecule: Molecule
    bulk_fraction: float = field(kw_only=True)

    def __post_init__(self):
        if not isinstance(self.molecule, Molecule):
            raise ValueError("species needs a Molecule")
        if not np.isfinite(self.bulk_fraction) or not 0 <= self.bulk_fraction <= 1:
            raise ValueError("bulk_fraction must be between zero and one")


@dataclass(frozen=True)
class FrozenSegment:
    """A segment type occupying entire, zero-based lattice sites."""
    segment: Segment
    sites: tuple[int, ...] = field(kw_only=True)

    def __post_init__(self):
        if not isinstance(self.segment, Segment):
            raise ValueError("frozen segment must be a Segment")
        sites = tuple(self.sites)
        if any(type(i) is not int or i < 0 for i in sites) or len(set(sites)) != len(sites):
            raise ValueError("Frozen sites must be distinct non-negative integers")
        object.__setattr__(self, "sites", sites)


@dataclass(frozen=True)
class System:
    lattice: PlanarLattice
    species: tuple[Species, ...] = field(kw_only=True)
    frozen: tuple[FrozenSegment, ...] = field(default=(), kw_only=True)

    def __post_init__(self):
        if not isinstance(self.lattice, PlanarLattice):
            raise ValueError("lattice must be a PlanarLattice")
        species, frozen = tuple(self.species), tuple(self.frozen)
        if not species or not all(isinstance(item, Species) for item in species):
            raise ValueError("system requires at least one Species")
        if not all(isinstance(item, FrozenSegment) for item in frozen):
            raise ValueError("frozen must contain FrozenSegment objects")
        if len({item.molecule.name for item in species}) != len(species):
            raise ValueError("Molecule names must be distinct")
        if not np.isclose(sum(item.bulk_fraction for item in species), 1, atol=1e-12, rtol=0):
            raise ValueError("Species bulk fractions must sum to one")
        sites = [site for item in frozen for site in item.sites]
        if any(site >= self.lattice.layers for site in sites) or len(set(sites)) != len(sites):
            raise ValueError("Frozen sites must be within the lattice and cannot overlap")
        if len(sites) == self.lattice.layers:
            raise ValueError("At least one lattice site must be accessible")
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "frozen", frozen)
        self.interactions()

    def segment_types(self) -> tuple[Segment, ...]:
        types = {}
        for segment in (segment for item in self.species for segment in item.molecule.segments):
            if segment.name in types and types[segment.name] != segment:
                raise ValueError("Conflicting segment definitions")
            types[segment.name] = segment
        for item in self.frozen:
            segment = item.segment
            if segment.name in types and types[segment.name] != segment:
                raise ValueError("Conflicting segment definitions")
            types[segment.name] = segment
        return tuple(types.values())

    def segment(self, name: str) -> Segment:
        for segment in self.segment_types():
            if segment.name == name:
                return segment
        raise KeyError(name)

    def molecule(self, name: str) -> Molecule:
        for item in self.species:
            if item.molecule.name == name:
                return item.molecule
        raise KeyError(name)

    def interactions(self) -> np.ndarray:
        segments = self.segment_types()
        indices = {segment.name: i for i, segment in enumerate(segments)}
        matrix = np.zeros((len(segments), len(segments)))
        for i, segment in enumerate(segments):
            for other, value in segment.chi.items():
                if other not in indices:
                    raise ValueError(f"Unknown chi target {other!r}")
                j = indices[other]
                reverse = segments[j].chi.get(segment.name)
                if reverse is not None and value != reverse:
                    raise ValueError(f"Conflicting chi values for {segment.name!r} and {other!r}")
                matrix[i, j] = matrix[j, i] = value
        return matrix



