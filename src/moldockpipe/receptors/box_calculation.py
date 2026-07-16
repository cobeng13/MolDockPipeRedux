from __future__ import annotations

import math
from collections.abc import Iterable

from .models import AtomRecord


def _heavy(atoms: Iterable[AtomRecord]) -> tuple[AtomRecord, ...]:
    result = tuple(atom for atom in atoms if atom.element.upper() not in {"H", "D"})
    if not result:
        raise ValueError("At least one heavy atom is required")
    return result


def center_from_atoms(atoms: Iterable[AtomRecord]) -> tuple[float, float, float]:
    values = _heavy(atoms)
    count = len(values)
    return tuple(sum(getattr(atom, axis) for atom in values) / count for axis in ("x", "y", "z"))  # type: ignore[return-value]


def envelope_box(atoms: Iterable[AtomRecord], padding: float = 8.0) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    if padding < 0 or not math.isfinite(padding):
        raise ValueError("Padding must be finite and non-negative")
    values = _heavy(atoms)
    center = center_from_atoms(values)
    size = tuple(max(getattr(a, axis) for a in values) - min(getattr(a, axis) for a in values) + 2 * padding
                 for axis in ("x", "y", "z"))
    return center, size  # type: ignore[return-value]


def radius_of_gyration_box(atoms: Iterable[AtomRecord], multiplier: float = 2.857) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    values = _heavy(atoms)
    center = center_from_atoms(values)
    rg = math.sqrt(sum((a.x-center[0])**2 + (a.y-center[1])**2 + (a.z-center[2])**2 for a in values) / len(values))
    edge = multiplier * rg
    if edge <= 0:
        raise ValueError("Radius of gyration did not produce a positive box")
    return center, (edge, edge, edge)
