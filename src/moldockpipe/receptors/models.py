from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any


@dataclass(frozen=True, order=True)
class ResidueKey:
    chain: str
    number: int
    insertion_code: str = ""
    name: str = ""

    def label(self) -> str:
        return f"{self.name} {self.chain}:{self.number}{self.insertion_code}".strip()


class ComponentRole(str, Enum):
    REFERENCE_LIGAND = "reference_ligand"
    RETAINED_COFACTOR = "retained_cofactor"
    RETAINED_ION = "retained_ion"
    RECEPTOR_COMPONENT = "receptor_component"
    REMOVE = "remove"
    UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class AtomRecord:
    name: str
    element: str
    x: float
    y: float
    z: float
    occupancy: float
    altloc: str = ""


@dataclass(frozen=True)
class Component:
    key: ResidueKey
    atoms: tuple[AtomRecord, ...]
    suggested_role: ComponentRole
    category: str
    reason: str
    covalently_linked: bool = False

    @property
    def heavy_atoms(self) -> tuple[AtomRecord, ...]:
        return tuple(atom for atom in self.atoms if atom.element.upper() not in {"H", "D"})


@dataclass(frozen=True)
class StructureAnalysis:
    source_path: Path
    models: tuple[str, ...]
    chains_by_model: dict[int, tuple[str, ...]]
    components: tuple[Component, ...]
    protein_residue_count: int
    water_count: int
    alternate_locations: tuple[str, ...]
    insertion_codes: tuple[str, ...]
    zero_occupancy_atoms: tuple[str, ...]
    connections: tuple[str, ...]
    protein_residue_issues: tuple["ProteinResidueIssue", ...] = ()
    protein_residues: tuple[ResidueKey, ...] = ()


@dataclass(frozen=True)
class ProteinResidueIssue:
    """A polymer residue which cannot be prepared without an explicit choice."""
    key: ResidueKey
    missing_atoms: tuple[str, ...]
    alternate_locations: tuple[str, ...]
    recommended_altloc: str = ""


@dataclass(frozen=True)
class ReceptorPreparationPlan:
    profile_id: str
    profile_name: str
    source_path: Path
    selected_model: int
    included_chains: tuple[str, ...]
    reference_ligand: ResidueKey | None
    removed_residues: tuple[ResidueKey, ...]
    retained_components: tuple[ResidueKey, ...]
    box_center: tuple[float, float, float]
    box_size: tuple[float, float, float]
    box_method: str
    box_parameters: dict[str, Any] = field(default_factory=dict)
    center_method: str = "manual"
    center_parameters: dict[str, Any] = field(default_factory=dict)
    altloc_choices: dict[ResidueKey, str] = field(default_factory=dict)
    excluded_receptor_residues: tuple[ResidueKey, ...] = ()
    preserve_hydrogens: bool = False
    chemistry_template_path: Path | None = None

    def as_record(self) -> dict[str, Any]:
        value = asdict(self)
        value["source_path"] = str(self.source_path)
        value["chemistry_template_path"] = str(self.chemistry_template_path) if self.chemistry_template_path else None
        value["altloc_choices"] = {key.label(): alt for key, alt in self.altloc_choices.items()}
        value["excluded_receptor_residues"] = [key.label() for key in self.excluded_receptor_residues]
        return value
