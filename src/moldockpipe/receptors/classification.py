from __future__ import annotations

from .models import AtomRecord, ComponentRole

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O"}
ION_NAMES = {"ZN", "MG", "CA", "NA", "CL", "K", "MN", "FE", "CU", "CO", "NI", "CD", "HG"}
COMMON_ADDITIVES = {
    "GOL", "EDO", "PEG", "PGE", "PG4", "MPD", "DMS", "ACT", "ACE", "FMT", "EOH",
    "IPA", "SO4", "PO4", "NO3", "TRS", "MES", "HEP", "BME", "DTT", "NAG",
}
COMMON_COFACTORS = {"HEM", "HEC", "FAD", "FMN", "NAD", "NAP", "SAM", "SAH", "PLP", "COA"}


def classify_component(name: str, atoms: tuple[AtomRecord, ...], *, polymer_like: bool = False,
                       covalently_linked: bool = False, near_protein: bool = True) -> tuple[str, ComponentRole, str]:
    name = name.upper()
    heavy = [atom for atom in atoms if atom.element.upper() not in {"H", "D"}]
    elements = {atom.element.upper() for atom in heavy}
    if name in WATER_NAMES:
        return "Water", ComponentRole.REMOVE, "recognized water"
    if len(heavy) == 1 and (name in ION_NAMES or next(iter(elements), "") in ION_NAMES):
        return "Metal / ion", ComponentRole.RETAINED_ION, "single-atom ion"
    if polymer_like:
        return "Modified polymer residue", ComponentRole.RECEPTOR_COMPONENT, "polymer residue classification"
    if name in COMMON_COFACTORS:
        return "Cofactor", ComponentRole.RETAINED_COFACTOR, "recognized cofactor"
    if name in COMMON_ADDITIVES:
        return "Crystallization additive", ComponentRole.REMOVE, "recognized solvent or buffer additive"
    if covalently_linked:
        return "Covalent component", ComponentRole.RECEPTOR_COMPONENT, "covalently connected to receptor"
    if "C" in elements and len(heavy) >= 5 and near_protein:
        return "Candidate ligand", ComponentRole.REFERENCE_LIGAND, "carbon-containing non-polymer near protein"
    return "Unclassified", ComponentRole.UNRESOLVED, "manual review required"
