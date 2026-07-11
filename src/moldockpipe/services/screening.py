from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from ..models import ScreeningPolicy


@dataclass(frozen=True)
class ScreeningResult:
    canonical_smiles: str | None
    inchikey: str | None
    descriptors: dict[str, float | int]
    rules: dict[str, bool]
    decision: str
    reason: str


def screen_smiles(smiles: str, enabled_rules: dict[str, bool], policy: ScreeningPolicy) -> ScreeningResult:
    """Compute explicit physicochemical rules; no rule is a full ADMET prediction."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError as exc:
        raise RuntimeError("RDKit is required for physicochemical screening") from exc

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ScreeningResult(None, None, {}, {}, "fail", "RDKit could not parse or sanitize SMILES")
    descriptors: dict[str, float | int] = {
        "mw": round(Descriptors.MolWt(molecule), 2), "alogp": round(Crippen.MolLogP(molecule), 2),
        "tpsa": round(rdMolDescriptors.CalcTPSA(molecule), 2), "hbd": Lipinski.NumHDonors(molecule),
        "hba": Lipinski.NumHAcceptors(molecule), "rotatable_bonds": Lipinski.NumRotatableBonds(molecule),
    }
    rules = {
        "lipinski": descriptors["mw"] <= 500 and descriptors["alogp"] <= 5 and descriptors["hbd"] <= 5 and descriptors["hba"] <= 10,
        "veber": descriptors["tpsa"] <= 140 and descriptors["rotatable_bonds"] <= 10,
        "egan": descriptors["alogp"] <= 5.88 and descriptors["tpsa"] <= 131,
        "ghose": 160 <= descriptors["mw"] <= 480 and -0.4 <= descriptors["alogp"] <= 5.6,
    }
    selected = [name for name, active in enabled_rules.items() if active]
    failed = [name for name in selected if not rules[name]]
    if policy is ScreeningPolicy.ANNOTATE_ONLY:
        decision = "warning" if failed else "pass"
    elif policy is ScreeningPolicy.EXCLUDE_FAILING_ALL:
        decision = "fail" if selected and len(failed) == len(selected) else "pass"
    elif policy is ScreeningPolicy.EXCLUDE_FAILING_ANY:
        decision = "fail" if failed else "pass"
    else:
        decision = "manual_review" if failed else "pass"
    reason = "All selected rules passed" if not failed else f"Rule flags: {', '.join(failed)}"
    return ScreeningResult(Chem.MolToSmiles(molecule, isomericSmiles=True), Chem.MolToInchiKey(molecule), descriptors, rules, decision, reason)
