from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class GeneratedState:
    state_smiles: str
    state_isomeric_smiles: str
    state_inchikey: str | None
    formal_charge: int
    structure_hash: str
    tautomer_index: int
    protomer_index: int
    molecule: object


class MolScrubService:
    def generate(self, smiles: str, *, ph: float = 7.4, enumerate_states: bool = True, max_states: int = 32) -> tuple[list[GeneratedState], bool]:
        """Return deduplicated chemical states. SDF writing belongs to the caller."""
        try:
            from rdkit import Chem
            from molscrub import Scrub
        except ImportError as exc:
            raise RuntimeError("MolScrub and RDKit are required for ligand preparation") from exc
        source = Chem.MolFromSmiles(smiles)
        if source is None:
            raise ValueError("Invalid or unsanitizable SMILES")
        if len(Chem.GetMolFrags(source)) > 1:
            raise ValueError("Disconnected structure requires an explicit fragment policy")
        scrubber = Scrub(ph_low=ph, ph_high=ph)
        molecules: Iterable[object] = scrubber(source) if enumerate_states else [source]
        states: list[GeneratedState] = []
        seen: set[str] = set()
        for molecule in molecules:
            isomeric = Chem.MolToSmiles(molecule, isomericSmiles=True)
            charge = Chem.GetFormalCharge(molecule)
            graph_hash = hashlib.sha256(f"{isomeric}|{charge}".encode()).hexdigest()
            if graph_hash in seen:
                continue
            seen.add(graph_hash)
            states.append(GeneratedState(Chem.MolToSmiles(molecule), isomeric, Chem.MolToInchiKey(molecule), charge, graph_hash, len(states), 0, molecule))
            if len(states) >= max_states:
                return states, True
        return states, False
