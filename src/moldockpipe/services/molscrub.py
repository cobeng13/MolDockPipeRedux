from __future__ import annotations

import hashlib
import base64
import json
import os
import subprocess
import sys
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
    def generate(self, smiles: str, *, ph: float = 7.4, enumerate_states: bool = True, max_states: int = 32,
                 timeout_seconds: float = 15) -> tuple[list[GeneratedState], bool]:
        """Bound slow chemistry in a process that can be stopped on timeout."""
        from rdkit import Chem

        # Reject invalid inputs before paying the cost of starting a worker.
        self._source_molecule(smiles)
        if not enumerate_states:
            return self._generate(smiles, ph=ph, enumerate_states=False, max_states=max_states)
        if timeout_seconds <= 0:
            raise ValueError("MolScrub time limit must be greater than zero")
        request = json.dumps({"smiles": smiles, "ph": ph, "max_states": max_states})
        try:
            result = subprocess.run(
                [sys.executable, "-m", "moldockpipe.services.molscrub_worker"],
                input=request, capture_output=True, text=True, encoding="utf-8",
                timeout=timeout_seconds,
                creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
            )
        except subprocess.TimeoutExpired as exc:
            raise TimeoutError(f"MolScrub exceeded {timeout_seconds:g}s; molecule skipped") from exc
        if result.returncode:
            raise RuntimeError(result.stderr.strip() or f"MolScrub worker exited with code {result.returncode}")
        response = json.loads(result.stdout)
        if "error" in response:
            raise RuntimeError(response["error"] + ("\n" + result.stderr.strip() if result.stderr.strip() else ""))
        states = []
        for record in response["states"]:
            record["molecule"] = Chem.Mol(base64.b64decode(record["molecule"]))
            states.append(GeneratedState(**record))
        return states, response["truncated"]

    @staticmethod
    def _source_molecule(smiles: str):
        from rdkit import Chem

        source = Chem.MolFromSmiles(smiles)
        if source is None:
            raise ValueError("Invalid or unsanitizable SMILES")
        if len(Chem.GetMolFrags(source)) > 1:
            raise ValueError("Disconnected structure requires an explicit fragment policy")
        return source

    def _generate(self, smiles: str, *, ph: float = 7.4, enumerate_states: bool = True, max_states: int = 32) -> tuple[list[GeneratedState], bool]:
        """Return deduplicated chemical states. SDF writing belongs to the caller."""
        try:
            from rdkit import Chem
            from molscrub import Scrub
        except ImportError as exc:
            raise RuntimeError(
                f"MolScrub/RDKit could not load in {sys.executable}: {type(exc).__name__}: {exc}. "
                "Open Tools > Troubleshooter for the dependency report."
            ) from exc
        source = self._source_molecule(smiles)
        scrubber = Scrub(ph_low=ph, ph_high=ph, num_etkdg_attempts=1, use_random_coords=False)
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
