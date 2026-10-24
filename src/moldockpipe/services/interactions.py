"""Read-only, geometric protein–ligand contact annotations (not energy scores).

Ligand chemistry comes from RDKit features; protein chemistry uses conservative
standard-residue templates. No hydrogens or coordinates are generated here.
See README for criteria and limitations; this is not a PLIP implementation.
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable
from functools import lru_cache
import math
from pathlib import Path


TYPES = {
    "hbond": ("H-bonds", "#0072B2"),
    "pi_stack": ("π–π", "#7B3294"),
    "cation_pi": ("Cation–π", "#D55E00"),
    "salt_bridge": ("Salt bridges", "#009E73"),
    "hydrophobic": ("Hydrophobic", "#777777"),
}
NOTES = (
    "Geometric predictions, not validated bonds or interaction energies. "
    "H-bond candidates use heavy atoms when donor hydrogens are absent. "
    "Standard amino acids only; waters, metals and cofactors are not classified. "
    "Charges assume ASP/GLU negative and LYS/ARG positive; histidine roles require explicit hydrogens."
)
STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
DONORS = {
    "ARG": {"NE", "NH1", "NH2"}, "ASN": {"ND2"}, "GLN": {"NE2"},
    "LYS": {"NZ"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"},
    "TRP": {"NE1"}, "CYS": {"SG"},
}
ACCEPTORS = {
    "ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"},
    "GLN": {"OE1"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}, "MET": {"SD"},
}
RINGS = {
    "PHE": [("CG", "CD1", "CE1", "CZ", "CE2", "CD2")],
    "TYR": [("CG", "CD1", "CE1", "CZ", "CE2", "CD2")],
    "HIS": [("CG", "ND1", "CE1", "NE2", "CD2")],
    "TRP": [("CG", "CD1", "NE1", "CE2", "CD2"), ("CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3")],
}


@lru_cache(maxsize=1)
def _factory():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def analyze_interactions(pdb: str, sdf: str, *, cancelled: Callable[[], bool] | None = None) -> tuple[list[dict], str]:
    """Return JSON-safe contacts with coordinates in the original reference frame."""
    import numpy as np
    from scipy.spatial import cKDTree
    from rdkit import Chem

    def check():
        if cancelled and cancelled():
            raise InterruptedError("Interaction analysis cancelled")

    check()
    mol = Chem.MolFromMolBlock(sdf.split("$$$$", 1)[0], removeHs=False)
    if mol is None or not mol.GetNumConformers():
        raise ValueError("Ligand chemistry could not be read for interaction analysis.")
    features = _factory().GetFeaturesForMol(mol)
    donors = {i for f in features if f.GetFamily() == "Donor" for i in f.GetAtomIds()}
    # RDKit's general feature set also accepts halogens; omit them for H-bonds.
    acceptors = {i for f in features if f.GetFamily() == "Acceptor" for i in f.GetAtomIds()
                 if mol.GetAtomWithIdx(i).GetSymbol() in {"N", "O", "S"}}
    ligand = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        ligand.append(dict(xyz=np.array(mol.GetConformer().GetAtomPosition(i)),
                           element=atom.GetSymbol(), name=f"{atom.GetSymbol()}{i + 1}",
                           donor=i in donors, acceptor=i in acceptors, serial=i + 1))
    for atom in mol.GetAtoms():
        ligand[atom.GetIdx()]["neighbors"] = [ligand[n.GetIdx()] for n in atom.GetNeighbors()]

    residues = defaultdict(dict)
    saw_model = False
    for line in pdb.splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or line[16:17] not in {" ", "A"}:
            continue
        resn, name = line[17:20].strip(), line[12:16].strip()
        if resn not in STANDARD:
            continue
        key = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), resn)
        if name in residues[key]:
            continue
        element = line[76:78].strip().title() or name.lstrip("0123456789")[0]
        if element == "D":
            element = "H"
        xyz = np.array([float(line[p:p + 8]) for p in (30, 38, 46)])
        if not np.isfinite(xyz).all():
            raise ValueError("Non-finite receptor coordinates in interaction analysis.")
        residues[key][name] = dict(xyz=xyz, element=element, name=name, serial=int(line[6:11]),
                                  residue=key, donor=False, acceptor=False)
    ligand_xyz = np.array([a["xyz"] for a in ligand])
    if not np.isfinite(ligand_xyz).all():
        raise ValueError("Non-finite ligand coordinates in interaction analysis.")
    tree = cKDTree(ligand_xyz)
    protein, protein_rings, protein_charges = [], [], []

    def ring(atoms):
        coords = np.array([a["xyz"] for a in atoms])
        center = coords.mean(axis=0)
        # Ordered ring vertices give an area-weighted normal without a LAPACK
        # dependency (important for the Conda desktop runtime on Windows).
        normal = np.cross(coords - center, np.roll(coords - center, -1, axis=0)).sum(axis=0)
        magnitude = np.linalg.norm(normal)
        if magnitude < 0.1:
            return None
        normal /= magnitude
        if max(abs((coords - center) @ normal)) > 0.15:
            return None
        return dict(xyz=center, normal=normal, atoms=atoms)

    for key, names in residues.items():
        check()
        # Restrict expensive chemistry/geometry to residues near this pose.
        if not any(tree.query(a["xyz"])[0] <= 8 for a in names.values()):
            continue
        resn = key[3]
        atoms = list(names.values())
        for atom in atoms:
            atom["neighbors"] = [a for a in atoms if a is not atom and
                                  0.4 < np.linalg.norm(a["xyz"] - atom["xyz"]) <=
                                  (1.35 if "H" in {a["element"], atom["element"]} else 1.95)]
            name = atom["name"]
            has_h = any(a["element"] == "H" for a in atom["neighbors"])
            atom["donor"] = (name == "N" and resn != "PRO") or name in DONORS.get(resn, set())
            atom["acceptor"] = name in {"O", "OXT"} or name in ACCEPTORS.get(resn, set())
            if resn == "HIS" and name in {"ND1", "NE2"}:
                atom["donor"] = has_h
                # Only assign the other N as acceptor if this tautomer is known.
                other = names.get("NE2" if name == "ND1" else "ND1")
                atom["acceptor"] = not has_h and other is not None and any(
                    a["element"] == "H" and np.linalg.norm(a["xyz"] - other["xyz"]) < 1.35 for a in atoms)
            if resn in {"ASP", "GLU"} and has_h and atom["element"] == "O":
                atom["acceptor"], atom["donor"] = False, True
        protein.extend(a for a in atoms if a["element"] != "H")
        for ring_names in RINGS.get(resn, []):
            if all(n in names for n in ring_names):
                value = ring([names[n] for n in ring_names])
                if value:
                    protein_rings.append(value)
        charged = {"ASP": (-1, ("OD1", "OD2")), "GLU": (-1, ("OE1", "OE2")),
                   "LYS": (1, ("NZ",)), "ARG": (1, ("CZ", "NE", "NH1", "NH2"))}.get(resn)
        if charged and all(n in names for n in charged[1]):
            group = [names[n] for n in charged[1]]
            # Explicit protonation of an acid overrides its default charge.
            if charged[0] < 0 and any(a["element"] == "H" for g in group for a in g["neighbors"]):
                continue
            protein_charges.append(dict(xyz=np.mean([a["xyz"] for a in group], axis=0),
                                        atoms=group, charge=charged[0]))
    ligand_rings = []
    for indices in mol.GetRingInfo().AtomRings():
        if len(indices) in {5, 6} and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in indices):
            value = ring([ligand[i] for i in indices])
            if value:
                ligand_rings.append(value)
    ligand_charges, assigned = [], set()
    for feature in features:
        if feature.GetFamily() not in {"PosIonizable", "NegIonizable"}:
            continue
        ids = feature.GetAtomIds()
        charge = sum(mol.GetAtomWithIdx(i).GetFormalCharge() for i in ids)
        if not charge or assigned.intersection(ids):
            continue
        assigned.update(ids)
        ligand_charges.append(dict(xyz=np.mean([ligand[i]["xyz"] for i in ids], axis=0), charge=charge))
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() and atom.GetIdx() not in assigned:
            ligand_charges.append(dict(xyz=ligand[atom.GetIdx()]["xyz"], charge=atom.GetFormalCharge()))

    contacts = []

    def add(kind, start, end, atoms, detail):
        key = atoms[0]["residue"]
        distance = float(np.linalg.norm(start - end))
        contacts.append(dict(kind=kind, start=dict(zip("xyz", map(float, start))),
                             end=dict(zip("xyz", map(float, end))), distance=round(distance, 2),
                             residue=f"{key[3]} {key[0] or '∅'}:{key[1]}{key[2]}",
                             protein_serials=[a["serial"] for a in atoms], detail=detail))

    def angle(a, vertex, b):
        u, v = a - vertex, b - vertex
        norm = np.linalg.norm(u) * np.linalg.norm(v)
        return math.degrees(math.acos(float(np.clip(np.dot(u, v) / norm, -1, 1)))) if norm > 1e-8 else 0

    def hbond(donor, acceptor):
        hydrogens = [a for a in donor["neighbors"] if a["element"] == "H"]
        if hydrogens:
            return "H-bond" if any(np.linalg.norm(h["xyz"] - acceptor["xyz"]) <= 2.6 and
                                   angle(donor["xyz"], h["xyz"], acceptor["xyz"]) >= 120
                                   for h in hydrogens) else None
        heavy = [a for a in donor["neighbors"] if a["element"] != "H"]
        if heavy and all(angle(a["xyz"], donor["xyz"], acceptor["xyz"]) >= 90 for a in heavy):
            return "H-bond candidate (donor H absent)"
        return None

    def nonpolar(atom):
        return atom["element"] == "C" and all(a["element"] in {"C", "H", "S"} for a in atom["neighbors"])

    hydrophobic = {}
    for atom in protein:
        check()
        for i in tree.query_ball_point(atom["xyz"], 4.0):
            other = ligand[i]
            distance = np.linalg.norm(atom["xyz"] - other["xyz"])
            if distance < 0.5 or other["element"] == "H":
                continue
            if distance <= 3.5:
                evidence = None
                if atom["donor"] and other["acceptor"]:
                    evidence = hbond(atom, other)
                if other["donor"] and atom["acceptor"]:
                    evidence = evidence or hbond(other, atom)
                if evidence:
                    add("hbond", atom["xyz"], other["xyz"], [atom],
                        f"{evidence}: {atom['name']}–{other['name']}")
            if nonpolar(atom) and nonpolar(other):
                old = hydrophobic.get(atom["residue"])
                if old is None or distance < old[0]:
                    hydrophobic[atom["residue"]] = (distance, atom, other)
    for _, atom, other in hydrophobic.values():
        add("hydrophobic", atom["xyz"], other["xyz"], [atom], f"Nonpolar C···C: {atom['name']}–{other['name']}")

    def offset(vector, normal):
        return float(np.linalg.norm(vector - np.dot(vector, normal) * normal))

    for pr in protein_rings:
        check()
        for lr in ligand_rings:
            delta = lr["xyz"] - pr["xyz"]
            distance = np.linalg.norm(delta)
            theta = math.degrees(math.acos(float(np.clip(abs(np.dot(pr["normal"], lr["normal"])), 0, 1))))
            lateral = min(offset(delta, pr["normal"]), offset(delta, lr["normal"]))
            if 0.5 <= distance <= 5.5 and lateral <= 2.0 and (theta <= 30 or theta >= 60):
                add("pi_stack", pr["xyz"], lr["xyz"], pr["atoms"],
                    f"π–π {'parallel' if theta <= 30 else 'T-shaped'}; plane angle {theta:.1f}°, offset {lateral:.1f} Å")
        for charge in ligand_charges:
            delta = charge["xyz"] - pr["xyz"]
            if charge["charge"] > 0 and 0.5 <= np.linalg.norm(delta) <= 6 and offset(delta, pr["normal"]) <= 2:
                add("cation_pi", pr["xyz"], charge["xyz"], pr["atoms"], "Receptor π ring–ligand cation")
    for pc in protein_charges:
        check()
        for lc in ligand_charges:
            if pc["charge"] * lc["charge"] < 0 and 0.5 <= np.linalg.norm(pc["xyz"] - lc["xyz"]) <= 4:
                add("salt_bridge", pc["xyz"], lc["xyz"], pc["atoms"], "Oppositely charged group centers")
        if pc["charge"] > 0:
            for lr in ligand_rings:
                delta = pc["xyz"] - lr["xyz"]
                if 0.5 <= np.linalg.norm(delta) <= 6 and offset(delta, lr["normal"]) <= 2:
                    add("cation_pi", pc["xyz"], lr["xyz"], pc["atoms"], "Receptor cation–ligand π ring")
    check()
    return sorted(contacts, key=lambda c: (c["kind"], c["residue"], c["distance"])), NOTES
