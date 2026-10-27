"""Versioned, coordinate-preserving ProLIF analysis of existing postdock exports.

All chemistry dependencies are lazy imports. Failed records never become empty
fingerprints. The viewer's older geometric annotations remain a separate API.
"""
from __future__ import annotations

from contextlib import ExitStack
import csv
import hashlib
import importlib.metadata
import io
import json
import math
from pathlib import Path
import shutil
import tempfile
import uuid

from ..fingerprints import file_sha256, fingerprint
from ..models import ProgressEvent
from ..project import utc_now
from ..execution_lock import project_execution_lock
from .postdock import _pdbqt_element, receptor_export_filename


DETECTORS = ["HBDonor", "HBAcceptor", "Hydrophobic", "PiStacking", "CationPi",
             "PiCation", "Anionic", "Cationic", "VdWContact", "XBDonor",
             "XBAcceptor", "MetalDonor", "MetalAcceptor"]
# PDB/PDBQT round to 0.001 A; SDF rounds to 0.0001 A. Compare each
# coordinate component with a conservative 0.001 A serialization allowance.
COORDINATE_TOLERANCE = 0.001
METALS = {"ZN", "MG", "CA", "MN", "FE", "CU", "CO", "NI", "CD", "NA", "K"}
IDENTITY = ["analysis_id", "receptor_profile_id", "run_id", "ligand_id", "state_id",
            "pose_id", "source_record_id", "source_sdf", "record_index"]
RESIDUE = ["protein_chain", "protein_residue", "residue_number", "insertion_code"]
METRICS = ["total_interaction_pairs", "conserved_interaction_pairs", "lost_reference_pairs",
           "new_interaction_pairs", "reference_retention_percent", "jaccard_similarity",
           "conserved_reference_residues", "total_reference_residues", "residue_retention_percent"]
SUMMARY = IDENTITY + ["source_ligand_name", "source_pose_id", "affinity", "scoring_protocol",
                      "source_properties_json", "mapping_status"] + METRICS + ["analysis_status", "error_message"]
DETAILED = IDENTITY + ["is_reference"] + RESIDUE + ["interaction_type", "distance_angstrom",
    "angle_name", "angle_degrees", "geometry_json", "ligand_atom_indices", "protein_atom_indices",
    "conserved_with_reference", "analysis_status"]
COMPARISON = IDENTITY + RESIDUE + ["interaction_type", "reference_present", "candidate_present",
                                  "comparison", "analysis_status"]
REPORTS = ("interactions_detailed.csv", "interaction_summary.csv", "interaction_comparison.csv")


def compare_fingerprints(candidate: set, reference: set) -> dict:
    conserved = candidate & reference
    ref_residues = {pair[:4] for pair in reference}
    residues = {pair[:4] for pair in candidate}
    union = candidate | reference
    return dict(zip(METRICS, [len(candidate), len(conserved), len(reference - candidate),
        len(candidate - reference), 100 * len(conserved) / len(reference) if reference else None,
        len(conserved) / len(union) if union else None, len(residues & ref_residues),
        len(ref_residues), 100 * len(residues & ref_residues) / len(ref_residues) if ref_residues else None]))


def sdf_records(path: Path):
    """Yield one record at a time, including invalid and unterminated records."""
    from rdkit import Chem
    block = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        index = 0
        for line in handle:
            block.append(line)
            if line.strip() == "$$$$":
                index += 1
                text = "".join(block)
                yield index, next(iter(Chem.ForwardSDMolSupplier(io.BytesIO(text.encode()), removeHs=False)), None)
                block = []
        if any(line.strip() for line in block):
            index += 1
            yield index, next(iter(Chem.ForwardSDMolSupplier(io.BytesIO("".join(block).encode()), removeHs=False)), None)
        elif not index:
            yield 1, None


def validate_molecule(mol):
    import numpy as np
    if mol is None or not mol.GetNumAtoms():
        raise ValueError("Invalid or unreadable SDF molecule")
    if mol.GetNumConformers() != 1:
        raise ValueError("A pose must have exactly one molecular conformer")
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError("Non-finite molecular coordinates")


def heavy_coordinates(mol):
    validate_molecule(mol)
    conf = mol.GetConformer()
    return [(atom.GetSymbol().upper(), tuple(conf.GetAtomPosition(atom.GetIdx())))
            for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]


def pdb_coordinates(text, *, pdbqt=False):
    result = []
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            element = (_pdbqt_element(line) if pdbqt else line[76:78].strip()).upper()
            if element not in {"H", "D", "TZ"}:
                xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
                if not all(math.isfinite(x) for x in xyz):
                    raise ValueError("Non-finite PDB coordinates")
                result.append((element, xyz))
    return result


def same_coordinates(left, right):
    """Element-preserving bijection, no fitting, invariant to file atom order."""
    import numpy as np
    from scipy.spatial import cKDTree
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching
    if not left or len(left) != len(right):
        return False
    for element in {item[0] for item in left + right}:
        a = np.array([xyz for symbol, xyz in left if symbol == element])
        b = np.array([xyz for symbol, xyz in right if symbol == element])
        if len(a) != len(b):
            return False
        if not len(a):
            continue
        nearby = cKDTree(b).query_ball_point(a, COORDINATE_TOLERANCE, p=float("inf"))
        if any(not indices for indices in nearby):
            return False
        rows = [i for i, indices in enumerate(nearby) for _ in indices]
        cols = [j for indices in nearby for j in indices]
        graph = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(a), len(b)))
        if (maximum_bipartite_matching(graph, perm_type="column") < 0).any():
            return False
    return True


def _raw_modes(path):
    """Keep only one docking run in memory; never the entire campaign."""
    modes, current, lines = {}, 1, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL"):
            current, lines = int(line.split()[1]), []
        elif line.startswith("ENDMDL"):
            modes[current] = pdb_coordinates("\n".join(lines), pdbqt=True)
            lines = []
        else:
            lines.append(line)
    if any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
        modes[current] = pdb_coordinates("\n".join(lines), pdbqt=True)
    return modes


def pose_mapping(repository, run, mol, modes):
    coordinates = heavy_coordinates(mol)
    matched = [mode for mode, atoms in modes.items() if same_coordinates(coordinates, atoms)]
    if not matched:
        raise ValueError("SDF heavy-atom coordinates do not match any source docking pose (no alignment performed)")
    with repository.connection() as conn:
        rows = [dict(row) for row in conn.execute(
            "SELECT pose_id,mode_index,affinity FROM docking_poses WHERE run_id=?", (run["run_id"],))]
    poses = [row for row in rows if row["mode_index"] in matched]
    if len(poses) == 1:
        return poses[0], "verified"
    return None, "ambiguous: multiple or missing database poses; affinity unavailable"


def register_export(repository, run, path, *, cancelled=lambda: False):
    """Index an existing export without converting or changing molecular files."""
    digest, relative = file_sha256(path), path.relative_to(repository.root).as_posix()
    modes = _raw_modes(repository.root / run["raw_path"])
    for index, mol in sdf_records(path):
        if cancelled():
            raise InterruptedError("Interaction analysis cancelled")
        try:
            pose, status = pose_mapping(repository, run, mol, modes)
        except ValueError as exc:
            pose, status = None, str(exc)
        with repository.connection() as conn:
            conn.execute("INSERT OR REPLACE INTO interaction_exports VALUES (?,?,?,?,?,?)",
                (relative, index, digest, run["run_id"], pose["pose_id"] if pose else None, status))


def coordinate_subset(atoms, pool):
    from scipy.spatial import cKDTree
    trees = {element: cKDTree([xyz for symbol, xyz in pool if symbol == element])
             for element in {symbol for symbol, _ in pool}}
    return bool(atoms) and all(symbol in trees and trees[symbol].query(xyz, p=float("inf"))[0] <= COORDINATE_TOLERANCE
                              for symbol, xyz in atoms)


def prepared_analysis_receptor(repository, profile, runs=()):
    receptor = repository.root / profile["receptor"]
    prepared = receptor.parent / "receptor_prepared.pdb"
    for path in (receptor, prepared):
        if not path.is_file():
            raise FileNotFoundError(f"Required existing analysis input missing: {path}")
    digest = file_sha256(receptor)
    if any(run["receptor_hash"] != digest for run in runs):
        raise ValueError("Docking run receptor hash differs from the current prepared receptor")
    if not same_coordinates(pdb_coordinates(prepared.read_text()), pdb_coordinates(receptor.read_text(), pdbqt=True)):
        raise ValueError("Prepared PDB heavy atoms do not match the docking receptor")
    exported = repository.root / "For_PostDocking" / str(profile["id"]) / receptor_export_filename(str(profile.get("name", profile["id"])))
    if exported.is_file() and file_sha256(exported) != file_sha256(prepared):
        raise ValueError("Post-docking receptor PDB differs from the Meeko prepared PDB")
    return receptor, exported if exported.is_file() else prepared, prepared


def single_reference(path):
    records = iter(sdf_records(path))
    _, molecule = next(records)
    validate_molecule(molecule)
    if next(records, None) is not None:
        raise ValueError("Crystallographic reference must contain exactly one SDF record")
    return molecule


def attach_crystallographic_reference(repository, profile_id, source_sdf, *, confirm_same_frame=False):
    """Attach a supplied crystal pose for interactions only; never prepare/redock.

    Without immutable source coordinates, an SDF alone cannot prove its frame or
    crystallographic origin. Require and preserve an explicit user declaration.
    """
    with project_execution_lock(repository.root):
        profiles = repository.get_receptor_profiles(include_archived=True)
        profile = next((p for p in profiles if str(p["id"]) == profile_id), None)
        if profile is None:
            raise ValueError(f"Unknown receptor profile: {profile_id}")
        receptor, _, prepared = prepared_analysis_receptor(repository, profile)
        source_sdf = Path(source_sdf)
        original_digest = file_sha256(source_sdf)
        reference = single_reference(source_sdf)
        evidence = "user_confirmed_same_frame"
        evidence_paths = {}
        report_path = receptor.parent / "preparation_report.json"
        sources = sorted(receptor.parent.glob("source.*"))
        if report_path.is_file() and sources:
            report = json.loads(report_path.read_text(encoding="utf-8"))
            source = next((p for p in sources if file_sha256(p) == report.get("source_sha256")), None)
            if source is None:
                raise ValueError("Crystallographic source provenance changed; cannot verify the supplied reference")
            from ..receptors.extraction import structure_as_pdb
            source_atoms = pdb_coordinates(structure_as_pdb(source))
            if report.get("prepared_receptor_sha256") != file_sha256(receptor):
                raise ValueError("Docking receptor no longer matches the crystallographic preparation report")
            if not coordinate_subset(pdb_coordinates(prepared.read_text()), source_atoms):
                raise ValueError("Prepared receptor coordinates differ from the crystallographic source; no alignment performed")
            if not coordinate_subset(heavy_coordinates(reference), source_atoms):
                raise ValueError("Supplied SDF does not preserve crystallographic source coordinates; no alignment performed")
            evidence = "source_coordinates_verified"
            evidence_paths = {p.relative_to(repository.root).as_posix(): file_sha256(p) for p in (source, report_path)}
        elif not confirm_same_frame:
            raise ValueError("Original crystal coordinates are unavailable. Confirm this SDF is the co-crystallized pose in the prepared receptor's coordinate frame")
        parent = repository.root / "inputs" / "receptors" / profile_id / "interaction_references"
        parent.mkdir(parents=True, exist_ok=True)
        destination = parent / uuid.uuid4().hex
        with tempfile.TemporaryDirectory(prefix=".pending-", dir=parent) as temporary:
            staging = Path(temporary) / "reference"
            staging.mkdir()
            sdf = staging / "reference.sdf"
            shutil.copyfile(source_sdf, sdf)
            if file_sha256(sdf) != original_digest:
                raise ValueError("Reference SDF changed while being attached")
            # Validate the copied bytes too, before publishing metadata.
            single_reference(sdf)
            record = {"identity": reference.GetProp("_Name") if reference.HasProp("_Name") and reference.GetProp("_Name") else source_sdf.stem,
                "sdf": (destination / "reference.sdf").relative_to(repository.root).as_posix(),
                "sha256": file_sha256(sdf), "source_kind": "manual_crystallographic_sdf",
                "coordinate_validation": evidence, "declared_crystallographic": True,
                "receptor_sha256": file_sha256(receptor), "prepared_pdb_sha256": file_sha256(prepared),
                "evidence_sha256": evidence_paths, "imported_at": utc_now(), "original_filename": source_sdf.name}
            (staging / "attachment.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
            staging.rename(destination)
        profile["interaction_reference"] = record
        repository.save_receptor_profiles(profiles)
        repository.add_artifact(destination / "reference.sdf", "interaction_reference_sdf", "interactions")
        repository.add_artifact(destination / "attachment.json", "interaction_reference_metadata", "interactions")
        repository.record_provenance_event(event_type="interaction_reference_attached", stage_name="interactions",
            receptor_profile_id=profile_id, data=record)
        return record


def resolve_manual_inputs(repository, profile, runs):
    record = profile["interaction_reference"]
    receptor, chosen, prepared = prepared_analysis_receptor(repository, profile, runs)
    sdf = repository.root / record["sdf"]
    if record.get("source_kind") != "manual_crystallographic_sdf" or not record.get("declared_crystallographic"):
        raise ValueError("Manual reference is missing its crystallographic declaration")
    if file_sha256(sdf) != record["sha256"]:
        raise ValueError("Attached crystallographic reference SDF has changed; attach it again")
    if file_sha256(receptor) != record["receptor_sha256"] or file_sha256(prepared) != record["prepared_pdb_sha256"]:
        raise ValueError("Receptor changed since reference attachment; attach the compatible reference again")
    paths = {"receptor": chosen, "prepared_receptor": prepared, "docking_receptor": receptor, "reference": sdf}
    for name, digest in record.get("evidence_sha256", {}).items():
        path = repository.root / name
        if file_sha256(path) != digest:
            raise ValueError("Manual reference coordinate evidence changed")
        paths[name] = path
    return chosen, single_reference(sdf), paths


def resolve_inputs(repository, profile, runs):
    """Require provenance that establishes a shared frame, not spatial proximity."""
    if profile.get("interaction_reference"):
        return resolve_manual_inputs(repository, profile, runs)
    root = repository.root
    receptor = root / profile["receptor"]
    prepared = receptor.parent / "receptor_prepared.pdb"
    exported = root / "For_PostDocking" / str(profile["id"]) / receptor_export_filename(str(profile.get("name", profile["id"])))
    reference = profile.get("reference_ligand") or {}
    if not reference.get("sdf") or not reference.get("pdb"):
        raise ValueError("Profile needs existing crystallographic reference SDF and PDB metadata")
    sdf, crystal = root / reference["sdf"], root / reference["pdb"]
    report = receptor.parent / "preparation_report.json"
    cleaned = receptor.parent / "receptor_cleaned.pdb"
    for path in (receptor, prepared, sdf, crystal, report, cleaned):
        if not path.is_file():
            raise FileNotFoundError(f"Required existing analysis input missing: {path}")
    provenance = json.loads(report.read_text(encoding="utf-8"))
    receptor_hash = file_sha256(receptor)
    if provenance.get("prepared_receptor_sha256") != receptor_hash:
        raise ValueError("Prepared docking receptor no longer matches its preparation report")
    if provenance.get("cleaned_receptor_sha256") != file_sha256(cleaned):
        raise ValueError("Cleaned receptor no longer matches its preparation report")
    if any(run["receptor_hash"] != receptor_hash for run in runs):
        raise ValueError("Docking run receptor hash differs from the current prepared receptor")
    source = receptor.parent / next((p.name for p in receptor.parent.glob("source.*")), "source.pdb")
    if not source.is_file() or file_sha256(source) != provenance.get("source_sha256"):
        raise ValueError("Crystallographic source provenance is missing or changed")
    from ..receptors.extraction import structure_as_pdb
    source_atoms = pdb_coordinates(structure_as_pdb(source))
    # Cleaned receptor and extracted reference must be unchanged subsets of the
    # same source. A reference in a different frame cannot pass this test.
    clean_atoms = pdb_coordinates(cleaned.read_text())
    crystal_atoms = pdb_coordinates(crystal.read_text())
    if not clean_atoms or not crystal_atoms or not coordinate_subset(clean_atoms, source_atoms) or not coordinate_subset(crystal_atoms, source_atoms):
        raise ValueError("Receptor/reference extraction coordinates differ from the crystallographic source")
    if not same_coordinates(clean_atoms, pdb_coordinates(prepared.read_text())):
        raise ValueError("Prepared receptor heavy atoms differ from cleaned receptor; frame cannot be verified")
    chosen = exported if exported.is_file() else prepared
    if chosen != prepared and file_sha256(chosen) != file_sha256(prepared):
        raise ValueError("Post-docking receptor PDB differs from the Meeko prepared PDB")
    records = iter(sdf_records(sdf))
    _, ref = next(records)
    validate_molecule(ref)
    if next(records, None) is not None:
        raise ValueError("Crystallographic reference must contain exactly one SDF record")
    if not same_coordinates(heavy_coordinates(ref), crystal_atoms):
        raise ValueError("Reference SDF does not preserve extracted crystallographic coordinates")
    return chosen, ref, {"receptor": chosen, "docking_receptor": receptor, "reference": sdf,
        "reference_crystal": crystal, "preparation_report": report, "cleaned_receptor": cleaned, "source": source}


class ProlifEngine:
    def __init__(self, receptor):
        import MDAnalysis as mda
        import gemmi
        import numpy as np
        import prolif as plf
        from rdkit import Chem
        self.plf, self.Chem = plf, Chem
        self.fp = plf.Fingerprint(DETECTORS, count=True)
        self.residues, self.excluded = {}, []
        protein_lines, metal_lines, assigned, original_atoms = [], [], {}, {}
        # Renumber only analysis copies, before ProLIF can merge insertion codes.
        for line in receptor.read_text().splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            key = (line[21:22].strip(), line[17:20].strip(), int(line[22:26]), line[26:27].strip())
            name, element = key[1], line[76:78].strip().upper()
            amino = gemmi.find_tabulated_residue(name).is_amino_acid()
            metal = element in METALS and name.upper() in METALS
            if not amino and not metal:
                if line.startswith("ATOM"):
                    raise ValueError(f"Unsupported protein residue chemistry: {key}")
                if key not in self.excluded:
                    self.excluded.append(key)
                continue
            if key not in assigned:
                assigned[key] = len(assigned) + 1
                if assigned[key] > 9999:
                    raise ValueError("Analysis supports at most 9999 receptor residues")
                self.residues[assigned[key]] = key
            replacement = line[:22] + f"{assigned[key]:4d}" + " " + line[27:]
            atom_key = (assigned[key], line[12:16].strip())
            if atom_key in original_atoms:
                raise ValueError(f"Ambiguous duplicate receptor atom identity: {key} {atom_key[1]}")
            original_atoms[atom_key] = {"pdb_serial": int(line[6:11]), "atom_name": atom_key[1], "residue": key}
            (protein_lines if amino else metal_lines).append(replacement)
        if not protein_lines:
            raise ValueError("Prepared receptor contains no supported protein residues")
        protein_pdb = "\n".join(protein_lines) + "\nEND\n"
        universe = mda.Universe(io.StringIO(protein_pdb), format="PDB")
        # The prepared PDB must already provide explicit protein hydrogens.
        if not any(atom.element == "H" for atom in universe.atoms):
            raise ValueError("Prepared receptor lacks explicit hydrogens required for ProLIF")
        # Supply PDB-reader connectivity to MDAnalysis instead of its generic
        # vdW bond guesser, which can join donor H to a nearby acceptor. MDAnalysis
        # still assigns bond orders/charges; this does not alter any coordinates.
        topology = Chem.MolFromPDBBlock(protein_pdb, removeHs=False, sanitize=False)
        if topology is None or topology.GetNumAtoms() != len(universe.atoms):
            raise ValueError("Prepared receptor connectivity could not be mapped")
        universe.add_TopologyAttr("bonds", [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in topology.GetBonds()])
        protein = universe.atoms.convert_to("RDKIT", implicit_hydrogens=False)
        Chem.SanitizeMol(protein)
        if protein.GetNumAtoms() != len(universe.atoms) or not np.allclose(
                protein.GetConformer().GetPositions(), universe.atoms.positions, atol=1e-5, rtol=0):
            raise ValueError("Receptor conversion changed atom mapping or coordinates")
        # Metals are disconnected physical ions, never docking TZ pseudoatoms.
        # Do not infer covalent bonds to their coordinating protein atoms.
        for line in metal_lines:
            atom = Chem.Atom(line[76:78].strip().title())
            info = Chem.AtomPDBResidueInfo(line[12:16], int(line[6:11]), "", line[17:20],
                int(line[22:26]), line[21:22], "", 1.0, 0.0, True)
            atom.SetMonomerInfo(info)
            edit = Chem.RWMol()
            edit.AddAtom(atom)
            ion = edit.GetMol()
            conf = Chem.Conformer(1)
            conf.SetAtomPosition(0, tuple(float(line[i:i+8]) for i in (30, 38, 46)))
            ion.AddConformer(conf)
            protein = Chem.CombineMols(protein, ion)
        Chem.SanitizeMol(protein)
        self.atom_mapping = []
        for atom in protein.GetAtoms():
            info = atom.GetPDBResidueInfo()
            source = original_atoms[(info.GetResidueNumber(), info.GetName().strip())]
            self.atom_mapping.append({"analysis_atom_index": atom.GetIdx(), **source})
        self.protein = plf.Molecule(protein)

    def analyze(self, molecule):
        import numpy as np
        validate_molecule(molecule)
        before = heavy_coordinates(molecule)
        ligand = self.Chem.AddHs(self.Chem.Mol(molecule), addCoords=True)
        if before != heavy_coordinates(ligand):
            raise ValueError("Hydrogen addition changed heavy-atom coordinates")
        metadata = self.fp.generate(self.plf.Molecule(ligand), self.protein, metadata=True)
        contacts = []
        for (_, resid), interactions in metadata.items():
            key = self.residues[resid.number]
            for kind, records in interactions.items():
                for record in records:
                    geometry = {name: float(value) for name, value in record.items()
                                if isinstance(value, (int, float, np.number))}
                    angles = {name: value for name, value in geometry.items() if "angle" in name.lower()}
                    angle_name = next(iter(angles)) if len(angles) == 1 else None
                    indices = record.get("parent_indices", record.get("indices", {}))
                    contacts.append(dict(zip(RESIDUE, key), interaction_type=kind,
                        distance_angstrom=geometry.get("distance"), angle_name=angle_name,
                        angle_degrees=angles.get(angle_name), geometry_json=json.dumps(geometry, sort_keys=True),
                        ligand_atom_indices=json.dumps(indices.get("ligand", ())),
                        protein_atom_indices=json.dumps(indices.get("protein", ()))))
        return contacts


def contact_pair(contact):
    return tuple(contact[field] for field in RESIDUE + ["interaction_type"])


def _writer(stack, directory, name, fields):
    handle = stack.enter_context((directory / name).open("w", newline="", encoding="utf-8-sig"))
    writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    def write(row):
        writer.writerow({key: "NA" if row.get(key) is None else row[key] for key in fields})
    return write


def reusable_bundle(folder, signature):
    try:
        manifest = json.loads((folder / "analysis_manifest.json").read_text(encoding="utf-8"))
        return (manifest["fingerprint"] == signature
            and set(manifest["report_sha256"]) == set(REPORTS)
            and all(file_sha256(folder / name) == manifest["report_sha256"][name] for name in REPORTS))
    except (OSError, ValueError, KeyError, TypeError):
        return False


class InteractionAnalysisService:
    def __init__(self, repository, *, progress=None, cancelled=None, engine_factory=ProlifEngine):
        self.repository, self.progress = repository, progress
        self.cancelled = cancelled or (lambda: False)
        self.engine_factory = engine_factory

    def check(self):
        if self.cancelled():
            raise InterruptedError("Interaction analysis cancelled")

    def emit(self, event, **kwargs):
        if self.progress:
            self.progress(ProgressEvent(event, "interactions", **kwargs))

    def run(self, *, profile_id=None, run_id=None, force=False):
        repo = self.repository
        self.check()
        profiles = repo.get_receptor_profiles(include_archived=True)
        profiles = [p for p in profiles if (str(p["id"]) == profile_id if profile_id else not p.get("archived"))]
        if profile_id and not profiles:
            raise ValueError(f"Unknown receptor profile: {profile_id}")
        total_success = total_failed = 0
        found_inputs = False
        self.emit("stage_started")
        for profile in profiles:
            self.check()
            with repo.connection() as conn:
                runs = [dict(row) for row in conn.execute("""SELECT d.*,s.parent_id,s.state_id,a.relative_path raw_path
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    JOIN artifacts a ON a.artifact_id=d.raw_output_artifact_id
                    WHERE d.receptor_profile_id=? AND d.status='completed' AND d.is_current=1 AND s.active=1
                    ORDER BY d.run_id""", (str(profile["id"]),)) if not run_id or row["run_id"] == run_id]
            inputs = [(run, path) for run in runs for path in sorted((repo.root / "For_PostDocking" /
                str(profile["id"]) / "SDF" / run["parent_id"] / run["state_id"] / run["run_id"]).glob("*.sdf"))]
            if not inputs:
                if profile_id or run_id:
                    raise ValueError("No existing SDF exports for this scope; use the existing SDF Export action first")
                continue
            selected_runs = [run for run in runs if any(item[0]["run_id"] == run["run_id"] for item in inputs)]
            found_inputs = True
            try:
                success, failed = self._profile(profile, selected_runs, inputs, force, run_id)
            except InterruptedError:
                raise
            except Exception as exc:
                total_failed += 1
                self.emit("item_failed", receptor_profile_id=str(profile["id"]), failed=total_failed, message=str(exc))
                if profile_id or run_id:
                    raise
                continue
            total_success += success
            total_failed += failed
        if not found_inputs:
            raise ValueError("No existing SDF exports; use the existing SDF Export action first")
        self.emit("stage_completed", succeeded=total_success, failed=total_failed, total=total_success + total_failed)
        return total_success, total_failed

    def _profile(self, profile, runs, inputs, force, scope_run_id=None):
        repo, pid = self.repository, str(profile["id"])
        workflow = repo.start_workflow_run("interactions", {"profile_id": pid,
            "requested_run_id": scope_run_id, "run_ids": [r["run_id"] for r in runs]})
        aid = uuid.uuid4().hex
        with repo.connection() as conn:
            conn.execute("INSERT INTO interaction_runs(analysis_id,receptor_profile_id,workflow_run_id,fingerprint,status,started_at) VALUES (?,?,?,?,?,?)",
                         (aid, pid, workflow, "pending", "running", utc_now()))
        repo.update_workflow_stage(workflow, "interactions", "running")
        try:
            prepared, reference, paths = resolve_inputs(repo, profile, runs)
            versions = {name: importlib.metadata.version(name) for name in ("prolif", "rdkit", "MDAnalysis", "pandas", "gemmi")}
            settings = {"schema": 1, "detectors": DETECTORS, "count": True,
                "coordinate_tolerance_angstrom_per_axis": COORDINATE_TOLERANCE,
                "hydrogens": "preserve protein; RDKit AddHs(addCoords=True) on ligand copies; no optimization",
                "components": "protein and physical metals; waters/cofactors excluded", "direction": "ligand perspective"}
            hashes = {path.relative_to(repo.root).as_posix(): file_sha256(path) for path in paths.values()}
            input_errors = {}
            for run, path in inputs:
                for name in (path.relative_to(repo.root).as_posix(), run["raw_path"]):
                    if name in hashes:
                        continue
                    try:
                        hashes[name] = file_sha256(repo.root / name)
                    except OSError as exc:
                        hashes[name] = "unreadable"
                        input_errors[name] = str(exc)
            with repo.connection() as conn:
                pose_metadata = [dict(row) for run in runs for row in conn.execute(
                    "SELECT pose_id,run_id,mode_index,affinity FROM docking_poses WHERE run_id=? ORDER BY mode_index",
                    (run["run_id"],))]
            signature = fingerprint(settings={**settings, "versions": versions,
                "profile": profile, "poses": pose_metadata, "requested_run_id": scope_run_id,
                "run_ids": [r["run_id"] for r in runs]}, inputs=hashes, tool_version="prolif-analysis-v1")
            with repo.connection() as conn:
                conn.execute("UPDATE interaction_runs SET fingerprint=? WHERE analysis_id=?", (signature, aid))
                previous = conn.execute("""SELECT a.* FROM interaction_runs a JOIN workflow_runs w
                    ON w.workflow_run_id=a.workflow_run_id WHERE a.receptor_profile_id=?
                    AND a.status IN ('completed','completed_with_errors')
                    AND json_extract(w.settings_json,'$.requested_run_id') IS ?
                    ORDER BY a.finished_at DESC,a.rowid DESC LIMIT 1""", (pid, scope_run_id)).fetchall()
            if not force:
                for prior in previous:
                    folder = repo.root / prior["output_path"]
                    if prior["status"] == "completed" and prior["fingerprint"] == signature and reusable_bundle(folder, signature):
                        self._finish(aid, workflow, "skipped", succeeded=prior["succeeded"], output=prior["output_path"])
                        self.emit("item_skipped", receptor_profile_id=pid, message="Unchanged interaction analysis reused")
                        return prior["succeeded"], 0
            engine = self.engine_factory(prepared)
            self.check()
            ref_contacts = engine.analyze(reference)
            ref_pairs = {contact_pair(c) for c in ref_contacts}
            parent = repo.root / "For_PostDocking" / pid / "Interactions"
            parent.mkdir(parents=True, exist_ok=True)
            destination = parent / aid
            succeeded = failed = 0
            with tempfile.TemporaryDirectory(prefix=".pending-", dir=parent) as temporary:
                staging = Path(temporary) / "reports"
                staging.mkdir()
                with ExitStack() as stack:
                    detail = _writer(stack, staging, REPORTS[0], DETAILED)
                    summary = _writer(stack, staging, REPORTS[1], SUMMARY)
                    matrix = _writer(stack, staging, REPORTS[2], COMPARISON)
                    ref_identity = {"analysis_id": aid, "receptor_profile_id": pid,
                        "ligand_id": str((profile.get("interaction_reference") or profile.get("reference_ligand") or {}).get("identity", "reference")),
                        "pose_id": "crystallographic_reference", "source_sdf": str(paths["reference"].relative_to(repo.root)), "record_index": 1}
                    for contact in ref_contacts:
                        detail({**ref_identity, **contact, "is_reference": True,
                                "conserved_with_reference": True, "analysis_status": "success"})
                    last_run, modes, mode_error = None, {}, None
                    for run, path in inputs:
                        self.check()
                        if last_run != run["run_id"]:
                            try:
                                modes = _raw_modes(repo.root / run["raw_path"])
                                mode_error = None
                            except (OSError, ValueError) as exc:
                                modes, mode_error = {}, str(exc)
                            last_run = run["run_id"]
                        relative = path.relative_to(repo.root).as_posix()
                        records = [(1, None)] if relative in input_errors else sdf_records(path)
                        for index, mol in records:
                            self.check()
                            record_id = hashlib.sha256(f"{relative}:{index}".encode()).hexdigest()[:24]
                            row = dict.fromkeys(SUMMARY)
                            row.update(analysis_id=aid, receptor_profile_id=pid, run_id=run["run_id"],
                                ligand_id=run["parent_id"], state_id=run["state_id"], pose_id=record_id,
                                source_record_id=record_id, source_sdf=relative, record_index=index,
                                scoring_protocol=json.loads(run["command_json"] or "{}").get("settings", {}).get("protocol", "vina"))
                            try:
                                if relative in input_errors:
                                    raise ValueError(f"Unreadable SDF: {input_errors[relative]}")
                                validate_molecule(mol)
                                row["source_ligand_name"] = mol.GetProp("_Name") if mol.HasProp("_Name") else None
                                row["source_properties_json"] = json.dumps({name: mol.GetProp(name) for name in mol.GetPropNames()}, sort_keys=True)
                                row["source_pose_id"] = next((mol.GetProp(k) for k in ("POSE_ID", "pose_id", "POSE_RANK") if mol.HasProp(k)), None)
                                if mode_error:
                                    raise ValueError(f"Source docking pose unavailable for coordinate verification: {mode_error}")
                                pose, mapping = pose_mapping(repo, run, mol, modes)
                                row["mapping_status"] = mapping
                                if pose:
                                    row.update(pose_id=pose["pose_id"], affinity=pose["affinity"])
                                contacts = engine.analyze(mol)
                                pairs = {contact_pair(c) for c in contacts}
                                row.update(compare_fingerprints(pairs, ref_pairs), analysis_status="success", error_message=None)
                            except InterruptedError:
                                raise
                            except Exception as exc:
                                contacts, pairs = [], None
                                row.update(analysis_status="failed", error_message=str(exc))
                                failed += 1
                            else:
                                succeeded += 1
                            # Commit output only after calculation completes for this record.
                            for contact in contacts:
                                detail({**row, **contact, "is_reference": False,
                                    "conserved_with_reference": contact_pair(contact) in ref_pairs})
                            for pair in sorted(ref_pairs | (pairs or set())):
                                present = pair in pairs if pairs is not None else None
                                matrix({**row, **dict(zip(RESIDUE + ["interaction_type"], pair)),
                                    "reference_present": pair in ref_pairs, "candidate_present": present,
                                    "comparison": "unknown" if present is None else "conserved" if present and pair in ref_pairs else "new" if present else "lost"})
                            summary(row)
                            with repo.connection() as conn:
                                conn.execute("INSERT OR REPLACE INTO interaction_exports VALUES (?,?,?,?,?,?)",
                                    (relative, index, hashes[relative], run["run_id"], pose["pose_id"] if row["analysis_status"] == "success" and pose else None,
                                     row["mapping_status"] or row["error_message"]))
                                conn.execute("INSERT INTO interaction_pose_summaries VALUES (?,?,?)", (aid, record_id, json.dumps(row)))
                            self.emit("item_succeeded" if pairs is not None else "item_failed", receptor_profile_id=pid,
                                item_id=record_id, index=succeeded+failed, succeeded=succeeded, failed=failed,
                                message=f"{succeeded} poses analyzed; {failed} failed")
                self.check()
                with (staging / REPORTS[1]).open(encoding="utf-8-sig", newline="") as handle:
                    reader = csv.DictReader(handle)
                    if reader.fieldnames != SUMMARY or sum(1 for _ in reader) != succeeded + failed:
                        raise ValueError("Interaction summary failed publication validation")
                # Detect source changes during a long-running analysis before publication.
                if any(file_sha256(repo.root / name) != digest for name, digest in hashes.items() if name not in input_errors):
                    raise ValueError("Analysis inputs changed during execution; reports were not published")
                manifest = dict(analysis_id=aid, fingerprint=signature, settings=settings, versions=versions,
                    input_sha256=hashes, reference_pairs=[list(p) for p in sorted(ref_pairs)],
                    input_errors=input_errors,
                    reference_origin=profile.get("interaction_reference") or {"source_kind": "pipeline_crystallographic_reference"},
                    residue_mapping=engine.residues, excluded_components=engine.excluded,
                    protein_atom_mapping=getattr(engine, "atom_mapping", []),
                    coordinate_checks=("prepared PDB vs docking receptor; reference evidence: "
                        + profile["interaction_reference"]["coordinate_validation"]
                        + "; each candidate vs raw pose; no alignment" if profile.get("interaction_reference")
                        else "source extraction, prepared receptor, reference SDF, each candidate vs raw pose; no alignment"),
                    succeeded=succeeded, failed=failed, status="completed" if not failed else "completed_with_errors",
                    report_sha256={name: file_sha256(staging / name) for name in REPORTS})
                (staging / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
                staging.rename(destination)
            relative_output = destination.relative_to(repo.root).as_posix()
            for path in destination.iterdir():
                repo.add_artifact(path, "interaction_report", "interactions")
            status = "completed" if not failed else "completed_with_errors"
            self._finish(aid, workflow, status, succeeded=succeeded, failed=failed, output=relative_output)
            if succeeded or not failed:
                pointer = parent / "latest.json.tmp"
                pointer.write_text(json.dumps({"analysis_id": aid, "output_path": relative_output}), encoding="utf-8")
                pointer.replace(parent / "latest.json")
            repo.record_provenance_event(event_type="interactions_published", workflow_run_id=workflow,
                stage_name="interactions", receptor_profile_id=pid, data={"analysis_id": aid, "output_path": relative_output})
            return succeeded, failed
        except BaseException as exc:
            self._finish(aid, workflow, "cancelled" if isinstance(exc, InterruptedError) else "failed", error=str(exc))
            raise

    def _finish(self, aid, workflow, status, *, succeeded=0, failed=0, output=None, error=None):
        repo = self.repository
        with repo.connection() as conn:
            conn.execute("UPDATE interaction_runs SET status=?,succeeded=?,failed=?,output_path=?,finished_at=?,error_message=? WHERE analysis_id=?",
                         (status, succeeded, failed, output, utc_now(), error, aid))
        workflow_status = "completed" if status in {"completed", "completed_with_errors", "skipped"} else status
        repo.update_workflow_stage(workflow, "interactions", workflow_status, summary={"succeeded": succeeded, "failed": failed, "output_path": output}, error=error)
        repo.finish_workflow_run(workflow, workflow_status, error=error)
