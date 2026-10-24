"""Read-only docking visualization inputs; deliberately independent of Qt."""
from __future__ import annotations

from collections import Counter, OrderedDict
from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re
import sqlite3
import tempfile
from collections.abc import Callable

from ..project import ProjectRepository
from .postdock import PostDockService, _pdbqt_element


@dataclass(frozen=True)
class ViewerStructure:
    receptor_pdb: str
    ligand_sdf: str
    receptor_name: str
    parent_id: str
    state_id: str
    run_id: str
    pose_id: str
    mode_index: int
    affinity: float
    interactions: list[dict] = field(default_factory=list)
    interaction_notes: str = ""


def extract_pose(text: str, mode_index: int) -> str:
    """Select a MODEL by its number, retaining global and model-local remarks."""
    lines = text.splitlines()
    if not any(line.startswith("MODEL") for line in lines):
        if mode_index != 1 or not any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
            raise ValueError(f"Docking output does not contain mode {mode_index}.")
        return text
    global_remarks, selected = [], []
    current = None
    matches = 0
    for line in lines:
        if line.startswith("MODEL"):
            if current is not None:
                raise ValueError("Malformed docking output: nested MODEL records.")
            current = int(line.split()[1])
            if current == mode_index:
                matches += 1
        elif line.startswith("ENDMDL"):
            current = None
        elif current == mode_index:
            selected.append(line)
        elif current is None and line.startswith("REMARK"):
            global_remarks.append(line)
    if matches != 1 or not any(line.startswith(("ATOM  ", "HETATM")) for line in selected):
        raise ValueError(f"Docking output does not contain one valid mode {mode_index}.")
    return "\n".join(global_remarks + [f"MODEL {mode_index}"] + selected + ["ENDMDL", ""])


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _heavy_atoms(text: str, *, pdbqt: bool) -> Counter:
    atoms = Counter()
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = _pdbqt_element(line) if pdbqt else line[76:78].strip()
        if not element:
            element = line[12:16].strip().lstrip("0123456789")[:1]
        if element.upper() in {"H", "D", "TZ"}:
            continue
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError("Receptor contains non-finite coordinates.")
        atoms[(line[12:27], xyz, element.upper())] += 1
    return atoms


class StructureViewerService:
    def __init__(self, repository: ProjectRepository, *, cache_size: int = 12) -> None:
        self.repository = repository
        self.postdock = PostDockService()
        self.cache_size = cache_size
        self._cache: OrderedDict[tuple, str] = OrderedDict()
        self._interaction_cache: OrderedDict[tuple, tuple[list[dict], str]] = OrderedDict()

    def _connection(self) -> sqlite3.Connection:
        # Viewing must not migrate, change journal mode, or commit project state.
        conn = sqlite3.connect(self.repository.database_path.as_uri() + "?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        return conn

    def poses(self, run_id: str) -> list[dict]:
        conn = self._connection()
        try:
            return [dict(row) for row in conn.execute(
                "SELECT pose_id, mode_index, affinity FROM docking_poses WHERE run_id=? ORDER BY mode_index", (run_id,))]
        finally:
            conn.close()

    def _raw_artifact(self, relative_path: str, checksum: str,
                      cancelled: Callable[[], bool] | None = None) -> bytes:
        """Read the recorded artifact, including exact paths moved into RunN folders.

        Keep the database authoritative when its path exists. An archived copy
        is usable only if its complete relative path and recorded digest match.
        """
        root = self.repository.root.resolve()
        relative = Path(relative_path)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("The docking artifact is outside this project.")
        original = (root / relative).resolve()
        if not original.is_relative_to(root):
            raise ValueError("The docking artifact is outside this project.")
        if original.exists():
            data = original.read_bytes()
            if _digest(data) != checksum:
                raise ValueError("The raw docking artifact checksum does not match the recorded result.")
            return data
        found = False
        for folder in sorted(root.iterdir()):
            if not re.fullmatch(r"Run[0-9]+", folder.name) or not folder.is_dir():
                continue
            if cancelled and cancelled():
                raise InterruptedError("Pose loading cancelled")
            candidate = (folder / relative).resolve()
            if not candidate.is_relative_to(root) or not candidate.is_file():
                continue
            found = True
            data = candidate.read_bytes()
            if _digest(data) == checksum:
                return data
        if found:
            raise ValueError("Archived docking artifacts were found, but their checksums do not match the recorded result.")
        raise FileNotFoundError(
            "The recorded docking output is missing from the project and its Run folders. "
            f"Restore the original artifact: {relative_path}")

    def load(self, run_id: str, pose_id: str, *, cancelled: Callable[[], bool] | None = None) -> ViewerStructure:
        if cancelled and cancelled():
            raise InterruptedError("Pose loading cancelled")
        conn = self._connection()
        try:
            row = conn.execute("""SELECT d.run_id, d.receptor_profile_id, d.receptor_hash,
                    d.receptor_profile_name, s.parent_id, s.state_id,
                    p.pose_id, p.mode_index, p.affinity, a.relative_path, a.sha256
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                JOIN docking_poses p ON p.run_id=d.run_id
                JOIN artifacts a ON a.artifact_id=d.raw_output_artifact_id
                WHERE d.run_id=? AND p.pose_id=? AND d.status='completed'""", (run_id, pose_id)).fetchone()
        finally:
            conn.close()
        if row is None:
            raise ValueError("This completed docking pose or its raw artifact is unavailable.")
        profiles = self.repository.get_receptor_profiles(include_archived=True)
        profile = next((p for p in profiles if str(p["id"]) == row["receptor_profile_id"]), None)
        if profile is None or not profile.get("receptor"):
            raise ValueError("The receptor profile for this docking run is unavailable.")
        receptor = self.repository.root / str(profile["receptor"])
        receptor_bytes = receptor.read_bytes()
        if not row["receptor_hash"] or _digest(receptor_bytes) != row["receptor_hash"]:
            raise ValueError("The receptor has changed since docking; this complex cannot be matched safely.")
        raw_bytes = self._raw_artifact(row["relative_path"], row["sha256"], cancelled)
        receptor_text = receptor_bytes.decode("utf-8")
        prepared = receptor.parent / "receptor_prepared.pdb"
        # Compare physical heavy atoms as well: older projects lack a PDB checksum,
        # and a neighboring PDB must not silently belong to a different receptor.
        if prepared.is_file():
            pdb = prepared.read_text(encoding="utf-8")
            physical = _heavy_atoms(receptor_text, pdbqt=True)
            if not physical or physical != _heavy_atoms(pdb, pdbqt=False):
                raise ValueError("The prepared receptor PDB does not match the docking receptor coordinates and atoms.")
        else:
            with tempfile.TemporaryDirectory(prefix="mdp_view_receptor_") as directory:
                output = Path(directory) / "receptor.pdb"
                self.postdock.export_receptor_pdb(receptor, output)
                pdb = output.read_text(encoding="utf-8")
        # AD4Zn TZ sites are scoring pseudoatoms, never displayed as real atoms.
        pdb = "\n".join(line for line in pdb.splitlines()
                        if not (line.startswith(("ATOM  ", "HETATM")) and line[76:78].strip().upper() == "TZ")) + "\n"
        key = (run_id, pose_id, row["sha256"])
        sdf = self._cache.get(key)
        if sdf is None:
            selected = extract_pose(raw_bytes.decode("utf-8"), row["mode_index"])
            with tempfile.TemporaryDirectory(prefix="mdp_view_pose_") as directory:
                work = Path(directory)
                source = work / "pose.pdbqt"
                source.write_text(selected, encoding="utf-8")
                output = work / "pose.sdf"
                options = {"cancelled": cancelled} if cancelled else {}
                self.postdock.export_sdf(source, output, work / "conversion.log", **options)
                sdf = output.read_text(encoding="utf-8")
                if "M  END" not in sdf or sdf.count("$$$$") != 1:
                    raise ValueError("Pose conversion did not produce exactly one SDF molecule.")
            self._cache[key] = sdf
            while len(self._cache) > self.cache_size:
                self._cache.popitem(last=False)
        self._cache.move_to_end(key)
        analysis_key = (*key, _digest(pdb.encode("utf-8")))
        analysis = self._interaction_cache.get(analysis_key)
        if analysis is None:
            from .interactions import analyze_interactions
            try:
                analysis = analyze_interactions(pdb, sdf, cancelled=cancelled)
            except (ImportError, ValueError, RuntimeError) as exc:
                # A chemistry/annotation failure must not hide a valid complex.
                analysis = ([], f"Interaction analysis unavailable: {exc}")
            else:
                self._interaction_cache[analysis_key] = analysis
                while len(self._interaction_cache) > self.cache_size:
                    self._interaction_cache.popitem(last=False)
        else:
            self._interaction_cache.move_to_end(analysis_key)
        return ViewerStructure(pdb, sdf, str(profile.get("name", row["receptor_profile_name"])),
                               row["parent_id"], row["state_id"], run_id, pose_id,
                               row["mode_index"], row["affinity"], *analysis)
