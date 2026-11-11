"""Portable, symmetry-corrected redocking results."""

import csv
import json
import logging
import math
import os
import subprocess
import sys
import tempfile
from importlib.metadata import version
from pathlib import Path

from ..fingerprints import file_sha256

LOG = logging.getLogger(__name__)
REFERENCE = "reference_ligand_heavy.mol2"


def validate_mol2(path: Path) -> None:
    if not path.is_file():
        raise ValueError(f"Missing MOL2: {path.name}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if "@<TRIPOS>ATOM" not in lines or "@<TRIPOS>BOND" not in lines:
        raise ValueError(f"Invalid MOL2 sections: {path.name}")
    counts = lines[lines.index("@<TRIPOS>MOLECULE") + 2].split()
    if len(counts) < 2 or int(counts[0]) <= 0 or int(counts[1]) <= 0:
        raise ValueError(f"MOL2 contains no atoms or bonds: {path.name}")


def validation_runs(repository, profile):
    profile_id = str(profile["id"])
    folder = repository.root / "inputs" / "receptors" / profile_id / "redocking"
    reference = profile.get("reference_ligand") or {}
    receptor_hash = _hash(repository.root / str(profile.get("receptor") or ""))
    reference_hash = _hash(repository.root / str(reference.get("sdf") or ""))
    with repository.connection() as conn:
        local = {row["run_id"]: dict(row) for row in conn.execute(
            """SELECT r.run_id, r.status, r.receptor_sha256, r.reference_sha256, r.finished_at,
                      s.artifacts_json AS mol2_artifacts
               FROM redocking_runs r LEFT JOIN redocking_stages s
                   ON s.run_id=r.run_id AND s.stage_name='mol2_generation'
               WHERE r.receptor_profile_id=?""", (profile_id,))}
    runs = []
    for run in folder.iterdir() if folder.is_dir() else ():
        if not run.is_dir():
            continue
        result = read_results(run)
        provenance = local.get(run.name)
        if provenance is not None:
            if provenance.pop("status") != "ARTIFACTS_READY":
                continue
            provenance.pop("run_id")
            try:
                provenance["validated_reference_mol2_sha256"] = json.loads(
                    provenance.pop("mol2_artifacts") or "{}").get(REFERENCE)
            except (ValueError, AttributeError):
                continue
            if result and result.get("provenance") != provenance:
                result = None
        elif result:
            provenance = result.get("provenance")
        else:
            try:
                metadata = json.loads((run / "validation.json").read_text(encoding="utf-8"))
                if metadata.get("status") == "ARTIFACTS_READY":
                    provenance = metadata.get("provenance") or {}
            except (OSError, ValueError, AttributeError):
                continue
        if not isinstance(provenance, dict):
            continue
        current = bool(receptor_hash and reference_hash
                       and receptor_hash == provenance.get("receptor_sha256")
                       and reference_hash == provenance.get("reference_sha256"))
        runs.append((run, provenance, result, current))
    return sorted(runs, key=lambda item: (str(item[1].get("finished_at") or ""), item[0].name), reverse=True)


def _hash(path):
    try:
        return file_sha256(path) if path.is_file() else None
    except OSError:
        return None


def _pose_path(folder, rank):
    return folder / f"pose_{rank:03d}_heavy.mol2"


def _poses(run):
    with (run / "pose_manifest.csv").open(encoding="utf-8", newline="") as handle:
        rows = sorted((int(row["pose_rank"]), float(row["affinity"]), row["artifact_sha256"])
                      for row in csv.DictReader(handle))
    if not rows or [row[0] for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("The pose manifest must contain consecutive ranks starting at 1.")
    if any(not math.isfinite(row[1]) for row in rows):
        raise ValueError("The pose manifest contains an invalid affinity.")
    return rows


def read_results(run):
    folder = run / "dockrmsd"
    try:
        data = json.loads((folder / "results.json").read_text(encoding="utf-8"))
        if data["schema"] != 1 or data["reference_mol2_sha256"] != _hash(folder / REFERENCE):
            return None
        rows = data["poses"]
        validated_reference = data["provenance"].get("validated_reference_mol2_sha256")
        reference_verified = not validated_reference or validated_reference == data["reference_mol2_sha256"]
        if len(rows) != len(expected := _poses(run)):
            return None
        for row, (rank, affinity, expected_hash) in zip(rows, expected):
            if (row["pose_rank"] != rank or row["affinity"] != affinity
                    or row["mol2_sha256"] != _hash(_pose_path(folder, rank))):
                return None
            value = row["rmsd_angstrom"]
            if value is not None and (not reference_verified or not math.isfinite(value) or value < 0
                                      or row["mol2_sha256"] != expected_hash or not data["reference_mol2_sha256"]):
                return None
        return _summarize(data)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, csv.Error):
        return None


def _summarize(data):
    poses = data["poses"]
    successful = [row for row in poses if row["rmsd_angstrom"] is not None]
    best = min(successful, key=lambda row: (row["rmsd_angstrom"], row["pose_rank"]), default=None)
    data["rank_1"] = poses[0]
    data["lowest"] = best
    data["status"] = ("Calculated" if len(successful) == len(poses)
                      else "Incomplete" if successful else "Unavailable")
    return data


def _engine():
    return _measure, version("pydockrmsd")


def _measure(reference, pose):
    # The upstream extension parses MOL2 in native code. Isolate it so a
    # parser failure cannot terminate the desktop app or its docking queue.
    code = """
import json
import sys
from pydockrmsd.dockrmsd import PyDockRMSD
result = PyDockRMSD(sys.argv[1], sys.argv[2])
print(json.dumps([result.rmsd, result.error]))
"""
    folder = Path(reference).resolve().parent
    process = subprocess.run(
        [sys.executable, "-c", code, Path(reference).name, os.path.relpath(pose, folder)],
        cwd=folder, capture_output=True, text=True,
        encoding="utf-8", errors="replace", timeout=120,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
    )
    if process.returncode:
        raise RuntimeError(process.stderr.strip() or f"PyDockRMSD exited with code {process.returncode}")
    value, error = json.loads(process.stdout.strip().splitlines()[-1])
    if error:
        raise ValueError(error)
    return float(value)


def calculate_results(run, provenance):
    folder = run / "dockrmsd"
    reference = folder / REFERENCE
    rows = _poses(run)
    reference_hash = _hash(reference)
    error = None
    engine_version = None
    try:
        engine, engine_version = _engine()
        validate_mol2(reference)
        expected_reference = provenance.get("validated_reference_mol2_sha256")
        if expected_reference and reference_hash != expected_reference:
            raise ValueError("Checksum differs from the validated reference")
    except Exception as exc:
        error = f"{REFERENCE}: {exc}" if engine_version else f"PyDockRMSD unavailable: {exc}"
    cached = read_results(run)
    if (not error and cached and cached["status"] == "Calculated"
            and cached.get("engine_version") == engine_version):
        if cached.get("provenance") != provenance:
            cached["provenance"] = provenance
            _write_results(folder, cached)
        return cached
    poses = []
    for rank, affinity, expected in rows:
        path = _pose_path(folder, rank)
        digest = _hash(path)
        value = None
        pose_error = error
        try:
            if not digest or digest != expected:
                raise ValueError("Missing file or checksum differs from the validated pose manifest")
            if not error:
                validate_mol2(path)
                measured = engine(str(reference), str(path))
                if not math.isfinite(measured) or measured < 0:
                    raise ValueError("Engine returned an invalid RMSD")
                value = measured
        except Exception as exc:
            pose_error = f"{path.name}: {exc}"
        poses.append(dict(pose_rank=rank, affinity=affinity, mol2_sha256=digest,
                          rmsd_angstrom=value, error=pose_error))
    data = _summarize(dict(schema=1, tool="PyDockRMSD", engine_version=engine_version,
                           reference_mol2_sha256=reference_hash, provenance=provenance, poses=poses))
    _write_results(folder, data)
    return data


def _write_results(folder, data):
    target = folder / "results.json"
    with tempfile.NamedTemporaryFile(dir=folder, delete=False, mode="w", encoding="utf-8") as handle:
        temporary = Path(handle.name)
    try:
        temporary.write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")
        os.replace(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)


def format_summary(result):
    if result is None:
        return "RMSD not calculated"
    rank_1 = result["rank_1"]["rmsd_angstrom"]
    best = result["lowest"]
    rank_text = f"{rank_1:.3f} Å" if rank_1 is not None else "unavailable"
    best_text = f"{best['rmsd_angstrom']:.3f} Å (pose {best['pose_rank']})" if best else "unavailable"
    return f"{result['status']} · Rank 1 {rank_text} · Lowest {best_text}"
