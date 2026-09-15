from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .project import ProjectRepository


MANIFEST_HEADERS = [
    "receptor_id", "receptor_name", "id", "smiles", "inchikey", "admet_status", "admet_reason",
    "sdf_status", "sdf_path", "sdf_reason", "pdbqt_status", "pdbqt_path", "pdbqt_reason",
    "vina_status", "vina_score", "vina_pose", "vina_reason", "config_hash", "receptor_sha1",
    "tools_rdkit", "tools_meeko", "tools_vina", "created_at", "updated_at", "protocol",
]
LEADERBOARD_HEADERS = [
    "receptor_id", "receptor_name", "parent_id", "state_id", "mode_index", "affinity",
    "rmsd_lb", "rmsd_ub", "run_id", "pose_rank", "protocol",
]


def _write_csv(path: Path, headers: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader(); writer.writerows(rows)
    temporary.replace(path)


def export_manifest_csv(repository: ProjectRepository) -> list[Path]:
    profiles = repository.get_receptor_profiles()
    if not profiles:
        raise ValueError("No active receptor profiles are available for manifest export.")
    outputs: list[Path] = []
    combined: list[dict[str, Any]] = []
    with repository.connection() as conn:
        parents = conn.execute("""SELECT p.*, COALESCE(s.reason, '') screening_reason,
            COALESCE(s.status, '') screening_status FROM parent_ligands p
            LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
            WHERE p.active=1 ORDER BY p.parent_id""").fetchall()
        for profile in profiles:
            profile_id = str(profile["id"]); rows = []
            for parent in parents:
                states = conn.execute("""SELECT s.*, c.status conformer_status, a.relative_path sdf_path
                    FROM molecular_states s LEFT JOIN conformers c ON c.state_id=s.state_id
                    LEFT JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                    WHERE s.parent_id=? AND s.active=1 ORDER BY s.state_id""", (parent["parent_id"],)).fetchall()
                run = conn.execute("""SELECT d.*, raw.relative_path pose_path, p.affinity
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    LEFT JOIN docking_poses p ON p.run_id=d.run_id
                    LEFT JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id
                    WHERE s.parent_id=? AND s.active=1 AND d.receptor_profile_id=? AND d.is_current=1
                    ORDER BY (p.affinity IS NULL),p.affinity LIMIT 1""", (parent["parent_id"], profile_id)).fetchone()
                pdbqt_paths = [repository.root / "artifacts" / "pdbqt" / str(parent["parent_id"]) / str(state["state_id"]) / "ligand.pdbqt" for state in states]
                row = {
                    "receptor_id": profile_id, "receptor_name": profile.get("name", profile_id),
                    "id": parent["parent_id"], "smiles": parent["source_smiles"], "inchikey": parent["parent_inchikey"] or "",
                    "admet_status": parent["screening_status"], "admet_reason": parent["screening_reason"],
                    "sdf_status": "DONE" if states else "", "sdf_path": ";".join(str(s["sdf_path"] or "") for s in states), "sdf_reason": "",
                    "pdbqt_status": "DONE" if states and all(s["conformer_status"] == "pdbqt_ready" for s in states) else "",
                    "pdbqt_path": ";".join(path.relative_to(repository.root).as_posix() for path in pdbqt_paths if path.is_file()), "pdbqt_reason": "",
                    "vina_status": run["status"] if run else "", "vina_score": run["affinity"] if run and run["affinity"] is not None else "",
                    "vina_pose": run["pose_path"] if run and run["pose_path"] else "", "vina_reason": run["reason"] if run and run["reason"] else "",
                    "config_hash": run["settings_fingerprint"] if run else "", "receptor_sha1": run["receptor_hash"] if run else "",
                    "tools_rdkit": "RDKit", "tools_meeko": "Meeko", "tools_vina": "Vina",
                    "created_at": parent["created_at"], "updated_at": parent["created_at"],
                    "protocol": json.loads(run["command_json"]).get("settings", {}).get("protocol", "vina") if run else profile.get("protocol", "vina"),
                }
                rows.append(row); combined.append(row)
            output = repository.root / "exports" / profile_id / "manifest.csv"
            _write_csv(output, MANIFEST_HEADERS, rows); outputs.append(output)
    aggregate = repository.root / "exports" / "manifest.csv"
    _write_csv(aggregate, MANIFEST_HEADERS, combined); outputs.append(aggregate)
    for output in outputs:
        repository.add_artifact(output, "manifest_csv", "export")
    repository.record_provenance_event(event_type="manifest_exported", stage_name="export",
        data={"outputs": [path.relative_to(repository.root).as_posix() for path in outputs], "rows": len(combined)})
    return outputs


def export_leaderboard_csv(repository: ProjectRepository) -> list[Path]:
    profiles = repository.get_receptor_profiles()
    if not profiles:
        raise ValueError("No active receptor profiles are available for leaderboard export.")
    outputs: list[Path] = []
    combined: list[dict[str, Any]] = []
    with repository.connection() as conn:
        for profile in profiles:
            profile_id = str(profile["id"])
            query_rows = conn.execute("""SELECT parent_id,state_id,mode_index,affinity,rmsd_lb,rmsd_ub,run_id,pose_rank,command_json
                FROM (SELECT s.parent_id,s.state_id,p.mode_index,p.affinity,p.rmsd_lb,p.rmsd_ub,d.run_id,d.command_json,
                    ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity ASC) pose_rank
                    FROM docking_poses p JOIN docking_runs d ON d.run_id=p.run_id
                    JOIN molecular_states s ON s.state_id=d.state_id JOIN parent_ligands l ON l.parent_id=s.parent_id
                    WHERE l.active=1 AND s.active=1 AND d.is_current=1 AND d.status='completed' AND d.receptor_profile_id=?)
                WHERE pose_rank<=3 ORDER BY affinity,parent_id,pose_rank""", (profile_id,)).fetchall()
            rows = [{"receptor_id": profile_id, "receptor_name": profile.get("name", profile_id), **dict(row)} for row in query_rows]
            for row in rows:
                row["protocol"] = json.loads(row.pop("command_json")).get("settings", {}).get("protocol", "vina")
            combined.extend(rows)
            output = repository.root / "exports" / profile_id / "leaderboard.csv"
            _write_csv(output, LEADERBOARD_HEADERS, rows); outputs.append(output)
    aggregate = repository.root / "exports" / "leaderboard.csv"
    _write_csv(aggregate, LEADERBOARD_HEADERS, combined); outputs.append(aggregate)
    for output in outputs:
        repository.add_artifact(output, "leaderboard_csv", "export")
    repository.record_provenance_event(event_type="leaderboard_exported", stage_name="export",
        data={"outputs": [path.relative_to(repository.root).as_posix() for path in outputs], "rows": len(combined)})
    return outputs
