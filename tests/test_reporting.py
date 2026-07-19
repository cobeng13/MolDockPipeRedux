from __future__ import annotations

import json
from pathlib import Path

from moldockpipe.project import ProjectRepository, utc_now
from moldockpipe.reporting import generate_project_report


def test_generate_project_report_uses_structured_provenance(tmp_path: Path) -> None:
    repo = ProjectRepository.create(tmp_path / "project", "Report Project")
    repo.save_receptor_profiles([{
        "id": "rec-a", "name": "Receptor A", "enabled": True, "archived": False,
        "receptor": "inputs/receptors/rec-a/receptor.pdbqt",
        "center_x": 1, "center_y": 2, "center_z": 3,
        "size_x": 20, "size_y": 21, "size_z": 22,
        "exhaustiveness": 16, "num_modes": 9, "energy_range": 3, "seed": 42, "cpu_count": 2,
    }])
    source = tmp_path / "source.pdb"; source.write_text("ATOM", encoding="ascii")
    prep_id = repo.start_receptor_preparation(
        profile_id="rec-a", receptor_name="Receptor A", source_path=source, source_sha256="source-hash",
        inventory={"selected_model": 0, "included_chains": ["A"],
                   "counts": {"protein_residues": 100, "waters": 4, "nonpolymer_components": 1},
                   "components": [{"identity": "LIG A:401", "category": "organic ligand",
                                   "suggested_role": "reference_ligand", "heavy_atom_count": 20,
                                   "classification_reason": "fixture"}]},
        decisions={"center_method": "reference_ligand_centroid", "box_center": [1, 2, 3],
                   "box_method": "ligand_envelope_padding", "box_size": [20, 21, 22],
                   "altloc_choices": {}, "excluded_receptor_residues": []},
    )
    repo.finish_receptor_preparation(prep_id, "completed")
    csv_path = tmp_path / "ligands.csv"; csv_path.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(csv_path)
    workflow_id = repo.start_workflow_run("docking", {"requested_dockings": 1})
    repo.update_workflow_stage(workflow_id, "vina", "completed", summary={"succeeded": 1})
    repo.finish_workflow_run(workflow_id, "completed")
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,state_structure_hash,
             generation_fingerprint,status,enumeration_truncated,created_at,active)
            VALUES ('state-1','lig1','CCO','CCO',0,'structure','generation','pdbqt_ready',0,?,1)""", (utc_now(),))
        conn.execute("""INSERT INTO docking_runs
            (run_id,state_id,receptor_hash,settings_fingerprint,status,command_json,return_code,
             started_at,ended_at,is_current,receptor_profile_id,receptor_profile_name,workflow_run_id)
            VALUES ('dock-1','state-1','receptor-hash','settings','completed','{}',0,?,?,1,'rec-a','Receptor A',?)""",
            (utc_now(), utc_now(), workflow_id))
        conn.execute("""INSERT INTO docking_poses
            (pose_id,run_id,mode_index,affinity,rmsd_lb,rmsd_ub,is_best_for_state)
            VALUES ('pose-1','dock-1',1,-8.5,0,0,1)""")

    html_path = generate_project_report(repo, output_root=repo.root / "exports" / "generated-report")
    data = json.loads((html_path.parent / "report_data.json").read_text(encoding="utf-8"))
    manifest = json.loads((html_path.parent / "source_run_manifest.json").read_text(encoding="utf-8"))
    receptor = data["receptors"][0]
    assert receptor["preparation"]["inventory"]["counts"]["protein_residues"] == 100
    assert receptor["docking"]["best_affinity_per_parent_statistics"]["best_minimum"] == -8.5
    assert receptor["docking"]["top_compounds"][0]["parent_id"] == "lig1"
    assert manifest["runs"]["current_docking_runs"] == ["dock-1"]
    assert "Receptor preparation" in html_path.read_text(encoding="utf-8")
    with repo.connection() as conn:
        snapshot = conn.execute("SELECT output_path,output_sha256 FROM report_snapshots").fetchone()
        report_artifacts = conn.execute("SELECT COUNT(*) FROM artifacts WHERE created_stage='report'").fetchone()[0]
    assert snapshot["output_path"].endswith("report.html") and len(snapshot["output_sha256"]) == 64
    assert report_artifacts == 3
