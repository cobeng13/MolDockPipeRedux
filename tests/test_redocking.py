from __future__ import annotations

import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from moldockpipe.project import ProjectRepository
from moldockpipe.receptors.ligand_chemistry import chemistry_summary, read_sdf, validate_same_heavy_graph, write_mol2, write_sdf
from moldockpipe.redocking.models import RedockingSettings
from moldockpipe.redocking.runner import RedockingRunner, validate_mol2, validate_redocking_prerequisites
from moldockpipe.services.meeko import MeekoResult
from moldockpipe.services.vina import VinaResult
from moldockpipe.services.validation import parse_vina_poses


def valid_pdbqt() -> str:
    atoms = [f"ATOM  {i:5d}  C{i:<2} LIG A   1    {i:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00     0.000 C" for i in range(1, 6)]
    return "ROOT\n" + "\n".join(atoms) + "\nENDROOT\nTORSDOF 0\n"


@pytest.fixture()
def redocking_project(tmp_path: Path):
    repo = ProjectRepository.create(tmp_path / "project"); profile_id = "receptor-a"
    root = repo.root / "inputs" / "receptors" / profile_id; reference = root / "reference_ligand"; reference.mkdir(parents=True)
    receptor = root / "receptor.pdbqt"; receptor.write_text(valid_pdbqt(), encoding="utf-8")
    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO")); AllChem.EmbedMolecule(molecule, randomSeed=7)
    sdf = reference / "reference_ligand.sdf"; mol2 = reference / "reference_ligand.mol2"
    write_sdf(molecule, sdf); write_mol2(molecule, mol2)
    (reference / "chemistry_mapping.json").write_text(json.dumps({"chemistry_source": "test_template", **chemistry_summary(molecule)}), encoding="utf-8")
    profile = {"id": profile_id, "name": "Test receptor", "enabled": True, "archived": False,
        "receptor": receptor.relative_to(repo.root).as_posix(), "center_x": 1, "center_y": 2, "center_z": 3,
        "size_x": 20, "size_y": 21, "size_z": 22, "exhaustiveness": 8, "num_modes": 9, "energy_range": 3,
        "seed": 42, "cpu_count": 1, "reference_ligand": {"identity": "LIG A:401", "sdf": sdf.relative_to(repo.root).as_posix(),
        "mol2": mol2.relative_to(repo.root).as_posix(), "mapping": (reference / "chemistry_mapping.json").relative_to(repo.root).as_posix()}}
    repo.save_receptor_profiles([profile]); return repo, profile, molecule


class FakeMeeko:
    def __init__(self, fail: bool = False): self.calls = 0; self.fail = fail
    def prepare_ligand(self, sdf, pdbqt, log):
        self.calls += 1
        if self.fail: raise RuntimeError("mock Meeko failure")
        pdbqt.write_text(valid_pdbqt(), encoding="utf-8"); log.write_text("Meeko OK", encoding="utf-8")
        return MeekoResult(("python", "mk_prepare_ligand", "-i", str(sdf), "-o", str(pdbqt)), 0, "OK", "")


class FakeVina:
    def __init__(self, fail: bool = False, cancel: bool = False): self.calls = 0; self.fail = fail; self.cancel = cancel
    def run(self, *, output, log, cancelled, **kwargs):
        self.calls += 1
        if self.cancel:
            kwargs.get("settings"); raise InterruptedError("Vina redocking was cancelled")
        if self.fail: return VinaResult(1, [], "", "mock Vina failure", ("vina",))
        output.write_text("MODEL 1\nREMARK VINA RESULT: -8.200 0.000 0.000\nENDMDL\nMODEL 2\nREMARK VINA RESULT: -7.100 1.000 2.000\nENDMDL\n", encoding="utf-8")
        log.write_text("Vina OK", encoding="utf-8"); return VinaResult(0, [(1, -8.2, 0, 0), (2, -7.1, 1, 2)], "OK", "", ("vina",))


class FakePostDock:
    def __init__(self, molecule, mismatch: bool = False): self.molecule = molecule; self.mismatch = mismatch; self.calls = 0
    def export_sdf(self, input_pdbqt, output_sdf, log):
        self.calls += 1; writer = Chem.SDWriter(str(output_sdf))
        for rank in range(2):
            molecule = Chem.Mol(self.molecule)
            if self.mismatch and rank == 1: molecule = Chem.RWMol(molecule); molecule.RemoveAtom(0); molecule = molecule.GetMol()
            writer.write(molecule)
        writer.close(); log.write_text("Export OK", encoding="utf-8")


def runner_for(repo, profile, molecule, **kwargs):
    return RedockingRunner(repo, profile, RedockingSettings(32, 2, 5, 123, 1), meeko=kwargs.get("meeko", FakeMeeko()),
        vina=kwargs.get("vina", FakeVina()), postdock=kwargs.get("postdock", FakePostDock(molecule)))


def test_prerequisites_report_actionable_missing(redocking_project) -> None:
    repo, profile, _ = redocking_project; assert validate_redocking_prerequisites(repo, profile) == []
    broken = dict(profile); broken["reference_ligand"] = {}
    assert "Reference ligand SDF" in validate_redocking_prerequisites(repo, broken)


def test_mol2_writer_and_graph_validation(redocking_project, tmp_path: Path) -> None:
    _, _, molecule = redocking_project; path = tmp_path / "ligand.mol2"; write_mol2(molecule, path); validate_mol2(path)
    text = path.read_text(); assert "@<TRIPOS>ATOM" in text and "@<TRIPOS>BOND" in text
    parsed = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    assert parsed is not None and parsed.GetNumAtoms() == molecule.GetNumHeavyAtoms()
    validate_same_heavy_graph(molecule, Chem.Mol(molecule))
    mismatch = Chem.RWMol(molecule); mismatch.RemoveAtom(0)
    with pytest.raises(ValueError, match="does not match"): validate_same_heavy_graph(molecule, mismatch.GetMol())


def test_pose_affinity_parsing(tmp_path: Path) -> None:
    output = tmp_path / "poses.pdbqt"; output.write_text("REMARK VINA RESULT: -9.100 0.000 0.000\nREMARK VINA RESULT: -8.200 1.000 2.000\n")
    assert parse_vina_poses(output) == [(1, -9.1, 0.0, 0.0), (2, -8.2, 1.0, 2.0)]


def test_successful_redocking_generates_all_artifacts(redocking_project) -> None:
    repo, profile, molecule = redocking_project; result = runner_for(repo, profile, molecule).run()
    assert result["status"] == "ARTIFACTS_READY" and result["poses"] == 2 and result["top_affinity"] == -8.2
    root = Path(result["run_root"]); assert (root / "dockrmsd/reference_ligand_heavy.mol2").is_file()
    assert (root / "dockrmsd/pose_001_heavy.mol2").is_file()
    assert (root / "dockrmsd/top_ranked_pose_heavy.mol2").read_bytes() == (root / "dockrmsd/pose_001_heavy.mol2").read_bytes()
    assert "DockRMSD" in (root / "dockrmsd/README.txt").read_text()
    with repo.connection() as conn:
        assert conn.execute("SELECT status FROM redocking_runs").fetchone()[0] == "ARTIFACTS_READY"
        assert conn.execute("SELECT COUNT(*) FROM redocking_poses").fetchone()[0] == 2


@pytest.mark.parametrize(("service", "message"), (("meeko", "mock Meeko failure"), ("vina", "mock Vina failure")))
def test_external_tool_failure_records_exact_stage(redocking_project, service: str, message: str) -> None:
    repo, profile, molecule = redocking_project
    kwargs = {service: FakeMeeko(True) if service == "meeko" else FakeVina(True)}
    with pytest.raises(RuntimeError, match=message): runner_for(repo, profile, molecule, **kwargs).run()
    with repo.connection() as conn:
        row = conn.execute("SELECT status,error_stage FROM redocking_runs").fetchone()
    assert row["status"] == "failed" and row["error_stage"] == ("ligand_preparation" if service == "meeko" else "vina_redocking")


def test_graph_mismatch_blocks_artifacts_ready(redocking_project) -> None:
    repo, profile, molecule = redocking_project
    with pytest.raises(RuntimeError, match="does not match"):
        runner_for(repo, profile, molecule, postdock=FakePostDock(molecule, True)).run()
    with repo.connection() as conn: assert conn.execute("SELECT status FROM redocking_runs").fetchone()[0] == "failed"


def test_cancelled_vina_is_not_success(redocking_project) -> None:
    repo, profile, molecule = redocking_project
    with pytest.raises(InterruptedError): runner_for(repo, profile, molecule, vina=FakeVina(cancel=True)).run()
    with repo.connection() as conn: assert conn.execute("SELECT status FROM redocking_runs").fetchone()[0] == "interrupted"


def test_resume_reuses_completed_stages_and_box_change_restarts_docking(redocking_project) -> None:
    repo, profile, molecule = redocking_project; meeko, vina, export = FakeMeeko(), FakeVina(), FakePostDock(molecule)
    first = runner_for(repo, profile, molecule, meeko=meeko, vina=vina, postdock=export).run(); run_id = first["run_id"]
    runner_for(repo, profile, molecule, meeko=meeko, vina=vina, postdock=export).run(run_id)
    assert (meeko.calls, vina.calls, export.calls) == (1, 1, 1)
    changed = dict(profile); changed["size_x"] = 30
    runner_for(repo, changed, molecule, meeko=meeko, vina=vina, postdock=export).run(run_id)
    # Only docking is invalidated; identical mocked multi-pose output safely reuses downstream export.
    assert meeko.calls == 1 and vina.calls == 2 and export.calls == 1


def test_invalid_reusable_artifact_is_reprocessed(redocking_project) -> None:
    repo, profile, molecule = redocking_project; export = FakePostDock(molecule)
    result = runner_for(repo, profile, molecule, postdock=export).run(); root = Path(result["run_root"])
    (root / "poses/pose_001.sdf").write_text("broken", encoding="utf-8")
    runner_for(repo, profile, molecule, postdock=export).run(str(result["run_id"]))
    assert export.calls == 2 and read_sdf(root / "poses/pose_001.sdf") is not None


def test_restart_marks_running_redocking_interrupted(redocking_project) -> None:
    repo, profile, _ = redocking_project
    with repo.connection() as conn:
        conn.execute("""INSERT INTO redocking_runs(run_id,receptor_profile_id,status,reference_ligand_id,receptor_path,receptor_sha256,
            reference_sdf_path,reference_mol2_path,reference_sha256,settings_json,started_at) VALUES
            ('run','receptor-a','running','LIG A:401','r','h','s','m','h','{}','now')""")
    assert repo.recover_interrupted_runs() == 1
    with repo.connection() as conn: assert conn.execute("SELECT status FROM redocking_runs").fetchone()[0] == "interrupted"
