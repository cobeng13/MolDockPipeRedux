from __future__ import annotations

import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from moldockpipe.project import ProjectRepository
from moldockpipe.receptors.ligand_chemistry import chemistry_summary, pose_in_reference_order, read_sdf, validate_same_heavy_graph, write_mol2, write_sdf
from moldockpipe.redocking.models import RedockingSettings
from moldockpipe.redocking import rmsd
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


def test_graph_validation_and_mol2_mapping_ignore_atom_order_and_stereo(redocking_project, tmp_path: Path) -> None:
    _, _, reference = redocking_project
    reordered = Chem.RenumberAtoms(reference, list(reversed(range(reference.GetNumAtoms()))))
    Chem.AssignStereochemistry(reordered, cleanIt=True, force=True)
    validate_same_heavy_graph(reference, reordered)
    aligned = pose_in_reference_order(reference, reordered)
    assert [atom.GetSymbol() for atom in aligned.GetAtoms()] == [atom.GetSymbol() for atom in Chem.RemoveHs(reference).GetAtoms()]
    path = tmp_path / "aligned_pose.mol2"; write_mol2(aligned, path)
    atom_lines = path.read_text(encoding="utf-8").split("@<TRIPOS>ATOM", 1)[1].split("@<TRIPOS>BOND", 1)[0]
    assert "MOLDOCK_REFERENCE_ATOM_INDEX" not in atom_lines  # metadata stays in-memory; labels are valid MOL2 fields


def test_pose_affinity_parsing(tmp_path: Path) -> None:
    output = tmp_path / "poses.pdbqt"; output.write_text("REMARK VINA RESULT: -9.100 0.000 0.000\nREMARK VINA RESULT: -8.200 1.000 2.000\n")
    assert parse_vina_poses(output) == [(1, -9.1, 0.0, 0.0), (2, -8.2, 1.0, 2.0)]


def test_successful_redocking_generates_all_artifacts(redocking_project, monkeypatch) -> None:
    monkeypatch.setattr(rmsd, "_engine", lambda: (
        lambda reference, pose: 2.5 if "001" in Path(pose).name else 0.8, "test-version"))
    repo, profile, molecule = redocking_project; result = runner_for(repo, profile, molecule).run()
    assert result["status"] == "ARTIFACTS_READY" and result["poses"] == 2 and result["top_affinity"] == -8.2
    root = Path(result["run_root"]); assert (root / "dockrmsd/reference_ligand_heavy.mol2").is_file()
    assert (root / "dockrmsd/pose_001_heavy.mol2").is_file()
    assert (root / "dockrmsd/top_ranked_pose_heavy.mol2").read_bytes() == (root / "dockrmsd/pose_001_heavy.mol2").read_bytes()
    assert "DockRMSD" in (root / "dockrmsd/README.txt").read_text()
    transformation = json.loads((root / "transformation.json").read_text(encoding="utf-8"))
    assert transformation["dockrmsd"]["results"] == "dockrmsd/results.json"
    saved = rmsd.read_results(root)
    assert saved["rank_1"]["rmsd_angstrom"] == 2.5
    assert saved["lowest"]["pose_rank"] == 2
    assert saved["lowest"]["rmsd_angstrom"] == 0.8
    assert rmsd.validation_runs(repo, profile)[0][2:] == (saved, True)
    assert transformation["operations"][-1]["atom_mapping"] == "pose graph mapped to reference atom order"
    with repo.connection() as conn:
        assert conn.execute("SELECT status FROM redocking_runs").fetchone()[0] == "ARTIFACTS_READY"
        assert conn.execute("SELECT COUNT(*) FROM redocking_poses").fetchone()[0] == 2


def test_rmsd_retries_missing_pose_and_reuses_results(redocking_project, monkeypatch):
    repo, profile, molecule = redocking_project
    meeko, vina = FakeMeeko(), FakeVina()
    monkeypatch.setattr(rmsd, "_engine", lambda: (lambda *_: 0.0, "test-version"))
    runner_for(repo, profile, molecule, meeko=meeko, vina=vina).run()
    run, provenance, _, _ = rmsd.validation_runs(repo, profile)[0]
    pose = run / "dockrmsd" / "pose_002_heavy.mol2"
    original = pose.read_bytes()
    pose.unlink()
    assert rmsd.read_results(run) is None
    partial = rmsd.calculate_results(run, provenance)
    assert partial["status"] == "Incomplete"
    assert pose.name in partial["poses"][1]["error"]
    pose.write_bytes(original)
    assert rmsd.calculate_results(run, provenance)["status"] == "Calculated"

    def unexpected(*_):
        pytest.fail("Matching results should be reused")

    monkeypatch.setattr(rmsd, "_engine", lambda: (unexpected, "test-version"))
    assert rmsd.calculate_results(run, provenance)["status"] == "Calculated"
    assert meeko.calls == vina.calls == 1


def test_rmsd_engine_failure_version_and_changed_inputs(redocking_project, monkeypatch):
    repo, profile, molecule = redocking_project

    def unavailable():
        raise ImportError("test missing engine")

    monkeypatch.setattr(rmsd, "_engine", unavailable)
    completed = runner_for(repo, profile, molecule).run()
    assert completed["status"] == "ARTIFACTS_READY"
    assert completed["dockrmsd"]["status"] == "Unavailable"
    run, provenance, _, _ = rmsd.validation_runs(repo, profile)[0]
    monkeypatch.setattr(rmsd, "_engine", lambda: (lambda *_: 0.25, "v1"))
    result = rmsd.calculate_results(run, provenance)
    assert result["status"] == "Calculated" and result["lowest"]["pose_rank"] == 1
    monkeypatch.setattr(rmsd, "_engine", lambda: (lambda *_: 1.25, "v2"))
    assert rmsd.calculate_results(run, provenance)["lowest"]["rmsd_angstrom"] == 1.25
    reference = run / "dockrmsd" / rmsd.REFERENCE
    reference.write_text(reference.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    assert rmsd.read_results(run) is None
    invalid = rmsd.calculate_results(run, provenance)
    assert invalid["status"] == "Unavailable"
    assert "Checksum" in invalid["poses"][0]["error"]


def test_real_pydockrmsd_on_generated_mol2(tmp_path):
    try:
        version("pydockrmsd")
    except PackageNotFoundError:
        pytest.skip("PyDockRMSD wheel is not installed")
    molecule = Chem.MolFromSmiles("CC")
    conformer = Chem.Conformer(2)
    conformer.SetAtomPosition(0, (0, 0, 0))
    conformer.SetAtomPosition(1, (1.5, 0, 0))
    molecule.AddConformer(conformer)
    reference = tmp_path / "reference.mol2"
    write_mol2(molecule, reference, heavy_only=True)
    assert rmsd._measure(str(reference), str(reference)) == pytest.approx(0)

    translated = Chem.Mol(molecule)
    for index in range(2):
        translated.GetConformer().SetAtomPosition(index, (index * 1.5, 0, 3))
    translated_path = tmp_path / "translated.mol2"
    write_mol2(translated, translated_path, heavy_only=True)
    assert rmsd._measure(str(reference), str(translated_path)) == pytest.approx(3)

    swapped = Chem.RenumberAtoms(molecule, [1, 0])
    swapped_path = tmp_path / "swapped.mol2"
    write_mol2(swapped, swapped_path, heavy_only=True)
    assert rmsd._measure(str(reference), str(swapped_path)) == pytest.approx(0)


def test_rmsd_package_and_report_round_trip(redocking_project, tmp_path, monkeypatch):
    pytest.importorskip("PyQt6.QtWidgets")
    from PyQt6.QtWidgets import QApplication, QFileDialog, QLabel, QMessageBox
    from moldockpipe.reporting import generate_project_report
    from moldockpipe.ui.receptor_manager import ReceptorManagerDialog
    from moldockpipe.ui.redocking_dialog import RedockingResultDialog

    app = QApplication.instance() or QApplication([])
    repo, profile, molecule = redocking_project
    monkeypatch.setattr(rmsd, "_engine", lambda: (
        lambda reference, pose: 2.5 if "001" in Path(pose).name else 0.8, "test-version"))
    completed = runner_for(repo, profile, molecule).run()
    result_dialog = RedockingResultDialog(repo, completed, profile, 123)
    assert any("Rank 1 2.500 Å" in label.text() for label in result_dialog.findChildren(QLabel))
    package = tmp_path / "receptor.moldock-receptor.zip"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(package), ""))
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(package), ""))
    monkeypatch.setattr(QMessageBox, "information", lambda *a: None)
    source = ReceptorManagerDialog(repo)
    source.table.selectRow(0)
    assert "Rank 1 2.500 Å" in source.rmsd_summary.text()
    source._export_profile()

    destination = ProjectRepository.create(tmp_path / "imported")
    target = ReceptorManagerDialog(destination)
    target._import_profile()
    target._save_locked()
    imported = next(item for item in destination.get_receptor_profiles() if item["name"] == profile["name"])
    run, _, result, current = rmsd.validation_runs(destination, imported)[0]
    assert current and result == completed["dockrmsd"]
    with destination.connection() as conn:
        assert conn.execute("SELECT COUNT(*) FROM redocking_runs").fetchone()[0] == 0
    report = generate_project_report(destination, output_root=destination.root / "exports" / "rmsd")
    data = json.loads(report.with_name("report_data.json").read_text(encoding="utf-8"))
    receptor = next(item for item in data["receptors"] if item["profile"]["id"] == imported["id"])
    assert receptor["redocking"]["dockrmsd"] == result
    assert "Lowest 0.800 Å (pose 2)" in report.read_text(encoding="utf-8")
    (destination.root / imported["receptor"]).write_text("different receptor", encoding="utf-8")
    assert rmsd.validation_runs(destination, imported)[0][3] is False
    (run / "dockrmsd" / "results.json").unlink()
    assert rmsd.validation_runs(destination, imported)[0][2] is None
    app.processEvents()


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


def test_ad4zn_redocking_uses_maps_and_invalidates_reuse(redocking_project, monkeypatch):
    from moldockpipe.receptors.ad4zn import AD4ZnEnvironment
    from test_ad4zn import pdbqt
    repo, profile, molecule = redocking_project
    profile['protocol'] = 'ad4zn'
    (repo.root / profile['receptor']).write_text(pdbqt('C', 'Zn'))
    monkeypatch.setattr('moldockpipe.redocking.runner.check_ad4zn_environment', lambda _: AD4ZnEnvironment({}))
    context = {'protocol': 'ad4zn', 'map_fingerprint': 'maps-a'}
    def ensure(repository, selected, ligands, executable, cancelled):
        assert selected['protocol'] == 'ad4zn'
        assert len(ligands) == 1 and ligands[0].is_file()
        return repo.root / 'maps/receptor_TZ', dict(context)
    monkeypatch.setattr('moldockpipe.redocking.runner.ensure_ad4zn_maps', ensure)
    class ZincVina(FakeVina):
        def run(self, **kwargs):
            assert kwargs['settings']['protocol'] == 'ad4zn'
            assert kwargs['maps'] == repo.root / 'maps/receptor_TZ'
            return super().run(**kwargs)
    vina = ZincVina()
    runner = runner_for(repo, profile, molecule, vina=vina)
    result = runner.run()
    assert vina.calls == 1
    runner.run(result['run_id'])
    assert vina.calls == 1
    context['map_fingerprint'] = 'maps-b'
    runner.run(result['run_id'])
    assert vina.calls == 2
    with repo.connection() as conn:
        row = conn.execute('SELECT settings_json FROM redocking_runs WHERE run_id=?', (result['run_id'],)).fetchone()
    assert json.loads(row[0])['ad4zn']['map_fingerprint'] == 'maps-b'
