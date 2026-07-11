import sqlite3

import pytest

from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


def test_repeated_chemistry_stages_reuse_artifacts(tmp_path) -> None:
    pytest.importorskip("rdkit")
    pytest.importorskip("molscrub")
    pytest.importorskip("meeko")
    repo = ProjectRepository.create(tmp_path / "project")
    csv_path = tmp_path / "input.csv"
    csv_path.write_text("id,smiles\nlig1,CCO\n", encoding="utf-8")
    repo.import_ligands_csv(csv_path)
    runner = PipelineRunner(repo)
    runner.run_screening()
    runner.run_molscrub()
    runner.run_meeko()
    runner.run_screening()
    runner.run_molscrub()
    runner.run_meeko()
    with sqlite3.connect(repo.database_path) as conn:
        assert conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM molecular_states").fetchone()[0] == 1
        assert conn.execute("SELECT COUNT(*) FROM artifacts").fetchone()[0] == 3
