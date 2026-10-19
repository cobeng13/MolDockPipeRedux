from copy import deepcopy
from pathlib import Path

import pytest
from moldockpipe.configuration import read_project_config
from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


@pytest.fixture
def receptor_project(tmp_path):
    repo = ProjectRepository.create(tmp_path / "project")
    profile = {"protocol": "vina", "id": "original", "name": "Original receptor", "enabled": True,
               "archived": False, "receptor": "inputs/receptors/original/receptor.pdbqt",
               "center_x": 1.25, "size_x": 20, "reference_ligand": {"name": "LIG"}}
    receptor = repo.root / profile["receptor"]
    receptor.parent.mkdir(parents=True)
    receptor.write_text("ORIGINAL RECEPTOR", encoding="ascii")
    repo.save_receptor_profiles([profile])
    return repo, profile


def test_repeated_create_cannot_reset_existing_receptors(receptor_project):
    repo, profile = receptor_project
    before = (repo.root / "project.yml").read_bytes()
    with pytest.raises(FileExistsError, match="Open it instead"):
        ProjectRepository.create(repo.root)
    assert (repo.root / "project.yml").read_bytes() == before
    assert repo.get_receptor_profiles() == [profile]
    assert (repo.root / profile["receptor"]).read_text() == "ORIGINAL RECEPTOR"


def test_settings_updates_preserve_profiles_and_unedited_settings(receptor_project):
    repo, profile = receptor_project
    before = repo.get_settings()
    repo.update_settings({"molscrub": {"workers": 2}, "meeko": {"workers": 0},
                          "postdock": {"selected_parents": ["lig1"]}})
    assert repo.get_receptor_profiles() == [profile]
    assert repo.get_settings()["molscrub"]["fragment_policy"] == before["molscrub"]["fragment_policy"]
    assert read_project_config(repo.root / "project.yml.bak") == before
    with pytest.raises(ValueError, match="Receptor profiles"):
        repo.update_settings({"vina": {"profiles": []}})
    assert repo.get_receptor_profiles() == [profile]


def test_stale_profile_snapshot_cannot_remove_newer_profiles(receptor_project):
    repo, profile = receptor_project
    snapshot = repo.get_receptor_profiles(include_archived=True)
    second = {**profile, "id": "second", "name": "Second"}
    repo.save_receptor_profiles([profile, second])
    with pytest.raises(RuntimeError, match="changed since"):
        repo.save_receptor_profiles([], expected_profiles=snapshot)
    assert repo.get_receptor_profiles() == [profile, second]
    with pytest.raises(RuntimeError, match="changed since"):
        repo.save_receptor_profiles([{**profile, "enabled": False}], expected_profiles=snapshot)
    assert repo.get_receptor_profiles() == [profile, second]


def test_failed_atomic_replace_preserves_configuration_and_backup(receptor_project, monkeypatch):
    import moldockpipe.configuration as configuration

    repo, profile = receptor_project
    path = repo.root / "project.yml"
    before = path.read_bytes()
    replace = configuration.os.replace

    def fail_configuration_write(source, destination):
        if Path(destination) == path:
            # The live file still contains the full prior configuration while
            # its replacement is being prepared. It was never truncated.
            assert repo.get_receptor_profiles() == [profile]
            raise OSError("simulated interrupted configuration replacement")
        return replace(source, destination)

    monkeypatch.setattr(configuration.os, "replace", fail_configuration_write)
    with pytest.raises(OSError, match="interrupted"):
        repo.save_receptor_profiles([])
    assert path.read_bytes() == before
    assert (repo.root / "project.yml.bak").read_bytes() == before
    assert list(repo.root.glob(".project.yml-*.tmp")) == []


@pytest.mark.parametrize("content", ["", "null\n", "- bad\n", "vina: [unterminated"])
def test_invalid_configuration_is_not_treated_as_missing_receptors(receptor_project, content):
    repo, _ = receptor_project
    backup = (repo.root / "project.yml.bak").read_bytes()
    (repo.root / "project.yml").write_text(content, encoding="utf-8")
    with pytest.raises(ValueError, match="configuration"):
        repo.get_receptor_profiles()
    with pytest.raises(ValueError, match="configuration"):
        repo.migrate()
    assert (repo.root / "project.yml.bak").read_bytes() == backup


def test_deliberately_disabled_receptors_are_not_automatically_reenabled(receptor_project):
    repo, profile = receptor_project
    disabled = deepcopy(profile)
    disabled["enabled"] = False
    repo.save_receptor_profiles([disabled])
    repo.migrate()
    assert repo.get_receptor_profiles() == [disabled]
    with pytest.raises(RuntimeError, match="Open Receptors"):
        PipelineRunner(repo).validate_docking_configuration()
