from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from moldockpipe import __main__
from moldockpipe.services import executables


def _make_executable(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("tool", encoding="ascii")
    if os.name != "nt":
        path.chmod(0o755)
    return path


def test_executable_override_is_platform_neutral(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    tool = _make_executable(tmp_path / executables.VINA.platform_names[0])
    monkeypatch.setenv("MOLDOCKPIPE_VINA", str(tool))

    assert executables.find_executable(executables.VINA, tmp_path / "project") == tool.resolve()


def test_bundled_executable_is_found_using_platform_name(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    application = tmp_path / "application"
    project = tmp_path / "project"
    tool = _make_executable(application / "tools" / "vina" / executables.VINA_SPLIT.platform_names[0])
    monkeypatch.delenv("MOLDOCKPIPE_VINA_SPLIT", raising=False)
    monkeypatch.setattr(executables, "application_root", lambda: application)
    monkeypatch.setattr(executables.shutil, "which", lambda _name: None)

    assert executables.find_executable(executables.VINA_SPLIT, project) == tool.resolve()


def test_invalid_explicit_override_does_not_silently_fall_back(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("MOLDOCKPIPE_VINA", str(tmp_path / "missing-vina"))

    with pytest.raises(FileNotFoundError, match="MOLDOCKPIPE_VINA"):
        executables.find_executable(executables.VINA, tmp_path / "project")


def test_headless_run_does_not_require_qt(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]) -> None:
    project = tmp_path / "project"
    project.mkdir()
    (project / "project.sqlite").touch()
    calls: list[str] = []

    class FakeRepository:
        def __init__(self, root: Path):
            self.root = root

        def recover_interrupted_runs(self) -> int:
            calls.append("recover")
            return 0

    class FakeRunner:
        def __init__(self, _repository: FakeRepository):
            pass

        def run_screening(self) -> tuple[int, int]:
            calls.append("screening")
            return 2, 0

        def run_molscrub(self) -> tuple[int, int]:
            raise AssertionError("unselected stage ran")

        run_meeko = run_molscrub
        run_vina = run_molscrub
        run_postdock = run_molscrub

    monkeypatch.setattr(__main__, "ProjectRepository", FakeRepository)
    monkeypatch.setattr("moldockpipe.pipeline.PipelineRunner", FakeRunner)
    monkeypatch.setattr("sys.argv", ["moldockpipe", "run", str(project), "--stages", "screening"])

    assert __main__.main() == 0
    assert calls == ["recover", "screening"]
    assert json.loads(capsys.readouterr().out) == {"screening": [2, 0]}


@pytest.mark.parametrize("system,machine,suffix", [
    ("Windows", "AMD64", "win.exe"),
    ("Linux", "x86_64", "linux_x86_64"),
    ("Darwin", "x86_64", "mac_x86_64"),
    ("Darwin", "arm64", "mac_aarch64"),
    ("Darwin", "aarch64", "mac_aarch64"),
])
@pytest.mark.parametrize("spec", [executables.VINA, executables.VINA_SPLIT])
def test_mixed_os_bundle_selects_matching_binary(tmp_path, monkeypatch, system, machine, suffix, spec):
    application = tmp_path / "application"
    prefix = spec.posix_names[0]
    names = [f"{prefix}_1.2.7_{value}" for value in
             ("linux_x86_64", "mac_x86_64", "mac_aarch64", "win.exe")]
    names.append("vina_split.exe")
    for name in names:
        _make_executable(application / "tools" / "vina" / name)
    monkeypatch.setattr(executables.platform, "system", lambda: system)
    monkeypatch.setattr(executables.platform, "machine", lambda: machine)
    monkeypatch.setattr(executables, "application_root", lambda: application)
    monkeypatch.setattr(executables.shutil, "which", lambda _: None)
    monkeypatch.delenv(spec.environment_variable, raising=False)
    expected = "vina_split.exe" if system == "Windows" and spec == executables.VINA_SPLIT else f"{prefix}_1.2.7_{suffix}"
    assert executables.find_executable(spec, tmp_path / "project").name == expected


@pytest.mark.parametrize("spec", [executables.VINA, executables.VINA_SPLIT])
def test_linux_arm_does_not_select_other_architectures(tmp_path, monkeypatch, spec):
    application = tmp_path / "application"
    for suffix in ("linux_x86_64", "mac_aarch64", "mac_x86_64"):
        _make_executable(application / "tools" / "vina" / f"{spec.posix_names[0]}_1.2.7_{suffix}")
    monkeypatch.setattr(executables.platform, "system", lambda: "Linux")
    monkeypatch.setattr(executables.platform, "machine", lambda: "aarch64")
    monkeypatch.setattr(executables, "application_root", lambda: application)
    monkeypatch.setattr(executables.shutil, "which", lambda _: None)
    monkeypatch.delenv(spec.environment_variable, raising=False)
    with pytest.raises(FileNotFoundError):
        executables.find_executable(spec, tmp_path / "project")
    generic = _make_executable(application / "tools" / "vina" / spec.posix_names[0])
    assert executables.find_executable(spec, tmp_path / "project") == generic.resolve()
