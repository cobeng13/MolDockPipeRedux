from pathlib import Path
import subprocess

import pytest

from moldockpipe.services.vinagpu import VinaGpuDockingBackend, find_vinagpu_executable, vinagpu_tool_hash


def inputs(tmp_path):
    executable = tmp_path / "engine" / "Vina-GPU+.exe"
    executable.parent.mkdir()
    executable.write_bytes(b"binary")
    for name in ("Kernel1_Opt.bin", "Kernel2_Opt.bin"):
        (executable.parent / name).write_bytes(b"kernel")
    atom = "ATOM      1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n"
    receptor, ligand = tmp_path / "receptor.pdbqt", tmp_path / "ligand.pdbqt"
    for path in (receptor, ligand):
        path.write_text(atom * 2, encoding="ascii")
    return dict(executable=executable, receptor=receptor, ligand=ligand,
                output=tmp_path / "run" / "output.pdbqt", log=tmp_path / "run" / "log.txt",
                settings={"center_x": 1, "center_y": 2, "center_z": 3, "size_x": 6, "size_y": 6,
                          "size_z": 6, "thread": 8000, "search_depth": 32, "seed": 42,
                          "num_modes": 9, "energy_range": 3, "exhaustiveness": 8, "cpu": 4})


def test_gpu_command_and_isolated_kernels(tmp_path, monkeypatch):
    kwargs = inputs(tmp_path)
    class Process:
        returncode = 0
        def __init__(self, command, **options):
            assert options["cwd"] == kwargs["output"].parent
            assert (options["cwd"] / "Kernel1_Opt.bin").read_bytes() == b"kernel"
            assert "--cpu" not in command and "--exhaustiveness" not in command and "--config" not in command
            assert command[command.index("--thread") + 1] == "8000"
            assert command[command.index("--search_depth") + 1] == "32"
            assert command[command.index("--seed") + 1] == "42"
            assert Path(command[command.index("--ligand") + 1]).is_absolute()
            kwargs["output"].write_text("REMARK VINA RESULT: -7.2 0.0 0.0\n" + kwargs["ligand"].read_text())
        def communicate(self, timeout=None):
            return "done", ""
    monkeypatch.setattr("moldockpipe.services.vinagpu.subprocess.Popen", Process)
    result = VinaGpuDockingBackend().run(**kwargs)
    assert result.poses == [(1, -7.2, 0.0, 0.0)]
    assert "Return code: 0" in kwargs["log"].read_text()
    assert result.command[0] == str(kwargs["executable"])


def test_missing_gpu_assets_fail_before_launch(tmp_path, monkeypatch):
    kwargs = inputs(tmp_path)
    monkeypatch.setenv("MOLDOCKPIPE_VINAGPU", str(kwargs["executable"]))
    assert find_vinagpu_executable(tmp_path) == kwargs["executable"]
    (kwargs["executable"].parent / "Kernel2_Opt.bin").unlink()
    with pytest.raises(FileNotFoundError, match="kernel missing"):
        VinaGpuDockingBackend().run(**kwargs)
    monkeypatch.setenv("MOLDOCKPIPE_VINAGPU", str(tmp_path / "missing.exe"))
    with pytest.raises(FileNotFoundError, match="does not point"):
        find_vinagpu_executable(tmp_path)


def test_gpu_cancellation_terminates_and_logs(tmp_path, monkeypatch):
    kwargs = inputs(tmp_path)
    checks = iter([False, True])
    class Process:
        returncode = -1
        stopped = False
        def __init__(self, *args, **options): pass
        def communicate(self, timeout=None):
            if not self.stopped:
                raise subprocess.TimeoutExpired("gpu", timeout)
            return "", "cancelled"
        def terminate(self): self.stopped = True
    monkeypatch.setattr("moldockpipe.services.vinagpu.subprocess.Popen", Process)
    with pytest.raises(InterruptedError, match="cancelled"):
        VinaGpuDockingBackend().run(**kwargs, cancelled=lambda: next(checks))
    assert "cancelled" in kwargs["log"].read_text()


def test_failed_gpu_output_is_not_accepted(tmp_path, monkeypatch):
    kwargs = inputs(tmp_path)
    class Process:
        returncode = 1
        def __init__(self, *args, **options):
            kwargs["output"].write_text("REMARK VINA RESULT: -7.2 0.0 0.0\n" + kwargs["ligand"].read_text())
        def communicate(self, timeout=None): return "", "OpenCL error"
    monkeypatch.setattr("moldockpipe.services.vinagpu.subprocess.Popen", Process)
    result = VinaGpuDockingBackend().run(**kwargs)
    assert not result.poses and result.return_code == 1


def test_gpu_hash_tracks_kernel_changes(tmp_path):
    kwargs = inputs(tmp_path)
    executable = kwargs["executable"]
    original = vinagpu_tool_hash(executable)
    (executable.parent / "Kernel1_Opt.bin").write_bytes(b"new kernel")
    assert vinagpu_tool_hash(executable) != original


def test_gpu_discovery_uses_dedicated_folder(tmp_path, monkeypatch):
    monkeypatch.delenv("MOLDOCKPIPE_VINAGPU", raising=False)
    app = tmp_path / "app"
    project = tmp_path / "project"
    monkeypatch.setattr("moldockpipe.services.vinagpu.application_root", lambda: app)
    wrong_folder = app / "tools" / "vina"
    wrong_folder.mkdir(parents=True)
    (wrong_folder / "Vina-GPU+.exe").write_bytes(b"wrong location")
    monkeypatch.setattr("moldockpipe.services.vinagpu.shutil.which", lambda _: str(wrong_folder / "Vina-GPU+.exe"))
    with pytest.raises(FileNotFoundError, match="tools/vinagpu"):
        find_vinagpu_executable(project)
    local = project / "tools" / "vinagpu" / "Vina-GPU+.exe"
    local.parent.mkdir(parents=True)
    local.write_bytes(b"project binary")
    local.chmod(0o755)
    assert find_vinagpu_executable(project) == local.resolve()
    shared = app / "tools" / "vinagpu" / "Vina-GPU+.exe"
    shared.parent.mkdir(parents=True)
    shared.write_bytes(b"application binary")
    shared.chmod(0o755)
    assert find_vinagpu_executable(project) == shared.resolve()
