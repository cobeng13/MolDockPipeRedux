import builtins
import subprocess

from moldockpipe import diagnostics
from moldockpipe.services.molscrub import MolScrubService


def test_failures_timeouts_and_native_crashes_do_not_stop_other_checks(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(diagnostics.platform, "platform", lambda: "Windows test")

    def run(command, **kwargs):
        name = command[3]
        calls.append(name)
        assert command[0] == diagnostics.sys.executable
        assert kwargs["timeout"] == diagnostics.TIMEOUT
        if name == "rdkit":
            raise subprocess.TimeoutExpired(command, 120, output=b"partial output", stderr=b"loading")
        if name == "molscrub":
            return subprocess.CompletedProcess(command, 1, "", "ImportError: missing dependency")
        if name == "meeko":
            return subprocess.CompletedProcess(command, -1073741819, "", "native crash")
        return subprocess.CompletedProcess(command, 0, "scored pose", "")

    monkeypatch.setattr(diagnostics.subprocess, "run", run)
    report = diagnostics.run_diagnostics(tmp_path)
    assert calls == list(diagnostics.CHECKS)
    assert [check["status"] for check in report["checks"]] == ["FAIL"] * 3 + ["PASS"] * 3
    rendered = diagnostics.format_report(report)
    assert "partial output" in rendered and "exceeded 120 seconds" in rendered
    assert "ImportError: missing dependency" in rendered and "native crash" in rendered
    assert diagnostics.sys.executable in rendered
    assert "pip check" in rendered


def test_molscrub_import_error_names_actual_missing_dependency(monkeypatch):
    original = builtins.__import__

    def broken_import(name, *args, **kwargs):
        if name == "molscrub":
            raise ImportError("cannot import name 'ExampleAPI' from 'rdkit.Chem'")
        # Keep this regression test runnable without a chemistry installation.
        if name == "rdkit":
            return type("RDKit", (), {"Chem": object()})
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", broken_import)
    import pytest
    with pytest.raises(RuntimeError, match="ExampleAPI") as error:
        MolScrubService()._generate("CCO")
    assert diagnostics.sys.executable in str(error.value)
    assert "Troubleshooter" in str(error.value)


def test_vina_probe_requires_scored_output(tmp_path, monkeypatch):
    from moldockpipe.services import vina

    monkeypatch.setattr(vina, "find_vina_executable", lambda root: tmp_path / "vina.exe")

    def run(command, **kwargs):
        # A successful exit alone does not establish working docking.
        return subprocess.CompletedProcess(command, 0, "help text", "")

    monkeypatch.setattr(diagnostics.subprocess, "run", run)
    import pytest
    with pytest.raises(RuntimeError, match="valid scored pose"):
        diagnostics._probe("vina", tmp_path, tmp_path)


def test_doctor_cli_saves_report_without_loading_qt(tmp_path, monkeypatch, capsys):
    from moldockpipe.__main__ import main

    report = {"checks": [{"name": "rdkit", "status": "FAIL"}]}
    monkeypatch.setattr(diagnostics, "run_diagnostics", lambda root: report)
    output = tmp_path / "report.json"
    monkeypatch.setattr("sys.argv", ["moldockpipe", "doctor", "--json", "--output", str(output)])
    assert main() == 1
    import json
    assert json.loads(output.read_text()) == report
    assert json.loads(capsys.readouterr().out) == report


def test_webengine_import_failure_is_reported_with_conda_guidance(tmp_path, monkeypatch):
    monkeypatch.setattr(diagnostics.platform, "platform", lambda: "Windows test")
    def run(command, **kwargs):
        if command[3] == "webengine":
            return subprocess.CompletedProcess(command, 1, "", "ImportError: QtWebEngineWidgets missing")
        return subprocess.CompletedProcess(command, 0, "passed", "")

    monkeypatch.setattr(diagnostics.subprocess, "run", run)
    report = diagnostics.run_diagnostics(tmp_path)
    check = next(check for check in report["checks"] if check["name"] == "webengine")
    assert check["status"] == "FAIL"
    assert "environment.yml" in check["suggestion"]
    assert "optional for headless" in check["suggestion"]
    assert all(check["status"] == "PASS" for check in report["checks"][:-1])


def test_webengine_probe_checks_assets_without_constructing_a_window(tmp_path, monkeypatch, capsys):
    from moldockpipe import __main__
    import pytest
    original = builtins.__import__

    class NoWindow:
        def __init__(self, *args, **kwargs):
            raise AssertionError("The import check must not create a graphics window")

    def fake_import(name, *args, **kwargs):
        if name == "PyQt6.QtWebEngineWidgets":
            return type("WebEngine", (), {"QWebEngineView": NoWindow})
        if name == "PyQt6.QtWebChannel":
            return type("WebChannel", (), {"QWebChannel": NoWindow})
        return original(name, *args, **kwargs)

    monkeypatch.setattr(__main__, "configure_qt_runtime", lambda: None)
    monkeypatch.setattr(builtins, "__import__", fake_import)
    monkeypatch.setattr(diagnostics, "__file__", str(tmp_path / "diagnostics.py"))
    assets = tmp_path / "ui/assets/molecular_viewer"
    assets.mkdir(parents=True)
    with pytest.raises(RuntimeError, match="asset missing"):
        diagnostics._probe("webengine", tmp_path, tmp_path)
    for filename in ("viewer.html", "viewer.js", "3Dmol-2.5.5.min.js"):
        (assets / filename).write_text("test asset")
    diagnostics._probe("webengine", tmp_path, tmp_path)
    assert "does not start a graphics window" in capsys.readouterr().out
