"""Exercise macOS discovery on every host without executing ADFR binaries."""
from pathlib import Path
from types import SimpleNamespace
import os

import pytest

from moldockpipe.receptors import ad4zn


@pytest.fixture
def mac_discovery(monkeypatch, tmp_path):
    monkeypatch.setattr(ad4zn.sys, "platform", "darwin")
    # Do not change os.name globally: pathlib must retain the test host's paths.
    monkeypatch.setattr(ad4zn, "os", SimpleNamespace(
        name="posix", environ=os.environ, access=os.access, X_OK=os.X_OK))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path / "home"))
    monkeypatch.setattr(ad4zn, "application_root", lambda: tmp_path / "application")
    monkeypatch.setattr(ad4zn.shutil, "which", lambda name: None)
    for name in ("PYTHON", "AUTOGRID"):
        monkeypatch.delenv("MOLDOCKPIPE_AD4ZN_" + name, raising=False)


def test_mac_search_locations(mac_discovery, tmp_path):
    roots = ad4zn._adfr_search_roots()
    home = tmp_path / "home"
    assert all(root in roots for root in (
        home, home / "Applications", Path("/Applications"), home / "Downloads",
        Path("/opt"), Path("/usr/local")))


@pytest.mark.parametrize("folder", [
    "home/ADFRsuite-1.0", "Applications/ADFRsuite-1.0",
    "home/Applications/ADFRsuite-1.0", "home/Downloads/ADFRsuite_x86_64Darwin_1.0",
    "opt/adfrsuite-1.1",
])
def test_mac_installation_discovers_both_tools(mac_discovery, monkeypatch, tmp_path, folder):
    installation = tmp_path / folder
    (installation / "bin").mkdir(parents=True)
    for name in ("pythonsh", "autogrid4"):
        tool = installation / "bin" / name
        tool.write_text("tool")
        tool.chmod(0o755)
    monkeypatch.setattr(ad4zn, "_adfr_search_roots", lambda: [installation.parent])

    environment = ad4zn.check_ad4zn_environment(tmp_path / "project")
    assert environment.components["python"] == (installation / "bin/pythonsh").resolve()
    assert environment.components["autogrid"] == (installation / "bin/autogrid4").resolve()
    monkeypatch.setenv("MOLDOCKPIPE_AD4ZN_PYTHON", str(tmp_path / "missing"))
    assert ad4zn.check_ad4zn_environment(tmp_path).components["python"] is None


def test_discovery_ignores_missing_roots_and_prefers_newer_versions(mac_discovery, monkeypatch, tmp_path):
    for name in ("ADFRsuite-1.9", "adfrsuite-1.10"):
        (tmp_path / name).mkdir()
    (tmp_path / "ADFRsuite-2.0").write_text("not a directory")
    monkeypatch.setattr(ad4zn, "_adfr_search_roots", lambda: [tmp_path / "missing", tmp_path, tmp_path])
    assert ad4zn.adfr_installations() == [tmp_path / "adfrsuite-1.10", tmp_path / "ADFRsuite-1.9"]


def test_mac_path_tools_are_still_discovered(mac_discovery, monkeypatch, tmp_path):
    monkeypatch.setattr(ad4zn, "_adfr_search_roots", lambda: [])
    tools = {}
    for name in ("pythonsh", "autogrid4"):
        tool = tmp_path / name
        tool.write_text("tool")
        tool.chmod(0o755)
        tools[name] = str(tool)
    monkeypatch.setattr(ad4zn.shutil, "which", tools.get)
    environment = ad4zn.check_ad4zn_environment(tmp_path)
    assert environment.components["python"] == (tmp_path / "pythonsh").resolve()
    assert environment.components["autogrid"] == (tmp_path / "autogrid4").resolve()
