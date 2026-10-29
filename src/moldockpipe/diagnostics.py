"""Bounded, isolated functional checks using the application's Python environment."""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import traceback
from datetime import datetime, timezone
from importlib import metadata, util
from pathlib import Path

CHECKS = ("rdkit", "molscrub", "meeko", "interactions", "vina", "vinagpu", "webengine")
TIMEOUT = 120


def _run(command: list[str], *, cwd: Path | None = None, timeout: int = 90) -> None:
    print("Command: " + subprocess.list2cmdline(command), flush=True)
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=timeout,
                            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    print(result.stdout)
    print(result.stderr, file=sys.stderr)
    if result.returncode:
        raise RuntimeError(f"Process exited with code {result.returncode}")


def _ethanol():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    if AllChem.EmbedMolecule(molecule, randomSeed=42) != 0:
        raise RuntimeError("RDKit could not generate ethanol coordinates")
    return molecule


def _probe(name: str, root: Path, work: Path) -> None:
    if name == "rdkit":
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from rdkit.Chem.Draw import rdMolDraw2D

        molecule = _ethanol()
        if not Chem.MolToInchiKey(molecule) or not 46 < Descriptors.MolWt(molecule) < 47:
            raise RuntimeError("RDKit molecule descriptors/InChI check failed")
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
        drawer.DrawMolecule(Chem.RemoveHs(molecule))
        drawer.FinishDrawing()
        if "<svg" not in drawer.GetDrawingText():
            raise RuntimeError("RDKit SVG drawing check failed")
        print("SMILES parsing, hydrogen addition, 3D embedding, descriptors, InChI and drawing passed.")
    elif name == "interactions":
        import MDAnalysis
        import gemmi
        import pandas
        import prolif as plf

        # Synthetic contact smoke test; this does not validate a study's results.
        probe = _ethanol()
        ligand = plf.Molecule.from_rdkit(probe, resname="LIG", resnumber=1)
        conformer = probe.GetConformer()
        for index in range(probe.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            conformer.SetAtomPosition(index, (point.x, point.y, point.z + 3.0))
        partner = plf.Molecule.from_rdkit(probe, resname="PRB", resnumber=2)
        fingerprint = plf.Fingerprint(["VdWContact"])
        if not fingerprint.generate(ligand, partner):
            raise RuntimeError("ProLIF synthetic contact check produced no interactions")
        print("ProLIF imports and synthetic contact detection passed (including RDKit drawing dependencies).")
    elif name == "molscrub":
        from .services.molscrub import MolScrubService

        states, _ = MolScrubService().generate("CCO", timeout_seconds=60)
        if not states or any(state.molecule.GetNumConformers() == 0 for state in states):
            raise RuntimeError("MolScrub returned no states with 3D coordinates")
        print(f"MolScrub worker generated {len(states)} ethanol state(s) with 3D coordinates.")
    elif name == "meeko":
        from rdkit import Chem
        from .services.validation import validate_pdbqt

        sdf, pdbqt = work / "ethanol.sdf", work / "ethanol.pdbqt"
        with Chem.SDWriter(str(sdf)) as writer:
            writer.write(_ethanol())
        _run([sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", str(sdf), "-o", str(pdbqt)], timeout=45)
        valid, reason = validate_pdbqt(pdbqt)
        if not valid:
            raise RuntimeError(reason)
        # Receptor preparation is a separate entry point used by the application.
        _run([sys.executable, "-m", "meeko.cli.mk_prepare_receptor", "--help"], timeout=20)
        _run([sys.executable, "-m", "meeko.cli.mk_export", "--help"], timeout=20)
        print("Ethanol PDBQT preparation passed; receptor and export commands load successfully.")
    elif name == "webengine":
        from .__main__ import configure_qt_runtime
        configure_qt_runtime()
        from PyQt6.QtWebEngineWidgets import QWebEngineView
        from PyQt6.QtWebChannel import QWebChannel
        assets = Path(__file__).parent / "ui" / "assets" / "molecular_viewer"
        for filename in ("viewer.html", "viewer.js", "3Dmol-2.5.5.min.js"):
            path = assets / filename
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(f"Offline viewer asset missing or empty: {path}")
        print("PyQt6 WebEngine/WebChannel bindings and offline viewer assets load successfully.")
        print("This check does not start a graphics window. Verify WebGL and PNG export in Results.")
    else:
        from .services.vina import find_vina_executable
        from .services.vinagpu import KERNEL_NAMES, find_vinagpu_executable, vinagpu_tool_hash
        from .services.validation import parse_vina_poses, validate_pdbqt

        gpu = name == "vinagpu"
        executable = find_vinagpu_executable(root) if gpu else find_vina_executable(root)
        print(f"Executable: {executable}", flush=True)
        if gpu:
            vinagpu_tool_hash(executable)
            for kernel in KERNEL_NAMES:
                shutil.copy2(executable.parent / kernel, work / kernel)
        # A rigid carbon probe avoids coupling engine checks to Python chemistry packages.
        atom = "ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n"
        receptor, ligand, output = (work / filename for filename in
                                     ("receptor.pdbqt", "ligand.pdbqt", "output.pdbqt"))
        receptor.write_text(atom, encoding="ascii")
        ligand.write_text("ROOT\n" + atom + "ENDROOT\nTORSDOF 0\n", encoding="ascii")
        command = [str(executable), "--receptor", str(receptor), "--ligand", str(ligand),
                   "--out", str(output), "--center_x", "0", "--center_y", "0", "--center_z", "0",
                   "--size_x", "10", "--size_y", "10", "--size_z", "10", "--num_modes", "1", "--seed", "42"]
        command += ["--thread", "1000", "--search_depth", "1"] if gpu else ["--cpu", "1", "--exhaustiveness", "1"]
        _run(command, cwd=work)
        valid, reason = validate_pdbqt(output)
        if not valid or not parse_vina_poses(output):
            raise RuntimeError(f"Docking did not produce a valid scored pose: {reason}")
        print("Small docking job produced a valid scored pose.")


def _package_info(name: str, module_name: str | None = None) -> dict[str, str | None]:
    info: dict[str, str | None] = {"version": None, "location": None}
    try:
        info["version"] = metadata.version(name)
    except metadata.PackageNotFoundError:
        pass
    try:
        spec = util.find_spec(module_name or name)
        info["location"] = spec.origin if spec else None
    except Exception as exc:
        info["location"] = str(exc)
    return info


def _advice(name: str, root: Path) -> str:
    python = subprocess.list2cmdline([sys.executable])
    if name == "webengine":
        return ("For the desktop viewer, create the Conda environment from environment.yml in the source checkout "
                "and launch the app from that environment. Keep Qt and WebEngine from the same distribution; "
                "do not replace Conda's Qt runtime with pip Qt wheels. Check README.md for first-time setup. "
                "WebEngine is optional for headless pipeline execution.")
    if name in {"vina", "vinagpu"}:
        return ("Check the executable path and runtime libraries. " +
                ("Check both kernel files beside Vina-GPU and the GPU vendor's OpenCL driver. "
                 "Vina-GPU is optional when using CPU docking." if name == "vinagpu" else
                 "Set MOLDOCKPIPE_VINA or place Vina in tools/vina."))
    if name in {"rdkit", "interactions"}:
        return ("Create a new environment from the current environment.yml and restart the app using that environment. "
                "Keep RDKit, Cairo and Qt together on conda-forge; mixing defaults and conda-forge can cause "
                "rdMolDraw2D DLL failures. Do not rename DLLs or overlay pip RDKit/Qt wheels. "
                f"Inspect Python dependencies with: {python} -m pip check")
    install = subprocess.list2cmdline([sys.executable, "-m", "pip", "install", "-e", str(root) + "[chemistry,dev]"])
    return (f"Use this app's Python to inspect dependencies: {python} -m pip check\n"
            f"From the source checkout, reinstall with: {install}\n"
            "Install RDKit into the same Conda environment as the app. Read the traceback for "
            "missing dependencies or DLL errors; an installed package can still fail to import. "
            "Restart the app after fixing the environment.")


def run_diagnostics(root: Path) -> dict:
    root = root.resolve()
    report = {"generated_at": datetime.now(timezone.utc).isoformat(),
              "python": sys.executable, "python_version": sys.version, "prefix": sys.prefix,
              "platform": platform.platform(), "project_root": str(root),
              "conda_prefix": os.environ.get("CONDA_PREFIX"),
              "virtual_env": os.environ.get("VIRTUAL_ENV"),
              "packages": {name: _package_info(name) for name in
                           ("rdkit", "molscrub", "meeko", "scipy", "gemmi", "prolif", "MDAnalysis", "pandas")},
              "checks": []}
    report["packages"]["PyQt6"] = _package_info("PyQt6")
    report["packages"]["PyQt6-WebEngine"] = _package_info("PyQt6-WebEngine", "PyQt6.QtWebEngineWidgets")
    for name in CHECKS:
        command = [sys.executable, "-m", "moldockpipe.diagnostics", name, str(root)]
        try:
            result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace",
                                    timeout=TIMEOUT,
                                    creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
            check = {"name": name, "status": "PASS" if result.returncode == 0 else "FAIL",
                     "return_code": result.returncode, "stdout": result.stdout, "stderr": result.stderr}
        except subprocess.TimeoutExpired as exc:
            def decoded(value):
                return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else value or ""
            check = {"name": name, "status": "FAIL", "stdout": decoded(exc.stdout),
                     "stderr": decoded(exc.stderr) + f"\nCheck exceeded {TIMEOUT} seconds."}
        except OSError as exc:
            check = {"name": name, "status": "FAIL", "stdout": "", "stderr": str(exc)}
        if check["status"] != "PASS":
            check["suggestion"] = _advice(name, root)
        report["checks"].append(check)
    return report


def format_report(report: dict) -> str:
    lines = ["MolDockPipe Troubleshooter", "",
             *(f"[{check['status']}] {check['name']}" for check in report['checks']), "",
             "Vina GPU is optional when using CPU docking. WebEngine is optional for headless execution.",
             "WebEngine checks bindings/assets; verify WebGL and PNG export in Results.", "", f"Python: {report['python']}",
             f"Version: {report['python_version']}", f"Environment: {report['prefix']}",
             f"Platform: {report['platform']}", f"Project/tool search root: {report['project_root']}",
             f"Conda: {report['conda_prefix']}", f"Virtual environment: {report['virtual_env']}",
             f"Generated (UTC): {report['generated_at']}", "", "Installed packages:"]
    lines.extend(f"  {name}: {info['version'] or 'unknown/not installed'} ({info['location']})"
                 for name, info in report["packages"].items())
    for check in report["checks"]:
        lines.extend(["", f"[{check['status']}] {check['name']}", check["stdout"].strip(), check["stderr"].strip()])
        if "suggestion" in check:
            lines.append("Suggested action: " + check["suggestion"])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    try:
        with tempfile.TemporaryDirectory(prefix="moldockpipe-doctor-") as directory:
            _probe(sys.argv[1], Path(sys.argv[2]), Path(directory))
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
