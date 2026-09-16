"""Optional AutoDock4Zn preparation; generated bundles are immutable and auditable."""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from ..fingerprints import file_sha256, fingerprint
from ..services.executables import application_root


def protocol_for(profile: dict) -> str:
    protocol = str(profile.get("protocol", "vina"))
    if protocol not in {"vina", "ad4zn"}:
        raise ValueError(f"Unsupported docking protocol: {protocol}")
    return protocol


@dataclass(frozen=True)
class AD4ZnEnvironment:
    components: dict[str, Path | None]

    @property
    def missing(self) -> list[str]:
        return [name for name, path in self.components.items() if path is None]

    @property
    def problems(self) -> list[str]:
        labels = {"python": "ADFR Python interpreter (install ADFR Suite)",
                  "autogrid": "AutoGrid4 (included with ADFR Suite)",
                  "zinc_pseudo": "zinc_pseudo.py", "prepare_gpf": "prepare_gpf4zn.py",
                  "parameters": "AD4Zn.dat"}
        problems = ["Missing: " + labels.get(name, name) for name in self.missing]
        parameter = self.components.get("parameters")
        if parameter:
            try:
                if not any(line.split()[:2] == ["atom_par", "TZ"] for line in parameter.read_text().splitlines()):
                    problems.append("Invalid AD4Zn.dat: download the actual parameter file, not the small symbolic-link placeholder. See docs/AD4Zn.md.")
            except (OSError, UnicodeError) as exc:
                problems.append(f"Cannot read AD4Zn.dat: {exc}")
        return problems

    def require(self) -> dict[str, Path]:
        if self.problems:
            raise FileNotFoundError("AD4Zn setup needs attention:\n" + "\n".join(self.problems))
        return dict(self.components)


def adfr_installations() -> list[Path]:
    """Find conventional ADFR installations without changing the host environment."""
    roots = [Path.home()]
    for variable in ("ProgramFiles", "ProgramFiles(x86)", "LOCALAPPDATA"):
        value = os.environ.get(variable)
        if value:
            roots.append(Path(value))
            if variable == "LOCALAPPDATA":
                roots.append(Path(value) / "Programs")
    found = []
    for root in roots:
        try:
            found.extend(path for path in root.glob("ADFRsuite*") if path.is_dir())
        except OSError:
            continue
    # A custom installer location may be registered instead of using a default folder.
    if os.name == "nt":
        import winreg
        for hive in (winreg.HKEY_CURRENT_USER, winreg.HKEY_LOCAL_MACHINE):
            for view in (winreg.KEY_WOW64_32KEY, winreg.KEY_WOW64_64KEY):
                try:
                    with winreg.OpenKey(hive, r"SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall", 0, winreg.KEY_READ | view) as entries:
                        for index in range(winreg.QueryInfoKey(entries)[0]):
                            try:
                                with winreg.OpenKey(entries, winreg.EnumKey(entries, index)) as entry:
                                    name = str(winreg.QueryValueEx(entry, "DisplayName")[0])
                                    if "adfr" in name.lower():
                                        location = str(winreg.QueryValueEx(entry, "InstallLocation")[0])
                                        if location and Path(location).is_dir():
                                            found.append(Path(location))
                            except OSError:
                                continue
                except OSError:
                    continue
    return sorted(set(found), key=lambda path: tuple(int(n) for n in re.findall(r"\d+", path.name)), reverse=True)


def check_ad4zn_environment(project_root: Path) -> AD4ZnEnvironment:
    """Discover only; explicit invalid overrides never fall back silently."""
    names = {"python": ("pythonsh.exe", "pythonsh") if os.name == "nt" else ("pythonsh",),
             "autogrid": ("autogrid4.exe",) if os.name == "nt" else ("autogrid4",),
             "zinc_pseudo": ("zinc_pseudo.py",), "prepare_gpf": ("prepare_gpf4zn.py",),
             "parameters": ("AD4Zn.dat",)}
    roots = [application_root() / "tools" / "ad4zn", Path(project_root) / "tools" / "ad4zn"]
    installations = adfr_installations()
    components = {}
    for key, candidates in names.items():
        executable = key in {"python", "autogrid"}
        override = os.environ.get("MOLDOCKPIPE_AD4ZN_" + key.upper(), "").strip()
        paths = [Path(override).expanduser()] if override else [root / name for root in roots for name in candidates]
        if override and executable:
            located = shutil.which(override)
            if located:
                paths.append(Path(located))
        if not override and key in {"python", "autogrid"}:
            for installation in installations:
                if key == "python":
                    paths += [installation / "python.exe"] if os.name == "nt" else [installation / "bin" / "pythonsh"]
                else:
                    paths += [installation / "bin" / name for name in candidates]
        if not override:
            paths += [Path(found) for name in candidates if (found := shutil.which(name))]
        components[key] = next((path.resolve() for path in paths if path.is_file()
            and path.stat().st_size and (not executable or os.name == "nt" or os.access(path, os.X_OK))), None)
    return AD4ZnEnvironment(components)


def atom_types(path: Path) -> tuple[str, ...]:
    values = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            value = line[77:].strip()
            if not re.fullmatch(r"[A-Za-z][A-Za-z0-9]*", value):
                raise ValueError(f"Invalid AutoDock atom type in {path.name}: {value!r}")
            try:
                coordinates = [float(line[start:start+8]) for start in (30, 38, 46)]
                if not all(math.isfinite(number) for number in coordinates):
                    raise ValueError()
            except ValueError as exc:
                raise ValueError(f"Invalid coordinates in {path.name}") from exc
            values.add(value)
    if not values:
        raise ValueError(f"No atom records in {path}")
    return tuple(sorted(values))


def collect_ligand_atom_types(paths: Iterable[Path]) -> tuple[str, ...]:
    values = set()
    for path in paths:
        values.update(atom_types(path))
    if not values:
        raise ValueError("AD4Zn map generation needs prepared ligand PDBQTs")
    if "TZ" in values:
        raise ValueError("TZ pseudoatoms belong to the receptor, not ligands")
    return tuple(sorted(values))


def require_zinc(receptor: Path) -> None:
    if not any(value.upper() == "ZN" for value in atom_types(receptor)):
        raise ValueError("AD4Zn requires retained Zn in the prepared receptor. Review receptor preparation.")


def _run(command: list[str], work: Path, records: list, cancelled: Callable[[], bool] | None = None) -> str:
    """Poll with communicate so large tool output cannot deadlock a pipe."""
    process = subprocess.Popen(command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    while True:
        try:
            stdout, stderr = process.communicate(timeout=0.2)
            break
        except subprocess.TimeoutExpired:
            if cancelled and cancelled():
                process.kill()
                stdout, stderr = process.communicate()
                records.append({"argv": command, "return_code": process.returncode, "stdout": stdout, "stderr": stderr})
                (work / "commands.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
                raise InterruptedError("AD4Zn preparation cancelled")
    records.append({"argv": command, "return_code": process.returncode, "stdout": stdout, "stderr": stderr})
    (work / "commands.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    if process.returncode:
        raise RuntimeError(f"AD4Zn tool failed ({Path(command[0]).name}): {stderr[-2000:] or stdout[-2000:]}. Logs: {work}")
    return stdout + "\n" + stderr


def probe_environment(tools: dict[str, Path], work: Path, records: list, cancelled=None) -> dict:
    runtime = _run([str(tools["python"]), "-c",
        "import sys, string, MolKit, AutoDockTools; assert hasattr(string, 'split'); print(sys.version)"], work, records, cancelled)
    version_text = _run([str(tools["autogrid"]), "--version"], work, records, cancelled)
    match = re.search(r"\b(4\.\d+\.\d+)(?:\.\d+)?", version_text)
    if not match or tuple(map(int, match[1].split("."))) < (4, 2, 7):
        raise ValueError("AD4Zn requires verified AutoGrid >= 4.2.7 (ADFR 4.2.7.x.2019-07-11 or newer). " + version_text[:300])
    return {"runtime": runtime.strip(), "autogrid_version": version_text.strip(),
            "tools": {key: {"path": str(path), "sha256": file_sha256(path)} for key, path in tools.items()}}


def prepare_zinc_receptor(base: Path, work: Path, tools: dict[str, Path], records: list, cancelled=None) -> Path:
    require_zinc(base)
    work.mkdir(parents=True, exist_ok=True)
    local = work / "receptor.pdbqt"
    shutil.copy2(base, local)
    output = work / "receptor_TZ.pdbqt"
    _run([str(tools["python"]), str(tools["zinc_pseudo"]), "-r", local.name, "-o", output.name], work, records, cancelled)
    require_zinc(output)
    if "TZ" not in atom_types(output):
        raise ValueError("Zinc helper produced no TZ pseudoatoms; review Zn coordination geometry. " + str(work))
    # Existing Zn must not disappear when generating pseudoatoms.
    def zinc_count(path):
        return sum(line[77:].strip().upper() == "ZN" for line in path.read_text().splitlines()
                   if line.startswith(("ATOM  ", "HETATM")))
    if zinc_count(base) != zinc_count(output):
        raise ValueError("Zinc helper changed the number of retained Zn atoms")
    return output


def grid_definition(profile: dict) -> dict:
    spacing = round(float(profile.get("ad4zn", {}).get("spacing", 0.375)), 3)
    center = [round(float(profile[f"center_{axis}"]), 3) for axis in "xyz"]
    sizes = [float(profile[f"size_{axis}"]) for axis in "xyz"]
    if not all(math.isfinite(v) for v in [spacing, *center, *sizes]) or spacing <= 0 or min(sizes) <= 0:
        raise ValueError("AD4Zn grid dimensions and spacing must be finite and positive")
    npts = [2 * math.ceil(size / (2 * spacing)) for size in sizes]
    if max(npts) > 126:
        raise ValueError("AD4Zn grid exceeds 126 intervals per axis; reduce the box or explicitly increase ad4zn.spacing")
    return {"spacing": spacing, "center": center, "npts": npts}


# Official prepare_gpf4zn.py zinc pair potentials. Validate values, not only keywords.
_ZINC_PAIRS = {("NA", "TZ"): (0.25, 23.2135, 12, 6), ("OA", "ZN"): (2.1, 3.8453, 12, 6),
               ("SA", "ZN"): (2.25, 7.5914, 12, 6), ("HD", "ZN"): (1, 0, 12, 6),
               ("NA", "ZN"): (2, 0.006, 12, 6), ("N", "ZN"): (2, 0.2966, 12, 6)}


def validate_gpf(path: Path, required: tuple[str, ...], grid: dict) -> None:
    entries: dict[str, list[list[str]]] = {}
    for line in path.read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if fields:
            entries.setdefault(fields[0], []).append(fields[1:])
    def one(key):
        values = entries.get(key, [])
        if len(values) != 1:
            raise ValueError(f"GPF needs exactly one {key}")
        return values[0]
    if one("parameter_file") != ["AD4Zn.dat"] or one("receptor") != ["receptor_TZ.pdbqt"]:
        raise ValueError("GPF must reference local AD4Zn.dat and receptor_TZ.pdbqt")
    if set(one("ligand_types")) != set(required):
        raise ValueError("GPF ligand types do not cover the requested library")
    if not {"ZN", "TZ"} <= {v.upper() for v in one("receptor_types")}:
        raise ValueError("GPF is missing Zn/TZ receptor types")
    if list(map(int, one("npts"))) != grid["npts"] or not math.isclose(float(one("spacing")[0]), grid["spacing"], abs_tol=1e-6):
        raise ValueError("GPF grid dimensions differ from the requested box")
    if any(not math.isclose(a, b, abs_tol=1e-4) for a, b in zip(map(float, one("gridcenter")), grid["center"])) or len(one("gridcenter")) != 3:
        raise ValueError("GPF grid center differs from the requested box")
    for key, expected in (("gridfld", "receptor_TZ.maps.fld"), ("elecmap", "receptor_TZ.e.map"), ("dsolvmap", "receptor_TZ.d.map")):
        if one(key) != [expected]:
            raise ValueError(f"Unexpected GPF {key}")
    if sorted(entries.get("map", [])) != sorted([[f"receptor_TZ.{t}.map"] for t in required]):
        raise ValueError("GPF affinity map list is incomplete")
    pairs = {(v[-2].upper(), v[-1].upper()): tuple(map(float, v[:4])) for v in entries.get("nbp_r_eps", [])}
    if any(pairs.get(pair) != values for pair, values in _ZINC_PAIRS.items()):
        raise ValueError("GPF is missing the official AD4Zn zinc pair potentials")


def validate_maps(work: Path, required: tuple[str, ...], grid: dict) -> list[Path]:
    paths = [work / "receptor_TZ.maps.fld"] + [work / f"receptor_TZ.{t}.map" for t in (*required, "e", "d")]
    for path in paths:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing/empty AD4Zn map artifact: {path.name}")
    field = paths[0].read_text()
    for path in paths[1:]:
        if path.name not in field:
            raise ValueError(f"Field descriptor does not reference {path.name}")
        with path.open() as stream:
            headers = [next(stream, "").split() for _ in range(6)]
            if [h[0] if h else "" for h in headers] != ["GRID_PARAMETER_FILE", "GRID_DATA_FILE", "MACROMOLECULE", "SPACING", "NELEMENTS", "CENTER"]:
                raise ValueError(f"Malformed AutoGrid header: {path.name}")
            if list(map(int, headers[4][1:])) != grid["npts"] or not math.isclose(float(headers[3][1]), grid["spacing"], abs_tol=1e-6):
                raise ValueError(f"Incorrect map dimensions: {path.name}")
            if len(headers[5]) != 4 or any(not math.isclose(float(a), b, abs_tol=1e-4) for a, b in zip(headers[5][1:], grid["center"])):
                raise ValueError(f"Incorrect map center: {path.name}")
            count = 0
            for line in stream:
                for value in line.split():
                    if not math.isfinite(float(value)):
                        raise ValueError(f"Non-finite map value: {path.name}")
                    count += 1
            if count != math.prod(n + 1 for n in grid["npts"]):
                raise ValueError(f"Truncated map: {path.name}")
    return paths


def ensure_ad4zn_maps(repository, profile: dict, ligands: Iterable[Path], executable: Path, cancelled=None) -> tuple[Path, dict]:
    """Build before docking any state; reuse only byte-verified complete bundles."""
    tools = check_ad4zn_environment(repository.root).require()
    receptor = (repository.root / str(profile["receptor"])).resolve()
    require_zinc(receptor)
    required = collect_ligand_atom_types(ligands)
    grid = grid_definition(profile)
    # Each profile owns its bundles, even when the receptor was imported/shared.
    profile_id = str(profile["id"])
    if not re.fullmatch(r"[A-Za-z0-9_-]+", profile_id):
        raise ValueError("Unsafe receptor profile identifier")
    root = repository.root / "inputs" / "receptors" / profile_id / "ad4zn"
    inputs = {"base_receptor": file_sha256(receptor), "vina": file_sha256(executable),
              **{key: file_sha256(path) for key, path in tools.items()}}
    parameters = {fields[1] for line in tools["parameters"].read_text().splitlines()
                  if len(fields := line.split()) > 1 and fields[0] == "atom_par"}
    if "TZ" not in parameters:
        raise ValueError("AD4Zn.dat must contain the TZ atom definition; a symlink placeholder is not a parameter file")
    unsupported = set(required) - parameters
    if unsupported:
        raise ValueError("AD4Zn.dat has no parameters for ligand atom types: " + ", ".join(sorted(unsupported)))
    request = fingerprint(settings={"protocol": "ad4zn", "grid": grid, "atom_types": required}, inputs=inputs, tool_version="ad4zn-v1")
    current = root / "current.json"
    if current.is_file():
        try:
            record = json.loads(current.read_text())
            folder = (repository.root / record["directory"]).resolve()
            folder.relative_to(root.resolve())
            if record["request"] == request and all((folder / name).is_file() and file_sha256(folder / name) == value
                                                     for name, value in record["hashes"].items()):
                require_zinc(folder / "receptor_TZ.pdbqt")
                if "TZ" not in atom_types(folder / "receptor_TZ.pdbqt"):
                    raise ValueError("Cached receptor has no TZ atoms")
                validate_gpf(folder / "receptor_TZ.gpf", required, grid)
                validate_maps(folder, required, grid)
                return folder / "receptor_TZ", record
        except (OSError, ValueError, KeyError, TypeError):
            pass
    work = root / ("bundle-" + uuid.uuid4().hex)
    work.mkdir(parents=True)
    records: list = []
    from ..project import utc_now
    try:
        versions = probe_environment(tools, work, records, cancelled)
        help_text = _run([str(executable), "--help"], work, records, cancelled)
        if not all(flag in help_text for flag in ("--maps", "--scoring", "ad4")):
            raise ValueError("Selected Vina does not advertise AD4 map scoring")
        versions["vina_version"] = _run([str(executable), "--version"], work, records, cancelled).strip()
        prepare_zinc_receptor(receptor, work, tools, records, cancelled)
        shutil.copy2(tools["parameters"], work / "AD4Zn.dat")
        _run([str(tools["python"]), str(tools["prepare_gpf"]), "-r", "receptor_TZ.pdbqt", "-o", "receptor_TZ.gpf", "-n",
              "-p", "ligand_types=" + ",".join(required), "-p", "npts=" + ",".join(map(str, grid["npts"])),
              "-p", "gridcenter=" + ",".join(map(str, grid["center"])), "-p", "spacing=" + str(grid["spacing"]),
              "-p", "parameter_file=AD4Zn.dat"], work, records, cancelled)
        validate_gpf(work / "receptor_TZ.gpf", required, grid)
        _run([str(tools["autogrid"]), "-p", "receptor_TZ.gpf", "-l", "autogrid.glg"], work, records, cancelled)
        paths = validate_maps(work, required, grid)
        paths += [work / name for name in ("AD4Zn.dat", "receptor_TZ.pdbqt", "receptor_TZ.gpf")]
        hashes = {path.name: file_sha256(path) for path in paths}
        record = {"protocol": "ad4zn", "request": request, "directory": work.relative_to(repository.root).as_posix(),
                  "maps_prefix": (work / "receptor_TZ").relative_to(repository.root).as_posix(),
                  "atom_types": list(required), "grid": grid, "inputs": inputs, "hashes": hashes,
                  "map_fingerprint": fingerprint(settings=grid, inputs=hashes, tool_version="ad4zn-v1"),
                  "versions": versions, "created_at": utc_now(), "warnings": ["Target-specific scientific validation is pending."]}
        (work / "manifest.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
        for path in work.iterdir():
            if path.is_file():
                repository.add_artifact(path, "ad4zn_preparation", "ad4zn")
        repository.record_provenance_event(event_type="ad4zn_maps_prepared", stage_name="ad4zn",
            receptor_profile_id=profile_id, data=record)
        temporary = root / ("current-" + uuid.uuid4().hex + ".tmp")
        temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
        temporary.replace(current)
        return work / "receptor_TZ", record
    except Exception as exc:
        (work / "failure.txt").write_text(str(exc), encoding="utf-8")
        repository.record_provenance_event(event_type="ad4zn_preparation_failed", stage_name="ad4zn", level="ERROR",
            receptor_profile_id=profile_id, message=str(exc), data={"artifacts": work.relative_to(repository.root).as_posix()})
        raise
