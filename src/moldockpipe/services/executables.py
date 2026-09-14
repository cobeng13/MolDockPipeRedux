from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ExecutableSpec:
    """Describe an external tool without coupling it to one operating system."""

    name: str
    environment_variable: str
    posix_names: tuple[str, ...]
    windows_names: tuple[str, ...]

    @property
    def platform_names(self) -> tuple[str, ...]:
        return self.windows_names if os.name == "nt" else self.posix_names


VINA = ExecutableSpec(
    name="AutoDock Vina",
    environment_variable="MOLDOCKPIPE_VINA",
    posix_names=("vina", "vina_1.2.7"),
    windows_names=("vina.exe", "vina_1.2.7_win.exe", "vina"),
)

VINA_SPLIT = ExecutableSpec(
    name="Vina Split",
    environment_variable="MOLDOCKPIPE_VINA_SPLIT",
    posix_names=("vina_split",),
    windows_names=("vina_split.exe", "vina_split"),
)


def application_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _usable_executable(path: Path) -> bool:
    if not path.is_file():
        return False
    return os.name == "nt" or os.access(path, os.X_OK)


def _explicit_executable(value: str) -> Path | None:
    path = Path(value).expanduser()
    if _usable_executable(path):
        return path.resolve()
    located = shutil.which(value)
    return Path(located).resolve() if located else None


def find_executable(spec: ExecutableSpec, project_root: Path) -> Path:
    """Resolve a tool from an override, bundled locations, or the active PATH."""

    override = os.environ.get(spec.environment_variable, "").strip()
    if override:
        executable = _explicit_executable(override)
        if executable is None:
            raise FileNotFoundError(
                f"{spec.name} configured by {spec.environment_variable} is not an executable: {override}"
            )
        return executable

    roots = (application_root(), Path(project_root).resolve())
    candidates = (
        root / relative / filename
        for root in roots
        for relative in (Path("tools") / "vina", Path("tools"))
        for filename in spec.platform_names
    )
    executable = next((candidate.resolve() for candidate in candidates if _usable_executable(candidate)), None)
    if executable is not None:
        return executable

    for filename in spec.platform_names:
        located = shutil.which(filename)
        if located:
            return Path(located).resolve()

    names = ", ".join(spec.platform_names)
    raise FileNotFoundError(
        f"{spec.name} executable not found. Set {spec.environment_variable}, place one of "
        f"[{names}] under tools/vina, or add it to PATH."
    )
