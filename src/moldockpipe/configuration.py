from __future__ import annotations

import os
import tempfile
from pathlib import Path

import yaml


def read_project_config(path: Path) -> dict:
    try:
        config = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ValueError(f"Project configuration is unreadable: {path}. Check {path.name}.bak before re-importing receptors.") from exc
    if not isinstance(config, dict):
        raise ValueError(f"Project configuration is empty or invalid: {path}. Check {path.name}.bak before re-importing receptors.")
    return config


def _replace_bytes(path: Path, content: bytes) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}-", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def write_project_config(path: Path, config: dict) -> None:
    """Replace a complete YAML file and keep its previous valid version."""
    content = yaml.safe_dump(config, sort_keys=False).encode("utf-8")
    if path.exists():
        read_project_config(path)
        _replace_bytes(path.with_suffix(path.suffix + ".bak"), path.read_bytes())
    _replace_bytes(path, content)
