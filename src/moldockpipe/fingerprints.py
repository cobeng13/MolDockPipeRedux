from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


PIPELINE_SCHEMA_VERSION = 1


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(*, settings: dict[str, Any], inputs: dict[str, str], tool_version: str) -> str:
    payload = {
        "schema": PIPELINE_SCHEMA_VERSION,
        "settings": settings,
        "inputs": inputs,
        "tool_version": tool_version,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()
