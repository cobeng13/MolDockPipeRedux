"""Isolated MolScrub generation, allowing the caller to stop slow chemistry."""
from __future__ import annotations

import base64
import contextlib
import json
import sys

from .molscrub import MolScrubService


def main() -> None:
    request = json.load(sys.stdin)
    try:
        # Keep third-party diagnostic output separate from the result protocol.
        with contextlib.redirect_stdout(sys.stderr):
            states, truncated = MolScrubService()._generate(**request)
        records = [dict(state.__dict__, molecule=base64.b64encode(state.molecule.ToBinary()).decode("ascii"))
                   for state in states]
        response = {"states": records, "truncated": truncated}
    except Exception as exc:
        response = {"error": str(exc)}
    sys.stdout.write(json.dumps(response))


if __name__ == "__main__":
    main()
