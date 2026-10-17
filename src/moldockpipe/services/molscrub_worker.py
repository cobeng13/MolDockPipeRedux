"""Isolated MolScrub generation, allowing the caller to stop slow chemistry."""
from __future__ import annotations

import base64
import contextlib
import json
import sys
import traceback

from .molscrub import MolScrubService


def main() -> None:
    if "--version" in sys.argv[1:]:
        with contextlib.redirect_stdout(sys.stderr):
            try:
                import molscrub
                version = getattr(molscrub, "__version__", "unknown")
            except ImportError:
                # Preserve the pipeline's previous fingerprint behavior when
                # the package cannot be imported; actual generation reports it.
                traceback.print_exc(file=sys.stderr)
                version = "unknown"
        sys.stdout.write(json.dumps({"version": version}))
        return
    request = json.load(sys.stdin)
    try:
        # Keep third-party diagnostic output separate from the result protocol.
        with contextlib.redirect_stdout(sys.stderr):
            states, truncated = MolScrubService()._generate(**request)
        records = [dict(state.__dict__, molecule=base64.b64encode(state.molecule.ToBinary()).decode("ascii"))
                   for state in states]
        response = {"states": records, "truncated": truncated}
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        response = {"error": str(exc)}
    sys.stdout.write(json.dumps(response))


if __name__ == "__main__":
    main()
