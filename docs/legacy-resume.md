# Resuming older MolScrub projects

When an older project stored MolScrub's tool version as `unknown` and the
installed version now has a version number, MolDockPipe stops preparation
before dispatching any new MolScrub jobs if it finds reusable historical
results. It explains how to retain the original campaign inputs. This prevents
an environment update from silently replacing structures and triggering
Meeko/docking again.

To continue the original campaign:

1. Open **Settings > MolScrub**.
2. Enable **Retain verified legacy results** and save.
3. Run the desired stages again.

The choice persists in the project's `project.yml`:

```yaml
resume:
  retain_legacy_molscrub: true
```

The same setting applies to GUI and command-line stage runs. It is outside
the chemical settings block, so enabling it cannot itself invalidate chemistry.
No application restart or database migration is required.

Saved states must match the configured generation settings under their
original fingerprint, and each SDF must match its registered SHA-256 value.
Projects predating MolScrub stage history can also reuse states when their
states and conformers have completed PDBQT preparation and the SDF artifacts
are active and valid. SDF files alone cannot establish completion. A newer failed,
cancelled, or interrupted MolScrub attempt prevents state reuse.

Retaining original structures does **not** declare an unknown historical
MolScrub build scientifically equivalent to the installed one. Fingerprints,
SDFs, conformer statuses, and docking identities stay unchanged. A provenance
event records the affected parents and the current and historical tool
identities. New molecules use the installed toolchain; use a separate project
when all molecules must be regenerated using one current toolchain.

Matching historical failures are also retained when this choice is enabled.
An explicit `run_molscrub(retry_failed=True)` retries them. Increasing their
time budget can also make them eligible for retry; merely saving an omitted
default 15-second limit does not.

Changed pH/state limits, unrecognized generation fingerprints, missing or
damaged SDFs, and incomplete legacy metadata remain ineligible for reuse.
Retaining structures cannot recover files that were already replaced or
deleted. Meeko and docking still independently verify their own reuse criteria.

New chemistry installations pin MolScrub to the tested `0.3.1` source revision
`0f8de09c7532d5b187c9376684196ab2c7755320`, rather than following Git HEAD.
Existing environments are not modified by this code patch.

## Patch verification (2026-11-03)

- Full suite: 358 passed, 3 skipped. The only warning was an upstream
  MDAnalysis deprecation.
- Final legacy/default-timeout and bulk-planning checks: 36 passed.
- An offscreen settings-dialog check verified the retention checkbox loads
  and saves both enabled and disabled values.
- Read-only planning for the 5,043-parent local project identified 1,425
  completed parents (2,045 intact SDF states) and 629 known failures for
  retention, leaving 2,989 parents to generate. The five older prepared
  projects without MolScrub history also remained eligible for legacy retention.

Regression tests cover the preflight stop before any scientific writes,
repeated fresh-runner resume through MolScrub/Meeko/Vina, preservation of
original fingerprints and docking records, incomplete-generation rejection,
damaged files, changed chemical settings, and explicit failure retries.
Student projects were not modified during verification.
