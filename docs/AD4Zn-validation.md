# AD4Zn implementation verification

Feature base: main at 99ca1a5. Branch: codex/ad4zn-support.
Implementation is isolated in the MolDockPipe-ad4zn worktree.

## Changes

Added optional AD4Zn tool discovery, runtime/version probes, Zn/TZ validation,
ligand atom-type union, GPF/map validation, immutable bundles, and provenance.
Receptor profiles and preparation UI expose explicit protocol selection.
Screening and reference redocking share the map preparation service. Existing
profiles default to Vina and retain legacy default fingerprint behavior.
Exports and reports identify the scoring protocol. No database schema migration
is required. CSV exports append a protocol column.

## Verification performed

From the feature worktree, using the existing project Conda environment:

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) 'src')
conda run -n moldockpipe-clean python -m pytest -q -p no:cacheprovider
```

Result: **90 passed, 1 skipped** in 12.80 seconds. The skipped test is the existing
interactive desktop window-construction test. `git diff --check` passed.

Tests include actual bundled Windows Vina 1.2.7 AD4 docking with synthetic
zero-energy maps, verifying execution and pose parsing. Helper/AutoGrid generation
uses controlled fake outputs. Actual subprocess failure and cancellation are tested.
Coverage includes legacy profiles, mixed campaigns, scientific-input changes,
map corruption, ligand-type coverage, invalid zinc potentials, missing TZ,
Meeko Zn loss, invalid parameters, reuse, redocking, and exports.

## Remaining external validation

Host preflight found these components missing: ADFR-compatible interpreter,
AutoGrid, zinc_pseudo.py, prepare_gpf4zn.py, and AD4Zn.dat. They were not silently
installed or bundled. Supply them following [the setup guide](AD4Zn.md).

A native end-to-end zinc-toolchain run has not been completed. macOS/Linux AD4Zn
execution and interactive UI behavior have not been validated. MMP-13 coordination
inspection, reference-ligand pose recovery, and external RMSD assessment remain
the requested manual acceptance step. Synthetic-map success is not scientific
validation. Keep this feature separate from main until that review.
