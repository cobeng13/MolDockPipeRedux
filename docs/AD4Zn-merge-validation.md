# AD4Zn integration into current main

Validation date: 2026-10-19. Integration branch: `codex/merge-ad4zn-main`.

## Recovery points and isolation

- Main parent: `4a0633a568faae0d07025956bf5f17bf90a00711`.
- Zinc parent: `ceb70636d947ecc1290249112248aa897340be1a` (`codex/ad4zn-support`).
- Both GitHub branch tips were verified with `git ls-remote` before integration.
- Integration uses a separate managed worktree. The original main checkout's
  uncommitted post-docking edits, deleted files, and untracked inputs are excluded.
- Both parent histories are retained. Main must remain at its original commit
  until the native acceptance gate below is satisfied.

## Conflict resolution and compatibility

The pipeline, repository, and receptor manager required manual conflict resolution.
The integration retains main's bounded planning reads, prepared-ligand hashes,
legacy Vina reuse, cancellation, logical CPU defaults, project locks, atomic
configuration writes, and stale-profile checks.

AD4Zn remains opt-in. Missing protocol values default to Vina. No database schema
migration is introduced. AD4Zn docking fingerprints include the map request and
validated map-content fingerprint, while ignoring incidental bundle timestamps
and locations. Switching between Vina and AD4Zn cannot reuse the other protocol's
scores. Changed prepared ligand bytes invalidate either protocol, preserving main's
newer behavior. Zinc preparation receives the pipeline cancellation callback and
propagates interruption instead of recording cancellation as failed docking.

Vina-GPU does not support the AD4Zn map-scoring interface. A zinc profile selected
with that backend fails with an actionable instruction to select Vina CPU; ordinary
profiles in the campaign still run. Exports retain protocol identity and use the
physical receptor instead of converting TZ pseudoatoms to physical atoms.

## Automated and runtime verification

- Final focused zinc/configuration/multidocking/redocking suite: **62 passed**.
- Native Vina AD4 smoke test: passed using synthetic maps; this checks command and
  pose parsing, not zinc scientific validity.
- Native standard Vina diagnostic: passed, producing a valid scored pose.
- PowerShell launcher syntax and `-CheckOnly -EnvironmentName moldockpiperedux`:
  passed; the launcher imported the integration checkout's source.
- Final full regression suite: **223 passed, 2 skipped** (255.06 seconds).
  The skipped checks require symlink privileges and an interactive native Qt
  desktop session. `git diff --check` also passed.

Run regression tests from the integration checkout with its source on PYTHONPATH:

```powershell
$env:PYTHONPATH = Join-Path (Get-Location) 'src'
conda run -n moldockpiperedux python -m pytest -q -p no:cacheprovider
```

Additional integration regressions cover protocol switches, map-provenance-only
changes, cancellation during map preparation, mixed Vina-GPU/zinc campaigns, and
post-docking protocol metadata. Configuration regression fixtures explicitly
include the now-normalized `protocol: vina` default.

## Native acceptance gate: BLOCKED

An actual map-preparation call was attempted with isolated copies of the existing
MMP13_Test 1GKC receptor and prepared reference ligand. The zinc helpers and valid
AD4Zn parameter file were found in the older checkout and supplied through explicit
process-local overrides. The attempt stopped before generating maps because:

- The ADFR-compatible Python interpreter is unavailable.
- AutoGrid4 is unavailable.

The former installation recorded in the historical reference manifests,
`C:/Program Files (x86)/ADFRsuite-1.0`, no longer contains those executables.
Conventional-installation discovery, PATH discovery, and the recorded executable
paths were checked. This integration does not install or bundle external tools.

The final full-suite log (`pytest-full-final.log`), focused log, native
preparation preflight JSON, copied reference
inputs, and the standard-Vina diagnostic are retained locally under the ignored
`Test_Bench/merge-validation/` directory in the integration worktree. Existing
reference projects were read without modification. Historical zinc run artifacts
are available but do not establish acceptance of this integrated implementation.

Before updating main, restore a compatible ADFR interpreter and AutoGrid >= 4.2.7,
supply the documented helper/parameter paths, run the official zinc tutorial and
a fresh MMP-13 reference redocking, inspect Zn/TZ geometry and metal coordination,
and assess pose recovery/RMSD against agreed target-specific acceptance criteria.
The existing guide specifies this review but no numerical RMSD acceptance cutoff;
do not silently invent one or treat synthetic-map tests as scientific acceptance.
Until these checks pass, retain the integration branch and original zinc branch
and leave main unchanged.
