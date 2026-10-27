# Interaction analysis implementation and validation

## Delivered behavior

The modular service is `src/moldockpipe/services/interaction_analysis.py`.
It discovers existing SDF exports, verifies a shared coordinate frame against
preparation provenance and raw docking poses, runs ProLIF once per pose against
one loaded receptor/reference, compares unique residue/type fingerprints, and
publishes three CSV reports and a reproducibility manifest atomically.

Integration uses `PipelineRunner.run_interactions`, existing project execution
locking and workflow/provenance history, schema migration 8, the existing Export
stage, and a dialog launched from Multidock Results. It supports manual receptor
batches, selected runs, optional post-export execution, cancellation, per-record
failure reporting, numeric sorting, pagination, and report downloads. No upstream
molecular preparation or docking workflow is introduced.

## Created files

- `src/moldockpipe/services/interaction_analysis.py`: scientific engine, input
  validation, pose mapping, comparison, streamed CSVs, reuse, safe publication.
- `src/moldockpipe/ui/interaction_dialog.py`: background worker and report UI.
- `tests/test_interaction_analysis.py`: chemistry, input, metrics, publication,
  workflow, CLI isolation, and migration regression tests.
- `tests/test_interaction_ui.py`: scope, numerical sorting, downloads,
  cancellation, close lifecycle, and opt-in settings tests.
- `docs/interaction-analysis.md`: usage, complete CSV schemas, scientific policy,
  reproduction, dependency installation, and limitations.
- `docs/interaction-analysis-implementation.md`: this delivery record.

## Modified files

- `src/moldockpipe/pipeline.py`: independent stage, export mappings, optional hook.
- `src/moldockpipe/project.py`: additive tables, versioning, crash recovery.
- `src/moldockpipe/__main__.py`: CLI stage and filters; analysis avoids input import.
- `src/moldockpipe/ui/results_dialog.py`: Analyze Interactions entry point.
- `src/moldockpipe/ui/receptor_manager.py`: rebase bundled reference paths when
  importing a receptor package with a new profile ID.
- `src/moldockpipe/ui/settings_dialog.py`: saved opt-in checkbox.
- `src/moldockpipe/ui/workers.py`: independent analysis-stage dispatch.
- `pyproject.toml`: ProLIF 2.2.2, MDAnalysis >=2.10,<3, pandas >=2.2,<4.
- `environment.yml`: defaults fontconfig >=2.18.2 for compatible native imports.
- `README.md`: feature entry point and documentation link.

## Executed validation

All tests ran in the existing `moldockpipe-clean` environment.

| Check | Result |
| --- | --- |
| Full regression: `python -m pytest -q --tb=short -p no:cacheprovider` | 321 passed, 3 skipped; 113 seconds |
| Final migration/export follow-up: `python -m pytest tests/test_interaction_analysis.py tests/test_project.py tests/test_idempotency.py tests/test_postdock.py -q --tb=short -p no:cacheprovider` | 75 passed; 70 seconds |
| Interaction chemistry and UI tests | 25 passed before the additional migration regression; included in the full/follow-up runs above |
| `python -m pip check` | No broken requirements |
| `git diff --check` | Passed |
| Real prepared receptor/reference, read-only | 530 residues, 72 reference contacts, 18 unique residue/type pairs |
| Public CLI on a copy of an existing MD_Test job | 3 poses analyzed, 0 failed; original pose IDs and affinities retained |

Chemical fixtures verify ligand donor and acceptor directionality, halogen bond
geometry and a near miss, metal coordination and a distance near miss, distinct
chains/insertion codes, and successful zero-contact calculations. Software tests
cover invalid/empty/multiple SDFs, missing conformers/files/raw poses, unreadable
files, zero denominators, duplicate contacts, stale and corrupted outputs,
mid-batch cancellation, publication failure, previous-report preservation,
completed-export reuse, and analysis without upstream work.

MDAnalysis emits a dependency deprecation warning about its old topology tables
import; the tests pass. The full suite's three optional tests remain skipped.

The real-job acceptance copy is under `build/interaction-acceptance-52f7584b`.
The original MD_Test project was not changed. On its three poses, retention was
55.56%, 61.11%, and 55.56%; Jaccard was 0.4545, 0.4231, and 0.4167. These values
demonstrate comparison/export behavior and are not affinity validation.

Validated package versions: Python 3.13.15, ProLIF 2.2.2, RDKit 2026.3.6,
MDAnalysis 2.11.0, pandas 3.0.3, Gemmi 0.7.5. The existing environment needed a
fontconfig/GLib repair to resolve an `intl-8.dll` conflict; the installation
constraint and repair instructions are recorded in the usage documentation.

## Use and limitations

Open Multidock Results → Analyze Interactions after normal SDF export. Enable
automatic execution in Settings → Post-docking if desired. CSVs are written under
`For_PostDocking/<profile>/Interactions/<analysis-id>/`; the dialog saves copies
to a selected location. Complete schemas and CLI examples are in
[the usage guide](interaction-analysis.md).

Profiles without compatible prepared receptor/reference/provenance artifacts
fail explicitly unless a reference is attached through the explicit manual route
described below. Waters and nonprotein cofactors are excluded and reported;
water bridges and ProLIF visualization are not implemented. Chemistry inference,
protonation, rigid receptor assumptions, and docking uncertainty remain scientific
limitations. The existing 3D viewer continues using its separate annotations.
One library call must finish before cancellation takes effect; analysis uses
sequential processing and supports at most 9999 included receptor residues.

## Manual reference attachment follow-up

**Load co-crystallized SDF…** in the analysis dialog attaches one immutable SDF to
the active receptor profile using `attach_crystallographic_reference`. The
attachment applies only to interaction analysis, preserves redocking metadata,
and does not invalidate completed docking runs. It requires no validation-run
database records or extracted reference PDB.

When available, original crystal-source evidence verifies both receptor and
ligand coordinates. Otherwise an explicit user declaration is required and
recorded as `user_confirmed_same_frame`. The manifest's `reference_origin` and
`coordinate_checks` distinguish this from source verification. The SDF and
prepared receptor hashes are checked on later runs. No alignment, extraction,
conversion, or receptor preparation is performed.

Receptor-package imports now rebase bundled reference paths and source-evidence
path keys to the new profile folder, fixing references to the exporting project.
Tests exercise actual package import, default and manual references with no
validation history, UI attachment, immutable copies, rejection of invalid or
multi-record SDFs, and both ligand and receptor frame mismatches.

Follow-up validation: full regression **327 passed, 3 skipped** (123 seconds);
final interaction/UI checks including receptor-frame rejection **32 passed**
(11 seconds). `pip check` and `git diff --check` passed.
