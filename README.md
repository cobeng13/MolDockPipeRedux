# MolDockPipe Redux

PyQt6 desktop workflow for traceable virtual screening projects. It keeps submitted
ligands, generated molecular states, docking runs, poses, and scientific artifacts
separate. The legacy scripts remain in `Reference/` as behavioral reference only.

## Installation

RDKit is best installed from conda-forge. Install the application into that environment:

```powershell
conda create -n moldockpipe-clean -c conda-forge python=3.13 rdkit pyqt6 pip
conda activate moldockpipe-clean
pip install -e ".[chemistry,dev]"
python -m pytest -q
moldockpipe
```

Keep RDKit and PyQt6 Conda-managed. The project installs MolScrub, Meeko, SciPy,
Gemmi, Joblib, and the test dependencies through pip. Do not install a pip PyQt6
wheel into the same environment.

For a headless Linux or macOS environment, PyQt is not required:

```bash
conda create -n moldockpipe -c conda-forge python=3.13 rdkit pip
conda activate moldockpipe
pip install -e ".[chemistry,dev]"
python -m pytest -q
moldockpipe run /path/to/project
```

The worker runs the existing pipeline and writes the same portable project
layout. Select individual stages with `--stages`, for example
`moldockpipe run /path/to/project --stages screening molscrub meeko`.

MolScrub generation has a 15-second time limit per molecule, adjustable in the
MolScrub settings. Slow generation is stopped, recorded as failed, and skipped
while the rest of the pipeline continues. Matching failures are reused on later
runs; changing MolScrub settings permits another attempt. Disconnected molecular
fragments fail before generation starts. Failed molecules remain visible in the
dashboard and do not contribute states to ligand preparation or docking.

## Project Inputs

In the UI, New Project asks for a project name and creates it under `Projects/`
in the application's current working folder. You can also create a portable
project at an explicit location with:

```powershell
moldockpipe create C:\path\to\project
```

Place the default ligand input at:

```text
project/inputs/input.csv
```

The CSV columns are:

```text
id,smiles,notes,params_json
```

`id` must be a stable, unique identifier. Imports are cumulative per project:
new compounds are added as pending and previously imported compounds and their
results remain in the library even when absent from the next CSV. A header-only
CSV is a no-op. Reimporting the same ID and SMILES preserves completed work;
notes and metadata updates do not reset processing. An existing ID with a
different SMILES is rejected: use a new ID for a different compound. CSV rows
with missing or duplicate IDs are rejected without changing the project.
Identity is based on ID and the submitted SMILES string; different IDs are
separate compounds, even if their SMILES match.

Use the ribbon **Receptors** manager to import one or more prepared receptor
PDBQT files. Imported files are copied into the portable project:

```text
project/inputs/receptors/profile_id/receptor.pdbqt
```

Opening a project, refreshing, and starting pipeline stages load the latest
`inputs/input.csv` into the cumulative library. The command-line run also imports
it before processing. Emptying the file never empties the library. Existing
projects upgrade archived compounds back into the library on opening.

## Workflow

The main window provides individual stage buttons and **Run All**:

```text
Screening -> MolScrub states -> Meeko PDBQT -> Vina docking -> post-docking export
```

Screening uses RDKit physicochemical rules. Its default policy is **Annotate only**;
the Settings dialog can instead exclude compounds failing any/all selected rules or
send them to manual review.

Run All and individual stages execute in the background. The interface shows:

- Overall checkpoint strip: red pending, yellow active, green complete.
- Current-stage item progress.
- Current ligand/state and explicit succeeded/skipped/failed summaries.

The project database is the processing ledger; exported CSVs are snapshots,
not required for resuming work. Screening checks its settings/input fingerprint
before calculating descriptors. MolScrub reuses matching generated states and
verified SDFs, Meeko skips verified prepared PDBQTs, and Vina reuses completed
runs per state and receptor when receptor/settings/tool fingerprints and output
checksums match. Changed scientific settings can require fresh work. Failed
preparation/docking or missing/damaged artifacts are retried. Closing during work marks unfinished runs as interrupted when
the project is reopened.

## Settings

The **Settings** button contains tabs for:

- Screening: policy and Lipinski, Veber, Egan, Ghose, and BOILED-Egg (BBB/Yolk) rules.
- MolScrub: pH, state enumeration, state limit, per-compound timeout, and parallel worker count.
- Meeko: parallel worker count.
- Post-docking: split/export mode, poses per compound, and successfully docked compound selection.
- Guardrails: purge generated workflow data and export data.

MolScrub and Meeko default to **Automatic**, using one worker per detected
logical CPU, limited by the number of compounds or states awaiting preparation.
In `project.yml`, `workers: 0` selects Automatic; positive values set a manual
limit. Older projects using Meeko's previous default of four workers switch to
Automatic once on upgrade. Other manual limits and choices made after upgrading
are preserved. Both stage Settings tabs offer Automatic and manual limits.

Both stages keep only one batch of jobs in flight. Stop prevents new submissions
and saves results from active jobs before exiting. MolScrub retains a separate
timeout for each compound. Changing worker count does not invalidate completed
states, prepared ligands, or docking results.

A project allows one pipeline run at a time across application windows and
terminal commands. Imports, recovery, and workflow resets cannot overlap that
run. Different projects can run independently. Closing an active pipeline asks
it to stop safely and waits for its worker thread before closing the window.

The separate **Receptors** manager contains each receptor's search box,
exhaustiveness, modes, energy range, seed, and CPU count. Docking runs every
enabled receptor sequentially against the same prepared ligand states. The
**Results** window displays one ranked tab per receptor.

## Vina and Post-Docking Tools

The development Vina executable is supplied in the repository under:

```text
tools/vina/vina.exe
```

`vina_1.2.7_win.exe` is also detected. Bundled Vina 1.2.7 binaries are selected
for Windows, Linux x86_64, macOS Intel, and macOS Apple Silicon automatically.
Linux ARM requires separately supplied compatible tools. Generic `vina` and
`vina_split` filenames remain supported on Linux and macOS. Tools can be placed in `tools/vina/`, supplied on `PATH`, or
selected explicitly with `MOLDOCKPIPE_VINA` and `MOLDOCKPIPE_VINA_SPLIT`.
Binaries under `tools/vina/` are tracked for fast development.

Optional Vina-GPU docking can be selected in **Settings > Docking**. CPU Vina
remains the default for existing and new projects, and receptor validation
redocking continues to use CPU Vina. Install a compatible Vina-GPU executable
and its matching `Kernel1_Opt.bin` and `Kernel2_Opt.bin` under `tools/vinagpu/`,
or set `MOLDOCKPIPE_VINAGPU` to the executable with both kernels beside it.
Automatic discovery checks only application/project `tools/vinagpu/` folders.
For the supplied Windows runtime, the executable is `Vina-GPU+.exe`; no
`VinaGPUConfig.txt`, kernel compiler executable, or AutoDock-GPU executable is needed.
The adapter targets the bundled legacy Vina-GPU command interface; other
builds must support the same flags. GPU-specific search controls default to
8000 computing lanes and search depth 32. Receptor boxes, seed, number of modes,
and energy range are shared; CPU count and exhaustiveness are not sent to GPU.
Each molecular state runs separately with isolated kernel copies, logs, and
normal pose artifacts. Executable, kernel, and GPU setting changes invalidate
GPU resume matches. GPU and CPU scores may differ. Compatible GPU drivers and
kernel binaries are required; no automatic CPU fallback occurs on GPU failure.

Post-docking supports:

- Split only: writes individual pose PDBQTs.
- Split and convert to SDF: splits selected poses, then runs `mk_export.py`.
- Convert multi-pose output to SDF: exports the Vina multi-model file directly.

The compound selector lists successfully docked compounds by most-negative score,
shows export status, and prevents selecting compounds already exported.

Outputs are organized as:

```text
project/For_PostDocking/profile_id/SDF/parent_id/state_id/run_id/
project/For_PostDocking/profile_id/PDBQTs/parent_id/state_id/run_id/
```

## Exports and Test Reset

The ribbon **Export Data** menu writes:

```text
project/exports/profile_id/manifest.csv
project/exports/profile_id/leaderboard.csv
```

Each manifest includes receptor identity and its reproducibility hashes. Each
leaderboard contains the top three poses per parent compound for that receptor,
sorted by most-negative affinity, across the cumulative project library.
Regenerating exports includes all imported batches. Workflow Maintenance can clear generated data
while preserving inputs and `project.yml`.
