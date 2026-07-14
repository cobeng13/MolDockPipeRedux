# MolDockPipe Redux

PyQt6 desktop workflow for traceable virtual screening projects. It keeps submitted
ligands, generated molecular states, docking runs, poses, and scientific artifacts
separate. The legacy scripts remain in `Reference/` as behavioral reference only.

## Installation

RDKit is best installed from conda-forge. Install the application into that environment:

```powershell
conda create -n moldockpipe-clean -c conda-forge python=3.13 rdkit pyqt=6 pip
conda activate moldockpipe-clean
pip install -e ".[chemistry,dev]"
python -m pytest -q
moldockpipe
```

Keep RDKit and PyQt6 Conda-managed. The project installs MolScrub, Meeko, SciPy,
Gemmi, Joblib, and the test dependencies through pip. Do not install a pip PyQt6
wheel into the same environment.

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

`id` must be a stable, unique identifier. Re-importing the CSV synchronizes the
project incrementally: unchanged compounds retain their results, new compounds
are added as pending, changed compounds are reprocessed from screening, and
compounds removed from the CSV are archived rather than deleted. CSV rows with
missing or duplicate IDs are rejected without changing the project.

Use the ribbon **Receptors** manager to import one or more prepared receptor
PDBQT files. Imported files are copied into the portable project:

```text
project/inputs/receptors/profile_id/receptor.pdbqt
```

Opening a project automatically loads `inputs/input.csv` when no active ligand
set is already present. Importing another CSV synchronizes it with the current
ligand set after confirmation when changes or archived compounds are detected.

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

Completed work is reused when fingerprints and artifacts match. Failed or missing
artifacts are retried. Closing during work marks unfinished runs as interrupted when
the project is reopened.

## Settings

The **Settings** button contains tabs for:

- Screening: policy and Lipinski, Veber, Egan, Ghose, and BOILED-Egg (BBB/Yolk) rules.
- MolScrub: pH, state enumeration, and state limit.
- Meeko: parallel worker count.
- Post-docking: split/export mode, poses per compound, and successfully docked compound selection.
- Guardrails: purge generated workflow data and export data.

The separate **Receptors** manager contains each receptor's search box,
exhaustiveness, modes, energy range, seed, and CPU count. Docking runs every
enabled receptor sequentially against the same prepared ligand states. The
**Results** window displays one ranked tab per receptor.

## Vina and Post-Docking Tools

The development Vina executable is supplied in the repository under:

```text
tools/vina/vina.exe
```

`vina_1.2.7_win.exe` is also detected. Place `vina_split.exe` in either
`tools/vina/` or `tools/`; binaries under `tools/vina/` are tracked for fast development.

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
sorted by most-negative affinity. Workflow Maintenance can clear generated data
while preserving inputs and `project.yml`.
