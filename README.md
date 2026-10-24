# MolDockPipe Redux

PyQt6 desktop workflow for traceable virtual screening projects. It keeps submitted
ligands, generated molecular states, docking runs, poses, and scientific artifacts
separate. The legacy scripts remain in `Reference/` as behavioral reference only.

## First-time installation (Windows desktop)

Install **64-bit Miniconda or Anaconda** and **Git for Windows** first. Open an
Anaconda/Miniconda PowerShell Prompt. Internet access is needed for installation;
the molecular viewer itself works offline afterward.

Clone the repository and run these commands from its root (the folder containing
`environment.yml` and `pyproject.toml`):

```powershell
git clone https://github.com/cobeng13/MolDockPipeRedux.git
cd MolDockPipeRedux
conda env create --file environment.yml
conda activate moldockpipe-clean
python -m pip check
python -m moldockpipe doctor --output diagnostics.txt
python -m moldockpipe ui
```

`environment.yml` installs Python 3.13, RDKit, PyQt 6.11, **PyQtWebEngine 6.11**,
and the application's chemistry/test dependencies. RDKit is explicitly sourced
from conda-forge so pip recognizes the installed package; all native Qt packages
come from the matched defaults-channel distribution. The environment name matches
both Windows launchers. After setup, double-click `Launch-MolDockPipe.cmd` to start
the application; the terminal command above also works.

There is no separate `requirements.txt`: `environment.yml` supplies the native
desktop packages and installs `.[chemistry,dev]` from `pyproject.toml` through
pip. Keep RDKit and Qt Conda-managed. Do not separately install pip `PyQt6` or
`PyQt6-WebEngine-Qt6` runtime wheels into this environment.

Review `diagnostics.txt` before continuing if required chemistry or CPU Vina
checks fail. An optional Vina-GPU failure is expected on unsupported GPUs. The
WebEngine check verifies bindings and bundled viewer files; to verify graphics,
open a project with completed docking results, select a compound in Results,
rotate its receptor/ligand view, and export a PNG. Working WebGL is required.

For development, `python -m pytest -q` runs the regression suite. CPU Vina and
Vina Split binaries are included in the source checkout. GPU docking and the
AD4Zn protocol have additional runtime requirements documented under `tools/`.

### Updating an existing installation

Close the application and run `git pull` from the repository root. For an
environment originally created using this `environment.yml`:

```powershell
conda env update --file environment.yml
conda activate moldockpipe-clean
python -m pip check
python -m moldockpipe ui
```

For an older or mixed-channel environment, create a separate clean environment
instead of mixing native Qt distributions. If `moldockpipe-clean` already exists,
choose a new name without deleting it:

```powershell
conda env create --name moldockpipe-desktop --file environment.yml
conda activate moldockpipe-desktop
python -m moldockpipe ui
```

To use that custom name with the double-click launcher, create a shortcut whose
target ends with `Launch-MolDockPipe.cmd -EnvironmentName moldockpipe-desktop`,
or run `./Launch-MolDockPipe.cmd -EnvironmentName moldockpipe-desktop` in
PowerShell. Existing project folders do not need to be recreated.

### Optional interactive Results viewer

Open **Multidock Results**, then select a completed compound. The panel beside
the receptor tabs displays its prepared receptor and best ligand pose. The mode
selector switches among poses from that same state and docking run. Drag to
rotate, scroll to zoom, choose receptor Cartoon/Sticks, toggle the receptor,
focus the ligand, or reset the initial camera. **Export PNG** saves only the
molecular canvas, at twice its displayed width and height with a white background.
Existing PDB/SDF exports remain available through the Export menu.

The viewer works offline: 3Dmol.js 2.5.5 and its license notices ship with the
application. Structures come directly from docking artifacts; running postdock
export first is unnecessary. Missing or changed structures produce an error
instead of displaying an unmatched complex.

Interaction toggles display colored dashed lines for **H-bonds** (blue),
**π–π** (purple), **cation–π** (orange), **salt bridges** (green), and
**hydrophobic contacts** (gray, initially hidden). Counts include all detected
contacts, including hidden types. Contacted residues appear as sticks on the
cartoon. Toggle **Contact labels** to show residue identity, interaction type,
and distance in Å; hiding the receptor also hides these annotations.
PNG exports include the visible contacts and labels. **Interaction criteria…**
explains the rules inside Results. No additional installation is required.

These are geometric predictions, not validated bonds, interaction energies,
or changes to docking scores. Analysis runs in the existing loading worker,
uses the original coordinates, and never generates hydrogens or changes
protonation. RDKit chemical features identify ligand donors/acceptors and
aromatic 5/6-member rings; ligand charge groups require actual SDF formal charge.
Protein roles use standard amino-acid templates. ASP/GLU are assumed negative
(unless an explicit acid hydrogen is present), LYS/ARG positive; histidine
H-bond roles require explicit ring hydrogens. Uncertain protonation can alter
assignments. Waters, metals, cofactors, modified residues, water bridges,
halogen bonds, and metal coordination are outside this first implementation.
Only the first receptor model and blank/A alternate locations are considered.

| Contact | Geometry used |
| --- | --- |
| H-bond | Donor–acceptor ≤3.5 Å. With explicit donor H: H–A ≤2.6 Å and D–H–A ≥120°. Without donor H: every heavy-neighbor–D–A angle ≥90°; explicitly labeled **candidate**. |
| π–π | Ring centers ≤5.5 Å; acute plane angle ≤30° (parallel) or ≥60° (T-shaped); minimum of the two ring-plane offsets ≤2 Å. Ring atoms must be planar within 0.15 Å. |
| Cation–π | Positive group–ring center ≤6 Å; in-plane offset ≤2 Å. |
| Salt bridge | Opposite charge group centers ≤4 Å. |
| Hydrophobic | Nonpolar C···C ≤4 Å; carbon neighbors limited to C/H/S; shortest contact per residue. |

Distances below 0.5 Å are excluded as coincident coordinates, not contacts.
Ring and charged-group distances are measured between centers, not nearest atoms.
This is a local geometric detector, **not PLIP**. Its ring/offset rules and
hydrophobic cutoff are informed by the
[PLIP documentation](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md);
H-bond and salt-bridge cutoffs are deliberately stricter than PLIP defaults.
Ligand feature definitions come from
[RDKit BaseFeatures](https://github.com/rdkit/rdkit/blob/master/Data/BaseFeatures.fdef).

WebEngine is included in the first-time desktop setup above. For older
installations, use the clean-environment instructions rather than installing
unmatched Qt packages into an existing environment. The native WebEngine library
alone does not supply the PyQt6 Python bindings.

Restart the application after installation. Working WebGL support is required.
If WebEngine/WebGL is unavailable, the ranked Results table still works; the
viewer displays setup guidance. Headless commands require neither WebEngine
nor PyQt.

The optional native rendering test can be run from a Windows desktop session:

```powershell
$env:MOLDOCKPIPE_WEBENGINE_TEST = "1"
$env:QT_QPA_PLATFORM = "windows"
python -m pytest tests/test_results_viewer.py -q
```

It verifies a structure payload over 2 MB and a nonblank PNG at exactly twice the
viewport dimensions. The regular test suite uses a fake web canvas for UI
selection/cancellation tests and skips this native graphics test by default.

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

Resume planning batches screening, MolScrub, Meeko and docking metadata reads
while retaining the existing fingerprints and full artifact-hash checks. Meeko
shares missing-directory resolution during preflight and checks cancellation
before dispatch. See [resume planning](docs/resume-planning.md) for behavior,
validation and the separate persistent-work architecture proposal.

## Troubleshooting student installations

Open **Tools → Troubleshooter…** (no project is required), click **Run checks**,
then copy or save the report. It records the app's Python executable, environment,
package versions and locations, full errors, and suggested actions. Each check
runs separately with a time limit, so an import crash does not stop other checks.
Reports include local file paths; review them before sharing.

The checks parse and embed a molecule with RDKit, generate 3D states through the
MolScrub worker, prepare an ethanol PDBQT with Meeko, load its receptor/export
commands, and run tiny docking jobs with Vina CPU and Vina GPU. GPU testing also
checks the kernel files and exercises the OpenCL runtime. They also import the
WebEngine/WebChannel bindings and verify the bundled offline viewer files; this
does not test the graphics driver or WebGL rendering. A missing GPU backend
is expected if only CPU docking is installed. These are smoke tests, not checks
of scientific docking accuracy or every receptor preparation scenario.

The same checks work without the desktop UI:

```powershell
python -m moldockpipe doctor --output diagnostics.txt
python -m moldockpipe doctor --project C:\path\to\project --json --output diagnostics.json
```

The command exits with code 1 if any check fails, including an unavailable optional
GPU backend or desktop viewer in a headless installation. It uses temporary
fixtures and does not change project data.

An installed package is not necessarily importable: a dependency, incompatible
API, or native DLL may fail at import time. MolScrub errors now retain that cause.
Compare the report's Python path with the environment where pip was run. Use
`python -m pip check` and `python -m pip install -e ".[chemistry,dev]"` with that
same environment activated, then restart the app. RDKit remains Conda-managed.

If a resumed run appears to pause after screening, the task display now names
the next stage before its setup begins. MolScrub dependency loading runs in a
separate process with a 30-second limit and responds to Stop; a stalled import
produces a troubleshooting error instead of waiting indefinitely. Completed
stages retain 100% progress. Active dashboard refreshes do not scan export
directories; existing export checks remain in the export workflow.

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
Enter **Start** and **End** docking scores (for example, -10 to -7 kcal/mol),
then click **Select score range** to replace the checked selection with matching
compounds. Both endpoints are included. Selection uses each compound's best
score across enabled receptors, for standard Vina or AD4Zn. Confirm with **OK**
to export; the existing export settings still determine which states and poses
of those compounds are written.

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

## Optional zinc-metalloprotein docking

The feature branch adds receptor-specific AutoDock4Zn preparation and Vina AD4
map scoring. Standard Vina remains the default. See [AD4Zn setup and validation](docs/AD4Zn.md)
for the required external tools and manual validation steps.

On Windows, after installing the environment, double-click `Launch-MolDockPipe.cmd`
to launch this checkout without manual activation or PYTHONPATH settings.
