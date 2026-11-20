# How MolDockPipe Redux works

MolDockPipe Redux is a desktop and command-line application that coordinates
ligand screening, molecular-state generation, ligand preparation, docking, and
post-docking export within a portable project. Its organizing principle is a
traceable relationship between a submitted compound, its generated states, each
receptor-specific docking attempt, and the resulting poses and files. It delegates
chemical preparation and docking to specialist tools rather than introducing a
new docking algorithm or scoring function.

This guide describes the implementation in this checkout. Installation and
operating instructions are in the [README](../README.md); the
[publication reporting guide](publication-reporting.md) explains what to retain
and report for a study. Default settings below describe new projects, not the
settings necessarily used in an existing campaign.

## Workflow and responsibilities

```text
Ligand CSV                         Receptor structure or prepared PDBQT
    |                                           |
RDKit descriptor screening         Receptor preparation/import and search box
    |                                           |
MolScrub molecular states                        |
    |                                           |
Meeko ligand PDBQT ------------------------------+
                                                |
                                 Dock each state against each enabled receptor
                                                |
                                 Store runs, scores, and pose identities
                                                |
                                 Split PDBQT and/or export SDF
                                                |
                                 Optional ProLIF/reference comparison
```

Receptor preparation is a separate prerequisite to campaign docking. The desktop
wizard analyzes the source structure, applies the user's component selections,
extracts an optional crystallographic reference ligand, and invokes Meeko receptor
preparation. It saves the structural inventory, preparation decisions, commands,
logs, and checked outputs. The structure and reference-template file pickers
include the Downloads folder in their sidebar when the system provides one.
Importing a prepared PDBQT provides a second route;
such an import does not by itself supply the full prepared PDB or reference
provenance needed for interaction comparison.

The default headless stage sequence is `screening molscrub meeko vina postdock`.
The `vina` stage name dispatches the configured docking backend; it does not
necessarily mean CPU execution. Interaction comparison is an explicit stage or an
optional action after SDF export. Selecting stages executes them in the order
supplied, without automatically adding missing upstream stages.

## Compound, state, run, and pose identity

| Entity | Meaning and relationship |
| --- | --- |
| Parent ligand | One imported compound identified by the submitted CSV `id`, with its original SMILES and annotations. |
| Molecular state | A generated chemical representation belonging to one parent, with isomeric SMILES, formal charge, structure hash, and generation fingerprint. |
| Conformer | A coordinate-bearing SDF artifact associated with a state. The current generation path stores one conformer per retained state. |
| Receptor profile | A named prepared receptor with its own identity, search box, protocol, parameters, and enabled status. |
| Docking run | One attempt for a state and receptor profile, with settings, input identity, command information, status, and output/log artifacts. |
| Docking pose | One mode in a run, identified separately from its score and mode index. |
| Artifact | A registered file with a project-relative location, SHA-256 digest, size, type, and producing stage. |

A parent is not a pose: a compound may generate several states, and each state
may yield several poses for several receptors. Counts and rankings must name the
entity being counted. State enumeration does not create independent submitted
compounds or experimental replicates.

Imports accumulate in the project library. Reimporting an existing ID with the
same submitted SMILES preserves its work; changing that ID's SMILES is rejected.
Compounds missing from a later CSV remain in the library. Different IDs remain
separate parents even when their SMILES match. Consequently, a later export can
include compounds from earlier import batches.

## What each stage computes

### 1. Descriptor screening

RDKit parses the submitted SMILES and computes molecular weight, Crippen MolLogP,
topological polar surface area (TPSA), hydrogen-bond donor/acceptor counts, and
rotatable bonds. Molecular weight, MolLogP, and TPSA are rounded to two decimal
places before applying the predicates below.

| Setting name | Predicate implemented in `screen_smiles` |
| --- | --- |
| `lipinski` | MW <= 500, MolLogP <= 5, donors <= 5, and acceptors <= 10. |
| `veber` | TPSA <= 140 and rotatable bonds <= 10. |
| `egan` | MolLogP <= 5.88 and TPSA <= 131. |
| `ghose` | 160 <= MW <= 480 and -0.4 <= MolLogP <= 5.6. |
| `boiled_egg` | YOLK flag: TPSA < 79 and 0.4 <= MolLogP <= 6.0. WHITE annotation: TPSA <= 130 and -0.4 <= MolLogP <= 5.6 when YOLK is false; otherwise GREY. |

These are the application's explicit descriptor predicates. In particular, the
Ghose setting includes only MW and MolLogP, and the BOILED-Egg labels use
rectangular thresholds. Do not describe these as complete reproductions of the
named published methods or as validated ADMET predictions.

New projects enable Lipinski and Veber with `annotate_only`: rule flags annotate
compounds rather than remove them. Other policies exclude compounds failing any
selected rule, exclude those failing all selected rules, or mark rule failures
`manual_review`. The current MolScrub eligibility check excludes decisions
`fail` and missing decisions under non-annotation policies; `manual_review` is
not an enforced approval gate. Invalid SMILES cannot be processed successfully
by later chemistry stages even when annotation-only policy permits an attempt.

### 2. Molecular-state generation

The MolScrub adapter passes the configured pH as both `ph_low` and `ph_high`.
New projects use pH 7.4, enumeration enabled, and a maximum of 32 retained states.
Generated states are deduplicated by isomeric SMILES and formal charge; their
coordinates are written to SDF. Reaching the configured state limit sets an
enumeration-truncation flag. A retained set is therefore bounded and must not be
described as exhaustive.

Disconnected structures are rejected before generation. With enumeration enabled,
generation runs in an isolated process with a per-compound timeout (15 seconds
when unspecified). With enumeration disabled, the adapter returns the parsed
source molecule in-process without running the enumerating/embedding call; this
route does not establish generated 3D coordinates.
Failures are recorded and supply no states to downstream preparation. Matching
known failures are normally reused on resume; explicit retry or relevant setting
changes can permit another attempt.

### 3. Ligand preparation

Meeko's ligand-preparation command converts each retained state's SDF into ligand
PDBQT. The adapter retains its log, and the pipeline validates the output before
marking the state ready for docking. Preparation failures are recorded per state.
MolScrub and Meeko can process independent jobs with bounded parallel workers;
worker count controls scheduling rather than molecular identity.

### 4. Receptor-specific docking

The pipeline docks prepared ligand states against enabled receptor profiles,
processing profiles and states sequentially. Each profile supplies its search-box
center and dimensions, modes, energy range, seed, and protocol-specific controls.
CPU Vina is the default. Optional Vina-GPU uses its own search controls and
executable/kernel identity; a GPU failure does not trigger automatic CPU fallback.
The optional AD4Zn protocol adds zinc-specific preparation and map scoring; see
[AD4Zn](AD4Zn.md) for toolchain requirements and validation boundaries.

A [redocking validation campaign](redocking-rmsd.md) can test a bounded series of
boxes against the crystallographic reference before production. Applying the
winner saves its box and validation search settings as production defaults for
that receptor. The validation record and each production run preserve the
settings used at execution time.

The raw multi-model output and log are retained, and parsed modes become pose
records. The database distinguishes current runs from historical attempts.
Results are ranked by ascending docking affinity (more negative first).
Leaderboards retain up to three poses per parent per receptor across active
states and current completed runs. This is a best-pose ranking, not an ensemble
free-energy calculation or a state-population-weighted estimate.

### 5. Post-docking export and interaction comparison

Export operates on current completed runs for active states and enabled profiles,
optionally restricted to selected parents. It can split multi-model PDBQT,
convert split poses to SDF, or convert the multi-model output directly to SDF.
The `poses_per_compound` setting limits the split poses selected for conversion
within each state/run; it is not a global cap across all states of a parent.
Direct multi-model SDF conversion does not apply that slice. The split tool can
leave additional pose PDBQTs in the output directory. Report actual exported SDF
record counts rather than deriving them from the setting name.

There are two distinct interaction features. The Results viewer draws local
geometric contacts for visual inspection. The optional ProLIF analysis compares
exported SDF records with a crystallographic reference in the prepared receptor's
coordinate frame and publishes detailed, summary, and comparison CSVs plus a
manifest. These features use different detectors and input requirements; their
counts are not interchangeable. See [interaction analysis](interaction-analysis.md)
for coordinate checks, detector roles, exclusions, and metric definitions.

## Storage and execution architecture

```text
project/
  project.yml                  Project configuration and receptor profiles
  project.sqlite               Processing ledger, identities, and provenance
  inputs/input.csv             Default cumulative-import source
  inputs/receptors/<profile>/  Receptor files and preparation/reference bundles
  artifacts/sdf/               Generated state coordinates
  artifacts/pdbqt/             Prepared ligand inputs
  artifacts/docking/           Receptor/parent/state/run docking outputs
  logs/                        Stage logs and preparation failure evidence
  exports/<profile>/           Manifest and leaderboard CSV snapshots
  For_PostDocking/<profile>/   Pose exports and interaction-analysis bundles
```

The SQLite database is the ledger used to find prior work; CSV reports are
snapshots and are not required for resuming. Registered artifact paths are
relative to the project root. Copy the complete project after stopping all
writers, rather than copying only exported CSVs or a live SQLite main file whose
transactions may still reside in its WAL. External tools, native dependencies,
and host-specific executable overrides must be installed separately on another
machine.

| Source module | Responsibility |
| --- | --- |
| [`project.py`](../src/moldockpipe/project.py) | Project creation, SQLite schema/migration, imports, artifacts, and provenance. |
| [`pipeline.py`](../src/moldockpipe/pipeline.py) | Stage orchestration, eligibility, reuse, and persistence of results. |
| [`services/`](../src/moldockpipe/services/) | Chemistry adapters, docking backends, export, and interaction analysis. |
| [`receptors/`](../src/moldockpipe/receptors/) | Structural inventory, component selection, reference chemistry, receptor preparation, and AD4Zn support. |
| [`redocking/`](../src/moldockpipe/redocking/) | Reference-ligand redocking, pose-level RMSD, and resumable box validation campaigns. |
| [`ui/`](../src/moldockpipe/ui/) | Desktop controls, background workers, and results visualization. |
| [`__main__.py`](../src/moldockpipe/__main__.py) | Headless command dispatch and desktop entry point. |
| [`reporting.py`](../src/moldockpipe/reporting.py) | Structured project summaries and report snapshots. |

The UI and CLI invoke the same `PipelineRunner` stage methods. A project-level
writer lock prevents overlapping pipeline work, imports, recovery, and resets
across processes. Background UI workers keep long operations off the GUI thread;
this does not mean that all workflow stages run simultaneously.

## Reuse, failure, and provenance boundaries

Reuse is stage-specific, not a universal promise of identical dependency tracking:

- Screening compares settings, submitted SMILES, and RDKit version fingerprints,
  including matching failed results.
- MolScrub checks generation compatibility and SDF hashes for successful states;
  known-failure eligibility also considers the attempted settings. Worker count
  is excluded, and a timeout change does not invalidate successful chemistry.
- Meeko currently reuses a prepared output when the conformer status is ready
  and the registered PDBQT hash matches. It does not have a complete independent
  signature for every preparation setting or installed tool version.
- Docking matches effective backend parameters, ligand/receptor hashes, and
  executable identity (including GPU kernels), then verifies the raw output hash.
  There is also an explicit compatibility path for older fingerprints.
- Post-docking reuse checks for expected output folders containing files of the
  required type. It does not establish complete requested-pose coverage or
  checksum integrity of every exported file.
- ProLIF snapshot reuse checks its analysis fingerprint and published report
  checksums; see the specialist guide for its input provenance requirements.

Missing or changed upstream artifacts can require regeneration. Failed and
interrupted attempts remain distinguishable from completed work. Stopping
parallel preparation prevents new submissions and allows active jobs to save
results; it need not terminate every active job immediately. Reopening a project
recovers unfinished records as interrupted.

The [persistent-signatures design](architecture/stage-signatures-and-pending-work.md)
is a proposal, not an implemented guarantee. Hashes establish file identity and
support traceability; they do not establish chemical correctness, experimentally
valid predictions, or bitwise reproduction across hardware and tool versions.

## What validation currently establishes

Automated tests cover project behavior, adapters, imports, resume rules, export,
and UI logic. Some tests substitute fake tools or fixtures. The `doctor` command
runs dependency and small execution probes, while the native graphics test is
opt-in. Neither constitutes a target-specific docking benchmark.

Reference redocking records `ARTIFACTS_READY` when its checked SDF/MOL2 files
are ready. With PyDockRMSD installed, the app calculates symmetry-aware,
heavy-atom RMSD in the receptor coordinate frame and stores separate pose-level
results. The rank-1 and lowest measured RMSDs answer different questions.
Vina's `rmsd_lb` and `rmsd_ub` compare poses within the docking output.
A campaign's rank-1 cutoff evaluates recovery of the reference pose for that
receptor and protocol. Report predictive performance using separate study data.
