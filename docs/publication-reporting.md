# Reporting a MolDockPipe study for publication

This guide supports methods writing and preparation of a reproducibility package.
It is not a completed manuscript or evidence of predictive performance. The
[workflow explanation](how-it-works.md) describes implemented behavior; the
study-specific values and results below must come from the campaign being
reported.

## Describe the software's contribution

A defensible implementation description is:

> MolDockPipe Redux coordinates descriptor screening, molecular-state generation,
> ligand preparation, receptor-specific docking, and pose export through a shared
> desktop and command-line workflow. A project-local SQLite ledger links submitted
> compounds to generated states, docking attempts, poses, and registered artifacts.
> Stage-specific reuse checks and retained settings, logs, and hashes support
> inspection and continuation of previous work.

Adapt this text to the released version and functions actually evaluated. Name
RDKit, MolScrub, Meeko, the docking engine, and ProLIF where used, and cite their
appropriate primary publications and exact software releases. The application
coordinates these tools; claims of a new preparation method, scoring function,
validated ADMET model, or improved predictive accuracy need separate evidence.

## Record the methods actually used

| Topic | Details to report and retain |
| --- | --- |
| Software identity | Repository URL, release/tag and full commit, local modifications, application version, and date of the experiment. Package version alone may identify several different development checkouts. |
| Execution environment | OS and architecture, Python and resolved package versions, CPU and logical cores, RAM, worker limits, and GPU/driver/OpenCL details when relevant. |
| Compound source | Library origin/version, access date, identifiers, submitted SMILES, preprocessing performed outside this application, duplicate policy, and all imported batches included in the cumulative project. |
| Screening | Enabled predicates, their implemented thresholds, RDKit version, policy, and numbers annotated or excluded. Specify the simplified Ghose and rectangular BOILED-Egg predicates if used; method names alone are insufficient. |
| State generation | MolScrub source revision/version, pH, enumeration enabled/disabled, state cap, timeout, retained-state counts, failures, and truncation flags. Do not infer tautomer/protomer populations from enumeration. |
| Ligand preparation | Meeko version, actual command/settings, failed states, and hashes of the SDF and prepared PDBQT inputs used for the reported runs. |
| Receptor preparation | Structural source/accession and model, chains and alternate locations, reference-ligand identity, retained/removed waters, ions and cofactors, any external repairs, component decisions, preparation commands, and final receptor hash. |
| Search space | Box center and dimensions with units and coordinate frame, how the box was chosen, and whether it was shared with reference redocking. |
| Docking | Backend/version and executable hash; seed, modes, energy range, CPU count and exhaustiveness for CPU Vina; GPU computing lanes/search depth and kernel hashes for Vina-GPU. Keep protocol-specific settings separate. |
| AD4Zn, if used | Protocol, physical zinc and pseudoatom handling, AutoGrid and helper versions/hashes, parameter file, map-generation settings and map identity, and target-specific validation. |
| Result selection | Receptor profiles and current run IDs, state/pose aggregation, ranking rule, score cutoffs, missing/failed runs, ties, and any manual exclusions. |
| Export | Mode, selected parent IDs, requested pose limit, actual exported records per state/run, and source run/pose mapping. The split/conversion limit is not a global parent-level pose cap. |
| Interaction analysis | Viewer contacts or ProLIF, detector version/settings, prepared protein and reference identity, coordinate-frame evidence, exclusions, candidate mapping failures, denominators, and exact snapshot/manifest. |

`project.yml` describes the current project configuration. Historical runs can
have different settings; use their recorded commands and fingerprints to explain
the results actually included. Review `unknown` version fields and imported
profiles with incomplete provenance rather than silently filling gaps from the
current environment.

Use **Export → Project report (HTML)** to capture `report.html`,
`report_data.json`, and `source_run_manifest.json`. The report shows recorded
redocking settings, recorded settings grouped across current completed production
runs for active states, and the receptor's saved production defaults. When
available, use the recorded settings for the results you cite. Applying a
validation winner updates the defaults for future runs. To report older
production results, use their original report snapshot or stored run commands.

## Account for compounds, states, and poses separately

Report a flow of submitted parents, descriptor-screening decisions, parents with
retained states, generated states, prepared states, attempted docking runs,
completed runs, retained poses, exported SDF records, and analyzed records.
Break docking and analysis counts down by receptor/protocol. Include failures,
known failures skipped on resume, and enumeration truncation.

Define the denominator of each percentage. A compound with several states and
poses contributes multiple computational observations, not multiple independent
experimental samples. Imports are cumulative, so the most recently supplied CSV
may not describe the entire exported library. A leaderboard contains up to three
poses per parent per receptor, not every generated mode.

## Separate engineering verification from scientific evaluation

| Evidence | What it can support | What it does not establish |
| --- | --- | --- |
| Unit/regression tests | The tested software contracts on the stated fixtures and environment. | Binding prediction accuracy or general operation of every external toolchain. |
| `doctor` report | Dependency imports and tiny chemistry/backend execution probes. | A validated target, exhaustive installation certification, or working WebGL rendering. |
| SHA-256 and reuse records | Identity of checked files and the recorded stage-specific reuse decision. | Chemical validity or universal detection of every dependency change. |
| Redocking `ARTIFACTS_READY` | Checked reference and pose files. Read the separate RMSD record to establish whether a measurement succeeded. | Prospective screening performance. |
| Docking score or interaction fingerprint | The engine's score or a geometric comparison under the declared protocol. | Measured affinity, experimentally confirmed contacts, or biological efficacy. |

For a redocking evaluation, record the RMSD tool and version, heavy-atom and
symmetry treatment, coordinate alignment policy, reference structure, acceptance
criterion, and pose ranking rule. The integrated PyDockRMSD calculation uses the
receptor coordinate frame. The validation campaign stops at a measured RMSD
strictly below 2 Å for affinity-ranked pose 1. If you report the minimum across
generated poses, give its rank separately. Vina's within-output RMSD bounds
describe a different comparison. See the [redocking guide](redocking-rmsd.md).

Choose additional scientific evaluations to match the manuscript's claims. For
example, a ranking-performance claim requires an appropriately described active/
inactive or decoy dataset, selection policy, metrics, and uncertainty estimates;
a runtime claim requires hardware, workload, backend settings, and repeated
measurements. Declare how failed jobs and resumed work enter those measurements.
CPU Vina, Vina-GPU, and AD4Zn results should retain their backend/protocol labels;
comparative claims require an explicit evaluation design. This repository's
integration tests do not supply those study results.

## Preserve a reproducibility package

After stopping all project writers, preserve the complete project directory,
including its database, configuration, inputs, raw docking outputs, preparation
logs, and analysis manifests. Preserve the exact code revision and local patch,
and capture the resolved environment rather than relying only on dependency
ranges. For example, from the repository root in the activated study environment:

```powershell
git rev-parse HEAD
git status --short
git diff --binary --output=study-code.patch
conda env export --no-builds > study-environment.yml
python -m pip freeze > study-python-packages.txt
python -m moldockpipe doctor --json --output study-diagnostics.json
```

These commands write evidence files into the working directory; retain them with
the study materials. The diff captures tracked-file changes, not untracked source
files: save any untracked files needed to reproduce the study separately.
A Conda export is useful evidence but does not guarantee a cross-platform or
permanently installable environment. MolScrub is installed from a Git dependency;
retain its resolved commit, not only the repository URL. Preserve external
engine/kernel/helper versions and redistribution information as applicable.

The diagnostic command returns a nonzero exit status if any probe fails,
including unavailable optional GPU or viewer components. Explain which probes
were required for the reported workflow and retain failures. An installation
probe is supplementary to the actual campaign logs and scientific evaluation.

Before sharing a package, review paths, source-data permissions, and third-party
licenses. Record a persistent archive identifier for the released code and data
where available. A manuscript's software/data availability statement should name
the exact archived version and explain any unavailable input data or external
tools; a development checkout alone is not a permanent archive.

## Methods paragraph template

Replace every bracketed field using the retained campaign evidence. Remove
optional clauses for features that were not used.

> Calculations were coordinated with MolDockPipe Redux [release and commit] on
> [OS, CPU/RAM, and GPU if used] using [resolved environment/archive]. The compound
> library [source/version] contained [N] submitted parent records. Descriptor
> screening used [implemented predicates and policy]. MolScrub [revision]
> generated states at pH [value], with enumeration [enabled/disabled], a cap of
> [value] states per parent, and a timeout of [value] seconds; [counts] were
> retained, failed, or truncated. Ligand PDBQT files were prepared using Meeko
> [version/settings]. Receptor [accession/model/chains] was prepared with
> [component decisions and method], and the search box was [center, dimensions,
> units, and selection procedure]. Docking used [engine/version/protocol] with
> [backend-specific parameters and seed]. Results were selected by [receptor,
> run scope, and parent/state/pose ranking rule]. [If performed: reference
> RMSD evaluation used tool/version, atom-mapping/alignment method, and criterion,
> with measured results reported in table/figure.] [If performed: ProLIF
> version/settings compared specified exported poses with reference identity,
> using coordinate provenance and exclusions recorded in manifest.] Code, inputs,
> run evidence, resolved environment, and outputs are available at [archive].
