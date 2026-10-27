# Protein–ligand interaction analysis

Open **Multidock Results**, choose a receptor tab, and click **Analyze
Interactions**. The dialog defaults to all existing exported poses from current
completed runs for that receptor. Choose **Selected docking run** to restrict it.
The numerical columns sort across the whole snapshot; pages contain 200 rows.
The Save buttons copy the three CSV reports to a chosen destination.

### Manually attach a co-crystallized reference

For imported receptors, choose **Load co-crystallized SDF…** in the same dialog.
Select an SDF containing exactly one crystallographic ligand pose, already in the
prepared receptor's coordinate frame. Confirm its identity and frame in the
dialog. MolDockPipe copies the SDF into the project, binds it to the receptor
hashes, and shows the selected reference above the summary table. Run **Analyze
Interactions** again after changing the reference. Existing reports retain the
reference recorded when they were generated.

No receptor-validation/redocking database run is required. The attachment is an
interaction-analysis override; existing redocking reference files and metadata
remain intact. New receptor-package imports also rebase bundled reference paths
to the imported profile's new directory.

When the original crystal source and preparation report are available, attachment
verifies the ligand's heavy-atom coordinates against that source and rejects a
mismatch. When those artifacts are unavailable, the SDF alone cannot establish
crystallographic origin or a shared frame. The explicit user declaration is then
recorded as `user_confirmed_same_frame`, rather than reported as independently
verified. The manifest records this distinction. No fitting or alignment occurs.
An existing full prepared receptor PDB matching the docking receptor is still
required; loading a ligand does not prepare or reconstruct the receptor.

Python callers can use
`attach_crystallographic_reference(repository, profile_id, sdf_path,
confirm_same_frame=True)` from `moldockpipe.services.interaction_analysis`.
Without source-coordinate evidence, omitting that declaration rejects attachment.

Use the existing Export action first if SDF poses have not been exported. No
conversion, extraction, receptor preparation, or docking is performed by analysis.
Completed runs can be exported and analyzed before the whole campaign completes.
The existing project writer lock prevents simultaneous docking and analysis.

In Settings → Post-docking, enable **Analyze interactions after SDF export** for
automatic execution. It is off by default. Existing reusable exports also qualify.
Split-only export does not supply SDFs and therefore skips automatic analysis.
An analysis error leaves docking results and export status available.

For headless use in the existing environment:

```powershell
conda activate moldockpipe-clean
python -m moldockpipe run "Projects/MyProject" --stages interactions
python -m moldockpipe run "Projects/MyProject" --stages interactions --profile-id PROFILE_ID
python -m moldockpipe run "Projects/MyProject" --stages interactions --profile-id PROFILE_ID --run-id RUN_ID --force-interactions
```

Analysis alone does not reimport the ligand input CSV or run upstream stages.
Unchanged successful analyses are reused after report checksum validation.
**Create a new snapshot** / `--force-interactions` requests a fresh calculation.

## Inputs and scientific policy

The profile's existing `receptor` and `reference_ligand` metadata determine input
locations. SDFs are discovered under the existing
`For_PostDocking/<profile>/SDF/<parent>/<state>/<run>/` hierarchy, restricted to
active states and current completed docking runs.

Analysis uses the existing post-Meeko `receptor_prepared.pdb`, preferring the
identical post-docking receptor export when present. It requires the preparation
report, unchanged docking receptor hash, cleaned PDB, original structural source,
and extracted reference PDB/SDF to verify provenance and coordinates. Older or
externally imported profiles without this evidence fail with a specific message;
the analyzer does not fabricate provenance or recreate structures. The explicit
manual-attachment route above supports imported profiles lacking crystal-source
provenance, recording the user's declaration and binding the reference to the
unchanged prepared receptor instead.

Reference heavy atoms must retain crystallographic coordinates. Every candidate
must match a source docking pose by an element-preserving coordinate bijection,
without fitting or alignment. The per-axis tolerance is 0.001 Å, accounting for
PDB/PDBQT (three decimal places) and SDF (four decimal places) serialization.
Spatial proximity alone is never used to certify a common frame.

MDAnalysis converts the prepared protein to RDKit, inferring bonds, bond orders,
and charges using connectivity from RDKit's PDB reader rather than the generic
MDAnalysis van der Waals bond guesser. Atom mapping and coordinates are checked. Existing explicit protein
hydrogens are required. RDKit adds missing ligand hydrogens only on in-memory
copies, without changing existing heavy atoms or formal charges. There is no
energy minimization, independent protonation workflow, or pose optimization.
Inference cannot resolve every nonstandard residue or ambiguous protonation state;
conversion failures block that receptor's analysis rather than omitting residues.

Protein amino-acid residues and retained physical metal ions are included.
Retained waters and other cofactors are explicitly listed as exclusions in the
manifest. Unsupported ATOM protein residues fail explicitly. AD4Zn TZ pseudoatoms
are never analyzed. No water-mediated interactions are implemented. Residues retain
chain, name, number, and insertion code. Temporary internal residue numbering
prevents ProLIF from merging insertion-code variants.

The enabled ProLIF 2.2.2 detectors are HBDonor, HBAcceptor, Hydrophobic,
PiStacking, CationPi, PiCation, Anionic, Cationic, VdWContact, XBDonor, XBAcceptor,
MetalDonor, and MetalAcceptor. They use their documented geometric defaults, with
`count=True` to retain distinct contacts. Labels describe the **ligand** role:
HBDonor means ligand donor/protein acceptor; CationPi means ligand cation/protein
ring; MetalAcceptor means a ligand coordinating a receptor metal. See
[ProLIF's docking guide](https://prolif.readthedocs.io/en/v2.2.2/notebooks/docking.html).

## Reports and schemas

Bundles are published under
`For_PostDocking/<profile>/Interactions/<analysis-id>/`. Temporary reports remain
hidden until all three CSVs and the manifest are written. Prior bundles remain
unchanged after cancellation or failure. `latest.json` points to the latest
published bundle with at least one successful candidate; partial failures are
explicitly counted. Database history records all attempts, including fatal errors.

Each CSV has analysis/profile/run identity and source identity as applicable:
`analysis_id`, `receptor_profile_id`, `run_id`, `ligand_id`, `state_id`, `pose_id`,
`source_record_id`, `source_sdf`, `record_index`. Record indices are one-based.
Original SDF names and pose properties are retained separately. Verified coordinate
mapping supplies the original database pose ID and affinity; ambiguous mappings
keep a stable source-record ID, explain ambiguity, and leave affinity `NA`.
Scores are reused, never recalculated. `scoring_protocol` distinguishes protocols.

Residue fields are `protein_chain`, `protein_residue`, `residue_number`, and
`insertion_code`. A blank chain or insertion code means the source has none.

| CSV | Additional columns and row meaning |
| --- | --- |
| `interactions_detailed.csv` | One detected contact per row, including reference contacts: `is_reference`, residue fields, `interaction_type`, `distance_angstrom`, `angle_name`, `angle_degrees`, `geometry_json`, `ligand_atom_indices`, `protein_atom_indices`, `conserved_with_reference`, `analysis_status` |
| `interaction_summary.csv` | One candidate SDF record per row: `source_ligand_name`, `source_pose_id`, `source_properties_json`, `affinity`, `scoring_protocol`, `mapping_status`, `total_interaction_pairs`, `conserved_interaction_pairs`, `lost_reference_pairs`, `new_interaction_pairs`, `reference_retention_percent`, `jaccard_similarity`, `conserved_reference_residues`, `total_reference_residues`, `residue_retention_percent`, `analysis_status`, `error_message` |
| `interaction_comparison.csv` | Long-form candidate/reference union matrix: residue fields, `interaction_type`, `reference_present`, `candidate_present`, `comparison`, `analysis_status` |

Atom indices are zero-based indices into analysis molecules, including hydrogens.
The manifest supplies receptor residue and atom mappings back to original PDB
serials and names. Ligand atom order is preserved, with added hydrogens appended.
Detector geometry is preserved
under its original named fields in `geometry_json`; a single exposed angle is
also copied to `angle_name`/`angle_degrees`. Multiple angles remain in JSON rather
than being ambiguously collapsed. Missing values are `NA`. CSVs use UTF-8 with a
BOM for Excel and are readable with pandas.

Fingerprint metrics count unique residue–interaction-type pairs, regardless of
how many contacts occur. Conserved pairs are the intersection, lost pairs are
reference minus candidate, and new pairs are candidate minus reference.
Reference retention is `100 * intersection / reference pair count`.
Jaccard is `intersection / union`. Residue retention ignores interaction type and
counts contacted reference residues using the complete residue identity.
Zero denominators yield `NA`, including two empty fingerprints for Jaccard.

A valid zero-contact candidate has success status and zero pair counts. Failed
records have failed status, an error, and `NA` metrics. Matrix presence for a
failed candidate is `NA`/unknown, never false/lost. A candidate/reference union
that is empty has no matrix rows; consult the summary for that candidate's status.

`analysis_manifest.json` records versions, settings, input hashes, report hashes,
reference fingerprint, residue mappings, exclusions, coordinate checks, and
success/failure counts. SQLite adds `interaction_exports`, `interaction_runs`,
and `interaction_pose_summaries` without replacing docking data.

## Dependencies and limits

Python dependency constraints are `prolif==2.2.2`, `MDAnalysis>=2.10,<3`, and
`pandas>=2.2,<4`. RDKit and Qt remain Conda-managed; update the existing environment
with the normal installation instructions. Imports are lazy so interaction
dependencies are unnecessary for viewing existing docking results.

The Windows environment also constrains defaults-channel `fontconfig>=2.18.2`
to avoid an `intl-8.dll` symbol conflict between Qt's GLib and RDKit's Cairo.
For an existing environment exhibiting `rdMolDraw2D` or Qt DLL import errors,
the validated repair is `conda install -n moldockpipe-clean defaults::fontconfig
defaults::libglib --force-reinstall`, followed by activation and `python -m pip
check`. Do not mix pip Qt runtime wheels into this environment.

Analysis streams one SDF record at a time and loads the receptor and reference
once per receptor analysis. Source docking coordinates for one run are retained
for mapping. Cancellation is checked between records; a single library operation
must return before cancellation can finish. The first version uses no parallel
chemistry workers. At most 9999 included receptor residues are supported by the
temporary PDB numbering adapter.

These are interaction-pattern similarities, not validated predictions of binding
affinity or biological activity. They depend on chemical typing, protonation,
hydrogen placement, detector geometry, the rigid receptor, and uncertain docking
poses. The existing viewer uses separate approximate geometric annotations and
does not display these ProLIF results.
