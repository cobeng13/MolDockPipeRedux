# AutoDock4Zn support (feature branch)

AD4Zn is an optional receptor-specific protocol. Existing profiles default to
standard Vina, and standard Vina does not require any AD4Zn tools. The ligand
Screening -> MolScrub -> Meeko pipeline is unchanged. This branch implements
software integration; it does not establish scientific validity for MMP-13.

## Windows setup

Install the ADFR suite separately from its official distribution. Its compatible
Python runtime must import MolKit and AutoDockTools; the legacy official grid
helper also requires Python 2's string.split. Do not replace the application's
Conda Python or install these legacy packages into its environment.

Supply these five components explicitly (the application never downloads them):

| Component | Environment override |
|---|---|
| ADFR-compatible interpreter / pythonsh | MOLDOCKPIPE_AD4ZN_PYTHON |
| AutoGrid >= 4.2.7; official recommended build 4.2.7.x.2019-07-11 or newer | MOLDOCKPIPE_AD4ZN_AUTOGRID |
| zinc_pseudo.py | MOLDOCKPIPE_AD4ZN_ZINC_PSEUDO |
| prepare_gpf4zn.py | MOLDOCKPIPE_AD4ZN_PREPARE_GPF |
| AD4Zn.dat | MOLDOCKPIPE_AD4ZN_PARAMETERS |

Each override takes the full path to one file. For example, in PowerShell:

```powershell
$env:MOLDOCKPIPE_AD4ZN_PYTHON = 'C:/path/to/ADFR/python.exe'
$env:MOLDOCKPIPE_AD4ZN_AUTOGRID = 'C:/path/to/ADFR/bin/autogrid4.exe'
$env:MOLDOCKPIPE_AD4ZN_ZINC_PSEUDO = 'C:/path/to/zinc_pseudo.py'
$env:MOLDOCKPIPE_AD4ZN_PREPARE_GPF = 'C:/path/to/prepare_gpf4zn.py'
$env:MOLDOCKPIPE_AD4ZN_PARAMETERS = 'C:/path/to/AD4Zn.dat'
# From the feature worktree, select its source without changing the stable
# checkout's editable installation:
$env:PYTHONPATH = (Join-Path (Get-Location) 'src')
python -m moldockpipe ui
```

Use the interpreter from the installed ADFR environment that passes the import
probe, not an arbitrary Python executable. Overrides are host configuration and
are not written into portable project profiles. Invalid overrides fail without
falling back. Without overrides, discovery checks application tools/ad4zn,
project tools/ad4zn, then PATH. Conventional names are pythonsh(.exe),
autogrid4(.exe), zinc_pseudo.py, prepare_gpf4zn.py, and AD4Zn.dat.
Other OSes can use compatible native tools through the same overrides; their
ADFR runtimes have not been validated by this Windows-first integration.

Official sources:

- [Zinc docking tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html)
- [ADFR downloads](https://ccsb.scripps.edu/adfr/downloads/)
- [Helper scripts](https://github.com/ccsb-scripps/AutoDock-Vina/tree/develop/example/autodock_scripts)
- [Actual AD4Zn.dat](https://github.com/ccsb-scripps/AutoDock-Vina/blob/develop/data/AD4Zn.dat)

The tutorial data directory's AD4Zn.dat can be a symbolic-link placeholder.
Use the actual parameter file above. Source versions and hashes are recorded
for every generated bundle. AutoGrid 4.2.6 is rejected because zinc pair-potential
handling changed in 4.2.7. A real zinc-toolchain run is still required after setup.

## Using the protocol

Choose **AutoDock4Zn / Vina AD4 scoring** in receptor preparation or edit an
existing prepared receptor through Receptors. Zn must be explicitly retained;
the app never switches protocols based only on metal detection. New preparations
run Meeko first and verify that retained Zn survives before adding TZ pseudoatoms.
Imported prepared PDBQTs are checked and processed when maps are generated.

The profile persists `protocol: ad4zn`; optional `ad4zn.spacing` defaults to
0.375 Angstrom. Existing search-box dimensions remain in Angstrom. AutoGrid
interval counts are rounded up to even values to cover that box, with a maximum
of 126 intervals per axis. Centers and spacing are rounded to 0.001 Angstrom,
matching map-header precision. Invalid or oversized grids fail before docking.

Maps are generated just before docking, after inspecting the complete active
prepared ligand-state library. Reference redocking generates maps for its prepared
reference ligand. Missing/unsupported ligand atom types fail before any AD4Zn
ligand is docked. Ligand chemistry is not changed to work around missing maps.

AD4Zn calls Vina with `--maps <prefix> --scoring ad4`; it never supplies
`--receptor` or box arguments. Standard Vina keeps its existing command. The
Receptors dialog shows missing tool files, Zn presence, and map status. Runtime,
version, GPF, zinc-potential, and full map checks happen in the worker before
AD4Zn docking starts. One failed AD4Zn profile does not stop standard profiles
in a screening campaign.

## Artifacts, reuse, and compatibility

No database schema migration is needed. JSON settings/provenance hold protocol
metadata. Generated bundles live at
`inputs/receptors/<id>/ad4zn/bundle-<id>/`; `current.json` points to the last
validated bundle. Each bundle includes the base/TZ receptors, parameter file,
GPF, maps/field descriptor, command output, versions, hashes, and manifest.
Failures retain their directory and logs and never replace a valid current bundle.

Reuse checks scientific input hashes and validates the complete cached outputs.
Changed tools, parameters, receptor, grid, required atom types, or map contents
invalidate the bundle. AD4Zn docking additionally includes the prepared ligand
PDBQT hash and search settings. Standard Vina fingerprint defaults remain
compatible with existing projects. Previous bundles remain available to audit
historical docking runs. Map generation uses immutable folders so concurrent
workers do not overwrite maps being consumed by another run.

CSV manifests and leaderboards append a `protocol` column; existing columns
remain intact. Consumers that require an exact column count must accept this
new column. Results identify the scoring protocol, and redocking SDFs carry a
DOCKING_PROTOCOL property. AD4Zn post-docking packages contain per-run protocol
JSON metadata; the exported physical receptor PDB does not turn TZ pseudoatoms
into real atoms. The full reproducible maps remain in the portable project.
Do not compare standard Vina and AD4Zn scores as one scoring scale.

## Acceptance and limitations

From this worktree, set PYTHONPATH to its src directory as shown above, then
run `python -m pytest -q -p no:cacheprovider` in the project Conda environment.
Tests cover failures, backward compatibility, mixed campaigns, atom-type union,
cache invalidation, redocking, and actual bundled Vina AD4 execution on synthetic
maps. Synthetic maps are a command/parsing test, not a force-field validation.

After installing the external tools, run the official zinc tutorial as a toolchain
check, then test the selected MMP-13 complex manually: retain the intended metals,
inspect TZ geometry, redock the reference ligand, inspect metal coordination and
pose recovery, and assess RMSD externally. MolDockPipe still does not compute
DockRMSD itself. Test standard Vina in a separate project as a release regression.
Keep this feature branch separate from main until that acceptance review.
