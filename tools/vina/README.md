# AutoDock Vina Binary

Place the Windows AutoDock Vina executable in this fixed application folder.
The application checks these names in order:

- `vina.exe`
- `vina_1.2.7_win.exe`

Place `vina_split.exe` here as well (the application also checks `tools/vina_split.exe`).

The receptor is expected at `inputs/receptor_prepared.pdbqt` in the selected project.
Each run writes a readable `vina.log.txt`; command/settings provenance is also
stored in the project database.

Meeko's `mk_export.py` is resolved from the active Python environment as
`python -m meeko.cli.mk_export`.
Development Vina binaries may be committed in this folder. Add compatible
replacement binaries here when testing a different Vina release.
