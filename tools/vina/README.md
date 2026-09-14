# AutoDock Vina Binary

Place the platform-compatible AutoDock Vina executable in this application
folder. The application recognizes:

- `vina.exe`
- `vina_1.2.7_win.exe`
- `vina` on Linux and macOS

Place `vina_split.exe` (Windows) or `vina_split` (Linux/macOS) here as well.
Project-local `tools/vina/` folders and executables on `PATH` are also supported.
Set `MOLDOCKPIPE_VINA` or `MOLDOCKPIPE_VINA_SPLIT` to select an explicit tool.

The receptor is expected at `inputs/receptor_prepared.pdbqt` in the selected project.
Each run writes a readable `vina.log.txt`; command/settings provenance is also
stored in the project database.

Meeko's `mk_export.py` is resolved from the active Python environment as
`python -m meeko.cli.mk_export`.
Development Vina binaries may be committed in this folder. Add compatible
replacement binaries here when testing a different Vina release.
