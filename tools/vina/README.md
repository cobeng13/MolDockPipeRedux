# AutoDock Vina Binary

Place the Windows AutoDock Vina executable in this fixed application folder.
The application checks these names in order:

- `vina.exe`
- `vina_1.2.7_win.exe`

The receptor is expected at `inputs/receptor_prepared.pdbqt` in the selected project.
Each run writes a readable `vina.log.txt`; command/settings provenance is also
stored in the project database.
Vina is user-supplied and is not committed to the repository.
