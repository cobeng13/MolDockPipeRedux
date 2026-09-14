# AutoDock Vina Binary

Place the platform-compatible AutoDock Vina executable in this application
folder. The application recognizes:

- `vina.exe`
- `vina_1.2.7_win.exe`
- `vina_1.2.7_linux_x86_64` on Linux x86_64
- `vina_1.2.7_mac_x86_64` on Intel macOS
- `vina_1.2.7_mac_aarch64` on Apple Silicon macOS
- `vina` on Linux and macOS as a generic fallback

The corresponding `vina_split_1.2.7_<OS>_<architecture>` binary is selected
using the same OS and CPU detection. Windows uses `vina_split.exe`; generic
`vina_split` remains supported on Linux/macOS. No Linux ARM binary is bundled.
Unix binaries are tracked with executable permissions. If copying files outside
Git strips those permissions, restore them with `chmod +x <binary>` on that host.
Native macOS/Linux execution still needs to be verified on those systems.
Project-local `tools/vina/` folders and executables on `PATH` are also supported.
Set `MOLDOCKPIPE_VINA` or `MOLDOCKPIPE_VINA_SPLIT` to select an explicit tool.

The receptor is expected at `inputs/receptor_prepared.pdbqt` in the selected project.
Each run writes a readable `vina.log.txt`; command/settings provenance is also
stored in the project database.

Meeko's `mk_export.py` is resolved from the active Python environment as
`python -m meeko.cli.mk_export`.
Development Vina binaries may be committed in this folder. Add compatible
replacement binaries here when testing a different Vina release.
