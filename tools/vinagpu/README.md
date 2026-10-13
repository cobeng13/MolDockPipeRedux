# Optional Vina-GPU runtime

Required files for the supplied Windows runtime:

```text
tools/vinagpu/
  Vina-GPU+.exe
  Kernel1_Opt.bin
  Kernel2_Opt.bin
```

Automatic discovery checks the application's `tools/vinagpu/` folder first,
then the project's `tools/vinagpu/` folder. It does not search `tools/vina/`
or PATH. An explicit `MOLDOCKPIPE_VINAGPU` override remains supported.
`VinaGPUConfig.txt` is unnecessary: MolDockPipe passes settings directly.
`Vina-GPU+_K.exe` and `AutoDock-GPU.exe` are unnecessary for this adapter.

The Windows `Vina-GPU+.exe` and matching OpenCL kernels in this folder are
copies of the legacy runtime in `Reference/`. They were verified with a
single-ligand smoke test on an NVIDIA GeForce RTX 3070 Ti. Kernel binaries
are device/driver dependent; this check does not establish compatibility
with other GPUs or operating systems.

Select **Settings > Docking > Vina-GPU** to use this engine for screening.
CPU Vina remains the default and remains the engine for receptor validation.
Use `MOLDOCKPIPE_VINAGPU` to select another compatible executable; supply
`Kernel1_Opt.bin` and `Kernel2_Opt.bin` beside it. The backend copies the kernels
to each run folder and supplies all inputs and outputs as absolute paths.
The kernel-compiling `_K` variant is not automatically selected.
