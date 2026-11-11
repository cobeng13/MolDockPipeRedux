# Redocking RMSD

MolDockPipe uses PyDockRMSD to compare each saved redocked pose with the crystallographic ligand. The calculation accounts for atom symmetry and uses heavy atoms in the receptor coordinate frame. The receptor summary shows the rank-1 RMSD and the smallest RMSD across the calculated poses, in Å.

Open **Redocking RMSD…** in Receptor Manager to view each pose or use **Calculate RMSD** on a saved run. Receptor packages include the result and its inputs. Runs whose receptor or reference ligand differs from the current profile appear as historical.

Each pose calculation runs in a separate process with a 120-second limit. An engine error or timeout is recorded against that pose. Retry with **Calculate RMSD** after resolving the problem. Successful saved results are reused when their input hashes and engine version match.

## Windows installation

Use Windows x64, the project's Python 3.13 Conda environment, Git, and Visual Studio Build Tools with the Desktop development with C++ workload. From the repository root, activate the desktop environment and run:

```powershell
./tools/build-pydockrmsd.ps1
python -m pip install munkres==1.1.4 ./build/wheels/pydockrmsd-1.0.0-cp313-cp313-win_amd64.whl
python -m pytest -q tests/test_redocking.py
```

The script builds the pinned upstream revision with pip's isolated build environment. PyDockRMSD imports `munkres` at startup, so install it alongside the wheel. The desktop app can still open projects and complete redocking if the wheel is unavailable. Install the wheel and use **Calculate RMSD** to recover the results from saved poses.
