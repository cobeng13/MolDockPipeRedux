# Redocking RMSD

MolDockPipe uses PyDockRMSD to compare each saved redocked pose with the crystallographic ligand. The calculation accounts for atom symmetry and uses heavy atoms in the receptor coordinate frame. The receptor summary shows the rank-1 RMSD and the smallest RMSD across the calculated poses, in Å.

Open **Redocking RMSD…** in Receptor Manager to view each pose or use **Calculate RMSD** on a saved run. Receptor packages include the result and its inputs. Runs whose receptor or reference ligand differs from the current profile appear as historical.

Each pose calculation runs in a separate process with a 120-second limit. An engine error or timeout is recorded against that pose. Retry with **Calculate RMSD** after resolving the problem. Successful saved results are reused when their input hashes and engine version match.

## Redocking campaigns

1. In Receptor Manager, select a prepared receptor with a crystallographic
   reference ligand and choose **Queue Redocking Campaign…**. Set the trial limit
   from 1 to 21 and choose the validation settings, including a positive seed
   shared by every trial. PyDockRMSD is required to evaluate the campaign.
2. The queue tries the saved box first, then shifts its center by 2 Å along one
   axis at a time. It repeats the sequence with all three box dimensions reduced
   by 2 Å and then increased by 2 Å. Invalid boxes and duplicate effective
   AutoGrid grids are skipped. The campaign stops at the first trial whose
   affinity-ranked pose 1 has a measured RMSD strictly below 2 Å.
3. In **View Validation Queue → View Details**, review the trial boxes, pose-1
   affinities and RMSDs, and output folders. If no trial meets the target, the
   details show the lowest measured pose-1 RMSD and the apply action is disabled.
   A failed campaign can be retried after resolving the error. On restart,
   interrupted work resumes from the saved trial records and reusable artifacts.

**Apply Winning Box and Settings** saves the successful box and the campaign's
exhaustiveness, modes, energy range, seed, and CPU count as the receptor's editable
production defaults. You can also use it on a campaign whose box was already
applied before settings were carried forward. A change to the saved defaults
makes earlier production results for that receptor historical. Run production
docking with the new defaults to generate current results. Each validation and
production run retains its recorded settings. Vina-GPU uses its own search depth
and computing lanes in place of CPU Vina's exhaustiveness and CPU count.

Campaigns preserve their trial records in the project database, including after
**Clear Completed Validations** clears ordinary queue entries. The queue checks
the receptor and reference input hashes, scoring protocol, and saved box before
resuming a campaign or applying its result. If these inputs changed, start a
new campaign. AD4Zn trials use their own map manifests so exploring trial boxes
leaves the receptor's current map-bundle pointer intact.

## Windows installation

Use Windows x64, the project's Python 3.13 Conda environment, Git, and Visual Studio Build Tools with the Desktop development with C++ workload. From the repository root, activate the desktop environment and run:

```powershell
./tools/build-pydockrmsd.ps1
python -m pip install munkres==1.1.4 ./build/wheels/pydockrmsd-1.0.0-cp313-cp313-win_amd64.whl
python -m pytest -q tests/test_redocking.py
```

The script builds the pinned upstream revision with pip's isolated build environment. PyDockRMSD imports `munkres` at startup, so install it alongside the wheel. The desktop app can still open projects and complete redocking if the wheel is unavailable. Install the wheel and use **Calculate RMSD** to recover the results from saved poses.
