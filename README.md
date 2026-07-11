# MolDockPipe Redux

PyQt6 desktop workflow for traceable virtual screening projects. It keeps submitted
ligands, generated molecular states, docking runs, poses, and scientific artifacts
separate. The legacy scripts remain in `Reference/` as behavioral reference only.

## Development

RDKit is best installed from conda-forge. Install the application into that environment:

```powershell
conda create -n moldockpipe -c conda-forge python=3.13 rdkit pyqt=6 pip
conda activate moldockpipe
pip install -e ".[chemistry,dev]"
pytest
moldockpipe
```

Create a project from the UI or with `moldockpipe create C:\path\to\project`.
The application accepts a CSV with `id,smiles,notes,params_json` and a prepared
receptor PDBQT. Vina is configured per project and is never bundled or downloaded.
