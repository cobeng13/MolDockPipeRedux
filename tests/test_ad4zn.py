"""AD4Zn contract tests use synthetic maps, never assert scientific validity."""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pytest

from moldockpipe.project import ProjectRepository, DEFAULT_VINA_PROFILE, utc_now
from moldockpipe.receptors import ad4zn
from moldockpipe.services.vina import build_vina_command, VinaDockingBackend


def pdbqt(*types):
    return "ROOT\n" + "\n".join(
        f"ATOM  {i:5d}  C   LIG A   1    {i * 0.1:8.3f}{0:8.3f}{0:8.3f}  1.00  0.00    {0:6.3f} {kind}"
        for i, kind in enumerate(types, 1)) + "\nENDROOT\nTORSDOF 0\n"


def write_gpf(work, required, grid):
    lines = ["parameter_file AD4Zn.dat", "receptor receptor_TZ.pdbqt", "receptor_types C ZN TZ",
             "ligand_types " + " ".join(required), "npts " + " ".join(map(str, grid['npts'])),
             "gridcenter " + " ".join(map(str, grid['center'])), "spacing " + str(grid['spacing']),
             "gridfld receptor_TZ.maps.fld", "elecmap receptor_TZ.e.map", "dsolvmap receptor_TZ.d.map"]
    lines += [f"map receptor_TZ.{t}.map" for t in required]
    # Independent expected values from the published zinc protocol.
    lines += ["nbp_r_eps 0.25 23.2135 12 6 NA TZ", "nbp_r_eps 2.1 3.8453 12 6 OA Zn",
              "nbp_r_eps 2.25 7.5914 12 6 SA Zn", "nbp_r_eps 1.0 0.0 12 6 HD Zn",
              "nbp_r_eps 2.0 0.0060 12 6 NA Zn", "nbp_r_eps 2.0 0.2966 12 6 N Zn"]
    (work / 'receptor_TZ.gpf').write_text('\n'.join(lines) + '\n')


def write_maps(work, required, grid):
    header = ('GRID_PARAMETER_FILE receptor_TZ.gpf\nGRID_DATA_FILE receptor_TZ.maps.fld\n'
              'MACROMOLECULE receptor_TZ.pdbqt\nSPACING ' + str(grid['spacing']) + '\nNELEMENTS '
              + ' '.join(map(str, grid['npts'])) + '\nCENTER ' + ' '.join(map(str, grid['center'])) + '\n')
    names = [f'receptor_TZ.{t}.map' for t in (*required, 'e', 'd')]
    (work / 'receptor_TZ.maps.fld').write_text('\n'.join(names))
    for name in names:
        (work / name).write_text(header + '0\n' * math.prod(n+1 for n in grid['npts']))


@pytest.fixture
def zinc_project(tmp_path, monkeypatch):
    repo = ProjectRepository.create(tmp_path / 'project')
    profile = {**DEFAULT_VINA_PROFILE, 'protocol': 'ad4zn', **{f'size_{a}': 0.75 for a in 'xyz'}}
    receptor = repo.root / profile['receptor']
    receptor.parent.mkdir(parents=True, exist_ok=True)
    receptor.write_text(pdbqt('C', 'Zn'))
    ligand = tmp_path / 'ligand.pdbqt'; ligand.write_text(pdbqt('C', 'OA'))
    tools = {}
    for name in ('python', 'autogrid', 'zinc_pseudo', 'prepare_gpf', 'parameters'):
        path = tmp_path / name; path.write_text(name); tools[name] = path
    tools['parameters'].write_text('\n'.join('atom_par ' + t + ' 0' for t in ('C','OA','N','NA','Br','Cl','A','HD','TZ','Zn')))
    vina = tmp_path / 'vina'; vina.write_text('vina')
    monkeypatch.setattr(ad4zn, 'check_ad4zn_environment', lambda _: ad4zn.AD4ZnEnvironment(tools))
    calls = []
    def run(command, work, records, cancelled=None):
        calls.append(command)
        if command[0] == str(tools['python']) and '-c' in command:
            return '2.7.18'
        if '--version' in command:
            return 'AutoGrid 4.2.7.x.2019-07-11' if command[0] == str(tools['autogrid']) else 'AutoDock Vina v1.2.7'
        if '--help' in command:
            return '--maps --scoring ad4'
        if str(tools['zinc_pseudo']) in command:
            (work / 'receptor_TZ.pdbqt').write_text(pdbqt('C', 'Zn', 'TZ'))
        elif str(tools['prepare_gpf']) in command:
            params = dict(command[i+1].split('=', 1) for i, v in enumerate(command) if v == '-p')
            grid = {'center': list(map(float, params['gridcenter'].split(','))),
                    'npts': list(map(int, params['npts'].split(','))), 'spacing': float(params['spacing'])}
            required = tuple(params['ligand_types'].split(','))
            write_gpf(work, required, grid)
        elif command[0] == str(tools['autogrid']):
            fields = (work / 'receptor_TZ.gpf').read_text().splitlines()
            entries = {line.split()[0]: line.split()[1:] for line in fields}
            grid = {'center': list(map(float, entries['gridcenter'])), 'npts': list(map(int, entries['npts'])),
                    'spacing': float(entries['spacing'][0])}
            write_maps(work, tuple(entries['ligand_types']), grid)
        return ''
    monkeypatch.setattr(ad4zn, '_run', run)
    return repo, profile, ligand, vina, tools, calls


def test_protocol_defaults_and_commands(tmp_path):
    repo = ProjectRepository.create(tmp_path / 'project')
    profile = {**DEFAULT_VINA_PROFILE}; profile.pop('protocol')
    repo.save_receptor_profiles([profile])
    assert repo.get_receptor_profiles()[0]['protocol'] == 'vina'
    with pytest.raises(ValueError, match='Unsupported'):
        repo.save_receptor_profiles([{**profile, 'protocol': 'unknown'}])
    args = dict(executable=Path('vina'), receptor=Path('receptor'), ligand=Path('ligand'), output=Path('out'))
    settings = {'center_x': 2, 'size_x': 20, 'cpu': 1}
    assert build_vina_command(**args, settings=settings) == ['vina', '--receptor', 'receptor', '--ligand', 'ligand', '--out', 'out', '--center_x', '2', '--size_x', '20', '--cpu', '1']
    command = build_vina_command(**args, settings={**settings, 'protocol': 'ad4zn'}, maps=Path('maps'))
    assert command == ['vina', '--maps', 'maps', '--scoring', 'ad4', '--ligand', 'ligand', '--out', 'out', '--cpu', '1']
    with pytest.raises(ValueError, match='requires'):
        build_vina_command(**args, settings={'protocol': 'ad4zn'})


def test_union_preserves_halogen_case_and_all_states(tmp_path):
    paths = [tmp_path / str(i) for i in range(3)]
    for path, types in zip(paths, [('A', 'C', 'OA', 'HD'), ('A', 'C', 'NA', 'Cl'), ('A', 'C', 'Br')]):
        path.write_text(pdbqt(*types))
    assert ad4zn.collect_ligand_atom_types(paths) == ('A', 'Br', 'C', 'Cl', 'HD', 'NA', 'OA')


def test_bundle_reuse_and_portable_metadata(zinc_project):
    repo, profile, ligand, vina, _, calls = zinc_project
    prefix, record = ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    count = len(calls)
    prefix2, record2 = ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    assert len(calls) == count and prefix == prefix2 and record == record2
    assert not Path(record['maps_prefix']).is_absolute()
    assert record['protocol'] == 'ad4zn'
    assert {'AD4Zn.dat', 'receptor_TZ.pdbqt', 'receptor_TZ.gpf', 'receptor_TZ.C.map'} <= record['hashes'].keys()
    assert (prefix.parent / 'manifest.json').is_file()


@pytest.mark.parametrize('change', ['parameters', 'zinc_pseudo', 'autogrid', 'vina', 'receptor', 'box', 'ligand_types', 'map'])
def test_scientific_changes_invalidate_bundle(zinc_project, change):
    repo, profile, ligand, vina, tools, calls = zinc_project
    prefix, before = ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    if change in tools:
        tools[change].write_text(tools[change].read_text() + '\n# changed')
    elif change == 'vina': vina.write_text('new vina')
    elif change == 'receptor': (repo.root / profile['receptor']).write_text(pdbqt('C', 'Zn', 'N'))
    elif change == 'box': profile['center_x'] = 1
    elif change == 'ligand_types': ligand.write_text(pdbqt('C', 'Br'))
    elif change == 'map': Path(str(prefix) + '.C.map').write_text('corrupt')
    prefix2, after = ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    assert prefix2 != prefix
    assert (prefix.parent / 'manifest.json').is_file()  # Keep old run artifacts.
    if change != 'map': assert after['request'] != before['request']


@pytest.mark.parametrize('failure', ['no_zinc', 'helper_fails', 'no_tz', 'gpf_parameter', 'gpf_potential', 'autogrid_fails', 'map_missing', 'map_truncated'])
def test_failures_do_not_publish_current_bundle(zinc_project, monkeypatch, failure):
    repo, profile, ligand, vina, tools, _ = zinc_project
    run = ad4zn._run
    if failure == 'no_zinc': (repo.root / profile['receptor']).write_text(pdbqt('C', 'N'))
    def broken(command, work, records, cancelled=None):
        result = run(command, work, records, cancelled)
        if str(tools['zinc_pseudo']) in command:
            if failure == 'helper_fails': raise RuntimeError('helper failed')
            if failure == 'no_tz': (work / 'receptor_TZ.pdbqt').write_text(pdbqt('C', 'Zn'))
        if str(tools['prepare_gpf']) in command and failure.startswith('gpf_'):
            path = work / 'receptor_TZ.gpf'
            text = path.read_text()
            text = text.replace('parameter_file AD4Zn.dat', 'parameter_file wrong.dat') if failure == 'gpf_parameter' else text.replace('23.2135', '1.0')
            path.write_text(text)
        if command[0] == str(tools['autogrid']) and '-p' in command:
            if failure == 'autogrid_fails': raise RuntimeError('AutoGrid failed')
            if failure == 'map_missing': (work / 'receptor_TZ.OA.map').unlink()
            if failure == 'map_truncated': (work / 'receptor_TZ.C.map').write_text('invalid')
        return result
    monkeypatch.setattr(ad4zn, '_run', broken)
    with pytest.raises((ValueError, RuntimeError)):
        ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    assert not (repo.root / 'inputs/receptors/default/ad4zn/current.json').exists()
    if failure != 'no_zinc': assert list(repo.root.glob('inputs/receptors/default/ad4zn/bundle-*/failure.txt'))


def test_missing_dependencies_and_invalid_override(tmp_path, monkeypatch):
    monkeypatch.setattr(ad4zn, 'application_root', lambda: tmp_path)
    monkeypatch.setattr(ad4zn.shutil, 'which', lambda _: None)
    for name in ('PYTHON', 'AUTOGRID', 'ZINC_PSEUDO', 'PREPARE_GPF', 'PARAMETERS'):
        monkeypatch.setenv('MOLDOCKPIPE_AD4ZN_' + name, str(tmp_path / 'missing'))
    environment = ad4zn.check_ad4zn_environment(tmp_path)
    assert len(environment.missing) == 5
    with pytest.raises(FileNotFoundError, match='AD4Zn.dat'):
        environment.require()


def test_old_autogrid_is_rejected(zinc_project, monkeypatch):
    repo, profile, ligand, vina, tools, _ = zinc_project
    run = ad4zn._run
    monkeypatch.setattr(ad4zn, '_run', lambda command, *args: 'AutoGrid 4.2.6' if '--version' in command else run(command, *args))
    with pytest.raises(ValueError, match='4.2.7'):
        ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)


def test_native_vina_ad4_smoke(tmp_path):
    """Actual bundled Vina command/parsing check against synthetic zero-energy maps."""
    from moldockpipe.services.vina import find_vina_executable
    try:
        executable = find_vina_executable(tmp_path)
    except FileNotFoundError:
        pytest.skip('No native Vina executable on this host')
    grid = {'npts': [40, 40, 40], 'center': [0, 0, 0], 'spacing': 0.375}
    write_maps(tmp_path, ('C',), grid)
    ligand = tmp_path / 'ligand.pdbqt'; ligand.write_text(pdbqt('C', 'C'))
    receptor = tmp_path / 'receptor.pdbqt'; receptor.write_text(pdbqt('C', 'Zn'))
    result = VinaDockingBackend().run(executable=executable, receptor=receptor, ligand=ligand,
        output=tmp_path / 'output.pdbqt', log=tmp_path / 'vina.log', maps=tmp_path / 'receptor_TZ',
        settings={'protocol': 'ad4zn', 'cpu': 1, 'exhaustiveness': 1, 'num_modes': 1, 'seed': 42})
    assert result.return_code == 0, result.stderr
    assert result.poses
    assert '--receptor' not in result.command


def test_mixed_campaign_reuse_and_ligand_hash(zinc_project, monkeypatch, tmp_path):
    from moldockpipe.pipeline import PipelineRunner
    from moldockpipe.services.vina import VinaResult
    from moldockpipe.csv_exports import export_manifest_csv, export_leaderboard_csv
    repo, zinc, _, vina, tools, _ = zinc_project
    csv_path = tmp_path / 'ligands.csv'; csv_path.write_text('id,smiles\nlig1,CCO\n')
    repo.import_ligands_csv(csv_path)
    sdf = repo.root / 'state.sdf'; sdf.write_text('SDF')
    artifact = repo.add_artifact(sdf, 'state_sdf', 'molscrub')
    with repo.connection() as conn:
        conn.execute("""INSERT INTO molecular_states
            (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,tautomer_index,protomer_index,
             state_structure_hash,generation_fingerprint,status,enumeration_truncated,created_at,active)
             VALUES ('state-1','lig1','CCO','CCO',0,0,0,'structure','generation','prepared',0,?,1)""", (utc_now(),))
        conn.execute("INSERT INTO conformers VALUES ('conf','state-1',0,?,'pdbqt_ready',?)", (artifact, utc_now()))
    ligand = repo.root / 'artifacts/pdbqt/lig1/state-1/ligand.pdbqt'
    ligand.parent.mkdir(parents=True); ligand.write_text(pdbqt('C', 'OA'))
    standard = {**zinc, 'id': 'standard', 'name': 'Standard', 'protocol': 'vina'}
    repo.save_receptor_profiles([standard, zinc])
    monkeypatch.setattr('moldockpipe.pipeline.find_vina_executable', lambda _: vina)
    docked = []
    def dock(self, **kwargs):
        docked.append(kwargs)
        assert ('maps' in kwargs) == (kwargs['settings'].get('protocol') == 'ad4zn')
        kwargs['output'].parent.mkdir(parents=True)
        kwargs['output'].write_text('REMARK VINA RESULT: -7.000 0.000 0.000\n')
        kwargs['log'].write_text('ok')
        return VinaResult(0, [(1, -7, 0, 0)], '', '')
    monkeypatch.setattr('moldockpipe.pipeline.VinaDockingBackend.run', dock)
    runner = PipelineRunner(repo)
    assert runner.run_vina() == (2, 0)
    assert len(docked) == 2
    assert runner.run_vina() == (2, 0)
    assert len(docked) == 2
    # Same atom-type union, different prepared ligand bytes: zinc docking must rerun.
    ligand.write_text(pdbqt('C', 'OA') + 'REMARK changed\n')
    assert runner.run_vina() == (2, 0)
    assert len(docked) == 3 and docked[-1]['settings']['protocol'] == 'ad4zn'
    export_manifest_csv(repo); export_leaderboard_csv(repo)
    assert 'ad4zn' in (repo.root / 'exports/leaderboard.csv').read_text()
    tools['parameters'].write_text(tools['parameters'].read_text() + '\n# new parameters')
    assert runner.run_vina() == (2, 0)
    assert len(docked) == 4
    # Missing zinc resources fails only that profile; standard results still reuse.
    monkeypatch.setattr(ad4zn, 'check_ad4zn_environment', lambda _: ad4zn.AD4ZnEnvironment({'autogrid': None}))
    assert runner.run_vina() == (1, 1)
    assert len(docked) == 4


@pytest.mark.parametrize('missing_zinc', [False, True])
def test_meeko_then_zinc_preparation_and_failure_retention(zinc_project, monkeypatch, tmp_path, missing_zinc):
    from test_receptors import PDB
    from moldockpipe.receptors.models import ReceptorPreparationPlan, ResidueKey
    from moldockpipe.receptors.preparation import prepare_receptor
    from types import SimpleNamespace
    repo, _, _, _, tools, _ = zinc_project
    monkeypatch.setattr('moldockpipe.receptors.preparation.check_ad4zn_environment', lambda _: ad4zn.AD4ZnEnvironment(tools))
    source = tmp_path / 'source.pdb'; source.write_text(PDB)
    plan = ReceptorPreparationPlan('prepared', 'Prepared', source, 0, ('A',), None,
        (ResidueKey('A', 501, '', 'HOH'), ResidueKey('A', 401, '', 'LIG')),
        (ResidueKey('A', 601, '', 'ZN'),), (15,16,17), (20,21,22), 'manual', protocol='ad4zn')
    def meeko(command, cwd, **kwargs):
        work = Path(cwd)
        (work / 'receptor.pdbqt').write_text(pdbqt('C', 'N') if missing_zinc else pdbqt('C', 'Zn'))
        return SimpleNamespace(returncode=0, stdout='prepared', stderr='')
    monkeypatch.setattr('moldockpipe.receptors.preparation.subprocess.run', meeko)
    if missing_zinc:
        with pytest.raises(RuntimeError, match='preserve all retained Zn'):
            prepare_receptor(repo.root, plan)
        assert not (repo.root / 'inputs/receptors/prepared').exists()
        assert list(repo.root.glob('logs/receptor_preparation/prepared/*/receptor.pdbqt'))
    else:
        folder = prepare_receptor(repo.root, plan)
        assert 'TZ' in ad4zn.atom_types(folder / 'ad4zn/receptor_preparation/receptor_TZ.pdbqt')
        report = json.loads((folder / 'preparation_report.json').read_text())
        assert report['protocol'] == 'ad4zn' and report['ad4zn']['tz_receptor_sha256']


@pytest.mark.parametrize('content', ['../../../data/AD4Zn.dat', 'atom_par C 0'])
def test_invalid_parameter_file_fails_before_tools_run(zinc_project, content):
    repo, profile, ligand, vina, tools, calls = zinc_project
    tools['parameters'].write_text(content)
    with pytest.raises(FileNotFoundError, match='Invalid AD4Zn.dat'):
        ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    assert not calls


def test_unknown_ligand_type_fails_before_map_generation(zinc_project):
    repo, profile, ligand, vina, _, calls = zinc_project
    ligand.write_text(pdbqt('C', 'G0'))
    with pytest.raises(ValueError, match='G0'):
        ad4zn.ensure_ad4zn_maps(repo, profile, [ligand], vina)
    assert not calls


def test_subprocess_failure_preserves_command_and_output(tmp_path):
    import sys
    with pytest.raises(RuntimeError, match='tool failed'):
        ad4zn._run([sys.executable, '-c', "import sys; print('diagnostic', file=sys.stderr); sys.exit(2)"], tmp_path, [])
    record = json.loads((tmp_path / 'commands.json').read_text())[0]
    assert record['return_code'] == 2 and 'diagnostic' in record['stderr']


def test_subprocess_cancellation_preserves_logs(tmp_path):
    import sys
    with pytest.raises(InterruptedError, match='cancelled'):
        ad4zn._run([sys.executable, '-c', 'import time; time.sleep(5)'], tmp_path, [], lambda: True)
    assert (tmp_path / 'commands.json').is_file()


def test_adfr_installer_layout_is_discovered(tmp_path, monkeypatch):
    import os
    installation = tmp_path / 'ADFRsuite-1.0'
    (installation / 'bin').mkdir(parents=True)
    interpreter = installation / ('python.exe' if os.name == 'nt' else 'bin/pythonsh')
    autogrid = installation / 'bin' / ('autogrid4.exe' if os.name == 'nt' else 'autogrid4')
    for path in (interpreter, autogrid):
        path.write_text('tool'); path.chmod(0o755)
    monkeypatch.setattr(ad4zn, 'adfr_installations', lambda: [installation])
    monkeypatch.setattr(ad4zn, 'application_root', lambda: tmp_path)
    monkeypatch.setattr(ad4zn.shutil, 'which', lambda _: None)
    for name in ('PYTHON', 'AUTOGRID'):
        monkeypatch.delenv('MOLDOCKPIPE_AD4ZN_' + name, raising=False)
    environment = ad4zn.check_ad4zn_environment(tmp_path)
    assert environment.components['python'] == interpreter.resolve()
    assert environment.components['autogrid'] == autogrid.resolve()
    monkeypatch.setenv('MOLDOCKPIPE_AD4ZN_PYTHON', str(tmp_path / 'missing'))
    assert ad4zn.check_ad4zn_environment(tmp_path).components['python'] is None


def test_parameter_placeholder_is_reported_in_ui_preflight(tmp_path):
    parameter = tmp_path / 'AD4Zn.dat'
    parameter.write_text('../../../data/AD4Zn.dat')
    environment = ad4zn.AD4ZnEnvironment({'parameters': parameter})
    assert 'symbolic-link placeholder' in environment.problems[0]
    parameter.write_text('atom_par TZ 0 0 0')
    assert environment.problems == []
