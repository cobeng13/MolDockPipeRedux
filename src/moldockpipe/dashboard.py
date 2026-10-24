"""Database and filesystem snapshot collection; contains no Qt widgets."""
from .services.postdock import receptor_export_filename


def _export_status(repo, profiles, settings, cancelled):
    """Stream export records so idle refreshes do not retain the whole campaign."""
    ids = [str(profile['id']) for profile in profiles]
    names = {str(profile['id']): str(profile.get('name', profile['id'])) for profile in profiles}
    root = repo.root / 'For_PostDocking'
    mode = settings.get('postdock', {}).get('mode', 'split_and_sdf')
    folders = ('PDBQTs',) if mode == 'split_only' else ('SDF',) if mode == 'sdf_only' else ('SDF', 'PDBQTs')
    seen = set()
    with repo.connection() as conn:
        conn.set_progress_handler(lambda: int(cancelled()), 5000)
        rows = conn.execute(f"""SELECT d.receptor_profile_id,s.parent_id,s.state_id,d.run_id
            FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
            JOIN parent_ligands p USING(parent_id)
            WHERE s.active=1 AND d.status='completed' AND d.is_current=1
            AND d.receptor_profile_id IN ({','.join('?' for _ in ids)})""", ids)
        for profile, parent, state, run in rows:
            if cancelled():
                raise InterruptedError('Dashboard refresh superseded')
            seen.add(profile)
            for folder in folders:
                target = root / profile / folder / parent / state / run
                suffix = '*.sdf' if folder == 'SDF' else '*.pdbqt'
                if not target.is_dir() or not any(target.glob(suffix)):
                    return 'Ready'
        validated = {row[0] for row in conn.execute(
            "SELECT DISTINCT receptor_profile_id FROM redocking_runs WHERE status='ARTIFACTS_READY'")}
    for profile in seen:
        if cancelled():
            raise InterruptedError('Dashboard refresh superseded')
        if not (root / profile / receptor_export_filename(names[profile])).is_file():
            return 'Ready'
        if profile in validated and not (root / profile / 'reference_ligand_heavy.mol2').is_file():
            return 'Ready'
    return 'Complete' if seen else 'Pending'


def collect_dashboard(repo, running=False, cancelled=lambda: False):
    settings = repo.get_settings()
    profiles = [profile for profile in repo.get_receptor_profiles() if profile.get("enabled")]
    profile_ids = [str(profile["id"]) for profile in profiles]
    placeholders = ",".join("?" for _ in profile_ids)
    structure_failures = repo.molscrub_failures(cancelled=cancelled)
    with repo.connection() as conn:
        conn.set_progress_handler(lambda: int(cancelled()), 5000)
        parents = conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0]
        screened = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE s.active=1 AND s.status IN ('completed','failed')").fetchone()[0]
        states = conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0]
        structure_passed = conn.execute("SELECT COUNT(DISTINCT parent_id) FROM molecular_states WHERE active=1").fetchone()[0]
        structure_total = conn.execute("""SELECT COUNT(*) FROM parent_ligands p
            LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
            WHERE s.decision NOT IN ('fail') OR ?='annotate_only'""",
            (settings.get("screening", {}).get("policy", "annotate_only"),)).fetchone()[0]
        prepared = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
        conformers = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1").fetchone()[0]
        meeko_terminal = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status IN ('pdbqt_ready','pdbqt_failed')").fetchone()[0]
        prepared_states = conn.execute("SELECT COUNT(DISTINCT s.state_id) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
        if profile_ids:
            docked = conn.execute(f"SELECT COUNT(*) FROM docking_runs d JOIN molecular_states s USING(state_id) WHERE s.active=1 AND d.status='completed' AND d.is_current=1 AND d.receptor_profile_id IN ({placeholders})", profile_ids).fetchone()[0]
            docking_terminal = conn.execute(f"SELECT COUNT(*) FROM docking_runs d JOIN molecular_states s USING(state_id) WHERE s.active=1 AND d.status IN ('completed','failed','interrupted') AND d.is_current=1 AND d.receptor_profile_id IN ({placeholders})", profile_ids).fetchone()[0]
        else:
            docked = docking_terminal = 0
    with repo.connection() as conn:
        conn.set_progress_handler(lambda: int(cancelled()), 5000)
        passed = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE s.active=1 AND s.decision != 'fail'").fetchone()[0]
        failed = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE s.active=1 AND s.decision = 'fail'").fetchone()[0]
        structure_failed = len(structure_failures)
        screening_failed_ids = {row[0] for row in conn.execute("SELECT parent_id FROM screening_results WHERE active=1 AND decision='fail'")}
        failed_ligands = len(screening_failed_ids | {row["parent_id"] for row in structure_failures})
        pdbqt_failed = conn.execute("SELECT COUNT(DISTINCT s.parent_id) FROM molecular_states s JOIN parent_ligands p USING(parent_id) WHERE s.active=1 AND s.status='pdbqt_failed'").fetchone()[0]
        pdbqt_passed = conn.execute("SELECT COUNT(DISTINCT s.parent_id) FROM molecular_states s JOIN parent_ligands p USING(parent_id) JOIN conformers c USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
        if profile_ids:
            docking_failed = conn.execute(f"""SELECT COUNT(DISTINCT s.parent_id || ':' || d.receptor_profile_id)
                FROM docking_runs d JOIN molecular_states s USING(state_id) JOIN parent_ligands p USING(parent_id)
                WHERE s.active=1 AND d.status IN ('failed','interrupted') AND d.is_current=1
                AND d.receptor_profile_id IN ({placeholders})""", profile_ids).fetchone()[0]
            docking_passed = conn.execute(f"""SELECT COUNT(DISTINCT s.parent_id || ':' || d.receptor_profile_id)
                FROM docking_runs d JOIN molecular_states s USING(state_id) JOIN parent_ligands p USING(parent_id)
                WHERE s.active=1 AND d.status='completed' AND d.is_current=1
                AND d.receptor_profile_id IN ({placeholders})""", profile_ids).fetchone()[0]
        else:
            docking_failed = docking_passed = 0
        stage_runs = conn.execute("SELECT COUNT(*) FROM stage_runs").fetchone()[0]
    has_docking_exports = bool(docking_passed)
    export_status = 'Ready' if has_docking_exports else 'Pending'
    if not running and has_docking_exports:
        export_status = _export_status(repo, profiles, settings, cancelled)
    expected_docking = prepared_states * len(profile_ids)
    return {
        "parents": parents,
        "screened": screened,
        "states": states,
        "prepared": prepared,
        "docked": docked,
        "structure_total": structure_total,
        "structure_passed": structure_passed,
        "structure_failed": structure_failed,
        "conformers": conformers,
        "meeko_terminal": meeko_terminal,
        "expected_docking": expected_docking,
        "docking_terminal": docking_terminal,
        "passed": passed,
        "failed": failed,
        "failed_ligands": failed_ligands,
        "pdbqt_failed": pdbqt_failed,
        "pdbqt_passed": pdbqt_passed,
        "docking_failed": docking_failed,
        "docking_passed": docking_passed,
        "stage_runs": stage_runs,
        "structure_failures": structure_failures[:6],
        "has_docking_exports": has_docking_exports,
        "export_status": export_status,
    }
