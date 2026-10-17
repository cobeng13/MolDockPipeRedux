"""Compare optimized planning with the previous per-item SQL/checks.

Reference decisions deliberately use separate queries and the original
predicates, rather than the optimized cache loaders.
"""
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

import moldockpipe.pipeline as pipeline
from moldockpipe.fingerprints import fingerprint
from moldockpipe.pipeline import PipelineRunner, _resolve_destination
from moldockpipe.project import ProjectRepository


@pytest.fixture
def cache_project(tmp_path):
    pytest.importorskip("rdkit")
    molscrub = pytest.importorskip("molscrub")
    repo = ProjectRepository.create(tmp_path / "project")
    config = repo.get_settings()
    config["molscrub"].update(workers=1, timeout_seconds=15)
    repo.update_settings({"molscrub": config["molscrub"]})
    sdf = repo.root / "input.sdf"
    sdf.write_text("unchanged input coordinates", encoding="ascii")
    artifact_id = repo.add_artifact(sdf, "state_sdf", "molscrub")
    return SimpleNamespace(repo=repo, sdf=sdf, artifact_id=artifact_id,
                           version=getattr(molscrub, "__version__", "unknown"),
                           settings=config["molscrub"])


def seed_parent(p, parent_id, *, previous="completed", state=True, generation="current",
                artifact=True, active=True, conformer=True):
    cache = {k: v for k, v in p.settings.items() if k != "workers"}
    generation_settings = {k: v for k, v in cache.items() if k != "timeout_seconds"}
    structure = f"structure-{parent_id}"
    state_fp = fingerprint(settings=cache if generation == "legacy" else generation_settings,
                           inputs={"state": structure}, tool_version=p.version)
    if generation == "stale":
        state_fp = "old-chemistry"
    state_id = parent_id + "-state"
    with p.repo.connection() as conn:
        conn.execute("""INSERT INTO parent_ligands(parent_id,source_smiles,parse_status,created_at)
            VALUES (?,'CCO','pending','now')""", (parent_id,))
        if previous:
            fp = fingerprint(settings=cache, inputs={"smiles": "CCO"}, tool_version=p.version)
            conn.execute("INSERT INTO stage_runs VALUES (?,?,?,?,?,?,?,?,?)",
                         ("run-" + parent_id, "molscrub", "parent", parent_id,
                          "stale" if previous == "stale_failed" else fp,
                          "failed" if previous == "stale_failed" else previous, "now", "now", "reason"))
        if state:
            conn.execute("""INSERT INTO molecular_states
                (state_id,parent_id,state_smiles,state_isomeric_smiles,formal_charge,
                 state_structure_hash,generation_fingerprint,status,created_at,active)
                VALUES (?,?,'CCO','CCO',0,?,?,'sdf_ready','now',?)""",
                         (state_id, parent_id, structure, state_fp, int(active)))
            if conformer:
                conn.execute("INSERT INTO conformers VALUES (?,?,1,?,'ready','now')",
                             ("conf-" + parent_id, state_id, p.artifact_id if artifact else None))
    return state_id


def reference_molscrub(p, retry_failed):
    """Original per-parent history/state reads and exact reuse conditions."""
    cache = {k: v for k, v in p.settings.items() if k != "workers"}
    generation = {k: v for k, v in cache.items() if k != "timeout_seconds"}
    reused, known_failed = set(), set()
    runner = PipelineRunner(p.repo)
    with p.repo.connection() as conn:
        parents = conn.execute("SELECT parent_id,source_smiles FROM parent_ligands").fetchall()
    for parent in parents:
        fp = fingerprint(settings=cache, inputs={"smiles": parent["source_smiles"]}, tool_version=p.version)
        with p.repo.connection() as conn:
            previous = conn.execute("""SELECT fingerprint,status,reason FROM stage_runs
                WHERE stage_name='molscrub' AND entity_type='parent' AND entity_id=?
                ORDER BY rowid DESC LIMIT 1""", (parent["parent_id"],)).fetchone()
            existing = conn.execute("""SELECT s.state_structure_hash,s.generation_fingerprint,
                a.relative_path,a.sha256 FROM molecular_states s LEFT JOIN conformers c ON c.state_id=s.state_id
                LEFT JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
                WHERE s.parent_id=? AND s.active=1""", (parent["parent_id"],)).fetchall()
        if not retry_failed and previous and previous["status"] == "failed" and previous["fingerprint"] == fp:
            known_failed.add(parent["parent_id"])
        elif previous and previous["status"] == "completed" and existing and all(
            state["generation_fingerprint"] in {
                fingerprint(settings=candidate, inputs={"state": state["state_structure_hash"]}, tool_version=p.version)
                for candidate in (generation, cache)}
            and runner._artifact_matches(state["relative_path"], state["sha256"]) for state in existing):
            reused.add(parent["parent_id"])
    return reused, known_failed


@pytest.mark.parametrize("retry_failed", [False, True])
@pytest.mark.parametrize("damage", [None, "missing", "corrupt"])
def test_molscrub_bulk_decisions_match_original(cache_project, monkeypatch, retry_failed, damage):
    p = cache_project
    for status in ("completed", "failed", "stale_failed", "running", "interrupted", "cancelled", None):
        seed_parent(p, str(status), previous=status)
    seed_parent(p, "no-states", state=False)
    seed_parent(p, "inactive", active=False)
    seed_parent(p, "no-artifact", artifact=False)
    seed_parent(p, "no-conformer", conformer=False)
    seed_parent(p, "legacy", generation="legacy")
    seed_parent(p, "stale-generation", generation="stale")
    # Latest history must be selected by rowid, not status or timestamp, and
    # a newer entry for another stage must not hide the MolScrub record.
    with p.repo.connection() as conn:
        conn.execute("INSERT INTO stage_runs VALUES ('unrelated','meeko','parent','completed','x','failed','later','later','x')")
        conn.execute("INSERT INTO stage_runs VALUES ('latest','molscrub','parent','legacy','x','completed','earlier','earlier',NULL)")
    if damage == "missing":
        p.sdf.unlink()
    elif damage == "corrupt":
        p.sdf.write_text("changed coordinates", encoding="ascii")
    expected_reused, expected_failed = reference_molscrub(p, retry_failed)
    events, pending = [], []

    def consume_jobs(function, jobs, workers, **kwargs):
        pending.extend(row["parent_id"] for row, run_id in jobs)
        return iter(())

    monkeypatch.setattr(pipeline, "PLANNING_BATCH_SIZE", 3)  # exercise batch boundaries
    monkeypatch.setattr(pipeline, "bounded_futures", consume_jobs)
    PipelineRunner(p.repo, progress=events.append).run_molscrub(retry_failed=retry_failed)
    reused = {event.item_id for event in events if event.event == "item_skipped" and "existing states" in event.message}
    failed = {event.item_id for event in events if event.event == "item_skipped" and "Known failure" in event.message}
    assert reused == expected_reused
    assert failed == expected_failed
    with p.repo.connection() as conn:
        all_ids = {r[0] for r in conn.execute("SELECT parent_id FROM parent_ligands")}
    assert set(pending) == all_ids - expected_reused - expected_failed


def test_screening_join_matches_original_fingerprint_and_status_checks(cache_project, monkeypatch):
    import rdkit
    p = cache_project
    settings = p.repo.get_settings()["screening"]
    enabled = {name: bool(settings.get(name, False)) for name in ("lipinski", "veber", "egan", "ghose", "boiled_egg")}
    fp = fingerprint(settings={"enabled": enabled, "policy": settings["policy"]},
                     inputs={"smiles": "CCO"}, tool_version=rdkit.__version__)
    cases = [("completed", fp, 1), ("failed", fp, 1), ("running", fp, 1),
             ("completed", "old-settings", 1), ("completed", fp, 0), (None, None, 0)]
    expected = set()
    for index, (status, prior_fp, active) in enumerate(cases):
        parent_id = str(index)
        seed_parent(p, parent_id, previous=None, state=False)
        if status:
            with p.repo.connection() as conn:
                conn.execute("INSERT INTO screening_results VALUES (?,?,?,'{}','{}','pass',NULL,'now',?)",
                             (parent_id, prior_fp, status, active))
        with p.repo.connection() as conn:
            previous = conn.execute("SELECT fingerprint,status FROM screening_results WHERE parent_id=? AND active=1", (parent_id,)).fetchone()
        if previous and previous["fingerprint"] == fp and previous["status"] in ("completed", "failed"):
            expected.add(parent_id)
    monkeypatch.setattr(pipeline, "screen_smiles", lambda *args: SimpleNamespace(
        canonical_smiles="CCO", inchikey="key", decision="pass", reason="", descriptors={}, rules={}))
    events = []
    PipelineRunner(p.repo, progress=events.append).run_screening()
    assert {event.item_id for event in events if event.event == "item_skipped"} == expected
    assert {event.item_id for event in events if event.event == "item_succeeded"} == set(map(str, range(len(cases)))) - expected


def reference_meeko_pending(repo):
    with repo.connection() as conn:
        rows = conn.execute("""SELECT s.state_id,s.parent_id,c.conformer_id FROM molecular_states s
            JOIN conformers c ON c.state_id=s.state_id JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id
            WHERE s.active=1 AND a.active=1 ORDER BY s.parent_id,s.state_id""").fetchall()
        pending, destinations = [], set()
        for row in rows:
            pdbqt = repo.root / "artifacts" / "pdbqt" / row["parent_id"] / row["state_id"] / "ligand.pdbqt"
            log = repo.root / "logs" / "meeko" / (row["state_id"] + ".log")
            for path in (pdbqt.resolve(), log.resolve()):
                if path in destinations:
                    raise ValueError("same output file")
                destinations.add(path)
            artifact = conn.execute("SELECT relative_path,sha256 FROM artifacts WHERE relative_path=? AND active=1",
                                    (pdbqt.relative_to(repo.root).as_posix(),)).fetchone()
            status = conn.execute("SELECT status FROM conformers WHERE conformer_id=?", (row["conformer_id"],)).fetchone()[0]
            if not (status == "pdbqt_ready" and artifact and PipelineRunner(repo)._artifact_matches(artifact["relative_path"], artifact["sha256"])):
                pending.append(row["conformer_id"])
    return pending


def test_meeko_bulk_and_resolved_paths_match_original(cache_project, monkeypatch):
    p = cache_project
    cases = ("valid", "missing", "corrupt", "failed", "inactive-artifact", "no-artifact", "inactive-state")
    for case in cases:
        state_id = seed_parent(p, case, active=case != "inactive-state")
        pdbqt = p.repo.root / "artifacts" / "pdbqt" / case / state_id / "ligand.pdbqt"
        pdbqt.parent.mkdir(parents=True)
        pdbqt.write_text("original PDBQT", encoding="ascii")
        artifact = None if case == "no-artifact" else p.repo.add_artifact(pdbqt, "ligand_pdbqt", "meeko")
        with p.repo.connection() as conn:
            conn.execute("UPDATE conformers SET status=? WHERE state_id=?",
                         ("pdbqt_failed" if case == "failed" else "pdbqt_ready", state_id))
            if case == "inactive-artifact":
                conn.execute("UPDATE artifacts SET active=0 WHERE artifact_id=?", (artifact,))
        if case == "missing":
            pdbqt.unlink()
        elif case == "corrupt":
            pdbqt.write_text("changed PDBQT", encoding="ascii")
    expected = reference_meeko_pending(p.repo)
    pending, events = [], []

    def consume_jobs(function, rows, workers, **kwargs):
        pending.extend(row["conformer_id"] for row in rows)
        return iter(())

    monkeypatch.setattr(pipeline, "PLANNING_BATCH_SIZE", 2)
    monkeypatch.setattr(pipeline, "bounded_futures", consume_jobs)
    PipelineRunner(p.repo, progress=events.append).run_meeko()
    assert pending == expected
    assert len([e for e in events if e.event == "stage_started"]) == 1


def test_cached_stages_use_bounded_connections_without_skipping_hashes(cache_project, monkeypatch):
    p = cache_project
    # A shared SDF is legal; every state's original integrity check must still run.
    for index in range(9):
        seed_parent(p, str(index))
    original_connection = p.repo.connection
    connections = []

    @contextmanager
    def counted():
        connections.append(1)
        with original_connection() as conn:
            yield conn

    monkeypatch.setattr(p.repo, "connection", counted)
    monkeypatch.setattr(pipeline, "PLANNING_BATCH_SIZE", 4)
    original_hash = pipeline.file_sha256
    hashes = []

    def counted_hash(path):
        hashes.append(path)
        return original_hash(path)

    monkeypatch.setattr(pipeline, "file_sha256", counted_hash)
    assert PipelineRunner(p.repo).run_molscrub() == (9, 0)
    assert len(connections) == 1 + 3  # library read + ceil(9 / 4) cache reads
    assert hashes == [p.sdf] * 9
    with original_connection() as conn:
        states = conn.execute("SELECT state_id,parent_id FROM molecular_states ORDER BY state_id").fetchall()
    connections.clear()
    assert len(list(PipelineRunner(p.repo)._docking_plan_rows(states, "profile"))) == 9
    assert len(connections) == 3
    connections.clear()
    assert len(list(PipelineRunner(p.repo)._meeko_plan_rows(states))) == 9
    assert len(connections) == 3


def test_cached_screening_reads_once_for_the_whole_library(cache_project, monkeypatch):
    import rdkit
    p = cache_project
    settings = p.repo.get_settings()["screening"]
    enabled = {name: bool(settings.get(name, False)) for name in ("lipinski", "veber", "egan", "ghose", "boiled_egg")}
    fp = fingerprint(settings={"enabled": enabled, "policy": settings["policy"]},
                     inputs={"smiles": "CCO"}, tool_version=rdkit.__version__)
    for index in range(12):
        seed_parent(p, str(index), state=False, previous=None)
    with p.repo.connection() as conn:
        conn.executemany("INSERT INTO screening_results VALUES (?,?,?,'{}','{}','pass',NULL,'now',1)",
                         [(str(index), fp, "completed") for index in range(12)])
    connections = []
    original = p.repo.connection

    @contextmanager
    def counted():
        connections.append(1)
        with original() as conn:
            yield conn

    monkeypatch.setattr(p.repo, "connection", counted)
    monkeypatch.setattr(pipeline, "screen_smiles", lambda *args: pytest.fail("cached screening called RDKit"))
    assert PipelineRunner(p.repo).run_screening() == (12, 0)
    assert len(connections) == 2  # one joined read, one tool provenance write


def test_destination_resolution_matches_pathlib_and_shares_missing_ancestors(tmp_path, monkeypatch):
    paths = [tmp_path / "missing" / "parent" / str(i) / "out.pdbqt" for i in range(12)]
    existing = tmp_path / "existing"
    existing.mkdir()
    (existing / "file").write_text("x")
    paths += [existing / "file", tmp_path / "missing" / ".." / "existing" / "file"]
    expected = [path.resolve() for path in paths]
    original = Path.resolve
    calls = []

    def counted(path, *args, **kwargs):
        calls.append(path)
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "resolve", counted)
    cache = {}
    assert [_resolve_destination(path, cache) for path in paths] == expected
    assert len(calls) < len(paths) // 2


def test_destination_resolution_retains_symlink_aliases(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("This host does not permit symlink creation")
    cache = {}
    for path in (alias / "new" / "out", target / "new" / "out", alias / ".." / "another"):
        assert _resolve_destination(path, cache) == path.resolve()
    assert _resolve_destination(alias / "new" / "out", cache) == _resolve_destination(target / "new" / "out", cache)


def test_cached_meeko_outputs_still_reject_physical_alias_collisions(cache_project, monkeypatch):
    p = cache_project
    for parent in ("first", "second"):
        state = seed_parent(p, parent)
        output = p.repo.root / "artifacts" / "pdbqt" / parent / state / "ligand.pdbqt"
        output.parent.mkdir(parents=True)
        output.write_text("unchanged")
        p.repo.add_artifact(output, "ligand_pdbqt", "meeko")
    with p.repo.connection() as conn:
        conn.execute("UPDATE conformers SET status='pdbqt_ready'")
    original = Path.resolve

    def aliased(path, *args, **kwargs):
        if path.name == "ligand.pdbqt":
            return p.repo.root / "same-physical-file"
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "resolve", aliased)
    monkeypatch.setattr(pipeline, "bounded_futures", lambda *args, **kwargs: pytest.fail("colliding jobs dispatched"))
    with pytest.raises(ValueError, match="same output file"):
        reference_meeko_pending(p.repo)
    with pytest.raises(ValueError, match="same output file"):
        PipelineRunner(p.repo).run_meeko()


def test_meeko_planning_can_be_cancelled_before_dispatch(cache_project, monkeypatch):
    p = cache_project
    seed_parent(p, "ligand")
    stop = []
    runner = PipelineRunner(p.repo, progress=lambda event: stop.append(True), cancelled=lambda: bool(stop))
    monkeypatch.setattr(pipeline, "bounded_futures", lambda *args, **kwargs: pytest.fail("cancelled planning dispatched jobs"))
    with pytest.raises(InterruptedError):
        runner.run_meeko()
