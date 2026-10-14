from contextlib import contextmanager
from threading import Barrier, Event, Lock, get_ident

import pytest
import yaml

from moldockpipe.parallel import bounded_futures, worker_count
from moldockpipe.pipeline import PipelineRunner
from moldockpipe.project import ProjectRepository


@pytest.mark.parametrize("cpus,requested,jobs,expected", [
    (12, 0, 100, 12), (64, 0, 100, 64), (12, 0, 3, 3),
    (None, 0, 100, 1), (12, 4, 100, 4), (12, 1, 100, 1), (12, 0, 0, 1),
])
def test_worker_limits_follow_cpu_count_or_explicit_override(monkeypatch, cpus, requested, jobs, expected):
    monkeypatch.setattr("moldockpipe.parallel.os.cpu_count", lambda: cpus)
    assert worker_count({"workers": requested}, jobs) == expected


def test_default_workers_are_automatic(tmp_path, monkeypatch):
    monkeypatch.setattr("moldockpipe.parallel.os.cpu_count", lambda: 12)
    repo = ProjectRepository.create(tmp_path / "project")
    for stage in ("molscrub", "meeko"):
        assert repo.get_settings()[stage]["workers"] == 0
        assert worker_count(repo.get_settings()[stage], 100) == 12
    assert worker_count({}, 100) == 12


@pytest.mark.parametrize("old_workers,expected", [(4, 0), (None, 0), (2, 2), (0, 0)])
def test_legacy_meeko_default_upgrades_once(tmp_path, monkeypatch, old_workers, expected):
    monkeypatch.setattr("moldockpipe.parallel.os.cpu_count", lambda: 12)
    repo = ProjectRepository.create(tmp_path / "project")
    config = repo.get_settings()
    config.pop("preparation_defaults_version")
    config["meeko"] = {} if old_workers is None else {"workers": old_workers}
    path = repo.root / "project.yml"
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    repo.migrate()
    assert repo.get_settings()["meeko"]["workers"] == expected
    if expected == 0:
        assert worker_count(repo.get_settings()["meeko"], 100) == 12
    # A manual choice made after upgrading must survive every later reopen.
    config = repo.get_settings()
    config["meeko"]["workers"] = 4
    path.write_text(yaml.safe_dump(config), encoding="utf-8")
    repo.migrate()
    assert repo.get_settings()["meeko"]["workers"] == 4


def test_parallel_queue_stays_bounded():
    submitted = consumed = 0

    def jobs():
        nonlocal submitted
        for index in range(5000):
            submitted += 1
            assert submitted - consumed <= 3
            yield index

    for item, future in bounded_futures(lambda value: value * 2, jobs(), 3):
        assert future.result() == item * 2
        consumed += 1
    assert consumed == 5000


def test_stop_during_job_selection_does_not_launch_selected_job():
    stop = Event()
    launched, discarded, completed = [], [], []

    def jobs():
        yield 0
        stop.set()
        yield 1

    with pytest.raises(InterruptedError):
        for item, future in bounded_futures(lambda value: launched.append(value), jobs(), 2,
                cancelled=stop.is_set, on_discard=discarded.append):
            future.result()
            completed.append(item)
    assert launched == [0]
    assert completed == [0]
    assert discarded == [1]


class ConcurrencyProbe:
    def __init__(self):
        self.lock = Lock()
        self.barrier = Barrier(2)
        self.calls = self.active = self.peak = 0

    @contextmanager
    def running(self):
        with self.lock:
            self.calls += 1
            ordinal = self.calls
            self.active += 1
            self.peak = max(self.peak, self.active)
        try:
            if ordinal <= 2:
                self.barrier.wait(timeout=10)
            yield
        finally:
            with self.lock:
                self.active -= 1


@pytest.fixture
def parallel_project(tmp_path, monkeypatch):
    Chem = pytest.importorskip("rdkit.Chem")
    from moldockpipe.services.molscrub import GeneratedState, MolScrubService
    from moldockpipe.services.meeko import MeekoService

    repo = ProjectRepository.create(tmp_path / "project")
    source = repo.root / "inputs" / "input.csv"
    source.write_text("id,smiles\n" + "".join(f"lig{i},CCO\n" for i in range(6)), encoding="utf-8")
    repo.import_project_input()
    config = repo.get_settings()
    for stage in ("molscrub", "meeko"):
        config[stage]["workers"] = 2
    (repo.root / "project.yml").write_text(yaml.safe_dump(config), encoding="utf-8")
    probes = {stage: ConcurrencyProbe() for stage in ("molscrub", "meeko")}

    def generate(self, smiles, **kwargs):
        with probes["molscrub"].running():
            return [GeneratedState(smiles, smiles, None, 0, "testhash", 0, 0,
                                   Chem.MolFromSmiles(smiles))], False

    def prepare(self, sdf, pdbqt, log):
        with probes["meeko"].running():
            pdbqt.parent.mkdir(parents=True, exist_ok=True)
            log.parent.mkdir(parents=True, exist_ok=True)
            pdbqt.write_text("REMARK parallel preparation\nROOT\n"
                "ATOM      1  C   LIG A   1       1.000   2.000   3.000  1.00  0.00     0.000 C\n"
                "ENDROOT\nTORSDOF 0\n", encoding="ascii")
            log.write_text("prepared", encoding="ascii")

    monkeypatch.setattr(MolScrubService, "generate", generate)
    monkeypatch.setattr(MeekoService, "prepare_ligand", prepare)
    return repo, PipelineRunner(repo), probes


def test_molscrub_runs_concurrently_but_persists_on_coordinator(parallel_project, monkeypatch):
    repo, runner, probes = parallel_project
    coordinator = get_ident()
    writers = []
    persist = runner._persist_states

    def record_writer(*args):
        writers.append(get_ident())
        return persist(*args)

    monkeypatch.setattr(runner, "_persist_states", record_writer)
    assert runner.run_molscrub() == (6, 0)
    assert probes["molscrub"].peak == 2
    assert writers == [coordinator] * 6
    with repo.connection() as conn:
        assert conn.execute("PRAGMA foreign_key_check").fetchall() == []
        assert conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0] == 6


@pytest.mark.parametrize("stage", ["molscrub", "meeko"])
def test_stop_drains_active_jobs_and_resume_reuses_completed_work(parallel_project, stage):
    repo, runner, probes = parallel_project
    if stage == "meeko":
        assert runner.run_molscrub() == (6, 0)
    stop = Event()

    def progress(event):
        if event.stage == stage and event.event == "item_succeeded":
            stop.set()

    runner.cancelled = stop.is_set
    runner.progress = progress
    method = getattr(runner, f"run_{stage}")
    with pytest.raises(InterruptedError, match="stopped safely"):
        method()
    assert probes[stage].calls == 2
    assert probes[stage].peak == 2
    with repo.connection() as conn:
        if stage == "molscrub":
            assert conn.execute("SELECT COUNT(*) FROM stage_runs WHERE status='completed'").fetchone()[0] == 2
            assert conn.execute("SELECT COUNT(*) FROM stage_runs WHERE status='running'").fetchone()[0] == 0
        else:
            assert conn.execute("SELECT COUNT(*) FROM conformers WHERE status='pdbqt_ready'").fetchone()[0] == 2
    runner.progress = None
    stop.clear()
    assert method() == ((6, 0) if stage == "molscrub" else (4, 0))
    assert probes[stage].calls == 6
    # Scheduling changes must never invalidate scientific results.
    config = repo.get_settings()
    config[stage]["workers"] = 3
    (repo.root / "project.yml").write_text(yaml.safe_dump(config), encoding="utf-8")
    assert method() == ((6, 0) if stage == "molscrub" else (0, 0))
    assert probes[stage].calls == 6


def test_stop_at_molscrub_start_leaves_no_running_or_failed_job(parallel_project):
    repo, runner, probes = parallel_project
    stop = Event()
    runner.cancelled = stop.is_set
    runner.progress = lambda event: stop.set() if event.event == "item_started" else None
    with pytest.raises(InterruptedError):
        runner.run_molscrub()
    assert probes["molscrub"].calls == 0
    with repo.connection() as conn:
        assert [row[0] for row in conn.execute("SELECT status FROM stage_runs")] == ["cancelled"]
        assert conn.execute("SELECT COUNT(*) FROM molecular_states").fetchone()[0] == 0


def test_concurrent_runner_and_crash_recovery_are_blocked(parallel_project, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from moldockpipe.services.molscrub import MolScrubService

    repo, runner, probes = parallel_project
    entered, release = Event(), Event()
    original = MolScrubService.generate

    def wait_for_release(self, *args, **kwargs):
        entered.set()
        assert release.wait(timeout=10)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(MolScrubService, "generate", wait_for_release)
    with ThreadPoolExecutor(max_workers=1) as pool:
        first = pool.submit(runner.run_molscrub)
        try:
            assert entered.wait(timeout=5)
            with pytest.raises(RuntimeError, match="already being processed"):
                PipelineRunner(repo).run_molscrub()
            with pytest.raises(RuntimeError, match="already being processed"):
                repo.recover_interrupted_runs()
            with pytest.raises(RuntimeError, match="already being processed"):
                repo.clear_workflow_data()
            with pytest.raises(RuntimeError, match="already being processed"):
                repo.import_project_input()
        finally:
            release.set()
        assert first.result(timeout=10) == (6, 0)
    assert probes["molscrub"].calls == 6


def test_meeko_rejects_colliding_jobs_before_any_worker_writes(parallel_project):
    repo, runner, probes = parallel_project
    assert runner.run_molscrub() == (6, 0)
    with repo.connection() as conn:
        row = conn.execute("SELECT state_id,sdf_artifact_id FROM conformers LIMIT 1").fetchone()
        conn.execute("INSERT INTO conformers VALUES ('duplicate',?,2,?,'ready','now')",
                     (row["state_id"], row["sdf_artifact_id"]))
    with pytest.raises(ValueError, match="same output file"):
        runner.run_meeko()
    assert probes["meeko"].calls == 0
