from __future__ import annotations

from threading import Event
from time import monotonic

from PyQt6.QtCore import QObject, QThread, pyqtSignal

from ..pipeline import PipelineRunner
from ..project import ProjectRepository
from ..execution_lock import project_execution_lock
from ..models import ProgressEvent


class TaskWorker(QThread):
    """Run repository work with connections created on the calling thread."""
    succeeded = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, function, parent=None):
        super().__init__(parent)
        self.function = function

    def run(self):
        try:
            self.succeeded.emit(self.function())
        except Exception as exc:
            self.failed.emit(str(exc))


class VinaWorker(QObject):
    finished = pyqtSignal(int, int)
    failed = pyqtSignal(str)

    def __init__(self, repository: ProjectRepository) -> None:
        super().__init__()
        self.repository = repository

    def run(self) -> None:
        try:
            completed, failed = PipelineRunner(self.repository).run_vina()
            self.finished.emit(completed, failed)
        except Exception as exc:
            self.failed.emit(str(exc))


class PipelineWorker(QObject):
    progress = pyqtSignal(object)
    finished = pyqtSignal(object)
    failed = pyqtSignal(str)
    stopped = pyqtSignal()

    def __init__(self, repository: ProjectRepository, stages: list[str] | None = None) -> None:
        super().__init__()
        self.repository = repository
        self.stages = stages
        self.cancel_event = Event()
        self._last_progress_at = float("-inf")

    def _report_progress(self, event: object) -> None:
        """Bound item updates, including failures already persisted in history."""
        now = monotonic()
        kind = getattr(event, "event", "")
        important = (kind in {"stage_preparing", "stage_started", "planning_completed", "stage_completed"}
                     or (getattr(event, "total", 0) > 0
                         and getattr(event, "index", 0) == getattr(event, "total", 0)))
        if important or now - self._last_progress_at >= 0.1:
            self._last_progress_at = now
            self.progress.emit(event)

    def stop(self) -> None:
        self.cancel_event.set()

    def run(self) -> None:
        try:
            with project_execution_lock(self.repository.root):
                summary = self._run_pipeline()
            if self.cancel_event.is_set():
                self.stopped.emit()
            else:
                self.finished.emit(summary)
        except InterruptedError:
            self.stopped.emit()
        except Exception as exc:
            self.failed.emit(str(exc))

    def _run_pipeline(self):
        runner = PipelineRunner(self.repository, progress=self._report_progress,
                                cancelled=self.cancel_event.is_set)
        if not self.stages or "vina" in self.stages:
            runner.validate_docking_configuration()
        if self.stages != ["interactions"]:
            self.repository.import_project_input()
        summary = {}
        for stage in ("screening", "molscrub", "meeko", "vina", "postdock", "interactions"):
            if (stage not in {"postdock", "interactions"} and not self.stages) or (self.stages and stage in self.stages):
                if self.cancel_event.is_set():
                    raise InterruptedError("Pipeline stopped safely")
                self._report_progress(ProgressEvent("stage_preparing", stage,
                    message=f"Preparing {stage}; loading dependencies and reading project inputs"))
                if self.cancel_event.is_set():
                    raise InterruptedError("Pipeline stopped safely")
                summary[stage] = getattr(runner, f"run_{stage}")()
        return summary
