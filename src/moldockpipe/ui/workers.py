from __future__ import annotations

from PyQt6.QtCore import QObject, pyqtSignal

from ..pipeline import PipelineRunner
from ..project import ProjectRepository


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

    def __init__(self, repository: ProjectRepository, stages: list[str] | None = None) -> None:
        super().__init__()
        self.repository = repository
        self.stages = stages

    def run(self) -> None:
        try:
            runner = PipelineRunner(self.repository, progress=self.progress.emit)
            summary = {}
            if self.stages is None or not self.stages or "screening" in self.stages:
                summary["screening"] = runner.run_screening()
            if self.stages is None or not self.stages or "molscrub" in self.stages:
                summary["molscrub"] = runner.run_molscrub()
            if self.stages is None or not self.stages or "meeko" in self.stages:
                summary["meeko"] = runner.run_meeko()
            if self.stages is None or not self.stages or "vina" in self.stages:
                summary["vina"] = runner.run_vina()
            if self.stages is not None and "postdock" in self.stages:
                summary["postdock"] = runner.run_postdock()
            self.finished.emit(summary)
        except Exception as exc:
            self.failed.emit(str(exc))
