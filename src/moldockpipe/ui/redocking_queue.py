from __future__ import annotations

from PyQt6.QtCore import QEventLoop, QObject, QThread, QTimer, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..project import ProjectRepository
from ..redocking.models import RedockingSettings
from ..redocking.runner import validate_redocking_prerequisites
from .redocking_dialog import RedockingWorker, STAGE_LABELS


class RedockingQueueController(QObject):
    """Runs a project's persistent receptor-validation queue one item at a time."""

    queue_changed = pyqtSignal()
    activity = pyqtSignal(str)
    validation_completed = pyqtSignal(object)

    def __init__(self, repository: ProjectRepository, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self.repository = repository
        self._thread: QThread | None = None
        self._worker: RedockingWorker | None = None
        self._item: dict[str, object] | None = None
        self._stopping = False

    @property
    def is_running(self) -> bool:
        return self._thread is not None and self._thread.isRunning()

    def start(self) -> None:
        self._stopping = False
        QTimer.singleShot(0, self._start_next)

    def shutdown(self, timeout_ms: int = 10000) -> bool:
        """Stop an active subprocess and leave its item resumable for the next launch."""
        self._stopping = True
        if not self.is_running:
            return True
        if self._worker:
            self._worker.cancel()
        assert self._thread
        loop = QEventLoop()
        self._thread.finished.connect(loop.quit)
        QTimer.singleShot(timeout_ms, loop.quit)
        loop.exec()
        return not self._thread or not self._thread.isRunning()

    def enqueue(self, profile: dict[str, object], settings: RedockingSettings) -> str:
        missing = validate_redocking_prerequisites(self.repository, profile)
        if missing:
            raise ValueError("Redocking cannot be queued. Missing:\n• " + "\n• ".join(missing))
        queue_id = self.repository.enqueue_redocking(profile, settings.as_dict())
        self.activity.emit(f"Queued redocking validation for {profile.get('name', profile.get('id'))}")
        self.queue_changed.emit()
        self.start()
        return queue_id

    def cancel(self, queue_id: str) -> bool:
        changed = self.repository.cancel_redocking_queue_item(queue_id)
        if changed and self._item and str(self._item["queue_id"]) == queue_id and self._worker:
            self._worker.cancel()
        if changed:
            self.queue_changed.emit()
        return changed

    def retry(self, queue_id: str) -> bool:
        changed = self.repository.retry_redocking_queue_item(queue_id)
        if changed:
            self.queue_changed.emit()
            self.start()
        return changed

    def _start_next(self) -> None:
        if self._stopping or self.is_running:
            return
        item = self.repository.claim_next_redocking()
        if not item:
            self.queue_changed.emit()
            return
        queue_id = str(item["queue_id"])
        profile = next((profile for profile in self.repository.get_receptor_profiles(include_archived=True)
                        if str(profile.get("id")) == str(item["receptor_profile_id"])), None)
        if not profile or profile.get("archived"):
            self.repository.finish_redocking_queue_item(queue_id, "failed", "Receptor profile is missing or archived")
            self.queue_changed.emit(); QTimer.singleShot(0, self._start_next); return
        missing = validate_redocking_prerequisites(self.repository, profile)
        if missing:
            self.repository.finish_redocking_queue_item(queue_id, "failed", "Missing: " + "; ".join(missing))
            self.queue_changed.emit(); QTimer.singleShot(0, self._start_next); return

        values = item["settings"]
        settings = RedockingSettings(
            int(values["exhaustiveness"]), int(values["num_modes"]), float(values["energy_range"]),
            int(values["seed"]), int(values["cpu_count"]),
        )
        thread = QThread(self)
        worker = RedockingWorker(self.repository, profile, settings, str(item["redocking_run_id"]))
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._progress)
        worker.succeeded.connect(self._succeeded)
        worker.failed.connect(self._failed)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(self._thread_finished)
        thread.finished.connect(thread.deleteLater)
        self._item, self._thread, self._worker = item, thread, worker
        self.activity.emit(f"Validation started: {item['receptor_profile_name']}")
        self.queue_changed.emit()
        thread.start()

    def _progress(self, event: dict[str, object]) -> None:
        if not self._item:
            return
        if event.get("event") == "stage_started":
            stage = str(event.get("stage", ""))
            self.repository.update_redocking_queue_stage(str(self._item["queue_id"]), stage)
            self.activity.emit(f"{self._item['receptor_profile_name']}: {STAGE_LABELS.get(stage, stage)}")
            self.queue_changed.emit()

    def _succeeded(self, result: dict[str, object]) -> None:
        if not self._item:
            return
        self.repository.finish_redocking_queue_item(str(self._item["queue_id"]), "completed")
        self.activity.emit(f"Validation artifacts ready: {self._item['receptor_profile_name']}")
        self.validation_completed.emit(result)
        self.queue_changed.emit()

    def _failed(self, message: str) -> None:
        if not self._item:
            return
        queue_id = str(self._item["queue_id"])
        if self._stopping:
            self.repository.requeue_running_redocking(queue_id, "Application closed; validation will resume from reusable stages")
            self.activity.emit(f"Validation paused for application shutdown: {self._item['receptor_profile_name']}")
        else:
            rows = {str(row["queue_id"]): row for row in self.repository.list_redocking_queue()}
            cancelled = bool(rows.get(queue_id, {}).get("cancel_requested"))
            self.repository.finish_redocking_queue_item(queue_id, "cancelled" if cancelled else "failed", message)
            self.activity.emit(f"Validation {'cancelled' if cancelled else 'failed'}: {self._item['receptor_profile_name']} — {message}")
        self.queue_changed.emit()

    def _thread_finished(self) -> None:
        self._thread = None
        self._worker = None
        self._item = None
        if not self._stopping:
            QTimer.singleShot(0, self._start_next)


class RedockingQueueDialog(QDialog):
    def __init__(self, controller: RedockingQueueController, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.controller = controller
        self.repository = controller.repository
        self.setWindowTitle("Redocking Validation Queue")
        self.resize(900, 460)
        description = QLabel(
            "Validations run one at a time in the background. You may close this window and leave MolDockPipe running. "
            "Queued work is retained if the application is restarted."
        )
        description.setWordWrap(True)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(("Receptor", "Status", "Current stage", "Queued", "Started", "Details"))
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        row = QHBoxLayout()
        cancel = QPushButton("Cancel Selected"); cancel.clicked.connect(self._cancel)
        retry = QPushButton("Retry Selected"); retry.clicked.connect(self._retry)
        output = QPushButton("Open Output Folder"); output.clicked.connect(self._open_output)
        clear = QPushButton("Clear Completed"); clear.clicked.connect(self._clear)
        row.addWidget(cancel); row.addWidget(retry); row.addWidget(output); row.addWidget(clear); row.addStretch()
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self); layout.addWidget(description); layout.addWidget(self.table); layout.addLayout(row); layout.addWidget(buttons)
        self.controller.queue_changed.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        selected = self._selected_id()
        items = self.repository.list_redocking_queue()
        self.table.setRowCount(0)
        for item in items:
            row = self.table.rowCount(); self.table.insertRow(row)
            stage = STAGE_LABELS.get(str(item.get("current_stage") or ""), str(item.get("current_stage") or "—"))
            details = item.get("error_message") or (f"Run {item['redocking_run_id']}" if item.get("redocking_run_id") else "")
            status_label = "Artifacts Ready" if item["status"] == "completed" else str(item["status"]).title()
            values = (item["receptor_profile_name"], status_label, stage,
                      item["created_at"], item.get("started_at") or "—", details)
            for column, value in enumerate(values):
                cell = QTableWidgetItem(str(value)); cell.setData(256, item["queue_id"]); self.table.setItem(row, column, cell)
            if selected == item["queue_id"]:
                self.table.selectRow(row)
        self.table.resizeColumnsToContents()

    def _selected_id(self) -> str | None:
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return None
        return str(self.table.item(rows[0].row(), 0).data(256))

    def _cancel(self) -> None:
        queue_id = self._selected_id()
        if queue_id and QMessageBox.question(self, "Cancel validation", "Cancel the selected queued or running validation?") == QMessageBox.StandardButton.Yes:
            self.controller.cancel(queue_id)

    def _retry(self) -> None:
        queue_id = self._selected_id()
        if queue_id:
            self.controller.retry(queue_id)

    def _clear(self) -> None:
        self.repository.clear_finished_redocking_queue(); self.refresh()

    def _open_output(self) -> None:
        queue_id = self._selected_id()
        item = next((row for row in self.repository.list_redocking_queue() if row["queue_id"] == queue_id), None)
        if not item or not item.get("redocking_run_id"):
            QMessageBox.information(self, "No output yet", "The selected validation has not created a run folder yet.")
            return
        folder = self.repository.root / "inputs" / "receptors" / str(item["receptor_profile_id"]) / "redocking" / str(item["redocking_run_id"])
        if not folder.is_dir():
            QMessageBox.warning(self, "Output unavailable", f"The run folder does not exist:\n{folder}")
            return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))
