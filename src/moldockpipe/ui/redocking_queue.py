from __future__ import annotations

from PyQt6.QtCore import QEventLoop, QObject, Qt, QThread, QTimer, QUrl, pyqtSignal
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QMessageBox, QPushButton,
    QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..project import ProjectRepository
from ..redocking.campaign import FINISHED_TRIALS, apply_winning_box, campaign_boxes, campaign_snapshot
from ..redocking.models import RedockingSettings
from ..redocking.runner import validate_redocking_prerequisites
from .redocking_dialog import RedockingWorker, STAGE_LABELS


class RedockingQueueController(QObject):
    """Runs a project's persistent receptor-validation queue one item at a time."""

    queue_changed = pyqtSignal()
    activity = pyqtSignal(str)
    validation_completed = pyqtSignal(object)
    box_applied = pyqtSignal(object)

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

    def enqueue(self, profile: dict[str, object], settings: RedockingSettings, *, trial_cap: int | None = None) -> str:
        missing = validate_redocking_prerequisites(self.repository, profile)
        if missing:
            raise ValueError("Redocking cannot be queued. Missing:\n• " + "\n• ".join(missing))
        values: dict[str, object] = {**settings.as_dict()}
        if trial_cap is None:
            queue_id = self.repository.enqueue_redocking(profile, values)
            description = "redocking validation"
        else:
            if settings.seed <= 0:
                raise ValueError("Use a positive seed so every trial starts with the same seed.")
            boxes = campaign_boxes(profile, trial_cap)
            values["campaign"] = {"snapshot": campaign_snapshot(self.repository, profile)}
            queue_id = self.repository.enqueue_redocking(profile, values, campaign_boxes=boxes)
            description = f"redocking campaign with up to {len(boxes)} trials"
        self.activity.emit(f"Queued {description} for {profile.get('name', profile.get('id'))}")
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
        if self._stopping or self._thread is not None:
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
        worker = RedockingWorker(self.repository, profile, settings, item["redocking_run_id"],
                                 campaign=item if "campaign" in values else None)
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
        event_type = event.get("event")
        if event_type == "campaign_trial":
            self._item["trial_label"] = f"Trial {event['trial']}/{event['total']}"
            stage = self._item["trial_label"]
        elif event_type == "stage_started":
            stage = str(event.get("stage", ""))
            if trial := self._item.get("trial_label"):
                stage = f"{trial} · {STAGE_LABELS.get(stage, stage)}"
        else:
            return
        self.repository.update_redocking_queue_stage(str(self._item["queue_id"]), stage)
        self.activity.emit(f"{self._item['receptor_profile_name']}: {STAGE_LABELS.get(stage, stage)}")
        self.queue_changed.emit()

    def _succeeded(self, result: dict[str, object]) -> None:
        if not self._item:
            return
        self.repository.finish_redocking_queue_item(str(self._item["queue_id"]), "completed")
        name = self._item["receptor_profile_name"]
        if "campaign" in self._item["settings"]:
            trials = self.repository.redocking_campaign_trials(str(self._item["queue_id"]))
            _, summary = campaign_summary(trials)
            self.activity.emit(f"Redocking campaign for {name}: {summary}")
        else:
            self.activity.emit(f"Validation artifacts ready: {name}")
        self.validation_completed.emit(result)
        self.queue_changed.emit()

    def _failed(self, message: str) -> None:
        if not self._item:
            return
        queue_id = str(self._item["queue_id"])
        cancelled = self.repository.redocking_cancel_requested(queue_id)
        if self._stopping and not cancelled:
            self.repository.requeue_running_redocking(queue_id, "Application closed; validation will resume from reusable stages")
            self.activity.emit(f"Validation paused for application shutdown: {self._item['receptor_profile_name']}")
        else:
            self.repository.finish_redocking_queue_item(queue_id, "cancelled" if cancelled else "failed", message)
            self.activity.emit(f"Validation {'cancelled' if cancelled else 'failed'}: {self._item['receptor_profile_name']} — {message}")
        self.queue_changed.emit()

    def _thread_finished(self) -> None:
        self._thread = None
        self._worker = None
        self._item = None
        if not self._stopping:
            QTimer.singleShot(0, self._start_next)


TRIAL_LABELS = {"hit": "Target met", "miss": "Measured", "rmsd_missing": "RMSD unavailable", "error": "Retry needed"}


def campaign_summary(trials: list[dict]) -> tuple[str, str]:
    winner = next((trial for trial in trials if trial["outcome"] == "hit"), None)
    if winner:
        return "Target met", f"Trial {winner['trial_number']} · Pose #1 RMSD {winner['rmsd']} Å"
    measured = [trial for trial in trials if trial["rmsd"] is not None]
    best = min(measured, key=lambda trial: trial["rmsd"], default=None)
    finished = sum(trial["outcome"] in FINISHED_TRIALS for trial in trials)
    detail = f"{finished}/{len(trials)} trials evaluated"
    if best:
        detail += f" · Best pose #1 RMSD {best['rmsd']} Å (trial {best['trial_number']})"
    missing = finished - len(measured)
    if missing:
        detail += f" · {missing} without pose #1 RMSD"
    return "Target not met", detail


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
        output = QPushButton("View Details"); output.clicked.connect(self._open_output)
        clear = QPushButton("Clear Completed Validations"); clear.clicked.connect(self._clear)
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
            name = item["receptor_profile_name"]
            if "campaign" in item["settings"]:
                name = f"{name} (campaign)"
                outcome, summary = campaign_summary(self.repository.redocking_campaign_trials(item["queue_id"]))
                if item["status"] == "completed":
                    status_label = outcome
                details = item.get("error_message") or summary
            values = (name, status_label, stage,
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
        if item and "campaign" in item["settings"]:
            RedockingCampaignDialog(self.controller, item, self).exec()
            return
        if not item or not item.get("redocking_run_id"):
            QMessageBox.information(self, "No output yet", "The selected validation has not created a run folder yet.")
            return
        open_trial_output(self, self.repository, item["receptor_profile_id"], item["redocking_run_id"])


def open_trial_output(parent: QWidget, repository: ProjectRepository, profile_id: str, run_id: str) -> None:
    folder = repository.root / "inputs" / "receptors" / profile_id / "redocking" / run_id
    if not folder.is_dir():
        QMessageBox.information(parent, "No output yet", "This trial has not created a run folder.")
        return
    QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))


class RedockingCampaignDialog(QDialog):
    def __init__(self, controller: RedockingQueueController, item: dict, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.controller = controller
        self.queue_id = item["queue_id"]
        self.profile_id = item["receptor_profile_id"]
        self.setWindowTitle(f"Redocking Campaign · {item['receptor_profile_name']}")
        self.resize(900, 420)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.TextFormat.PlainText)
        self.table = QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(("Trial", "Center (Å)", "Dimensions (Å)", "Pose #1 affinity",
                                              "Pose #1 RMSD (Å)", "Outcome", "Details"))
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        output = buttons.addButton("Open Trial Output", QDialogButtonBox.ButtonRole.ActionRole)
        output.clicked.connect(self._open_output)
        self.apply_button = buttons.addButton("Apply Winning Box and Settings", QDialogButtonBox.ButtonRole.ActionRole)
        self.apply_button.clicked.connect(self._apply)
        layout = QVBoxLayout(self)
        layout.addWidget(self.summary)
        layout.addWidget(self.table)
        layout.addWidget(buttons)
        controller.queue_changed.connect(self.refresh)
        self.refresh()

    def refresh(self) -> None:
        selected = self.table.currentRow()
        trials = self.controller.repository.redocking_campaign_trials(self.queue_id)
        item = next(row for row in self.controller.repository.list_redocking_queue() if row["queue_id"] == self.queue_id)
        outcome, summary = campaign_summary(trials)
        self.summary.setText(f"{outcome} · {summary}" if item["status"] == "completed" else summary)
        self.table.setRowCount(len(trials))
        best = min((trial for trial in trials if trial["rmsd"] is not None),
                   key=lambda trial: trial["rmsd"], default=None)
        focus_trial = best["trial_number"] - 1 if best else 0
        for row, trial in enumerate(trials):
            box = trial["box"]
            center = ", ".join(str(box[f"center_{axis}"]) for axis in "xyz")
            size = " × ".join(str(box[f"size_{axis}"]) for axis in "xyz")
            status = trial["outcome"]
            if status == "hit":
                focus_trial = row
            elif item["status"] == "running" and str(item["current_stage"]).startswith(f"Trial {trial['trial_number']}/"):
                status = "running"
            pending = "Not run" if item["status"] in {"completed", "cancelled"} else "Pending"
            values = (trial["trial_number"], center, size,
                      trial["affinity"] if trial["affinity"] is not None else "—",
                      trial["rmsd"] if trial["rmsd"] is not None else "—",
                      TRIAL_LABELS.get(status, "Running" if status == "running" else pending), trial["error"] or "")
            for column, value in enumerate(values):
                cell = QTableWidgetItem(str(value))
                cell.setToolTip(str(value))
                cell.setData(Qt.ItemDataRole.UserRole, trial["run_id"])
                self.table.setItem(row, column, cell)
        if trials:
            self.table.selectRow(selected if 0 <= selected < len(trials) else focus_trial)
        self.table.resizeColumnsToContents()
        self.table.setColumnWidth(6, min(self.table.columnWidth(6), 260))
        self.apply_button.setEnabled(item["status"] == "completed" and outcome == "Target met")

    def _open_output(self) -> None:
        selected = self.table.selectionModel().selectedRows()
        if selected:
            run_id = selected[0].data(Qt.ItemDataRole.UserRole)
            open_trial_output(self, self.controller.repository, self.profile_id, run_id)

    def _apply(self) -> None:
        try:
            profile = apply_winning_box(self.controller.repository, self.queue_id)
        except (KeyError, OSError, RuntimeError, ValueError) as exc:
            QMessageBox.warning(self, "Winning box and settings could not be applied", str(exc))
            return
        self.controller.box_applied.emit(profile)
        self.controller.activity.emit(f"Applied winning docking box and settings for {profile.get('name', self.profile_id)}")
        self.apply_button.setEnabled(False)
        QMessageBox.information(self, "Winning box and settings applied",
                                "The winning box and validation settings are now the production docking defaults for this receptor.")
