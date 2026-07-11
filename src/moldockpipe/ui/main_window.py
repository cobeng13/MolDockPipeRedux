from __future__ import annotations

import shutil
import csv
import time
from pathlib import Path

from PyQt6.QtCore import QObject, QThread, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QDoubleSpinBox,
    QFrame, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton,
    QProgressBar, QSpinBox, QSplitter, QStatusBar, QStyle, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QToolBar, QToolButton, QVBoxLayout, QWidget,
    QMenu,
)

import yaml

from ..project import ProjectRepository
from ..pipeline import PipelineRunner


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.repo: ProjectRepository | None = None
        self.stage_running = False
        self.stage_started_at = 0.0
        self.setWindowTitle("MolDockPipe Redux")
        self.resize(1280, 760)
        self._build_ui()

    def _build_ui(self) -> None:
        toolbar = QToolBar("Project")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        project_button = QToolButton()
        project_button.setText("Project")
        project_menu = QMenu(project_button)
        project_menu.addAction("New Project", self.create_project)
        project_menu.addAction("Open Project", self.open_project)
        project_button.setMenu(project_menu)
        project_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(project_button)
        toolbar.addSeparator()
        import_action = toolbar.addAction("Import")
        import_action.triggered.connect(self.import_csv)
        toolbar.addSeparator()
        self.run_all_action = toolbar.addAction("Run Pipeline")
        self.run_all_action.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.run_all_action.setEnabled(False)
        self.run_all_action.triggered.connect(self.run_all)
        self.resume_action = toolbar.addAction("Resume")
        self.resume_action.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        self.resume_action.setEnabled(False)
        self.resume_action.triggered.connect(self.run_all)
        export_button = QToolButton()
        export_button.setText("Export")
        export_button.setToolTip("Export manifest or leaderboard CSV")
        export_menu = QMenu(export_button)
        export_menu.addAction("Manifest CSV", self.export_manifest_csv)
        export_menu.addAction("Leaderboard CSV", self.export_leaderboard_csv)
        export_button.setMenu(export_menu)
        export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(export_button)
        toolbar.addSeparator()
        self.settings_action = toolbar.addAction("Settings")
        self.settings_action.setEnabled(False)
        self.settings_action.triggered.connect(self.edit_docking_settings)

        self.project_path = QLineEdit()
        self.project_path.setReadOnly(True)
        self.project_path.setPlaceholderText("Create or open a portable project")
        project_label = QLabel("Project")
        self.project_summary = QLabel("No project loaded")
        self.project_summary.setObjectName("projectSummary")
        self.dashboard_status = QLabel("Idle")
        self.dashboard_status.setObjectName("dashboardStatus")

        self.stage_buttons: dict[str, QPushButton] = {}
        self.stage_details: dict[str, QLabel] = {}
        self._dashboard_values = (0, 0, 0, 0, 0)
        stages_row = QHBoxLayout()
        for key, label, callback in (
            ("screening", "Screening", self.run_screening),
            ("molscrub", "States", self.run_molscrub),
            ("meeko", "Ligands", self.run_meeko),
            ("vina", "Docking", self.run_vina),
            ("postdock", "Export", self.run_postdock),
        ):
            button = QPushButton(label)
            button.setMinimumHeight(42)
            button.clicked.connect(callback)
            button.setEnabled(False)
            self.stage_buttons[key] = button
            stages_row.addWidget(button)
        self.stage_action = self.stage_buttons["screening"]
        self.molscrub_action = self.stage_buttons["molscrub"]
        self.meeko_action = self.stage_buttons["meeko"]
        self.vina_action = self.stage_buttons["vina"]
        self.postdock_action = self.stage_buttons["postdock"]

        self.stat_cards: dict[str, QLabel] = {}
        stats_row = QHBoxLayout()
        for key, title in (("ligands", "Ligands"), ("passed", "Passed"), ("failed", "Failed"), ("running", "Running")):
            card = QFrame()
            card.setFrameShape(QFrame.Shape.StyledPanel)
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 8, 12, 8)
            card_layout.addWidget(QLabel(title))
            value = QLabel("0")
            value.setObjectName("statValue")
            value.setAlignment(Qt.AlignmentFlag.AlignCenter)
            value.setWordWrap(True)
            card_layout.addWidget(value)
            self.stat_cards[key] = value
            stats_row.addWidget(card)

        self.overall_progress = CheckpointProgress()
        self.overall_progress.setVisible(False)
        self.stage_progress = QProgressBar()
        self.stage_progress.setRange(0, 100)
        self.stage_progress.setValue(0)
        self.stage_progress.setFormat("%p%")
        self.current_stage_label = QLabel("Current: Idle")
        self.current_item_label = QLabel("")
        self.current_task_label = QLabel("Current task: None")
        top = QWidget()
        top_layout = QVBoxLayout(top)
        top_layout.setContentsMargins(0, 0, 0, 8)
        header = QHBoxLayout()
        header.addWidget(project_label)
        header.addWidget(self.project_path, 1)
        header.addWidget(self.project_summary)
        header.addWidget(self.dashboard_status)
        top_layout.addLayout(header)
        top_layout.addLayout(stats_row)
        top_layout.addLayout(stages_row)
        top_layout.addWidget(self.current_task_label)
        top_layout.addWidget(self.stage_progress)
        self.failure_summary = QLabel("Failures\nNo failures recorded")
        self.failure_summary.setWordWrap(True)
        self.failure_summary.setObjectName("failureSummary")
        top_layout.addWidget(self.failure_summary)

        self.parent_table = QTableWidget(0, 6)
        self.parent_table.setHorizontalHeaderLabels(["Parent ID", "Source SMILES", "Canonical SMILES", "InChIKey", "Parse status", "Screening"])
        self.parent_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.state_table = QTableWidget(0, 6)
        self.state_table.setHorizontalHeaderLabels(["State ID", "Parent ID", "State SMILES", "Charge", "Status", "Truncated"])
        self.state_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.results_table = QTableWidget(0, 5)
        self.results_table.setHorizontalHeaderLabels(["Parent", "Best-state Vina affinity", "Best state", "Prepared", "Docked"])
        self.results_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)

        tabs = QTabWidget()
        tabs.addTab(self.parent_table, "Parent Ligands")
        tabs.addTab(self.state_table, "Molecular States")
        tabs.addTab(self.results_table, "Results")
        # Keep the data views available for later, but hide them while the
        # progress and execution workflow is being refined.
        self.data_tabs = tabs
        self.data_tabs.setVisible(False)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(180)
        self.log.setPlaceholderText("Structured pipeline events and execution logs appear here.")
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(self.log)
        splitter.setSizes([560])
        layout = QVBoxLayout()
        layout.addWidget(top)
        layout.addWidget(splitter)
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
        self.setStatusBar(QStatusBar())

    def create_project(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Choose new project folder")
        if not path:
            return
        self.repo = ProjectRepository.create(Path(path))
        self._activate_project("Created")

    def open_project(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Open project folder")
        if not path:
            return
        repo = ProjectRepository(Path(path))
        if not repo.database_path.exists():
            QMessageBox.warning(self, "Not a project", "This folder does not contain project.sqlite.")
            return
        repo.migrate()
        interrupted = repo.recover_interrupted_runs()
        self.repo = repo
        self._activate_project("Opened" + (f"; recovered {interrupted} interrupted stage runs" if interrupted else ""))

    def import_csv(self) -> None:
        if not self.repo:
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import ligand CSV", filter="CSV files (*.csv)")
        if not path:
            return
        with self.repo.connection() as conn:
            existing = conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0]
        if existing:
            answer = QMessageBox.question(self, "Replace ligand set?", "Replace the current ligand set and generated workflow records with this CSV?", QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
            self.repo.clear_workflow_data()
        try:
            count = self.repo.import_ligands_csv(Path(path))
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))
            return
        self._write_log(f"Imported {count} parent ligands from {path}")
        self.refresh_tables()

    def _activate_project(self, action: str) -> None:
        assert self.repo
        self.project_path.setText(str(self.repo.root))
        self.stage_action.setEnabled(True)
        self.molscrub_action.setEnabled(True)
        self.meeko_action.setEnabled(True)
        self.vina_action.setEnabled(True)
        self.postdock_action.setEnabled(True)
        self.settings_action.setEnabled(True)
        self.run_all_action.setEnabled(True)
        self.resume_action.setEnabled(True)
        self._auto_import_input()
        self._write_log(f"{action} project: {self.repo.root}")
        self.refresh_tables()
        self._refresh_checkpoint_state()

    def _auto_import_input(self) -> None:
        if not self.repo:
            return
        with self.repo.connection() as conn:
            has_ligands = conn.execute("SELECT EXISTS(SELECT 1 FROM parent_ligands)").fetchone()[0]
        if has_ligands:
            return
        for candidate in (self.repo.root / "inputs" / "input.csv", self.repo.root / "input" / "input.csv"):
            if candidate.is_file():
                count = self.repo.import_ligands_csv(candidate)
                self._write_log(f"Auto-loaded {count} ligands from {candidate}")
                return

    def refresh_tables(self) -> None:
        if not self.repo:
            return
        with self.repo.connection() as conn:
            parents = conn.execute("""SELECT p.parent_id, p.source_smiles, p.canonical_source_smiles, p.parent_inchikey,
                p.parse_status, COALESCE(s.decision, 'not screened') AS decision
                FROM parent_ligands p LEFT JOIN screening_results s USING(parent_id) ORDER BY p.parent_id""").fetchall()
            states = conn.execute("SELECT state_id, parent_id, state_isomeric_smiles, formal_charge, status, enumeration_truncated FROM molecular_states ORDER BY parent_id, state_id").fetchall()
            results = conn.execute("""SELECT s.parent_id, MIN(p.affinity) score, s.state_id,
                SUM(CASE WHEN s.status='prepared' THEN 1 ELSE 0 END) prepared, COUNT(p.pose_id) docked
                FROM molecular_states s LEFT JOIN docking_runs d ON d.state_id=s.state_id AND d.is_current=1
                LEFT JOIN docking_poses p ON p.run_id=d.run_id GROUP BY s.parent_id, s.state_id ORDER BY score""").fetchall()
        self._load_table(self.parent_table, parents)
        self._load_table(self.state_table, states)
        self._load_table(self.results_table, results)

    def _refresh_checkpoint_state(self) -> None:
        if not self.repo:
            return
        with self.repo.connection() as conn:
            parents = conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0]
            screened = conn.execute("SELECT COUNT(*) FROM screening_results WHERE status IN ('completed','failed')").fetchone()[0]
            states = conn.execute("SELECT COUNT(*) FROM molecular_states").fetchone()[0]
            prepared = conn.execute("SELECT COUNT(*) FROM conformers WHERE status='pdbqt_ready'").fetchone()[0]
            conformers = conn.execute("SELECT COUNT(*) FROM conformers").fetchone()[0]
            meeko_terminal = conn.execute("SELECT COUNT(*) FROM conformers WHERE status IN ('pdbqt_ready','pdbqt_failed')").fetchone()[0]
            docked = conn.execute("SELECT COUNT(DISTINCT state_id) FROM docking_runs WHERE status='completed' AND is_current=1").fetchone()[0]
            docking_terminal = conn.execute("SELECT COUNT(DISTINCT state_id) FROM docking_runs WHERE status IN ('completed','failed','interrupted')").fetchone()[0]
        self.overall_progress.reset()
        if parents and screened == parents:
            self.overall_progress.complete_stage("screening")
        elif screened:
            self.overall_progress.start_stage("screening")
        if states:
            self.overall_progress.complete_stage("molscrub")
        if conformers and meeko_terminal == conformers:
            self.overall_progress.complete_stage("meeko")
        elif prepared:
            self.overall_progress.start_stage("meeko")
        if states and docking_terminal == states:
            self.overall_progress.complete_stage("vina")
        elif docked:
            self.overall_progress.start_stage("vina")
        self._update_dashboard(parents, screened, states, prepared, docked)

    def _update_dashboard(self, parents: int, screened: int, states: int, prepared: int, docked: int) -> None:
        """Refresh the compact statistics and stage labels shown in the ribbon."""
        if not self.repo:
            return
        self._dashboard_values = (parents, screened, states, prepared, docked)
        with self.repo.connection() as conn:
            passed = conn.execute("SELECT COUNT(*) FROM screening_results WHERE decision != 'fail'").fetchone()[0]
            failed = conn.execute("SELECT COUNT(*) FROM screening_results WHERE decision = 'fail'").fetchone()[0]
            screening_failures = conn.execute("""SELECT p.parent_id, COALESCE(s.reason, 'Screening rule failure')
                FROM screening_results s JOIN parent_ligands p USING(parent_id)
                WHERE s.decision='fail' ORDER BY p.parent_id LIMIT 6""").fetchall()
            pdbqt_failures = conn.execute("""SELECT state_id, COALESCE(reason, 'PDBQT generation failed')
                FROM molecular_states WHERE status='pdbqt_failed' ORDER BY state_id LIMIT 6""").fetchall()
            docking_failures = conn.execute("""SELECT s.state_id, CASE
                WHEN COUNT(p.pose_id)=0 THEN 'No valid Vina poses'
                ELSE 'Docking failed' END
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                LEFT JOIN docking_poses p ON p.run_id=d.run_id
                WHERE d.status IN ('failed','interrupted') GROUP BY d.run_id, s.state_id ORDER BY s.state_id LIMIT 6""").fetchall()
        self.stat_cards["ligands"].setText(str(parents))
        self.stat_cards["passed"].setText(str(passed))
        self.stat_cards["failed"].setText(f"Drug-likeness  {failed}\nPDBQT  {len(pdbqt_failures)}\nDocking  {len(docking_failures)}")
        self.stat_cards["running"].setText(self.current_stage_label.text() if self.stage_running else "Idle")
        self.project_summary.setText(f"{parents} ligands | {states} states | {prepared} prepared | {docked} docked")
        failure_lines = []
        if screening_failures:
            failure_lines.append("Druglikeness Screening: " + "; ".join(f"{row[0]} ({row[1]})" for row in screening_failures))
        if pdbqt_failures:
            failure_lines.append("PDBQT Generation: " + "; ".join(f"{row[0]} ({row[1]})" for row in pdbqt_failures))
        if docking_failures:
            failure_lines.append("Docking: " + "; ".join(f"{row[0]} ({row[1]})" for row in docking_failures))
        self.failure_summary.setText("Failures\n" + ("\n".join(failure_lines) if failure_lines else "No failures recorded"))
        labels = {
            "screening": f"Screening\n{passed} passed; {failed} failed",
            "molscrub": f"States\n{states} generated",
            "meeko": f"Ligands\n{prepared} prepared",
            "vina": f"Docking\n{docked} states complete",
            "postdock": "Export\nReady" if docked else "Export\nPending",
        }
        for stage, button in self.stage_buttons.items():
            status = self.overall_progress.states.get(stage, "pending")
            icon = "✓" if status == "complete" else "⟳" if status == "active" else "○"
            button.setText(f"{icon} {labels[stage]}")
            button.setToolTip(labels[stage].replace("\n", "\n"))
            if status == "complete":
                button.setStyleSheet("QPushButton { color: #166534; background: #dcfce7; border: 1px solid #86efac; }")
            elif status == "active":
                button.setStyleSheet("QPushButton { color: #166534; background: #bbf7d0; border: 2px solid #22c55e; font-weight: 600; }")
            elif status == "warning":
                button.setStyleSheet("QPushButton { color: #92400e; background: #fef3c7; border: 1px solid #f59e0b; }")
            else:
                button.setStyleSheet("QPushButton { color: #374151; background: #f3f4f6; border: 1px solid #d1d5db; }")

    def _refresh_live_dashboard(self) -> None:
        if not self.repo:
            return
        with self.repo.connection() as conn:
            values = (
                conn.execute("SELECT COUNT(*) FROM parent_ligands").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM screening_results WHERE status IN ('completed','failed')").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM molecular_states").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM conformers WHERE status='pdbqt_ready'").fetchone()[0],
                conn.execute("SELECT COUNT(DISTINCT state_id) FROM docking_runs WHERE status='completed' AND is_current=1").fetchone()[0],
            )
        self._update_dashboard(*values)

    @staticmethod
    def _load_table(table: QTableWidget, rows: list[object]) -> None:
        table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            for column, value in enumerate(row):
                table.setItem(row_index, column, QTableWidgetItem("" if value is None else str(value)))
        table.resizeColumnsToContents()

    def _write_log(self, message: str) -> None:
        self.log.append(message)
        self.statusBar().showMessage(message, 5000)

    def run_screening(self) -> None:
        self._start_pipeline(["screening"], "Screening")

    def run_molscrub(self) -> None:
        self._start_pipeline(["molscrub"], "MolScrub")

    def run_meeko(self) -> None:
        self._start_pipeline(["meeko"], "Meeko")

    def _start_pipeline(self, stages: list[str], label: str) -> None:
        if not self.repo or self.stage_running:
            return
        self.stage_running = True
        self.stage_started_at = time.monotonic()
        self.pipeline_is_full = not stages
        self.overall_progress.start_stage(stages[0] if stages else "screening")
        self._set_stage_actions(False)
        self.run_all_action.setEnabled(False)
        self.current_stage_label.setText(f"{label} running")
        self.current_task_label.setText(f"Current task: {label}")
        self._write_log(f"{label} started in background...")
        self.pipeline_thread = QThread(self)
        self.pipeline_worker = PipelineWorker(self.repo, stages)
        self.pipeline_worker.moveToThread(self.pipeline_thread)
        self.pipeline_thread.started.connect(self.pipeline_worker.run)
        self.pipeline_worker.progress.connect(self._on_pipeline_progress)
        self.pipeline_worker.finished.connect(self._pipeline_finished)
        self.pipeline_worker.failed.connect(self._pipeline_failed)
        self.pipeline_worker.finished.connect(self.pipeline_thread.quit)
        self.pipeline_worker.failed.connect(self.pipeline_thread.quit)
        self.pipeline_thread.finished.connect(self.pipeline_worker.deleteLater)
        self.pipeline_thread.finished.connect(self.pipeline_thread.deleteLater)
        self.pipeline_thread.start()

    def run_vina(self) -> None:
        self._start_pipeline(["vina"], "Vina")

    def run_postdock(self) -> None:
        self._start_pipeline(["postdock"], "Post-docking")

    def run_all(self) -> None:
        self.overall_progress.reset()
        self.stage_progress.setValue(0)
        self._start_pipeline([], "Run All")

    def _set_stage_actions(self, enabled: bool) -> None:
        # Keep completed/active ribbon stages visually readable while a run is
        # active. Their handlers still no-op through _start_pipeline's guard.
        for action in (self.stage_action, self.molscrub_action, self.meeko_action, self.vina_action, self.postdock_action):
            action.setEnabled(True)
        self.settings_action.setEnabled(enabled)

    def _on_pipeline_progress(self, event: object) -> None:
        stage_index = {"screening": 0, "molscrub": 1, "meeko": 2, "vina": 3, "postdock": 3}
        stage = getattr(event, "stage", "")
        index = getattr(event, "index", 0)
        total = getattr(event, "total", 0)
        if total:
            self.stage_progress.setValue(int(index * 100 / total))
            self.stage_progress.setFormat(f"{index} / {total}   %p%")
        if getattr(event, "event", "") == "stage_started":
            self.overall_progress.start_stage(stage)
            self._update_dashboard(*self._dashboard_values)
        elif getattr(event, "event", "") == "stage_completed":
            if getattr(event, "failed", 0):
                self.overall_progress.states[stage] = "warning"
                self.overall_progress.update()
            else:
                self.overall_progress.complete_stage(stage)
        self.current_stage_label.setText(f"{stage}: {getattr(event, 'event', '')}")
        self.current_item_label.setText(str(getattr(event, "item_id", "") or ""))
        self.current_task_label.setText(f"Current task: {stage.title()} | {self.current_item_label.text()} ({index} / {total})")
        self.dashboard_status.setText(self.current_stage_label.text())
        self.stat_cards["running"].setText(stage.title() if stage else "Running")
        if total and index:
            elapsed = max(time.monotonic() - self.stage_started_at, 0.1)
            rate = index / elapsed * 60
            remaining = max(total - index, 0)
            eta = remaining / (index / elapsed) if index else 0
            self.dashboard_status.setText(f"{stage.title()} | {rate:.1f}/min | ETA {eta / 60:.1f} min")
        self._refresh_live_dashboard()
        item = getattr(event, "item_id", None)
        message = getattr(event, "message", "")
        self.statusBar().showMessage(f"{stage}: {item or ''} {message}")
        if getattr(event, "event", "") in {"item_failed", "stage_completed"}:
            self._write_log(self._event_summary(event))
        if getattr(event, "event", "") == "stage_completed":
            self.refresh_tables()
            self._refresh_checkpoint_state()

    def _pipeline_finished(self, summary: object) -> None:
        self.stage_running = False
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        if getattr(self, "pipeline_is_full", False):
            self.overall_progress.complete_all()
        else:
            stage = next(iter(summary), "") if isinstance(summary, dict) else ""
            self.overall_progress.complete_stage(stage)
        self.stage_progress.setValue(100)
        self.current_stage_label.setText("Idle")
        self.current_item_label.setText("")
        self.current_task_label.setText("Current task: None")
        self.dashboard_status.setText("Ready")
        self._write_log(f"Pipeline completed: {summary}")
        self.refresh_tables()
        self._refresh_checkpoint_state()

    def _pipeline_failed(self, message: str) -> None:
        self.stage_running = False
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        self._write_log(f"Run All stopped: {message}")
        QMessageBox.critical(self, "Run All failed", message)

    @staticmethod
    def _event_summary(event: object) -> str:
        stage = getattr(event, "stage", "")
        kind = getattr(event, "event", "")
        if kind == "stage_completed":
            return (f"{stage} completed: succeeded={getattr(event, 'succeeded', 0)}, "
                    f"skipped={getattr(event, 'skipped', 0)}, failed={getattr(event, 'failed', 0)}")
        return f"{stage} {kind}: {getattr(event, 'item_id', '')} {getattr(event, 'message', '')}".strip()

    def _vina_finished(self, completed: int, failed: int) -> None:
        self.stage_running = False
        self.vina_action.setEnabled(True)
        self._write_log(f"Vina docking completed: {completed}; failed: {failed}")
        self.refresh_tables()

    def _vina_failed(self, message: str) -> None:
        self.stage_running = False
        self.vina_action.setEnabled(True)
        self._write_log(f"Vina docking failed: {message}")
        QMessageBox.critical(self, "Vina docking failed", message)

    def edit_docking_settings(self) -> None:
        if not self.repo:
            return
        dialog = SettingsDialog(self.repo.get_settings(), self, self.repo)
        dialog.workflow_reset.connect(self._on_workflow_reset)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        config_path = self.repo.root / "project.yml"
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
            config.update(dialog.values())
            config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            self._write_log(f"Updated docking settings in {config_path}")
        except Exception as exc:
            QMessageBox.critical(self, "Settings update failed", str(exc))

    def _on_workflow_reset(self) -> None:
        self.refresh_tables()
        self.overall_progress.reset()
        self.stage_progress.setValue(0)
        self.current_stage_label.setText("Idle")
        self.current_item_label.setText("")
        self._write_log("Workflow reset: refreshed project tables and cleared generated state")
        self._refresh_checkpoint_state()

    def export_manifest_csv(self) -> None:
        if self.repo:
            SettingsDialog(self.repo.get_settings(), self, self.repo)._export_manifest()

    def export_leaderboard_csv(self) -> None:
        if self.repo:
            SettingsDialog(self.repo.get_settings(), self, self.repo)._export_leaderboard()


class CheckpointProgress(QWidget):
    """Workflow state model used by the ribbon buttons."""

    STAGES = ("screening", "molscrub", "meeko", "vina", "postdock")

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.states = {stage: "pending" for stage in self.STAGES}
        self.setMinimumHeight(68)

    def reset(self) -> None:
        self.states = {stage: "pending" for stage in self.STAGES}
        self.update()

    def start_stage(self, stage: str) -> None:
        if stage in self.states:
            self.states[stage] = "active"
            self.update()

    def complete_stage(self, stage: str) -> None:
        if stage in self.states:
            self.states[stage] = "complete"
            self.update()

    def complete_all(self) -> None:
        self.states = {stage: "complete" for stage in self.STAGES}
        self.update()

    def paintEvent(self, event: object) -> None:
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        width = self.width()
        segment = width / len(self.STAGES)
        colors = {"pending": "#d9534f", "active": "#f0ad4e", "complete": "#5cb85c"}
        labels = {"screening": "Screening", "molscrub": "MolScrub", "meeko": "Meeko", "vina": "Vina"}
        y = 18
        painter.setPen(Qt.PenStyle.NoPen)
        for index, stage in enumerate(self.STAGES):
            center_x = int(segment * index + segment / 2)
            if index < len(self.STAGES) - 1:
                painter.setBrush(QColor("#c8c8c8"))
                painter.drawRect(center_x + 9, y - 2, int(segment - 18), 4)
            painter.setBrush(QColor(colors[self.states[stage]]))
            painter.drawEllipse(center_x - 9, y - 9, 18, 18)
            painter.setPen(QColor("#303030"))
            painter.drawText(int(segment * index), 48, int(segment), 18, Qt.AlignmentFlag.AlignCenter, labels[stage])
            painter.setPen(Qt.PenStyle.NoPen)


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

class SettingsDialog(QDialog):
    """Tabbed editor for all project pipeline settings."""

    workflow_reset = pyqtSignal()

    def __init__(self, settings: dict[str, object], parent: QWidget | None = None, repository: ProjectRepository | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Pipeline Settings")
        self.resize(520, 500)
        self.settings = settings
        self.repository = repository
        self.fields: dict[str, object] = {}
        tabs = QTabWidget()
        tabs.addTab(self._screening_tab(settings.get("screening", {})), "Screening")
        tabs.addTab(self._molscrub_tab(settings.get("molscrub", {})), "MolScrub")
        tabs.addTab(self._meeko_tab(settings.get("meeko", {})), "Meeko")
        tabs.addTab(self._vina_tab(settings.get("vina", {})), "Docking")
        tabs.addTab(self._postdock_tab(settings.get("postdock", {})), "Post-docking")
        tabs.addTab(self._workflow_maintenance_tab(), "Workflow Maintenance")
        layout = QVBoxLayout(self)
        layout.addWidget(tabs)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def _workflow_maintenance_tab(self) -> QWidget:
        form = QFormLayout()
        description = QLabel("Reset generated workflow data while preserving source inputs and project settings. This removes derived states, preparations, docking outputs, logs, and reports.")
        description.setWordWrap(True)
        form.addRow(description)
        purge_button = QPushButton("Purge generated workflow data")
        purge_button.clicked.connect(self._purge_workflow)
        form.addRow("Destructive action", purge_button)
        widget = QWidget(); widget.setLayout(form); return widget

    def _purge_workflow(self) -> None:
        if not self.repository:
            return
        answer = QMessageBox.warning(self, "Purge workflow data", "This removes generated states, PDBQT files, docking outputs, logs, reports, and database workflow rows. Inputs and project settings are preserved. Continue?", QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.repository.clear_workflow_data()
        for name in ("artifacts", "logs", "exports"):
            target = self.repository.root / name
            if target.exists():
                shutil.rmtree(target)
            target.mkdir(parents=True, exist_ok=True)
        QMessageBox.information(self, "Purge complete", "Generated workflow data was removed. The input files and project.yml were preserved.")
        self.workflow_reset.emit()

    # Kept as dialog helpers for the Export Data toolbar menu. They are no
    # longer exposed as actions in the destructive-operations tab.
    def _export_manifest(self) -> None:
        if not self.repository:
            return
        output = self.repository.root / "exports" / "manifest.csv"
        output.parent.mkdir(parents=True, exist_ok=True)
        headers = ["id", "smiles", "inchikey", "admet_status", "admet_reason", "sdf_status", "sdf_path", "sdf_reason", "pdbqt_status", "pdbqt_path", "pdbqt_reason", "vina_status", "vina_score", "vina_pose", "vina_reason", "config_hash", "receptor_sha1", "tools_rdkit", "tools_meeko", "tools_vina", "created_at", "updated_at"]
        with self.repository.connection() as conn:
            parents = conn.execute("""SELECT p.*, COALESCE(s.reason, '') screening_reason,
                COALESCE(s.status, '') screening_status FROM parent_ligands p
                LEFT JOIN screening_results s USING(parent_id) ORDER BY p.parent_id""").fetchall()
            rows = []
            for parent in parents:
                states = conn.execute("""SELECT s.*, c.status conformer_status, a.relative_path sdf_path
                    FROM molecular_states s LEFT JOIN conformers c ON c.state_id=s.state_id
                    LEFT JOIN artifacts a ON a.artifact_id=c.sdf_artifact_id WHERE s.parent_id=?""", (parent["parent_id"],)).fetchall()
                pose = conn.execute("""SELECT p.affinity, d.receptor_hash, raw.relative_path pose_path, d.settings_fingerprint
                    FROM docking_poses p JOIN docking_runs d ON d.run_id=p.run_id
                    LEFT JOIN artifacts raw ON raw.artifact_id=d.raw_output_artifact_id
                    JOIN molecular_states s ON s.state_id=d.state_id
                    WHERE s.parent_id=? AND d.is_current=1 ORDER BY p.affinity LIMIT 1""", (parent["parent_id"],)).fetchone()
                rows.append({
                    "id": parent["parent_id"], "smiles": parent["source_smiles"], "inchikey": parent["parent_inchikey"] or "",
                    "admet_status": parent["screening_status"], "admet_reason": parent["screening_reason"],
                    "sdf_status": "DONE" if states else "", "sdf_path": ";".join(str(s["sdf_path"] or "") for s in states), "sdf_reason": "",
                    "pdbqt_status": "DONE" if states and all(s["conformer_status"] == "pdbqt_ready" for s in states) else "", "pdbqt_path": "", "pdbqt_reason": "",
                    "vina_status": "DONE" if pose else "", "vina_score": pose["affinity"] if pose else "", "vina_pose": pose["pose_path"] if pose else "", "vina_reason": "",
                    "config_hash": pose["settings_fingerprint"] if pose else "", "receptor_sha1": pose["receptor_hash"] if pose else "",
                    "tools_rdkit": "RDKit", "tools_meeko": "Meeko", "tools_vina": "Vina", "created_at": parent["created_at"], "updated_at": parent["created_at"],
                })
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader(); writer.writerows(rows)
        QMessageBox.information(self, "Export complete", f"Wrote {output}")

    def _export_leaderboard(self) -> None:
        if not self.repository:
            return
        output = self.repository.root / "exports" / "leaderboard.csv"
        output.parent.mkdir(parents=True, exist_ok=True)
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT parent_id, state_id, mode_index, affinity, rmsd_lb, rmsd_ub, run_id, pose_rank
                FROM (SELECT s.parent_id, s.state_id, p.mode_index, p.affinity, p.rmsd_lb, p.rmsd_ub, d.run_id,
                    ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity ASC) pose_rank
                    FROM docking_poses p JOIN docking_runs d ON d.run_id=p.run_id
                    JOIN molecular_states s ON s.state_id=d.state_id WHERE d.is_current=1)
                WHERE pose_rank <= 3 ORDER BY affinity ASC, parent_id, pose_rank""").fetchall()
        headers = ["parent_id", "state_id", "mode_index", "affinity", "rmsd_lb", "rmsd_ub", "run_id", "pose_rank"]
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader(); writer.writerows([dict(row) for row in rows])
        QMessageBox.information(self, "Export complete", f"Wrote {output}")

    def _screening_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        policy = QComboBox()
        policy.addItem("Annotate only", "annotate_only")
        policy.addItem("Exclude failing all selected rules", "exclude_failing_all")
        policy.addItem("Exclude failing any selected rule", "exclude_failing_any")
        policy.addItem("Manual review", "manual_review")
        policy.setCurrentIndex(max(0, policy.findData(values.get("policy", "annotate_only"))))
        self.fields["screening.policy"] = policy
        form.addRow("Policy", policy)
        for key, label in (("lipinski", "Lipinski"), ("veber", "Veber"), ("egan", "Egan"), ("ghose", "Ghose")):
            box = QCheckBox(label)
            box.setChecked(bool(values.get(key, key in ("lipinski", "veber"))))
            self.fields[f"screening.{key}"] = box
            form.addRow("Rule", box)
        widget = QWidget(); widget.setLayout(form); return widget

    def _molscrub_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        self._double(form, "pH", "molscrub.ph", values.get("ph", 7.4), 0, 14)
        enumerate_box = QCheckBox("Enumerate molecular states")
        enumerate_box.setChecked(bool(values.get("enumerate_states", True)))
        self.fields["molscrub.enumerate_states"] = enumerate_box
        form.addRow("", enumerate_box)
        self._integer(form, "Maximum states", "molscrub.max_states", values.get("max_states", 32), 1, 10000)
        widget = QWidget(); widget.setLayout(form); return widget

    def _meeko_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        self._integer(form, "Parallel workers", "meeko.workers", values.get("workers", 4), 1, 32)
        widget = QWidget(); widget.setLayout(form); return widget

    def _vina_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        self._line(form, "Receptor", "vina.receptor", values.get("receptor", "inputs/receptor_prepared.pdbqt"))
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
            self._double(form, key.replace("_", " ").title(), f"vina.{key}", values.get(key, 0.0 if key.startswith("center") else 20.0))
        self._integer(form, "Exhaustiveness", "vina.exhaustiveness", values.get("exhaustiveness", 8), 1, 10000)
        self._integer(form, "Number of modes", "vina.num_modes", values.get("num_modes", 9), 1, 100)
        self._double(form, "Energy range", "vina.energy_range", values.get("energy_range", 3), 0, 100)
        self._integer(form, "Random seed", "vina.seed", values.get("seed", 42), 0, 2147483647)
        self._integer(form, "CPU count", "vina.cpu_count", values.get("cpu_count", 1), 1, 256)
        widget = QWidget(); widget.setLayout(form); return widget

    def _postdock_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        mode = QComboBox()
        mode.addItem("Split only", "split_only")
        mode.addItem("Split and convert to SDF", "split_and_sdf")
        mode.addItem("Convert multi-pose output to SDF", "sdf_only")
        mode.setCurrentIndex(max(0, mode.findData(values.get("mode", "split_and_sdf"))))
        self.fields["postdock.mode"] = mode
        form.addRow("Operation", mode)
        self._integer(form, "Poses per compound", "postdock.poses_per_compound", values.get("poses_per_compound", 3), 1, 100)
        selector = QPushButton("Select successfully docked compounds")
        selector.clicked.connect(self._select_compounds)
        self.fields["postdock.selected_parents"] = list(values.get("selected_parents", []))
        form.addRow("Compounds", selector)
        self.postdock_selection_label = QLabel(self._selection_text())
        form.addRow("Selection", self.postdock_selection_label)
        widget = QWidget(); widget.setLayout(form); return widget

    def _selection_text(self) -> str:
        selected = self.fields.get("postdock.selected_parents", [])
        return "All successfully docked compounds" if not selected else f"{len(selected)} compounds selected"

    def _select_compounds(self) -> None:
        if not self.repository:
            return
        with self.repository.connection() as conn:
            rows = conn.execute("""SELECT s.parent_id, MIN(p.affinity) best_affinity FROM docking_runs d
                JOIN molecular_states s ON s.state_id=d.state_id JOIN docking_poses p ON p.run_id=d.run_id
                WHERE d.status='completed' AND d.is_current=1 GROUP BY s.parent_id ORDER BY best_affinity ASC, s.parent_id""").fetchall()
        dialog = CompoundSelectorDialog(rows, set(self.fields.get("postdock.selected_parents", [])), self.repository, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.fields["postdock.selected_parents"] = dialog.selected()
            self.postdock_selection_label.setText(self._selection_text())

    def _line(self, form: QFormLayout, label: str, key: str, value: object) -> None:
        field = QLineEdit(str(value))
        self.fields[key] = field
        form.addRow(label, field)

    def _double(self, form: QFormLayout, label: str, key: str, value: object, minimum: float = -100000, maximum: float = 100000) -> None:
        field = QDoubleSpinBox()
        field.setRange(minimum, maximum)
        field.setDecimals(4)
        field.setValue(float(value))
        self.fields[key] = field
        form.addRow(label, field)

    def _integer(self, form: QFormLayout, label: str, key: str, value: object, minimum: int, maximum: int) -> None:
        field = QSpinBox()
        field.setRange(minimum, maximum)
        field.setValue(int(value))
        self.fields[key] = field
        form.addRow(label, field)

    def values(self) -> dict[str, dict[str, object]]:
        values: dict[str, dict[str, object]] = {"screening": {}, "molscrub": {}, "meeko": {}, "vina": {}, "postdock": {}}
        for key, field in self.fields.items():
            if isinstance(field, QLineEdit):
                value = field.text().strip()
            elif isinstance(field, QDoubleSpinBox):
                value = field.value()
            elif isinstance(field, QCheckBox):
                value = field.isChecked()
            elif isinstance(field, QComboBox):
                value = field.currentData()
            else:
                value = field.value()
            section, name = key.split(".", 1)
            values[section][name] = value
        return values


class CompoundSelectorDialog(QDialog):
    def __init__(self, compounds: list[object], selected: set[str], repository: ProjectRepository, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Select Docked Compounds")
        self.table = QTableWidget(len(compounds), 4)
        self.table.setHorizontalHeaderLabels(["", "ID", "Vina score", "Exported Status"])
        self.table.setSortingEnabled(False)
        for row_index, compound in enumerate(compounds):
            parent_id = str(compound["parent_id"])
            check = QTableWidgetItem()
            check.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled)
            check.setCheckState(Qt.CheckState.Checked if parent_id in selected else Qt.CheckState.Unchecked)
            self.table.setItem(row_index, 1, QTableWidgetItem(parent_id))
            self.table.setItem(row_index, 2, QTableWidgetItem(f"{float(compound['best_affinity']):.3f}"))
            exported = any((repository.root / folder / parent_id).exists() for folder in ("For_PostDocking/SDF", "For_PostDocking/PDBQTs"))
            status = QTableWidgetItem("Exported" if exported else "Not exported")
            if exported:
                check.setFlags(Qt.ItemFlag.ItemIsEnabled)
                check.setCheckState(Qt.CheckState.Unchecked)
            self.table.setItem(row_index, 0, check)
            self.table.setItem(row_index, 3, status)
        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)
        layout = QVBoxLayout(self)
        layout.addWidget(self.table)
        if not compounds:
            layout.addWidget(QLabel("No successfully docked compounds are available."))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def selected(self) -> list[str]:
        selected = []
        for row in range(self.table.rowCount()):
            item = self.table.item(row, 0)
            if item and item.checkState() == Qt.CheckState.Checked:
                selected.append(self.table.item(row, 1).text())
        return selected
