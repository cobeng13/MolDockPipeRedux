from __future__ import annotations

import time
from pathlib import Path

from PyQt6.QtCore import QThread, QTimer, Qt, QUrl
from PyQt6.QtGui import QCloseEvent, QDesktopServices, QFont
from PyQt6.QtWidgets import (
    QDialog, QFileDialog, QFrame, QHBoxLayout, QInputDialog, QLabel, QLineEdit, QMainWindow, QMessageBox,
    QPushButton, QProgressBar, QProgressDialog, QSplitter, QStatusBar, QStyle, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QToolBar, QToolButton, QVBoxLayout, QWidget, QMenu,
)

from ..parallel import logical_cpu_count
from ..project import ProjectRepository
from ..execution_lock import project_execution_lock
from ..csv_exports import export_leaderboard_csv as write_leaderboards, export_manifest_csv as write_manifests
from ..reporting import generate_project_report
from .compound_selector import CompoundSelectorDialog
from .progress import CheckpointProgress
from .receptor_manager import ReceptorManagerDialog
from .redocking_queue import RedockingQueueController, RedockingQueueDialog
from .receptor_wizard import ReceptorPreparationWizard, ReceptorPreparationWorker
from .results_dialog import DockingResultsDialog
from .settings_dialog import SettingsDialog
from .workers import PipelineWorker, TaskWorker
from ..dashboard import collect_dashboard


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.repo: ProjectRepository | None = None
        self.stage_running = False
        self.pipeline_thread = None
        self.pipeline_worker = None
        self._close_after_pipeline = False
        self.stage_started_at = 0.0
        self.redocking_queue: RedockingQueueController | None = None
        self.setWindowTitle("MolDockPipe Redux")
        self.resize(1280, 760)
        self._closing = False
        self._project_task = None
        self._dashboard_worker = None
        self._dashboard_pending = False
        self._dashboard_timer = QTimer(self)
        self._dashboard_timer.setSingleShot(True)
        self._dashboard_timer.timeout.connect(self._start_dashboard_snapshot)
        self._gui_last_tick = time.monotonic()
        self._gui_max_gap = 0.0
        self._snapshot_seconds = 0.0
        self._heartbeat = QTimer(self)
        self._heartbeat.timeout.connect(self._record_gui_heartbeat)
        self._heartbeat.start(250)
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
        import_action = toolbar.addAction("Import Ligands")
        import_action.triggered.connect(self.import_csv)
        self.receptors_action = toolbar.addAction("Receptors")
        self.receptors_action.setEnabled(False)
        self.receptors_action.triggered.connect(self.manage_receptors)
        tools_button = QToolButton()
        tools_button.setText("Tools")
        tools_menu = QMenu(tools_button)
        tools_menu.addAction("Troubleshooter…", self.show_troubleshooter)
        tools_menu.addAction("Save performance report…", self.save_performance_report)
        self.prepare_receptor_action = tools_menu.addAction("Receptor Preparation…")
        self.prepare_receptor_action.setEnabled(False)
        self.prepare_receptor_action.triggered.connect(self.prepare_receptor)
        self.validation_queue_action = tools_menu.addAction("Redocking Validation Queue…")
        self.validation_queue_action.setEnabled(False)
        self.validation_queue_action.triggered.connect(self.show_validation_queue)
        tools_button.setMenu(tools_menu)
        tools_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(tools_button)
        self.results_action = toolbar.addAction("Results")
        self.results_action.setEnabled(False)
        self.results_action.triggered.connect(self.show_multidock_results)
        advanced_button = QToolButton()
        advanced_button.setText("Advanced")
        advanced_menu = QMenu(advanced_button)
        advanced_menu.addAction("Run Screening", self.run_screening)
        advanced_menu.addAction("Generate States", self.run_molscrub)
        advanced_menu.addAction("Prepare Ligands", self.run_meeko)
        advanced_menu.addAction("Run Docking", self.run_vina)
        advanced_button.setMenu(advanced_menu)
        advanced_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(advanced_button)
        export_button = QToolButton()
        export_button.setText("Export")
        export_button.setToolTip("Export manifest or leaderboard CSV")
        export_menu = QMenu(export_button)
        export_menu.addAction("Manifest CSV", self.export_manifest_csv)
        export_menu.addAction("Leaderboard CSV", self.export_leaderboard_csv)
        export_menu.addAction("Project report (HTML)", self.export_project_report)
        export_menu.addSeparator()
        self.docked_export_action = export_menu.addAction("Export docked compounds (PDBQT/SDF + receptor PDB)", self.select_compounds_for_export)
        export_button.setMenu(export_menu)
        export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        toolbar.addWidget(export_button)
        toolbar.addSeparator()
        self.settings_action = toolbar.addAction("Settings")
        self.settings_action.setEnabled(False)
        self.settings_action.triggered.connect(self.edit_docking_settings)

        pipeline_toolbar = QToolBar("Pipeline Controls")
        pipeline_toolbar.setMovable(False)
        self.addToolBarBreak(Qt.ToolBarArea.TopToolBarArea)
        self.addToolBar(pipeline_toolbar)
        self.run_all_action = pipeline_toolbar.addAction("Run All")
        self.run_all_action.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay))
        self.run_all_action.setEnabled(False)
        self.run_all_action.triggered.connect(self.run_all)
        self.refresh_action = pipeline_toolbar.addAction("Refresh")
        self.refresh_action.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_BrowserReload))
        self.refresh_action.setEnabled(False)
        self.refresh_action.triggered.connect(self.refresh_dashboard)
        self.stop_action = pipeline_toolbar.addAction("Stop")
        self.stop_action.setIcon(self.style().standardIcon(QStyle.StandardPixmap.SP_MediaStop))
        self.stop_action.setEnabled(False)
        self.stop_action.triggered.connect(self.stop_pipeline)

        self.project_path = QLineEdit()
        self.project_path.setReadOnly(True)
        self.project_path.setPlaceholderText("Create or open a portable project")
        project_label = QLabel("Project")
        self.project_summary = QLabel("Current Project: None")
        self.project_summary.setObjectName("projectSummary")
        project_font = QFont(self.project_summary.font())
        project_font.setPointSize(16)
        project_font.setBold(True)
        self.project_summary.setFont(project_font)
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
            button.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
            button.setEnabled(False)
            self.stage_buttons[key] = button
            stages_row.addWidget(button)
        self.stage_action = self.stage_buttons["screening"]
        self.molscrub_action = self.stage_buttons["molscrub"]
        self.meeko_action = self.stage_buttons["meeko"]
        self.vina_action = self.stage_buttons["vina"]
        self.postdock_action = self.stage_buttons["postdock"]

        self.stat_cards: dict[str, QLabel] = {}
        self.stat_details: dict[str, QLabel] = {}
        stats_row = QHBoxLayout()
        stats_row.setContentsMargins(0, 0, 0, 0)
        stats_row.setSpacing(8)
        for key, title in (("ligands", "Input Ligands"), ("passed", "Passed"), ("failed", "Failed"), ("running", "Pipeline Status")):
            card = QFrame()
            card.setFrameShape(QFrame.Shape.StyledPanel)
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(10, 0, 10, 6)
            title_label = QLabel(title)
            title_font = QFont(title_label.font())
            title_font.setPointSize(14)
            title_font.setBold(True)
            title_label.setFont(title_font)
            title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            card_layout.addWidget(title_label)
            value = QLabel("0")
            value.setObjectName("statValue")
            value.setAlignment(Qt.AlignmentFlag.AlignCenter)
            value.setWordWrap(True)
            value_font = QFont(value.font())
            value_font.setPointSize(30)
            value_font.setBold(True)
            value.setFont(value_font)
            if key == "passed":
                value.setVisible(False)
            card_layout.addWidget(value)
            detail = QLabel("")
            detail.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            detail.setWordWrap(True)
            detail_font = QFont(detail.font())
            detail_font.setPointSize(13 if key == "passed" else 9)
            detail.setFont(detail_font)
            if key == "passed":
                detail.setAlignment(Qt.AlignmentFlag.AlignCenter)
            card_layout.addWidget(detail)
            self.stat_cards[key] = value
            self.stat_details[key] = detail
            stats_row.addWidget(card)

        self.overall_progress = CheckpointProgress()
        self.overall_progress.setVisible(False)
        self.stage_progress = QProgressBar()
        self.stage_progress.setRange(0, 100)
        self.stage_progress.setValue(0)
        self.stage_progress.setFormat("%p%")
        self.stage_progress.setMinimumHeight(28)
        self.current_stage_label = QLabel("Current: Idle")
        self.current_item_label = QLabel("")
        self.current_task_label = QLabel("Current task: None")
        top = QWidget()
        top_layout = QVBoxLayout(top)
        top_layout.setContentsMargins(0, 0, 0, 8)
        top_layout.setSpacing(4)
        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.addWidget(self.project_summary)
        top_layout.addLayout(header)
        top_layout.addLayout(stats_row)
        top_layout.addLayout(stages_row)
        top_layout.addWidget(self.current_task_label)
        top_layout.addWidget(self.dashboard_status)
        top_layout.addWidget(self.stage_progress)

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
        self.log.document().setMaximumBlockCount(1000)
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(120)
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

    def closeEvent(self, event: QCloseEvent) -> None:
        self._closing = True
        timer = getattr(self, "_dashboard_timer", None)
        if timer is not None:
            timer.stop()
        heartbeat = getattr(self, "_heartbeat", None)
        if heartbeat is not None:
            heartbeat.stop()
        if getattr(self, "_dashboard_worker", None) is not None or getattr(self, "_project_task", None) is not None:
            if self._dashboard_worker is not None:
                self._dashboard_worker.requestInterruption()
            event.ignore()
            if self._pipeline_active():
                self._close_after_pipeline = True
                self.stop_pipeline()
            return
        if self._pipeline_active():
            event.ignore()
            self._close_after_pipeline = True
            self.stop_pipeline()
            return
        if self.redocking_queue and self.redocking_queue.is_running and not self.redocking_queue.shutdown():
            self._closing = False
            self._heartbeat.start(250)
            event.ignore()
            QMessageBox.warning(self, "Validation is stopping",
                "MolDockPipe is still stopping the active validation process. Try closing the application again in a few seconds.")
            return
        super().closeEvent(event)

    def _pipeline_active(self) -> bool:
        thread = getattr(self, "pipeline_thread", None)
        return self.stage_running or (thread is not None and thread.isRunning())

    def _pipeline_thread_finished(self) -> None:
        self.pipeline_thread = None
        self.pipeline_worker = None
        if self._close_after_pipeline:
            self._close_after_pipeline = False
            self.close()

    def _start_project_task(self, function, ready, label):
        if self._project_task is not None or self._closing:
            return
        worker = TaskWorker(function, self)
        self._project_task = worker
        self.statusBar().showMessage(label)
        self._set_stage_actions(False)
        self.run_all_action.setEnabled(False)
        self.refresh_action.setEnabled(False)
        worker.succeeded.connect(lambda result: ready(result) if not self._closing else None)
        worker.failed.connect(lambda message: QMessageBox.critical(self, label, message) if not self._closing else None)
        worker.finished.connect(self._project_task_finished)
        worker.start()

    def _project_task_finished(self):
        worker, self._project_task = self._project_task, None
        worker.deleteLater()
        self._set_stage_actions(not self._pipeline_active())
        self.run_all_action.setEnabled(self.repo is not None and not self._pipeline_active())
        self.refresh_action.setEnabled(self.repo is not None and not self._pipeline_active())
        self.statusBar().clearMessage()
        if self._closing:
            self.close()

    def _block_project_change(self) -> bool:
        if self._closing:
            return True
        if self._project_task is not None:
            self.statusBar().showMessage("Wait for the current project operation to finish.", 5000)
            return True
        if self._pipeline_active():
            QMessageBox.warning(self, "Pipeline running", "Finish or stop the pipeline before changing projects or importing ligands.")
            return True
        return False

    def create_project(self) -> None:
        if self._block_project_change():
            return
        if self.redocking_queue and self.redocking_queue.is_running:
            QMessageBox.warning(self, "Validation running", "Wait for or cancel the active redocking validation before switching projects.")
            return
        name, accepted = QInputDialog.getText(self, "New Project", "Project name:")
        name = name.strip()
        if not accepted or not name:
            return
        if name in {".", ".."} or Path(name).name != name or any(character in name for character in '<>:"/\\|?*'):
            QMessageBox.warning(self, "Invalid project name", "Choose a name without path separators or reserved filename characters.")
            return
        path = Path.cwd() / "Projects" / name
        if path.exists():
            QMessageBox.warning(self, "Project already exists", f"A project named '{name}' already exists in {path.parent}.")
            return
        def ready(repo):
            self.repo = repo
            self._activate_project("Created")
        self._start_project_task(lambda: ProjectRepository.create(path), ready, "Creating project...")

    def open_project(self) -> None:
        if self._block_project_change():
            return
        if self.redocking_queue and self.redocking_queue.is_running:
            QMessageBox.warning(self, "Validation running", "Wait for or cancel the active redocking validation before switching projects.")
            return
        path = QFileDialog.getExistingDirectory(self, "Open project folder")
        if not path:
            return
        repo = ProjectRepository(Path(path))
        if not repo.database_path.exists():
            QMessageBox.warning(self, "Not a project", "This folder does not contain project.sqlite.")
            return
        def load():
            with project_execution_lock(repo.root):
                repo.migrate()
                interrupted = repo.recover_interrupted_runs()
                repo.import_project_input()
            return interrupted

        def ready(interrupted):
            self.repo = repo
            self._activate_project("Opened" + (f"; recovered {interrupted} interrupted stage runs" if interrupted else ""))
        self._start_project_task(load, ready, "Opening project...")

    def import_csv(self) -> None:
        if self._block_project_change():
            return
        if not self.repo:
            return
        path, _ = QFileDialog.getOpenFileName(self, "Import ligand CSV", filter="CSV files (*.csv)")
        if not path:
            return
        csv_path = Path(path)
        repo = self.repo
        def ready(result):
            self._write_log(
                f"Imported into project library from {path}: added={result.added}, unchanged={result.unchanged}, "
                f"metadata updated={result.changed}"
            )
            self._refresh_checkpoint_state()
        self._start_project_task(lambda: repo.sync_ligands_csv(csv_path), ready, "Importing ligands...")

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
        self.refresh_action.setEnabled(True)
        self.receptors_action.setEnabled(True)
        self.prepare_receptor_action.setEnabled(True)
        self.validation_queue_action.setEnabled(True)
        self.results_action.setEnabled(True)
        if self.redocking_queue:
            self.redocking_queue.deleteLater()
        self.redocking_queue = RedockingQueueController(self.repo, self)
        self.redocking_queue.activity.connect(self._write_log)
        self.redocking_queue.activity.connect(lambda message: self.statusBar().showMessage(message, 10000))
        self.redocking_queue.validation_completed.connect(
            lambda result: self._write_log(f"Redocking artifacts ready: run {result.get('run_id')}")
        )
        self.redocking_queue.start()
        self._write_log(f"{action} project: {self.repo.root}")
        self.refresh_tables()
        self._refresh_checkpoint_state()

    def _auto_import_input(self) -> None:
        if not self.repo:
            return
        candidate = self._project_input_csv()
        if candidate and MainWindow._input_csv_has_ligand_rows(candidate):
            try:
                result = self.repo.import_project_input()
            except Exception as exc:
                QMessageBox.critical(self, "Input import failed", str(exc))
                return
            self._write_log(f"Imported project input: added={result.added}, unchanged={result.unchanged}")

    @staticmethod
    def _input_csv_has_ligand_rows(path: Path) -> bool:
        """Return whether a CSV contains data rows rather than only its header."""
        with path.open(encoding="utf-8-sig") as handle:
            return sum(1 for line in handle if line.strip()) > 1

    def _project_input_csv(self) -> Path | None:
        """Return the project-owned ligand input CSV, including the legacy path."""
        assert self.repo
        for candidate in (self.repo.root / "inputs" / "input.csv", self.repo.root / "input" / "input.csv"):
            if candidate.is_file():
                return candidate
        return None

    def refresh_tables(self) -> None:
        if not self.repo:
            return
        # These views are currently hidden. Populating thousands of invisible
        # table cells and sizing their columns still blocks the GUI thread.
        if not self.data_tabs.isVisible():
            return
        with self.repo.connection() as conn:
            parents = conn.execute("""SELECT p.parent_id, p.source_smiles, p.canonical_source_smiles, p.parent_inchikey,
                p.parse_status, COALESCE(s.decision, 'not screened') AS decision
                FROM parent_ligands p LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
                ORDER BY p.parent_id""").fetchall()
            states = conn.execute("SELECT state_id, parent_id, state_isomeric_smiles, formal_charge, status, enumeration_truncated FROM molecular_states ORDER BY parent_id, state_id").fetchall()
            results = conn.execute("""SELECT s.parent_id, MIN(p.affinity) score, s.state_id,
                SUM(CASE WHEN s.status='prepared' THEN 1 ELSE 0 END) prepared, COUNT(p.pose_id) docked
                FROM molecular_states s LEFT JOIN docking_runs d ON d.state_id=s.state_id AND d.is_current=1
                LEFT JOIN docking_poses p ON p.run_id=d.run_id WHERE s.active=1 GROUP BY s.parent_id, s.state_id ORDER BY score""").fetchall()
        self._load_table(self.parent_table, parents)
        self._load_table(self.state_table, states)
        self._load_table(self.results_table, results)

    def _refresh_checkpoint_state(self) -> None:
        self._request_dashboard()

    def _refresh_live_dashboard(self) -> None:
        self._request_dashboard()

    def _request_dashboard(self) -> None:
        if not self.repo or self._closing:
            return
        self._dashboard_pending = True
        if self._dashboard_worker is not None and self._dashboard_context != (self.repo, self.stage_running):
            self._dashboard_worker.requestInterruption()
        if self._dashboard_worker is not None or self._dashboard_timer.isActive():
            return
        self._start_dashboard_snapshot()

    def _start_dashboard_snapshot(self) -> None:
        if not self.repo or self._closing:
            return
        self._dashboard_pending = False
        repo, running = self.repo, self.stage_running
        self._snapshot_started = time.monotonic()
        worker = TaskWorker(lambda: collect_dashboard(repo, running, worker.isInterruptionRequested), self)
        self._dashboard_context = (repo, running)
        self._dashboard_worker = worker
        worker.succeeded.connect(lambda snapshot: self._dashboard_ready(repo, running, snapshot))
        worker.failed.connect(lambda message: self.statusBar().showMessage(f"Dashboard update failed: {message}", 5000)
                              if not self._closing and not worker.isInterruptionRequested() else None)
        worker.finished.connect(self._dashboard_finished)
        worker.start()

    def _dashboard_ready(self, repo, running, snapshot) -> None:
        if not self._closing and self.repo is repo and self.stage_running == running:
            self._apply_dashboard(snapshot, running)
        else:
            self._dashboard_pending = True

    def _dashboard_finished(self) -> None:
        self._snapshot_seconds = time.monotonic() - self._snapshot_started
        worker, self._dashboard_worker = self._dashboard_worker, None
        worker.deleteLater()
        # Delay from completion, so slow queries cannot create a refresh storm.
        if not self._closing and (self.stage_running or self._dashboard_pending):
            self._dashboard_timer.start(1000)
        if self._closing:
            self.close()

    def _apply_dashboard(self, snapshot, running) -> None:
        parents = snapshot["parents"]
        screened = snapshot["screened"]
        states = snapshot["states"]
        prepared = snapshot["prepared"]
        docked = snapshot["docked"]
        structure_total = snapshot["structure_total"]
        structure_passed = snapshot["structure_passed"]
        structure_failed = snapshot["structure_failed"]
        conformers = snapshot["conformers"]
        meeko_terminal = snapshot["meeko_terminal"]
        expected_docking = snapshot["expected_docking"]
        docking_terminal = snapshot["docking_terminal"]
        passed = snapshot["passed"]
        failed = snapshot["failed"]
        failed_ligands = snapshot["failed_ligands"]
        pdbqt_failed = snapshot["pdbqt_failed"]
        pdbqt_passed = snapshot["pdbqt_passed"]
        docking_failed = snapshot["docking_failed"]
        docking_passed = snapshot["docking_passed"]
        stage_runs = snapshot["stage_runs"]
        structure_failures = snapshot["structure_failures"]
        has_docking_exports = snapshot["has_docking_exports"]
        export_status = snapshot["export_status"]
        self._dashboard_values = (parents, screened, states, prepared, docked)
        if not running:
            self.overall_progress.reset()
            if parents and screened == parents:
                self.overall_progress.complete_stage("screening")
            elif screened:
                self.overall_progress.start_stage("screening")
            if structure_total and structure_passed + structure_failed >= structure_total:
                if structure_failed:
                    self.overall_progress.states["molscrub"] = "warning"
                    self.overall_progress.update()
                else:
                    self.overall_progress.complete_stage("molscrub")
            elif structure_passed or structure_failed:
                self.overall_progress.start_stage("molscrub")
            if conformers and meeko_terminal == conformers:
                self.overall_progress.complete_stage("meeko")
            elif prepared:
                self.overall_progress.start_stage("meeko")
            if expected_docking and docking_terminal == expected_docking:
                self.overall_progress.complete_stage("vina")
            elif docked:
                self.overall_progress.start_stage("vina")
        self.stat_cards["ligands"].setText(str(parents))
        self.stat_cards["passed"].setText(str(passed))
        self.stat_cards["failed"].setText(str(failed_ligands))
        self.stat_cards["running"].setText("Running" if self.stage_running else "Idle" if stage_runs else "Staging")
        screening_percent = (failed / parents * 100) if parents else 0
        self.stat_details["passed"].setText(
            f"Successfully Screened: {passed}\n"
            f"Successful 3D Gen: {structure_passed}\n"
            f"Successful PDBQT: {pdbqt_passed}\n"
            f"Successful Docked Pairs: {docking_passed}"
        )
        self.stat_details["failed"].setText(
            f"Failed Screening: {failed} ({screening_percent:.1f}% vs total ligands)\n"
            f"Failed 3D Structure Generation: {structure_failed}\n"
            f"Failed PDBQT Generation: {pdbqt_failed}\n"
            f"Failed Docking: {docking_failed}"
        )
        self.project_summary.setText(f"Current Project: {self.repo.root.name}")
        self.docked_export_action.setEnabled(has_docking_exports)
        labels = {
            "screening": f"Screening\n{passed} passed; {failed} failed",
            "molscrub": f"States\n{states} generated; {structure_failed} failed",
            "meeko": f"Ligands\n{prepared} prepared",
            "vina": f"Docking\n{docked} receptor-state runs complete",
            "postdock": f"Export\n{export_status}",
        }
        for stage, button in self.stage_buttons.items():
            status = self.overall_progress.states.get(stage, "pending")
            icon = "\u2713" if status == "complete" else "\u27f3" if status == "active" else "\u26a0" if status == "warning" else "\u25cb"
            button.setText(f"{icon} {labels[stage]}")
            button.setToolTip(labels[stage].replace("\n", "\n"))
            if stage == "molscrub" and structure_failures:
                button.setToolTip(labels[stage] + "\n" + "\n".join(
                    f"{row['parent_id']}: {row['reason']}" for row in structure_failures[:6]))
            if stage == "postdock":
                button.setEnabled(has_docking_exports)
            if status == "complete":
                button.setStyleSheet("QPushButton { color: #166534; background: #dcfce7; border: 1px solid #86efac; }")
            elif status == "active":
                button.setStyleSheet("QPushButton { color: #1d4ed8; background: #dbeafe; border: 2px solid #3b82f6; font-weight: 600; }")
            elif status == "warning":
                button.setStyleSheet("QPushButton { color: #92400e; background: #fef3c7; border: 1px solid #f59e0b; }")
            else:
                button.setStyleSheet("QPushButton { color: #374151; background: #f3f4f6; border: 1px solid #d1d5db; }")

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

    def show_troubleshooter(self) -> None:
        from .troubleshooter import TroubleshooterDialog
        from ..services.executables import application_root

        TroubleshooterDialog(self.repo.root if self.repo else application_root(), self).exec()

    def _record_gui_heartbeat(self):
        now = time.monotonic()
        self._gui_max_gap = max(self._gui_max_gap, now - self._gui_last_tick)
        self._gui_last_tick = now

    def save_performance_report(self):
        if not self.repo or self._project_task is not None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save performance report", "moldockpipe-performance.json", "JSON (*.json)")
        if not path:
            return
        from ..performance import save_performance_report
        repo = self.repo
        metrics = {"heartbeat_max_gap_seconds": self._gui_max_gap,
                   "last_dashboard_seconds": self._snapshot_seconds, "pipeline_running": self.stage_running}
        self._start_project_task(lambda: save_performance_report(repo, Path(path), metrics),
            lambda result: self._write_log(f"Saved performance report: {result}"), "Saving performance report...")

    def run_molscrub(self) -> None:
        self._start_pipeline(["molscrub"], "MolScrub")

    def run_meeko(self) -> None:
        self._start_pipeline(["meeko"], "Meeko")

    def _start_pipeline(self, stages: list[str], label: str) -> None:
        if not self.repo or self._pipeline_active() or self._closing or self._project_task is not None:
            return
        self.stage_running = True
        self.stage_started_at = time.monotonic()
        self._last_live_refresh_at = float("-inf")
        self.pipeline_is_full = not stages
        self.overall_progress.start_stage(stages[0] if stages else "screening")
        self._set_stage_actions(False)
        self.run_all_action.setEnabled(False)
        self.refresh_action.setEnabled(False)
        self.stop_action.setEnabled(True)
        self.current_stage_label.setText(f"{label} running")
        self.current_task_label.setText(f"Current task: {label}")
        self.dashboard_status.setText("Running")
        self._write_log(f"{label} started in background...")
        self.pipeline_thread = QThread(self)
        self.pipeline_worker = PipelineWorker(self.repo, stages)
        self.pipeline_worker.moveToThread(self.pipeline_thread)
        self.pipeline_thread.started.connect(self.pipeline_worker.run)
        self.pipeline_worker.progress.connect(self._on_pipeline_progress)
        self.pipeline_worker.finished.connect(self._pipeline_finished)
        self.pipeline_worker.failed.connect(self._pipeline_failed)
        self.pipeline_worker.stopped.connect(self._pipeline_stopped)
        self.pipeline_worker.finished.connect(self.pipeline_thread.quit)
        self.pipeline_worker.failed.connect(self.pipeline_thread.quit)
        self.pipeline_worker.stopped.connect(self.pipeline_thread.quit)
        self.pipeline_thread.finished.connect(self.pipeline_worker.deleteLater)
        self.pipeline_thread.finished.connect(self._pipeline_thread_finished)
        self.pipeline_thread.finished.connect(self.pipeline_thread.deleteLater)
        self.pipeline_thread.start()

    def stop_pipeline(self) -> None:
        if not self.stage_running or self.pipeline_worker is None:
            return
        self.stop_action.setEnabled(False)
        self.dashboard_status.setText("Stopping safely…")
        self.current_task_label.setText("Current task: stopping safely")
        self._write_log("Stop requested; waiting for active work to clean up...")
        self.pipeline_worker.stop()

    def _pipeline_stopped(self) -> None:
        self.stage_running = False
        self.stage_progress.setRange(0, 100)
        self.stage_progress.setFormat("Stopped")
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        self.refresh_action.setEnabled(True)
        self.stop_action.setEnabled(False)
        self.current_stage_label.setText("Idle")
        self.current_item_label.setText("")
        self.current_task_label.setText("Current task: None")
        self.dashboard_status.setText("Pipeline Stopped Safely")
        self.statusBar().showMessage("Pipeline Stopped Safely", 10000)
        self._write_log("Pipeline Stopped Safely")
        self.refresh_tables()
        self._refresh_checkpoint_state()

    def run_vina(self) -> None:
        self._start_pipeline(["vina"], "Vina")

    def run_postdock(self) -> None:
        if self.repo and not self.repo.get_settings().get("postdock", {}).get("selected_parents"):
            QMessageBox.warning(self, "Select compounds", "Choose compounds from the Export menu before starting docked-compound export.")
            return
        self._start_pipeline(["postdock"], "Post-docking")

    def select_compounds_for_export(self) -> None:
        if not self.repo or self.stage_running:
            return
        repo = self.repo
        def load():
            profiles = [str(profile['id']) for profile in repo.get_receptor_profiles() if profile.get('enabled')]
            with repo.connection() as conn:
                rows = [dict(row) for row in conn.execute(f"""SELECT s.parent_id, MIN(p.affinity) best_affinity
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    JOIN parent_ligands l ON l.parent_id=s.parent_id JOIN docking_poses p ON p.run_id=d.run_id
                    WHERE s.active=1 AND d.status='completed' AND d.is_current=1
                    AND d.receptor_profile_id IN ({','.join('?' for _ in profiles)})
                    GROUP BY s.parent_id ORDER BY best_affinity,s.parent_id""", profiles)] if profiles else []
            for row in rows:
                # Preserve the selector's legacy export markers. A partially
                # exported receptor campaign must remain selectable for resume.
                row['exported'] = any((repo.root / 'For_PostDocking' / folder / row['parent_id']).exists()
                                      for folder in ('SDF', 'PDBQTs'))
            return rows, repo.get_settings().get('postdock', {}).get('selected_parents', [])

        def ready(result):
            rows, current = result
            dialog = CompoundSelectorDialog(rows, set(current), repo, self)
            if dialog.exec() != QDialog.DialogCode.Accepted:
                return
            selected = dialog.selected()
            if not selected:
                QMessageBox.warning(self, "Select compounds", "Select at least one successfully docked compound to export.")
                return
            repo.update_settings({'postdock': {'selected_parents': selected}})
            self._write_log(f"Export selected compounds: {len(selected)}")
            QTimer.singleShot(0, self.run_postdock)
        self._start_project_task(load, ready, "Loading docked compounds...")

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
        self.receptors_action.setEnabled(enabled)

    def _on_pipeline_progress(self, event: object) -> None:
        kind = getattr(event, "event", "")
        stage = getattr(event, "stage", "")
        index = getattr(event, "index", 0)
        total = getattr(event, "total", 0)
        if kind == "stage_preparing":
            self.stage_started_at = time.monotonic()
            self.stage_progress.setRange(0, 0)
            self.stage_progress.setFormat("Preparing…")
        elif kind == "stage_completed":
            index = total
            self.stage_progress.setRange(0, 100)
            self.stage_progress.setValue(100)
            self.stage_progress.setFormat(f"{total} / {total}   %p%")
        elif total:
            self.stage_progress.setRange(0, 100)
            self.stage_progress.setValue(int(index * 100 / total))
            self.stage_progress.setFormat(f"{index} / {total}   %p%")
        if getattr(event, "event", "") == "stage_started":
            self.stage_progress.setRange(0, 100)
            self.stage_started_at = time.monotonic()
            self.overall_progress.start_stage(stage)
        elif getattr(event, "event", "") == "stage_completed":
            if getattr(event, "failed", 0):
                self.overall_progress.states[stage] = "warning"
                self.overall_progress.update()
            else:
                self.overall_progress.complete_stage(stage)
        stage_names = {"screening": "Drug-likeness Screening", "molscrub": "Molecular States", "meeko": "Ligand Preparation", "vina": "Docking", "postdock": "Export"}
        display_stage = stage_names.get(stage, stage.title())
        self.current_stage_label.setText(display_stage)
        profile_name = getattr(event, "receptor_profile_name", None)
        item_text = str(getattr(event, "item_id", "") or "")
        self.current_item_label.setText(f"{profile_name}: {item_text}" if profile_name else item_text)
        detail = getattr(event, "message", "") if kind in {"stage_preparing", "stage_started", "planning_completed"} else self.current_item_label.text()
        self.current_task_label.setText(f"{display_stage}\n{detail}\n{index} / {total}" if kind != "stage_preparing" else f"{display_stage}\n{detail}")
        self.dashboard_status.setText(self.current_stage_label.text())
        self.stat_cards["running"].setText("Running")
        if kind == "stage_preparing":
            self.statusBar().showMessage(f"{stage}: {detail}")
            self._write_log(f"{stage}: {detail}")
            # Show the handoff immediately, without blocking behind dashboard I/O.
            return
        if total and index:
            elapsed = max(time.monotonic() - self.stage_started_at, 0.1)
            rate = index / elapsed * 60
            remaining = max(total - index, 0)
            eta = remaining / (index / elapsed) if index else 0
            self.dashboard_status.setText(f"{stage.title()} | {rate:.1f}/min | ETA {eta / 60:.1f} min")
        # Item events arrive in bursts, especially when reusing completed work.
        # Full-project queries and export-file checks must not run per ligand.
        now = time.monotonic()
        if (getattr(event, "event", "") == "stage_completed"
                or now - getattr(self, "_last_live_refresh_at", float("-inf")) >= 1.0):
            self._last_live_refresh_at = now
            self._refresh_live_dashboard()
        item = getattr(event, "item_id", None)
        message = getattr(event, "message", "")
        self.statusBar().showMessage(f"{stage}: {item or ''} {message}")
        if (getattr(event, "event", "") in {"item_failed", "stage_completed"}
                or (getattr(event, "event", "") == "item_skipped" and getattr(event, "failed", 0))):
            self._write_log(self._event_summary(event))
        if getattr(event, "event", "") == "stage_completed":
            self.refresh_tables()
            if not getattr(self, "stage_running", False):
                self._refresh_checkpoint_state()

    def _pipeline_finished(self, summary: object) -> None:
        self.stage_running = False
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        self.refresh_action.setEnabled(True)
        self.stop_action.setEnabled(False)
        if getattr(self, "pipeline_is_full", False):
            self.overall_progress.complete_all()
        else:
            stage = next(iter(summary), "") if isinstance(summary, dict) else ""
            self.overall_progress.complete_stage(stage)
        self.stage_progress.setRange(0, 100)
        self.stage_progress.setValue(100)
        self.current_stage_label.setText("Idle")
        self.current_item_label.setText("")
        self.current_task_label.setText("Current task: None")
        self.dashboard_status.setText("Ready")
        self._write_log(f"Pipeline completed: {summary}")
        self.refresh_tables()
        self._refresh_checkpoint_state()
        if isinstance(summary, dict) and "vina" in summary:
            QMessageBox.information(
                self,
                "Docking complete",
                "Docking has finished. Export is ready.\n\nUse the Export menu, then choose 'Export docked compounds (PDBQT/SDF + receptor PDB)'."
            )

    def _pipeline_failed(self, message: str) -> None:
        self.stage_running = False
        self.stage_progress.setRange(0, 100)
        self.stage_progress.setValue(0)
        self.stage_progress.setFormat("Stopped")
        self.current_stage_label.setText("Failed")
        self.current_item_label.setText("")
        self.current_task_label.setText(f"Stopped: {message}")
        self.dashboard_status.setText("Pipeline failed")
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        self.refresh_action.setEnabled(True)
        self.stop_action.setEnabled(False)
        self._write_log(f"Run All stopped: {message}")
        self._refresh_checkpoint_state()
        QMessageBox.critical(self, "Run All failed", message)

    def refresh_dashboard(self) -> None:
        if not self.repo or self.stage_running:
            return
        if self._project_task is not None:
            return
        repo = self.repo
        def ready(result):
            self._write_log(f"Reloaded ligands: added={result.added}, unchanged={result.unchanged}, metadata updated={result.changed}")
            self._refresh_checkpoint_state()
        self._start_project_task(repo.import_project_input, ready, "Refreshing project inputs...")

    @staticmethod
    def _event_summary(event: object) -> str:
        stage = getattr(event, "stage", "")
        kind = getattr(event, "event", "")
        if kind == "stage_completed":
            if stage == "vina":
                reused = getattr(event, "skipped", 0)
                newly_docked = max(0, getattr(event, "succeeded", 0) - reused)
                return (f"Docking completed: newly docked={newly_docked}, "
                        f"reused={reused}, failed={getattr(event, 'failed', 0)}")
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
            self.repo.update_settings(dialog.values())
            self._write_log(f"Updated docking settings in {config_path}")
        except Exception as exc:
            QMessageBox.critical(self, "Settings update failed", str(exc))

    def manage_receptors(self) -> None:
        if not self.repo or self.stage_running:
            return
        assert self.redocking_queue
        dialog = ReceptorManagerDialog(self.repo, self, self.redocking_queue)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            enabled = sum(1 for profile in self.repo.get_receptor_profiles() if profile.get("enabled"))
            self._write_log(f"Updated receptor profiles: {enabled} enabled")
            self._refresh_checkpoint_state()

    def show_validation_queue(self) -> None:
        if self.redocking_queue:
            RedockingQueueDialog(self.redocking_queue, self).exec()

    def prepare_receptor(self) -> None:
        if not self.repo or self.stage_running:
            return
        wizard = ReceptorPreparationWizard(self)
        if wizard.exec() != QDialog.DialogCode.Accepted:
            return
        plan = wizard.plan()
        if any(str(profile.get("name", "")).casefold() == plan.profile_name.casefold()
               for profile in self.repo.get_receptor_profiles() if not profile.get("archived")):
            QMessageBox.warning(self, "Duplicate receptor name", "Choose a unique receptor profile name.")
            return
        activate = wizard.finish_page.activate.isChecked()
        progress = QProgressDialog("Preparing receptor with Meeko…", "", 0, 0, self)
        progress.setWindowTitle("Receptor Preparation")
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        thread = QThread(self)
        worker = ReceptorPreparationWorker(self.repo.root, plan)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)

        def succeeded(folder: str) -> None:
            assert self.repo
            profiles = self.repo.get_receptor_profiles(include_archived=True)
            expected_profiles = list(profiles)
            profiles.append({
                "id": plan.profile_id, "name": plan.profile_name, "enabled": activate, "archived": False, "protocol": plan.protocol,
                "receptor": (Path(folder) / "receptor.pdbqt").relative_to(self.repo.root).as_posix(),
                "center_x": plan.box_center[0], "center_y": plan.box_center[1], "center_z": plan.box_center[2],
                "size_x": plan.box_size[0], "size_y": plan.box_size[1], "size_z": plan.box_size[2],
                "exhaustiveness": 8, "num_modes": 9, "energy_range": 3, "seed": 42, "cpu_count": logical_cpu_count(),
                "reference_ligand": {
                    "identity": plan.reference_ligand.label() if plan.reference_ligand else None,
                    "name": plan.reference_ligand.name if plan.reference_ligand else None,
                    "chain": plan.reference_ligand.chain if plan.reference_ligand else None,
                    "residue": plan.reference_ligand.number if plan.reference_ligand else None,
                    "pdb": f"inputs/receptors/{plan.profile_id}/reference_ligand/reference_ligand.pdb" if plan.reference_ligand else None,
                    "sdf": f"inputs/receptors/{plan.profile_id}/reference_ligand/reference_ligand.sdf" if plan.reference_ligand else None,
                    "mol2": f"inputs/receptors/{plan.profile_id}/reference_ligand/reference_ligand.mol2" if plan.reference_ligand else None,
                    "mapping": f"inputs/receptors/{plan.profile_id}/reference_ligand/chemistry_mapping.json" if plan.reference_ligand else None,
                },
            })
            try:
                self.repo.save_receptor_profiles(profiles, expected_profiles=expected_profiles)
            except (RuntimeError, ValueError, OSError) as exc:
                QMessageBox.warning(self, "Receptor profile was not registered",
                    f"{exc}\n\nPrepared files are saved in {folder}.")
                return
            self._write_log(f"Created receptor profile '{plan.profile_name}' in {folder}")
            self._refresh_checkpoint_state()
            QMessageBox.information(self, "Receptor profile created", f"{plan.profile_name}\n\nBox: {plan.box_size[0]:.1f} × {plan.box_size[1]:.1f} × {plan.box_size[2]:.1f} Å")

        worker.succeeded.connect(succeeded)
        worker.failed.connect(lambda message: QMessageBox.critical(self, "Receptor preparation failed", message))
        worker.finished.connect(progress.close)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(lambda: self.prepare_receptor_action.setEnabled(True))
        self._receptor_preparation_thread = thread
        self._receptor_preparation_worker = worker
        self.prepare_receptor_action.setEnabled(False)
        thread.start()

    def show_multidock_results(self) -> None:
        if self.repo:
            dialog = DockingResultsDialog(self.repo, self)
            try:
                dialog.exec()
            finally:
                dialog.deleteLater()

    def _on_workflow_reset(self) -> None:
        if not self.repo:
            return
        project_root = self.repo.root
        self.repo = ProjectRepository(project_root)
        self.repo.migrate()
        self.repo.recover_interrupted_runs()
        self._activate_project("Reloaded")

    def export_manifest_csv(self) -> None:
        if not self.repo:
            return
        try:
            outputs = write_manifests(self.repo)
        except Exception as exc:
            QMessageBox.critical(self, "Manifest export failed", str(exc)); return
        self._write_log(f"Exported manifest CSV files: {len(outputs)}")
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.repo.root / "exports")))
        QMessageBox.information(self, "Manifest export complete",
            f"Combined manifest:\n{self.repo.root / 'exports' / 'manifest.csv'}\n\n"
            f"Also wrote {len(outputs) - 1} receptor-specific manifest files.")

    def export_leaderboard_csv(self) -> None:
        if not self.repo:
            return
        try:
            outputs = write_leaderboards(self.repo)
        except Exception as exc:
            QMessageBox.critical(self, "Leaderboard export failed", str(exc)); return
        self._write_log(f"Exported leaderboard CSV files: {len(outputs)}")
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.repo.root / "exports")))
        QMessageBox.information(self, "Leaderboard export complete",
            f"Combined leaderboard:\n{self.repo.root / 'exports' / 'leaderboard.csv'}\n\n"
            f"Also wrote {len(outputs) - 1} receptor-specific leaderboard files.")

    def export_project_report(self) -> None:
        if not self.repo or self.stage_running:
            return
        try:
            output = generate_project_report(self.repo)
        except Exception as exc:
            QMessageBox.critical(self, "Report generation failed", str(exc))
            return
        self._write_log(f"Generated project report: {output}")
        QMessageBox.information(self, "Report generated",
                                f"The HTML report and its provenance files were written to:\n\n{output.parent}")
