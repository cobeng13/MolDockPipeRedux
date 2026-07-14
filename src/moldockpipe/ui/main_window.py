from __future__ import annotations

import time
from pathlib import Path

from PyQt6.QtCore import QThread, Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QDialog, QFileDialog, QFrame, QHBoxLayout, QInputDialog, QLabel, QLineEdit, QMainWindow, QMessageBox,
    QPushButton, QProgressBar, QSplitter, QStatusBar, QStyle, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QToolBar, QToolButton, QVBoxLayout, QWidget, QMenu,
)

import yaml

from ..project import ProjectRepository
from .compound_selector import CompoundSelectorDialog
from .progress import CheckpointProgress
from .receptor_manager import ReceptorManagerDialog
from .results_dialog import DockingResultsDialog
from .settings_dialog import SettingsDialog
from .workers import PipelineWorker


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
        import_action = toolbar.addAction("Import Ligands")
        import_action.triggered.connect(self.import_csv)
        self.receptors_action = toolbar.addAction("Receptors")
        self.receptors_action.setEnabled(False)
        self.receptors_action.triggered.connect(self.manage_receptors)
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
        export_menu.addSeparator()
        self.docked_export_action = export_menu.addAction("Export docked compounds (PDBQT/SDF)", self.select_compounds_for_export)
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

    def create_project(self) -> None:
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
        self.repo = ProjectRepository.create(path)
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
        csv_path = Path(path)
        try:
            preview = self.repo.preview_ligand_sync(csv_path)
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))
            return
        if preview.changed or preview.archived:
            answer = QMessageBox.question(
                self,
                "Synchronize ligand set?",
                f"Add {preview.added}, change {preview.changed}, and archive {preview.archived} compounds?\n\nUnchanged compounds and their results will be preserved.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return
        try:
            result = self.repo.sync_ligands_csv(csv_path)
        except Exception as exc:
            QMessageBox.critical(self, "Import failed", str(exc))
            return
        self._write_log(
            f"Synchronized ligands from {path}: added={result.added}, unchanged={result.unchanged}, "
            f"changed={result.changed}, archived={result.archived}"
        )
        self.refresh_tables()
        self._refresh_checkpoint_state()

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
        self.results_action.setEnabled(True)
        self._auto_import_input()
        self._write_log(f"{action} project: {self.repo.root}")
        self.refresh_tables()
        self._refresh_checkpoint_state()

    def _auto_import_input(self) -> None:
        if not self.repo:
            return
        with self.repo.connection() as conn:
            has_ligands = conn.execute("SELECT EXISTS(SELECT 1 FROM parent_ligands WHERE active=1)").fetchone()[0]
        if has_ligands:
            return
        candidate = self._project_input_csv()
        if candidate:
            count = self.repo.import_ligands_csv(candidate)
            self._write_log(f"Auto-loaded {count} ligands from {candidate}")

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
        with self.repo.connection() as conn:
            parents = conn.execute("""SELECT p.parent_id, p.source_smiles, p.canonical_source_smiles, p.parent_inchikey,
                p.parse_status, COALESCE(s.decision, 'not screened') AS decision
                FROM parent_ligands p LEFT JOIN screening_results s ON s.parent_id=p.parent_id AND s.active=1
                WHERE p.active=1 ORDER BY p.parent_id""").fetchall()
            states = conn.execute("SELECT state_id, parent_id, state_isomeric_smiles, formal_charge, status, enumeration_truncated FROM molecular_states WHERE active=1 ORDER BY parent_id, state_id").fetchall()
            results = conn.execute("""SELECT s.parent_id, MIN(p.affinity) score, s.state_id,
                SUM(CASE WHEN s.status='prepared' THEN 1 ELSE 0 END) prepared, COUNT(p.pose_id) docked
                FROM molecular_states s LEFT JOIN docking_runs d ON d.state_id=s.state_id AND d.is_current=1
                LEFT JOIN docking_poses p ON p.run_id=d.run_id WHERE s.active=1 GROUP BY s.parent_id, s.state_id ORDER BY score""").fetchall()
        self._load_table(self.parent_table, parents)
        self._load_table(self.state_table, states)
        self._load_table(self.results_table, results)

    def _refresh_checkpoint_state(self) -> None:
        if not self.repo:
            return
        profile_ids = [str(profile["id"]) for profile in self.repo.get_receptor_profiles() if profile.get("enabled")]
        placeholders = ",".join("?" for _ in profile_ids)
        with self.repo.connection() as conn:
            parents = conn.execute("SELECT COUNT(*) FROM parent_ligands WHERE active=1").fetchone()[0]
            screened = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE p.active=1 AND s.active=1 AND s.status IN ('completed','failed')").fetchone()[0]
            states = conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0]
            prepared = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
            conformers = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1").fetchone()[0]
            meeko_terminal = conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status IN ('pdbqt_ready','pdbqt_failed')").fetchone()[0]
            prepared_states = conn.execute("SELECT COUNT(DISTINCT s.state_id) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
            if profile_ids:
                docked = conn.execute(f"SELECT COUNT(*) FROM docking_runs d JOIN molecular_states s USING(state_id) WHERE s.active=1 AND d.status='completed' AND d.is_current=1 AND d.receptor_profile_id IN ({placeholders})", profile_ids).fetchone()[0]
                docking_terminal = conn.execute(f"SELECT COUNT(*) FROM docking_runs d JOIN molecular_states s USING(state_id) WHERE s.active=1 AND d.status IN ('completed','failed','interrupted') AND d.is_current=1 AND d.receptor_profile_id IN ({placeholders})", profile_ids).fetchone()[0]
            else:
                docked = docking_terminal = 0
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
        expected_docking = prepared_states * len(profile_ids)
        if expected_docking and docking_terminal == expected_docking:
            self.overall_progress.complete_stage("vina")
        elif docked:
            self.overall_progress.start_stage("vina")
        self._update_dashboard(parents, screened, states, prepared, docked)

    def _update_dashboard(self, parents: int, screened: int, states: int, prepared: int, docked: int) -> None:
        """Refresh the compact statistics and stage labels shown in the ribbon."""
        if not self.repo:
            return
        self._dashboard_values = (parents, screened, states, prepared, docked)
        profile_ids = [str(profile["id"]) for profile in self.repo.get_receptor_profiles() if profile.get("enabled")]
        placeholders = ",".join("?" for _ in profile_ids)
        with self.repo.connection() as conn:
            passed = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE p.active=1 AND s.active=1 AND s.decision != 'fail'").fetchone()[0]
            failed = conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE p.active=1 AND s.active=1 AND s.decision = 'fail'").fetchone()[0]
            structure_failed = conn.execute("SELECT COUNT(*) FROM parent_ligands WHERE active=1 AND parse_reason LIKE 'MolScrub failed:%'").fetchone()[0]
            structure_passed = conn.execute("SELECT COUNT(DISTINCT parent_id) FROM molecular_states WHERE active=1").fetchone()[0]
            pdbqt_failed = conn.execute("SELECT COUNT(DISTINCT s.parent_id) FROM molecular_states s JOIN parent_ligands p USING(parent_id) WHERE p.active=1 AND s.active=1 AND s.status='pdbqt_failed'").fetchone()[0]
            pdbqt_passed = conn.execute("SELECT COUNT(DISTINCT s.parent_id) FROM molecular_states s JOIN parent_ligands p USING(parent_id) JOIN conformers c USING(state_id) WHERE p.active=1 AND s.active=1 AND c.status='pdbqt_ready'").fetchone()[0]
            if profile_ids:
                docking_failed = conn.execute(f"""SELECT COUNT(DISTINCT s.parent_id || ':' || d.receptor_profile_id)
                    FROM docking_runs d JOIN molecular_states s USING(state_id) JOIN parent_ligands p USING(parent_id)
                    WHERE p.active=1 AND s.active=1 AND d.status IN ('failed','interrupted') AND d.is_current=1
                    AND d.receptor_profile_id IN ({placeholders})""", profile_ids).fetchone()[0]
                docking_passed = conn.execute(f"""SELECT COUNT(DISTINCT s.parent_id || ':' || d.receptor_profile_id)
                    FROM docking_runs d JOIN molecular_states s USING(state_id) JOIN parent_ligands p USING(parent_id)
                    WHERE p.active=1 AND s.active=1 AND d.status='completed' AND d.is_current=1
                    AND d.receptor_profile_id IN ({placeholders})""", profile_ids).fetchone()[0]
            else:
                docking_failed = docking_passed = 0
            stage_runs = conn.execute("SELECT COUNT(*) FROM stage_runs").fetchone()[0]
            screening_failures = conn.execute("""SELECT p.parent_id, COALESCE(s.reason, 'Screening rule failure')
                FROM screening_results s JOIN parent_ligands p USING(parent_id)
                WHERE p.active=1 AND s.active=1 AND s.decision='fail' ORDER BY p.parent_id LIMIT 6""").fetchall()
            pdbqt_failures = conn.execute("""SELECT state_id, COALESCE(reason, 'PDBQT generation failed')
                FROM molecular_states WHERE active=1 AND status='pdbqt_failed' ORDER BY state_id LIMIT 6""").fetchall()
            docking_failures = conn.execute("""SELECT s.state_id, CASE
                WHEN COUNT(p.pose_id)=0 THEN 'No valid Vina poses'
                ELSE 'Docking failed' END
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                LEFT JOIN docking_poses p ON p.run_id=d.run_id
                WHERE s.active=1 AND d.status IN ('failed','interrupted') GROUP BY d.run_id, s.state_id ORDER BY s.state_id LIMIT 6""").fetchall()
            docked_pairs = [(row[0], row[1]) for row in conn.execute(f"""SELECT DISTINCT d.receptor_profile_id, s.parent_id
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id JOIN parent_ligands p USING(parent_id)
                WHERE p.active=1 AND s.active=1 AND d.status='completed' AND d.is_current=1
                AND d.receptor_profile_id IN ({placeholders})""", profile_ids).fetchall()] if profile_ids else []
        self.stat_cards["ligands"].setText(str(parents))
        self.stat_cards["passed"].setText(str(passed))
        self.stat_cards["failed"].setText(str(failed))
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
        exported_root = self.repo.root / "For_PostDocking"
        ready_for_export = [(profile_id, parent_id) for profile_id, parent_id in docked_pairs
                            if not any((exported_root / profile_id / folder / parent_id).exists() for folder in ("SDF", "PDBQTs"))]
        export_status = "Ready" if ready_for_export else "Pending"
        self.docked_export_action.setEnabled(bool(ready_for_export))
        labels = {
            "screening": f"Screening\n{passed} passed; {failed} failed",
            "molscrub": f"States\n{states} generated",
            "meeko": f"Ligands\n{prepared} prepared",
            "vina": f"Docking\n{docked} receptor-state runs complete",
            "postdock": f"Export\n{export_status}",
        }
        for stage, button in self.stage_buttons.items():
            status = self.overall_progress.states.get(stage, "pending")
            icon = "\u2713" if status == "complete" else "\u27f3" if status == "active" else "\u25cb"
            button.setText(f"{icon} {labels[stage]}")
            button.setToolTip(labels[stage].replace("\n", "\n"))
            if stage == "postdock":
                button.setEnabled(bool(ready_for_export))
            if status == "complete":
                button.setStyleSheet("QPushButton { color: #166534; background: #dcfce7; border: 1px solid #86efac; }")
            elif status == "active":
                button.setStyleSheet("QPushButton { color: #1d4ed8; background: #dbeafe; border: 2px solid #3b82f6; font-weight: 600; }")
            elif status == "warning":
                button.setStyleSheet("QPushButton { color: #92400e; background: #fef3c7; border: 1px solid #f59e0b; }")
            else:
                button.setStyleSheet("QPushButton { color: #374151; background: #f3f4f6; border: 1px solid #d1d5db; }")

    def _refresh_live_dashboard(self) -> None:
        if not self.repo:
            return
        profile_ids = [str(profile["id"]) for profile in self.repo.get_receptor_profiles() if profile.get("enabled")]
        placeholders = ",".join("?" for _ in profile_ids)
        with self.repo.connection() as conn:
            docked = conn.execute(f"SELECT COUNT(*) FROM docking_runs d JOIN molecular_states s USING(state_id) WHERE s.active=1 AND d.status='completed' AND d.is_current=1 AND d.receptor_profile_id IN ({placeholders})", profile_ids).fetchone()[0] if profile_ids else 0
            values = (
                conn.execute("SELECT COUNT(*) FROM parent_ligands WHERE active=1").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM screening_results s JOIN parent_ligands p USING(parent_id) WHERE p.active=1 AND s.active=1 AND s.status IN ('completed','failed')").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0],
                conn.execute("SELECT COUNT(*) FROM conformers c JOIN molecular_states s USING(state_id) WHERE s.active=1 AND c.status='pdbqt_ready'").fetchone()[0],
                docked,
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
        self.refresh_action.setEnabled(False)
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
        self.pipeline_worker.finished.connect(self.pipeline_thread.quit)
        self.pipeline_worker.failed.connect(self.pipeline_thread.quit)
        self.pipeline_thread.finished.connect(self.pipeline_worker.deleteLater)
        self.pipeline_thread.finished.connect(self.pipeline_thread.deleteLater)
        self.pipeline_thread.start()

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
        profile_ids = [str(profile["id"]) for profile in self.repo.get_receptor_profiles() if profile.get("enabled")]
        placeholders = ",".join("?" for _ in profile_ids)
        with self.repo.connection() as conn:
            rows = conn.execute(f"""SELECT s.parent_id, MIN(p.affinity) best_affinity
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                JOIN parent_ligands l ON l.parent_id=s.parent_id
                JOIN docking_poses p ON p.run_id=d.run_id
                WHERE l.active=1 AND s.active=1 AND d.status='completed' AND d.is_current=1
                AND d.receptor_profile_id IN ({placeholders})
                GROUP BY s.parent_id ORDER BY best_affinity ASC, s.parent_id""", profile_ids).fetchall() if profile_ids else []
        current = self.repo.get_settings().get("postdock", {}).get("selected_parents", [])
        dialog = CompoundSelectorDialog(rows, set(current), self.repo, self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        selected = dialog.selected()
        if not selected:
            QMessageBox.warning(self, "Select compounds", "Select at least one successfully docked compound to export.")
            return
        config_path = self.repo.root / "project.yml"
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        config.setdefault("postdock", {})["selected_parents"] = selected
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        self._write_log(f"Export selected compounds: {len(selected)}")
        self.run_postdock()

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
        stage_names = {"screening": "Drug-likeness Screening", "molscrub": "Molecular States", "meeko": "Ligand Preparation", "vina": "Docking", "postdock": "Export"}
        display_stage = stage_names.get(stage, stage.title())
        self.current_stage_label.setText(display_stage)
        profile_name = getattr(event, "receptor_profile_name", None)
        item_text = str(getattr(event, "item_id", "") or "")
        self.current_item_label.setText(f"{profile_name}: {item_text}" if profile_name else item_text)
        self.current_task_label.setText(f"{display_stage}\n{self.current_item_label.text()}\n{index} / {total}")
        self.dashboard_status.setText(self.current_stage_label.text())
        self.stat_cards["running"].setText("Running")
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
        self.refresh_action.setEnabled(True)
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
        if isinstance(summary, dict) and "vina" in summary:
            QMessageBox.information(
                self,
                "Docking complete",
                "Docking has finished. Export is ready.\n\nUse the Export menu, then choose 'Export docked compounds (PDBQT/SDF)'."
            )

    def _pipeline_failed(self, message: str) -> None:
        self.stage_running = False
        self._set_stage_actions(True)
        self.run_all_action.setEnabled(True)
        self.refresh_action.setEnabled(True)
        self._write_log(f"Run All stopped: {message}")
        self._refresh_checkpoint_state()
        QMessageBox.critical(self, "Run All failed", message)

    def refresh_dashboard(self) -> None:
        if not self.repo or self.stage_running:
            return
        candidate = self._project_input_csv()
        if candidate:
            try:
                result = self.repo.sync_ligands_csv(candidate)
            except Exception as exc:
                QMessageBox.critical(self, "Input reload failed", str(exc))
                return
            self._write_log(
                f"Reloaded ligands from {candidate}: added={result.added}, unchanged={result.unchanged}, "
                f"changed={result.changed}, archived={result.archived}"
            )
        self.refresh_tables()
        self._refresh_checkpoint_state()
        self._write_log("Dashboard refreshed")

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

    def manage_receptors(self) -> None:
        if not self.repo or self.stage_running:
            return
        dialog = ReceptorManagerDialog(self.repo, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            enabled = sum(1 for profile in self.repo.get_receptor_profiles() if profile.get("enabled"))
            self._write_log(f"Updated receptor profiles: {enabled} enabled")
            self._refresh_checkpoint_state()

    def show_multidock_results(self) -> None:
        if self.repo:
            DockingResultsDialog(self.repo, self).exec()

    def _on_workflow_reset(self) -> None:
        if not self.repo:
            return
        project_root = self.repo.root
        self.repo = ProjectRepository(project_root)
        self.repo.migrate()
        self.repo.recover_interrupted_runs()
        self._activate_project("Reloaded")

    def export_manifest_csv(self) -> None:
        if self.repo:
            SettingsDialog(self.repo.get_settings(), self, self.repo)._export_manifest()

    def export_leaderboard_csv(self) -> None:
        if self.repo:
            SettingsDialog(self.repo.get_settings(), self, self.repo)._export_leaderboard()
