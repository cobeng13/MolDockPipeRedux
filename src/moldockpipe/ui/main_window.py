from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import QObject, QThread, Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QDoubleSpinBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton,
    QSpinBox, QSplitter, QStatusBar, QTableWidget, QTableWidgetItem,
    QTabWidget, QTextEdit, QToolBar, QVBoxLayout, QWidget,
)

import yaml

from ..project import ProjectRepository
from ..pipeline import PipelineRunner


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.repo: ProjectRepository | None = None
        self.setWindowTitle("MolDockPipe Redux")
        self.resize(1280, 760)
        self._build_ui()

    def _build_ui(self) -> None:
        toolbar = QToolBar("Project")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)
        new_action = toolbar.addAction("New Project")
        new_action.triggered.connect(self.create_project)
        open_action = toolbar.addAction("Open Project")
        open_action.triggered.connect(self.open_project)
        import_action = toolbar.addAction("Import CSV")
        import_action.triggered.connect(self.import_csv)
        self.stage_action = toolbar.addAction("Run Screening")
        self.stage_action.setEnabled(False)
        self.stage_action.triggered.connect(self.run_screening)
        self.molscrub_action = toolbar.addAction("Generate States")
        self.molscrub_action.setEnabled(False)
        self.molscrub_action.triggered.connect(self.run_molscrub)
        self.meeko_action = toolbar.addAction("Prepare Ligands")
        self.meeko_action.setEnabled(False)
        self.meeko_action.triggered.connect(self.run_meeko)
        self.vina_action = toolbar.addAction("Run Vina")
        self.vina_action.setEnabled(False)
        self.vina_action.triggered.connect(self.run_vina)
        settings_action = toolbar.addAction("Docking Settings")
        settings_action.setEnabled(False)
        settings_action.triggered.connect(self.edit_docking_settings)
        self.settings_action = settings_action

        self.project_path = QLineEdit()
        self.project_path.setReadOnly(True)
        self.project_path.setPlaceholderText("Create or open a portable project")
        top = QWidget()
        top_layout = QFormLayout(top)
        top_layout.addRow("Project", self.project_path)

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
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setPlaceholderText("Structured pipeline events and execution logs appear here.")
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(tabs)
        splitter.addWidget(self.log)
        splitter.setSizes([530, 180])
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
        self.settings_action.setEnabled(True)
        self._auto_import_input()
        self._write_log(f"{action} project: {self.repo.root}")
        self.refresh_tables()

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
        if not self.repo:
            return
        try:
            completed, failed = PipelineRunner(self.repo).run_screening()
        except Exception as exc:
            QMessageBox.critical(self, "Screening failed", str(exc))
            return
        self._write_log(f"Screening completed: {completed}; failed: {failed}")
        self.refresh_tables()

    def run_molscrub(self) -> None:
        if not self.repo:
            return
        try:
            created, failed = PipelineRunner(self.repo).run_molscrub()
        except Exception as exc:
            QMessageBox.critical(self, "State generation failed", str(exc))
            return
        self._write_log(f"Generated molecular states: {created}; failed parents: {failed}")
        self.refresh_tables()

    def run_meeko(self) -> None:
        if not self.repo:
            return
        try:
            prepared, failed = PipelineRunner(self.repo).run_meeko()
        except Exception as exc:
            QMessageBox.critical(self, "Meeko preparation failed", str(exc))
            return
        self._write_log(f"Meeko preparation completed: {prepared}; failed: {failed}")
        self.refresh_tables()

    def run_vina(self) -> None:
        if not self.repo:
            return
        self.vina_action.setEnabled(False)
        self._write_log("Vina docking started in background...")
        self.vina_thread = QThread(self)
        self.vina_worker = VinaWorker(self.repo)
        self.vina_worker.moveToThread(self.vina_thread)
        self.vina_thread.started.connect(self.vina_worker.run)
        self.vina_worker.finished.connect(self._vina_finished)
        self.vina_worker.failed.connect(self._vina_failed)
        self.vina_worker.finished.connect(self.vina_thread.quit)
        self.vina_worker.failed.connect(self.vina_thread.quit)
        self.vina_thread.finished.connect(self.vina_worker.deleteLater)
        self.vina_thread.finished.connect(self.vina_thread.deleteLater)
        self.vina_thread.start()

    def _vina_finished(self, completed: int, failed: int) -> None:
        self.vina_action.setEnabled(True)
        self._write_log(f"Vina docking completed: {completed}; failed: {failed}")
        self.refresh_tables()

    def _vina_failed(self, message: str) -> None:
        self.vina_action.setEnabled(True)
        self._write_log(f"Vina docking failed: {message}")
        QMessageBox.critical(self, "Vina docking failed", message)

    def edit_docking_settings(self) -> None:
        if not self.repo:
            return
        dialog = DockingSettingsDialog(self.repo.get_settings().get("vina", {}), self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        config_path = self.repo.root / "project.yml"
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
            config["vina"] = dialog.values()
            config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            self._write_log(f"Updated docking settings in {config_path}")
        except Exception as exc:
            QMessageBox.critical(self, "Settings update failed", str(exc))


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

class DockingSettingsDialog(QDialog):
    """Temporary editor for the Vina section of project.yml."""

    def __init__(self, settings: dict[str, object], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Docking Settings")
        self.resize(440, 420)
        form = QFormLayout(self)
        self.fields: dict[str, QLineEdit | QSpinBox | QDoubleSpinBox] = {}
        self._line(form, "Receptor", "receptor", settings.get("receptor", "inputs/receptor_prepared.pdbqt"))
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
            self._double(form, key.replace("_", " ").title(), key, settings.get(key, 0.0 if key.startswith("center") else 20.0))
        self._integer(form, "Exhaustiveness", "exhaustiveness", settings.get("exhaustiveness", 8), 1, 10_000)
        self._integer(form, "Number of modes", "num_modes", settings.get("num_modes", 9), 1, 100)
        self._double(form, "Energy range", "energy_range", settings.get("energy_range", 3), 0, 100)
        self._integer(form, "Random seed", "seed", settings.get("seed", 42), 0, 2_147_483_647)
        self._integer(form, "CPU count", "cpu_count", settings.get("cpu_count", 1), 1, 256)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)

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

    def values(self) -> dict[str, object]:
        values: dict[str, object] = {}
        for key, field in self.fields.items():
            if isinstance(field, QLineEdit):
                values[key] = field.text().strip()
            elif isinstance(field, QDoubleSpinBox):
                values[key] = field.value()
            else:
                values[key] = field.value()
        return values
