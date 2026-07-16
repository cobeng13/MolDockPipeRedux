from __future__ import annotations

import json
import threading
from pathlib import Path

from PyQt6.QtCore import QObject, QThread, QTimer, QUrl, pyqtSignal, pyqtSlot
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QFormLayout, QHBoxLayout, QLabel, QMessageBox, QPushButton, QSpinBox,
    QDoubleSpinBox, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..project import ProjectRepository
from ..redocking.models import RedockingSettings
from ..redocking.runner import RedockingRunner, STAGES, validate_redocking_prerequisites

STAGE_LABELS = {
    "reference_validation": "Reference ligand validation", "ligand_preparation": "Ligand preparation",
    "vina_redocking": "Vina redocking", "pose_export": "Pose export", "mol2_generation": "MOL2 generation",
    "final_validation": "Final validation",
}


class RedockingWorker(QObject):
    progress = pyqtSignal(object); succeeded = pyqtSignal(object); failed = pyqtSignal(str); finished = pyqtSignal()

    def __init__(self, repository: ProjectRepository, profile: dict[str, object], settings: RedockingSettings,
                 run_id: str | None = None) -> None:
        super().__init__(); self.cancel_event = threading.Event(); self.runner = RedockingRunner(repository, profile, settings,
            progress=self.progress.emit, cancel_event=self.cancel_event); self.run_id = run_id

    @pyqtSlot()
    def run(self) -> None:
        try: self.succeeded.emit(self.runner.run(self.run_id))
        except Exception as exc: self.failed.emit(str(exc))
        finally: self.finished.emit()

    def cancel(self) -> None: self.cancel_event.set()


class RedockingSetupDialog(QDialog):
    def __init__(self, repository: ProjectRepository, profile: dict[str, object], parent: QWidget | None = None) -> None:
        super().__init__(parent); self.repository = repository; self.profile = profile; self.setWindowTitle("Redocking Validation")
        reference = profile.get("reference_ligand") if isinstance(profile.get("reference_ligand"), dict) else {}
        form = QFormLayout(); form.addRow("Receptor profile", QLabel(str(profile.get("name", profile.get("id", "")))))
        form.addRow("Reference ligand", QLabel(str(reference.get("identity", "Missing"))))
        center = ", ".join(f"{float(profile.get('center_'+axis, 0)):.3f}" for axis in "xyz")
        size = " × ".join(f"{float(profile.get('size_'+axis, 0)):.1f}" for axis in "xyz")
        form.addRow("Docking box", QLabel(f"Center: {center}\nSize: {size} Å"))
        self.values: dict[str, QSpinBox | QDoubleSpinBox] = {}
        defaults = {"exhaustiveness": max(32, int(profile.get("exhaustiveness", 32))), "num_modes": max(20, int(profile.get("num_modes", 20))),
                    "energy_range": max(5.0, float(profile.get("energy_range", 5))), "seed": int(profile.get("seed", 123456)),
                    "cpu_count": int(profile.get("cpu_count", 1))}
        for key, label, minimum, maximum in (("exhaustiveness", "Exhaustiveness", 1, 10000), ("num_modes", "Modes", 1, 100),
            ("energy_range", "Energy range", 0, 100), ("seed", "Seed", 0, 2147483647), ("cpu_count", "CPU count", 1, 256)):
            spin = QDoubleSpinBox() if key == "energy_range" else QSpinBox(); spin.setRange(minimum, maximum); spin.setValue(defaults[key])
            self.values[key] = spin; form.addRow(label, spin)
        mapping_path = Path(str(reference.get("mapping", ""))); mapping_path = mapping_path if mapping_path.is_absolute() else repository.root / mapping_path
        warning = QLabel("")
        if mapping_path.is_file():
            mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
            if mapping.get("chemistry_source") == "rdkit_pdb_inference": warning.setText("Warning: connectivity was inferred from crystallographic coordinates. Review the chemistry mapping before use.")
            if int(mapping.get("heavy_atom_count", 0)) > 80: warning.setText((warning.text() + "\n" if warning.text() else "") + "Warning: the reference ligand is very large for routine redocking.")
        warning.setWordWrap(True); form.addRow(warning)
        missing = validate_redocking_prerequisites(repository, profile)
        if missing: form.addRow(QLabel("Redocking cannot start. Missing:\n• " + "\n• ".join(missing)))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel); self.run_button = buttons.addButton("Run Redocking", QDialogButtonBox.ButtonRole.AcceptRole)
        self.run_button.setEnabled(not missing); buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self); layout.addLayout(form); layout.addWidget(buttons)

    def settings(self) -> RedockingSettings:
        return RedockingSettings(int(self.values["exhaustiveness"].value()), int(self.values["num_modes"].value()),
            float(self.values["energy_range"].value()), int(self.values["seed"].value()), int(self.values["cpu_count"].value()))


class RedockingProgressDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent); self.setWindowTitle("Redocking Validation"); self.setModal(True); self.resize(480, 360)
        self.labels = {}; layout = QVBoxLayout(self); self.active = QLabel("Starting…"); layout.addWidget(self.active)
        for stage in STAGES:
            label = QLabel(f"○ {STAGE_LABELS[stage]}"); self.labels[stage] = label; layout.addWidget(label)
        self.elapsed = QLabel("Elapsed: 0 s"); layout.addWidget(self.elapsed); layout.addStretch()
        self.cancel_button = QPushButton("Cancel"); layout.addWidget(self.cancel_button)
        self.seconds = 0; self.timer = QTimer(self); self.timer.timeout.connect(self._tick); self.timer.start(1000)

    def _tick(self) -> None: self.seconds += 1; self.elapsed.setText(f"Elapsed: {self.seconds} s")

    def update_event(self, event: dict[str, object]) -> None:
        stage = str(event.get("stage", "")); label = self.labels.get(stage)
        if not label: return
        if event.get("event") == "stage_started": self.active.setText(f"Active: {STAGE_LABELS[stage]}"); label.setText(f"▶ {STAGE_LABELS[stage]}")
        elif event.get("event") in {"stage_completed", "stage_reused"}: label.setText(f"✓ {STAGE_LABELS[stage]}"); label.setStyleSheet("color: #238636")


class RedockingResultDialog(QDialog):
    activate_requested = pyqtSignal()
    run_again_requested = pyqtSignal()

    def __init__(self, repository: ProjectRepository, result: dict[str, object], profile: dict[str, object], seed: int, parent: QWidget | None = None) -> None:
        super().__init__(parent); self.repository = repository; self.result = result; self.setWindowTitle("Redocking Artifacts Ready")
        reference = profile.get("reference_ligand", {})
        text = QLabel(f"Redocking Artifacts Ready\n\nReceptor: {profile.get('name')}\nReference: {reference.get('identity')}\n"
            f"Vina poses: {result.get('poses')}\nTop affinity: {result.get('top_affinity')} kcal/mol\nSeed: {seed}\n\n"
            "DockRMSD inputs\n✓ Reference ligand MOL2\n✓ Redocked pose MOL2 files\n✓ Top-ranked pose MOL2\n✓ Command examples\n\nRMSD status\nNot calculated")
        open_button = QPushButton("Open Output Folder"); open_button.clicked.connect(lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(result["run_root"]))))
        poses = QPushButton("View Pose Table"); poses.clicked.connect(self._show_poses)
        again = QPushButton("Run Again"); again.clicked.connect(lambda: (self.run_again_requested.emit(), self.accept()))
        activate = QPushButton("Activate Receptor"); activate.clicked.connect(self.activate_requested.emit)
        close = QPushButton("Close"); close.clicked.connect(self.accept)
        row = QHBoxLayout(); row.addWidget(open_button); row.addWidget(poses); row.addWidget(again); row.addWidget(activate); row.addWidget(close)
        layout = QVBoxLayout(self); layout.addWidget(text); layout.addLayout(row)

    def _show_poses(self) -> None:
        dialog = QDialog(self); dialog.setWindowTitle("Redocked Poses"); table = QTableWidget(0, 4)
        table.setHorizontalHeaderLabels(("Rank", "Affinity", "SDF", "MOL2"))
        with self.repository.connection() as conn:
            rows = conn.execute("SELECT pose_rank,affinity,sdf_path,mol2_path FROM redocking_poses WHERE run_id=? ORDER BY pose_rank", (str(self.result["run_id"]),)).fetchall()
        for row_data in rows:
            row = table.rowCount(); table.insertRow(row)
            for column, value in enumerate(row_data): table.setItem(row, column, QTableWidgetItem(str(value)))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); buttons.rejected.connect(dialog.reject)
        layout = QVBoxLayout(dialog); layout.addWidget(table); layout.addWidget(buttons); dialog.resize(760, 400); dialog.exec()
