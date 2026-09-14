from __future__ import annotations

import shutil

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout,
    QLabel, QLineEdit, QMessageBox, QPushButton, QSpinBox, QTabWidget, QVBoxLayout, QWidget,
)

from ..project import ProjectRepository
from ..csv_exports import export_leaderboard_csv, export_manifest_csv
from .compound_selector import CompoundSelectorDialog


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
        for name in ("artifacts", "logs", "exports", "For_PostDocking"):
            target = self.repository.root / name
            if target.exists():
                shutil.rmtree(target)
            target.mkdir(parents=True, exist_ok=True)
        QMessageBox.information(self, "Purge complete", "Generated workflow data was removed. The input files and project.yml were preserved.")
        self.workflow_reset.emit()

    def _export_manifest(self) -> None:
        if not self.repository:
            return
        outputs = export_manifest_csv(self.repository)
        QMessageBox.information(self, "Export complete", f"Wrote {len(outputs)} manifest CSV files")

    def _export_leaderboard(self) -> None:
        if not self.repository:
            return
        outputs = export_leaderboard_csv(self.repository)
        QMessageBox.information(self, "Export complete", f"Wrote {len(outputs)} leaderboard CSV files")

    def _screening_tab(self, values: dict[str, object]) -> QWidget:
        form = QFormLayout()
        policy = QComboBox()
        for label, value in (("Annotate only", "annotate_only"), ("Exclude failing all selected rules", "exclude_failing_all"), ("Exclude failing any selected rule", "exclude_failing_any"), ("Manual review", "manual_review")):
            policy.addItem(label, value)
        policy.setCurrentIndex(max(0, policy.findData(values.get("policy", "annotate_only"))))
        self.fields["screening.policy"] = policy
        form.addRow("Policy", policy)
        for key, label in (("lipinski", "Lipinski"), ("veber", "Veber"), ("egan", "Egan"), ("ghose", "Ghose"), ("boiled_egg", "BOILED-Egg (BBB/Yolk)")):
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
        widget = QWidget(); widget.setLayout(form); return widget

    def _selection_text(self) -> str:
        selected = self.fields.get("postdock.selected_parents", [])
        return "All successfully docked compounds" if not selected else f"{len(selected)} compounds selected"

    def _select_compounds(self) -> None:
        if not self.repository:
            return
        profile_ids = [str(profile["id"]) for profile in self.repository.get_receptor_profiles() if profile.get("enabled")]
        placeholders = ",".join("?" for _ in profile_ids)
        with self.repository.connection() as conn:
            rows = conn.execute(f"""SELECT s.parent_id, MIN(p.affinity) best_affinity FROM docking_runs d
                JOIN molecular_states s ON s.state_id=d.state_id JOIN docking_poses p ON p.run_id=d.run_id
                JOIN parent_ligands l ON l.parent_id=s.parent_id
                WHERE l.active=1 AND s.active=1 AND d.status='completed' AND d.is_current=1
                AND d.receptor_profile_id IN ({placeholders}) GROUP BY s.parent_id ORDER BY best_affinity ASC, s.parent_id""", profile_ids).fetchall() if profile_ids else []
        dialog = CompoundSelectorDialog(rows, set(self.fields.get("postdock.selected_parents", [])), self.repository, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.fields["postdock.selected_parents"] = dialog.selected()

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
        values: dict[str, dict[str, object]] = {"screening": {}, "molscrub": {}, "meeko": {}, "postdock": {}}
        values["postdock"]["selected_parents"] = list(self.settings.get("postdock", {}).get("selected_parents", []))
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
