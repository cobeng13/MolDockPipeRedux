from __future__ import annotations

import shutil
import uuid
from pathlib import Path

from PyQt6.QtWidgets import (
    QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QSpinBox, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..project import DEFAULT_VINA_PROFILE, ProjectRepository
from .redocking_dialog import RedockingSetupDialog
from .redocking_queue import RedockingQueueController, RedockingQueueDialog
from ..redocking.runner import validate_redocking_prerequisites


class ReceptorProfileDialog(QDialog):
    def __init__(self, repository: ProjectRepository, profile: dict[str, object], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.repository = repository
        self.profile = dict(profile)
        self.setWindowTitle("Receptor Profile")
        form = QFormLayout()
        self.name = QLineEdit(str(profile.get("name", "")))
        self.enabled = QCheckBox("Use this receptor during Docking and Run All")
        self.enabled.setChecked(bool(profile.get("enabled", True)))
        self.receptor = QLineEdit(str(profile.get("receptor", "")))
        browse = QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        receptor_row = QHBoxLayout(); receptor_row.addWidget(self.receptor); receptor_row.addWidget(browse)
        form.addRow("Name", self.name)
        form.addRow("", self.enabled)
        form.addRow("Receptor PDBQT", receptor_row)
        self.numbers: dict[str, QDoubleSpinBox | QSpinBox] = {}
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
            field = QDoubleSpinBox(); field.setRange(-100000, 100000); field.setDecimals(4)
            field.setValue(float(profile.get(key, 0 if key.startswith("center") else 20)))
            self.numbers[key] = field; form.addRow(key.replace("_", " ").title(), field)
        for key, label, minimum, maximum, default in (
            ("exhaustiveness", "Exhaustiveness", 1, 10000, 8), ("num_modes", "Number of modes", 1, 100, 9),
            ("energy_range", "Energy range", 0, 100, 3), ("seed", "Random seed", 0, 2147483647, 42),
            ("cpu_count", "CPU count", 1, 256, 1),
        ):
            field = QSpinBox(); field.setRange(minimum, maximum); field.setValue(int(profile.get(key, default)))
            self.numbers[key] = field; form.addRow(label, field)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._validate); buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self); layout.addLayout(form); layout.addWidget(buttons)

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select prepared receptor", filter="PDBQT files (*.pdbqt);;All files (*)")
        if path:
            self.receptor.setText(path)

    def _validate(self) -> None:
        name = self.name.text().strip()
        receptor_text = self.receptor.text().strip()
        receptor = Path(receptor_text)
        resolved = receptor if receptor.is_absolute() else self.repository.root / receptor
        if not name:
            QMessageBox.warning(self, "Missing name", "Enter a receptor profile name.")
            return
        if not receptor_text or not resolved.is_file():
            QMessageBox.warning(self, "Missing receptor", "Select an existing prepared receptor PDBQT file.")
            return
        if any(key.startswith("size_") and field.value() <= 0 for key, field in self.numbers.items()):
            QMessageBox.warning(self, "Invalid search box", "All search-box sizes must be greater than zero.")
            return
        self.accept()

    def values(self) -> tuple[dict[str, object], Path]:
        profile = dict(self.profile)
        profile["name"] = self.name.text().strip()
        profile["enabled"] = self.enabled.isChecked()
        profile["archived"] = False
        profile["receptor"] = self.receptor.text().strip()
        for key, field in self.numbers.items():
            profile[key] = field.value()
        text = self.receptor.text().strip()
        source = Path(text) if Path(text).is_absolute() else self.repository.root / text
        return profile, source


class ReceptorManagerDialog(QDialog):
    def __init__(self, repository: ProjectRepository, parent: QWidget | None = None,
                 queue_controller: RedockingQueueController | None = None) -> None:
        super().__init__(parent)
        self.repository = repository
        self.profiles = repository.get_receptor_profiles(include_archived=True)
        self.original_profiles = {str(profile["id"]): dict(profile) for profile in self.profiles}
        self.sources: dict[str, Path] = {}
        self.deleted_profile_ids: list[str] = []
        self.queue_controller = queue_controller or RedockingQueueController(repository, self)
        self.setWindowTitle("Receptor Manager")
        self.resize(760, 440)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Name", "Enabled", "Status", "Receptor"])
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        controls = QHBoxLayout()
        for label, callback in (("Add", self._add), ("Edit", self._edit), ("Duplicate", self._duplicate),
                                ("Enable / Disable", self._toggle), ("Archive", self._archive), ("Delete", self._delete)):
            button = QPushButton(label); button.clicked.connect(callback); controls.addWidget(button)
        controls.addStretch()
        validation_controls = QHBoxLayout()
        self.redock_button = QPushButton("Add to Validation Queue"); self.redock_button.clicked.connect(self._redock); validation_controls.addWidget(self.redock_button)
        queue_button = QPushButton("View Validation Queue"); queue_button.clicked.connect(self._show_queue); validation_controls.addWidget(queue_button)
        validation_controls.addStretch()
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self._save); buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        description = QLabel("Each enabled profile docks the same prepared ligand states against its own receptor and search settings.")
        description.setWordWrap(True)
        layout.addWidget(description); layout.addWidget(self.table); layout.addLayout(controls); layout.addLayout(validation_controls); layout.addWidget(buttons)
        self._refresh()
        self.table.itemSelectionChanged.connect(self._update_redock_button)

    def _refresh(self) -> None:
        self.table.setRowCount(0)
        for profile in self.profiles:
            row = self.table.rowCount(); self.table.insertRow(row)
            values = (profile.get("name", ""), "Yes" if profile.get("enabled") else "No",
                      "Archived" if profile.get("archived") else "Active", profile.get("receptor", ""))
            for column, value in enumerate(values):
                self.table.setItem(row, column, QTableWidgetItem(str(value)))
        self._update_redock_button()

    def _update_redock_button(self) -> None:
        index = self._selected()
        eligible = index is not None and not self.profiles[index].get("archived") and not validate_redocking_prerequisites(self.repository, self.profiles[index])
        self.redock_button.setEnabled(bool(eligible))

    def _redock(self) -> None:
        index = self._selected()
        if index is None: return
        profile = self.profiles[index]; missing = validate_redocking_prerequisites(self.repository, profile)
        if missing:
            QMessageBox.warning(self, "Redocking cannot start", "Missing:\n• " + "\n• ".join(missing) + "\n\nReturn to Receptor Preparation to complete the profile."); return
        setup = RedockingSetupDialog(self.repository, profile, self, action_label="Add to Queue")
        if setup.exec() != QDialog.DialogCode.Accepted: return
        self.repository.save_receptor_profiles(self.profiles)
        try:
            self.queue_controller.enqueue(profile, setup.settings())
        except ValueError as exc:
            QMessageBox.warning(self, "Validation not queued", str(exc)); return
        QMessageBox.information(self, "Validation queued",
            f"{profile.get('name')} was added to the validation queue.\n\n"
            "Validations run one at a time. You can close the queue window and leave MolDockPipe running.")
        self._show_queue()

    def _show_queue(self) -> None:
        RedockingQueueDialog(self.queue_controller, self).exec()

    def _activate_index(self, index: int) -> None:
        self.profiles[index]["enabled"] = True; self.repository.save_receptor_profiles(self.profiles); self._refresh(); self.table.selectRow(index)

    def _selected(self) -> int | None:
        rows = self.table.selectionModel().selectedRows()
        return rows[0].row() if rows else None

    def _edit_profile(self, index: int) -> bool:
        dialog = ReceptorProfileDialog(self.repository, self.profiles[index], self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return False
        profile, source = dialog.values()
        if any(other != index and not current.get("archived") and str(current.get("name", "")).casefold() == str(profile["name"]).casefold()
               for other, current in enumerate(self.profiles)):
            QMessageBox.warning(self, "Duplicate name", "Receptor profile names must be unique.")
            return False
        self.profiles[index] = profile
        self.sources[str(profile["id"])] = source
        self._refresh(); self.table.selectRow(index)
        return True

    def _add(self) -> None:
        profile = dict(DEFAULT_VINA_PROFILE)
        profile["id"] = uuid.uuid4().hex
        profile["name"] = "New receptor"
        profile["receptor"] = ""
        self.profiles.append(profile)
        index = len(self.profiles) - 1
        if not self._edit_profile(index):
            self.profiles.pop(index); self._refresh()

    def _edit(self) -> None:
        index = self._selected()
        if index is not None and not self.profiles[index].get("archived"):
            self._edit_profile(index)

    def _duplicate(self) -> None:
        index = self._selected()
        if index is None:
            return
        original = self.profiles[index]
        duplicate = dict(original)
        duplicate["id"] = uuid.uuid4().hex
        duplicate["name"] = f"{original.get('name', 'Receptor')} Copy"
        duplicate["enabled"] = False
        duplicate["archived"] = False
        source = self.repository.root / str(original.get("receptor", ""))
        duplicate["receptor"] = str(source)
        self.profiles.append(duplicate)
        new_index = len(self.profiles) - 1
        if not self._edit_profile(new_index):
            self.profiles.pop(new_index); self._refresh()

    def _toggle(self) -> None:
        index = self._selected()
        if index is not None and not self.profiles[index].get("archived"):
            self.profiles[index]["enabled"] = not bool(self.profiles[index].get("enabled"))
            self._refresh(); self.table.selectRow(index)

    def _archive(self) -> None:
        index = self._selected()
        if index is None or self.profiles[index].get("archived"):
            return
        answer = QMessageBox.question(self, "Archive receptor", "Archive this receptor profile? Historical docking results will be retained.")
        if answer == QMessageBox.StandardButton.Yes:
            self.profiles[index]["archived"] = True
            self.profiles[index]["enabled"] = False
            self._refresh()

    def _delete(self) -> None:
        index = self._selected()
        if index is None:
            return
        profile = self.profiles[index]
        with self.repository.connection() as conn:
            running = conn.execute("""SELECT EXISTS(SELECT 1 FROM redocking_queue_items
                WHERE receptor_profile_id=? AND status='running')""", (str(profile["id"]),)).fetchone()[0]
        if running:
            QMessageBox.warning(self, "Validation running",
                "Cancel this receptor's active validation and wait for it to stop before deleting the receptor.")
            return
        name = str(profile.get("name", "this receptor"))
        dialog = QMessageBox(QMessageBox.Icon.Warning, "Delete receptor profile",
            f"Delete '{name}' from this project?\n\n"
            "When you click Save, its prepared receptor and redocking files will be deleted. "
            "Ligand inputs and other receptor profiles are not affected. Historical database records are retained.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel, self)
        dialog.button(QMessageBox.StandardButton.Yes).setText("Delete")
        dialog.setDefaultButton(QMessageBox.StandardButton.Cancel)
        if dialog.exec() != QMessageBox.StandardButton.Yes:
            return
        profile_id = str(profile["id"])
        self.profiles.pop(index)
        self.sources.pop(profile_id, None)
        # Profiles added in this unsaved dialog have no persisted directory.
        if profile_id in self.original_profiles:
            self.deleted_profile_ids.append(profile_id)
        self._refresh()

    def _save(self) -> None:
        docking_keys = {"receptor", "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
                        "exhaustiveness", "num_modes", "energy_range", "seed", "cpu_count"}
        for profile in self.profiles:
            profile_id = str(profile["id"])
            source = self.sources.get(profile_id)
            original = self.original_profiles.get(profile_id)
            original_source = self.repository.root / str(original.get("receptor", "")) if original else None
            source_changed = bool(source and (original_source is None or source.resolve() != original_source.resolve()))
            if source:
                destination = self.repository.root / "inputs" / "receptors" / profile_id / "receptor.pdbqt"
                destination.parent.mkdir(parents=True, exist_ok=True)
                if source.resolve() != destination.resolve():
                    shutil.copy2(source, destination)
                profile["receptor"] = destination.relative_to(self.repository.root).as_posix()
            if original and (source_changed or any(original.get(key) != profile.get(key) for key in docking_keys)):
                with self.repository.connection() as conn:
                    conn.execute("UPDATE docking_runs SET is_current=0 WHERE receptor_profile_id=?", (profile_id,))
        for profile_id in self.deleted_profile_ids:
            self.repository.remove_receptor_profile(profile_id)
        self.repository.save_receptor_profiles(self.profiles)
        self.accept()
