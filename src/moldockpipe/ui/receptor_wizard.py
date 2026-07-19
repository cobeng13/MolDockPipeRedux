from __future__ import annotations

import re
import uuid
from pathlib import Path

from PyQt6.QtCore import QObject, Qt, pyqtSignal, pyqtSlot
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
    QDialog, QFrame, QMessageBox, QPushButton, QRadioButton, QStackedWidget, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

from ..receptors.analysis import analyze_structure, atoms_for_residues
from ..receptors.box_calculation import center_from_atoms, envelope_box, radius_of_gyration_box
from ..receptors.models import ComponentRole, ProteinResidueIssue, ReceptorPreparationPlan, ResidueKey, StructureAnalysis
from ..receptors.preparation import prepare_receptor


ROLE_LABELS = {
    ComponentRole.REFERENCE_LIGAND: "Reference ligand",
    ComponentRole.RETAINED_COFACTOR: "Retained cofactor",
    ComponentRole.RETAINED_ION: "Retained metal / ion",
    ComponentRole.RECEPTOR_COMPONENT: "Receptor component",
    ComponentRole.REMOVE: "Remove",
    ComponentRole.UNRESOLVED: "Choose an action…",
}


def _page_surface(page: QWidget) -> QVBoxLayout:
    """Give every wizard page an opaque, full-height surface.

    Qt's Windows wizard style paints unused page space with its own light
    palette, even when the page itself is styled.  An expanding frame avoids
    that split light/dark page and gives every step a consistent content area.
    """
    outer = QVBoxLayout(page)
    outer.setContentsMargins(0, 0, 0, 0)
    surface = QFrame()
    surface.setObjectName("wizardPageSurface")
    outer.addWidget(surface, 1)
    content = QVBoxLayout(surface)
    content.setContentsMargins(22, 18, 22, 18)
    content.setSpacing(12)
    return content


def _page_heading(layout: QVBoxLayout, title: str, detail: str | None = None) -> None:
    heading = QLabel(title)
    heading.setObjectName("wizardPageHeading")
    layout.addWidget(heading)
    if detail:
        subtitle = QLabel(detail)
        subtitle.setObjectName("wizardPageSubtitle")
        subtitle.setWordWrap(True)
        layout.addWidget(subtitle)


class SourcePage(QWidget):
    def __init__(self, wizard: "ReceptorPreparationWizard") -> None:
        super().__init__(); self.owner = wizard; self.setObjectName("wizardStep")
        self.name = QLineEdit(); self.source = QLineEdit(); self.model = QComboBox(); self.model.addItem("Model 1", 0)
        self.chains = QLineEdit(); self.chains.setPlaceholderText("All chains, or comma-separated: A, B")
        browse = QPushButton("Browse…"); browse.clicked.connect(self._browse)
        row = QHBoxLayout(); row.addWidget(self.source); row.addWidget(browse)
        content = _page_surface(self)
        _page_heading(content, "Select structure", "Choose the source structure and the receptor chains to prepare.")
        form = QFormLayout(); form.setContentsMargins(0, 0, 0, 0); form.setHorizontalSpacing(14); form.setVerticalSpacing(10)
        form.addRow("Receptor profile name", self.name); form.addRow("Input structure", row)
        form.addRow("Structural model", self.model); form.addRow("Chains to include", self.chains)
        note = QLabel("The source is copied into the portable project. The original file is never modified.")
        note.setWordWrap(True); note.setProperty("secondary", True)
        content.addLayout(form); content.addWidget(note); content.addStretch()

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select receptor structure", filter="Structures (*.pdb *.cif *.mmcif)")
        if path:
            self.source.setText(path)
            if not self.name.text(): self.name.setText(Path(path).stem)
            try:
                overview = analyze_structure(Path(path))
                self.model.clear()
                for index, name in enumerate(overview.models): self.model.addItem(f"Model {name}", index)
                self.chains.setPlaceholderText(", ".join(overview.chains_by_model.get(0, ())))
            except Exception as exc:
                QMessageBox.warning(self, "Structure analysis failed", str(exc))

    def validatePage(self) -> bool:
        name, source = self.name.text().strip(), Path(self.source.text().strip())
        if not name or not source.is_file():
            QMessageBox.warning(self, "Missing structure", "Enter a profile name and select an existing structure file."); return False
        model = int(self.model.currentData() or 0)
        chains = tuple(value.strip() for value in self.chains.text().split(",") if value.strip())
        try:
            analysis = analyze_structure(source, model, chains)
        except Exception as exc:
            QMessageBox.critical(self, "Structure analysis failed", str(exc)); return False
        available = set(analysis.chains_by_model[model])
        if chains and not set(chains) <= available:
            QMessageBox.warning(self, "Unknown chain", f"Available chains: {', '.join(sorted(available))}"); return False
        self.owner.analysis = analysis
        self.owner.included_chains = chains or tuple(analysis.chains_by_model[model])
        self.owner.components_page.load_analysis(analysis)
        self.owner.integrity_page.load_analysis(analysis)
        return True


class ComponentsPage(QWidget):
    def __init__(self, wizard: "ReceptorPreparationWizard") -> None:
        super().__init__(); self.owner = wizard; self.setObjectName("wizardStep")
        self.summary = QLabel(); self.summary.setWordWrap(True)
        self.table = QTableWidget(0, 6); self.table.setHorizontalHeaderLabels(["Component", "Chain", "Residue", "Heavy atoms", "Suggested role", "Action"])
        self._inference_confirmed = False
        self.template = QLineEdit(); self.template.setPlaceholderText("Optional chemically complete SDF for the selected reference ligand")
        browse = QPushButton("Browse SDF…"); browse.clicked.connect(self._browse_template)
        template_row = QHBoxLayout(); template_row.addWidget(self.template); template_row.addWidget(browse)
        self.selectors: list[QComboBox] = []
        layout = _page_surface(self); _page_heading(layout, "Review receptor components", "Confirm what remains in the receptor and choose an optional ligand chemistry template.")
        layout.addWidget(self.summary); layout.addWidget(self.table, 1); layout.addLayout(template_row)

    def _browse_template(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Select reference-ligand chemical template", filter="SDF files (*.sdf)")
        if path: self.template.setText(path)

    def load_analysis(self, analysis: StructureAnalysis) -> None:
        self.summary.setText(f"{analysis.protein_residue_count} protein residues; {analysis.water_count} waters; "
                             f"{len(analysis.alternate_locations)} alternate-location atoms; {len(analysis.zero_occupancy_atoms)} zero-occupancy atoms")
        self.table.setRowCount(0); self.selectors.clear(); reference_assigned = False
        for component in analysis.components:
            row = self.table.rowCount(); self.table.insertRow(row)
            values = (component.key.name, component.key.chain, f"{component.key.number}{component.key.insertion_code}",
                      len(component.heavy_atoms), f"{component.category}: {component.reason}")
            for column, value in enumerate(values): self.table.setItem(row, column, QTableWidgetItem(str(value)))
            selector = QComboBox()
            for role, label in ROLE_LABELS.items(): selector.addItem(label, role.value)
            suggested = component.suggested_role
            if suggested == ComponentRole.REFERENCE_LIGAND:
                if reference_assigned: suggested = ComponentRole.UNRESOLVED
                else: reference_assigned = True
            selector.setCurrentIndex(selector.findData(suggested.value)); self.table.setCellWidget(row, 5, selector); self.selectors.append(selector)

    def validatePage(self) -> bool:
        roles = [ComponentRole(selector.currentData()) for selector in self.selectors]
        if ComponentRole.UNRESOLVED in roles:
            QMessageBox.warning(self, "Unresolved components", "Choose an action for every unclassified component."); return False
        if roles.count(ComponentRole.REFERENCE_LIGAND) > 1:
            QMessageBox.warning(self, "Multiple reference ligands", "Choose at most one reference ligand."); return False
        if self.template.text().strip() and not Path(self.template.text().strip()).is_file():
            QMessageBox.warning(self, "Missing chemical template", "The selected reference-ligand SDF does not exist."); return False
        if ComponentRole.REFERENCE_LIGAND in roles and not self.template.text().strip() and not self._inference_confirmed:
            answer = QMessageBox.warning(self, "Use RDKit connectivity inference?",
                "No chemically complete SDF was supplied. RDKit will infer bonds from crystallographic coordinates. "
                "Preparation will stop if a complete mapping cannot be established. Continue with this explicit fallback?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes: return False
            self._inference_confirmed = True
        return True


class IntegrityPage(QWidget):
    """Make incomplete polymer residues an informed user decision, not a Meeko surprise."""
    def __init__(self, wizard: "ReceptorPreparationWizard") -> None:
        super().__init__(); self.owner = wizard; self.setObjectName("wizardStep")
        self.summary = QLabel(); self.summary.setWordWrap(True)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Residue", "Missing atoms", "Alt. locations", "Recommended", "Action"])
        self.issues: tuple[ProteinResidueIssue, ...] = ()
        self.actions: list[QComboBox] = []
        self.altlocs: list[QComboBox] = []
        layout = _page_surface(self)
        _page_heading(layout, "Review receptor integrity", "Incomplete residues can make Meeko infer invalid bonds. Choose whether to exclude each one from this prepared receptor.")
        layout.addWidget(self.summary); layout.addWidget(self.table, 1)

    def load_analysis(self, analysis: StructureAnalysis) -> None:
        self.issues = analysis.protein_residue_issues
        self.actions.clear(); self.altlocs.clear(); self.table.setRowCount(0)
        if not self.issues:
            self.summary.setText("No incomplete standard amino-acid residues were detected in the selected chains.")
            return
        self.summary.setText(f"{len(self.issues)} residue(s) are incomplete after alternate-location selection. "
                             "Excluding a residue removes only that residue from the prepared receptor; it never changes the source file.")
        for issue in self.issues:
            row = self.table.rowCount(); self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(issue.key.label()))
            self.table.setItem(row, 1, QTableWidgetItem(", ".join(issue.missing_atoms)))
            self.table.setItem(row, 2, QTableWidgetItem(", ".join(issue.alternate_locations) or "—"))
            altloc = QComboBox(); altloc.addItem("No alternate location", "")
            for label in issue.alternate_locations:
                altloc.addItem(label, label)
            altloc.setCurrentIndex(max(0, altloc.findData(issue.recommended_altloc)))
            action = QComboBox(); action.addItem("Choose an action…", "unresolved")
            action.addItem("Exclude from prepared receptor", "exclude")
            action.addItem("Keep and let Meeko attempt it", "keep")
            self.table.setCellWidget(row, 3, altloc); self.table.setCellWidget(row, 4, action)
            self.altlocs.append(altloc); self.actions.append(action)
        self.table.resizeColumnsToContents()

    def validatePage(self) -> bool:
        if any(combo.currentData() == "unresolved" for combo in self.actions):
            QMessageBox.warning(self, "Resolve receptor integrity", "Choose an action for every incomplete receptor residue before continuing.")
            return False
        kept = [issue.key.label() for issue, action in zip(self.issues, self.actions) if action.currentData() == "keep"]
        if kept:
            answer = QMessageBox.warning(self, "Attempt preparation with incomplete residues?",
                "Meeko may still fail for: " + ", ".join(kept) + ". Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return False
        return True

    def selected_altlocs(self) -> dict[ResidueKey, str]:
        return {issue.key: str(choice.currentData()) for issue, choice in zip(self.issues, self.altlocs) if choice.currentData()}

    def excluded_residues(self) -> tuple[ResidueKey, ...]:
        return tuple(issue.key for issue, action in zip(self.issues, self.actions) if action.currentData() == "exclude")


class BindingSitePage(QWidget):
    def __init__(self, wizard: "ReceptorPreparationWizard") -> None:
        super().__init__(); self.owner = wizard; self.setObjectName("wizardStep")
        self.ligand = QRadioButton("Reference-ligand centroid"); self.ligand.setChecked(True)
        self.residues = QRadioButton("Selected residues"); self.residue_text = QLineEdit(); self.residue_text.setPlaceholderText("A:120, A:122, A:198")
        self.manual = QRadioButton("Manual coordinates")
        self.center = [QDoubleSpinBox() for _ in range(3)]
        for spin in self.center: spin.setRange(-100000, 100000); spin.setDecimals(4)
        content = _page_surface(self)
        _page_heading(content, "Choose the binding-site center")
        form = QFormLayout(); form.setContentsMargins(0, 0, 0, 0); form.setVerticalSpacing(10)
        form.addRow(self.ligand); form.addRow(self.residues, self.residue_text); form.addRow(self.manual)
        for axis, spin in zip("XYZ", self.center): form.addRow(f"Center {axis}", spin)
        content.addLayout(form); content.addStretch()

    def validatePage(self) -> bool:
        try:
            if self.ligand.isChecked():
                component = self.owner.reference_component()
                if component is None: raise ValueError("Select a reference ligand or use another center method")
                center = center_from_atoms(component.atoms)
                self.owner.center_method = "reference_ligand_centroid"
                self.owner.center_parameters = {"reference_ligand": component.key.label(),
                                                "heavy_atom_count": len(component.heavy_atoms)}
            elif self.residues.isChecked():
                keys = self._parse_residues(self.residue_text.text())
                center = center_from_atoms(atoms_for_residues(Path(self.owner.source_page.source.text()), int(self.owner.source_page.model.currentData()), keys))
                self.owner.center_method = "selected_residue_centroid"
                self.owner.center_parameters = {"residues": [key.label() for key in keys]}
            else:
                center = tuple(spin.value() for spin in self.center)
                self.owner.center_method = "manual"
                self.owner.center_parameters = {}
        except Exception as exc:
            QMessageBox.warning(self, "Invalid binding site", str(exc)); return False
        self.owner.box_center = center
        for spin, value in zip(self.center, center): spin.setValue(value)
        return True

    @staticmethod
    def _parse_residues(text: str) -> tuple[ResidueKey, ...]:
        result = []
        for token in text.split(","):
            match = re.fullmatch(r"\s*([^:\s]+):(-?\d+)([A-Za-z]?)\s*", token)
            if not match: raise ValueError("Use chain:number notation, for example A:120, A:122")
            result.append(ResidueKey(match.group(1), int(match.group(2)), match.group(3)))
        if not result: raise ValueError("Enter at least one binding-site residue")
        return tuple(result)


class BoxPage(QWidget):
    def __init__(self, wizard: "ReceptorPreparationWizard") -> None:
        super().__init__(); self.owner = wizard; self.setObjectName("wizardStep")
        self.method = QComboBox(); self.method.addItem("Ligand envelope plus padding (recommended)", "ligand_envelope_padding")
        self.method.addItem("Literature-based cubic box (2.857 × reference-ligand Rg)", "radius_of_gyration")
        self.method.addItem("Manual dimensions", "manual")
        self.padding = QDoubleSpinBox(); self.padding.setRange(0, 100); self.padding.setValue(8); self.padding.setSuffix(" Å")
        self.size = [QDoubleSpinBox() for _ in range(3)]
        for spin in self.size: spin.setRange(0.1, 100000); spin.setValue(22); spin.setDecimals(3); spin.setSuffix(" Å")
        self.method.currentIndexChanged.connect(self.refresh_dimensions)
        self.padding.valueChanged.connect(self.refresh_dimensions)
        content = _page_surface(self)
        _page_heading(content, "Choose docking-box dimensions")
        form = QFormLayout(); form.setContentsMargins(0, 0, 0, 0); form.setVerticalSpacing(10)
        form.addRow("Method", self.method); form.addRow("Padding", self.padding)
        for axis, spin in zip("XYZ", self.size): form.addRow(f"Size {axis}", spin)
        content.addLayout(form); content.addStretch()

    def refresh_dimensions(self, *_args, show_error: bool = False) -> bool:
        """Immediately reflect the selected calculation method in the size fields."""
        method = str(self.method.currentData())
        manual = method == "manual"
        self.padding.setEnabled(method == "ligand_envelope_padding")
        for spin in self.size:
            spin.setReadOnly(not manual)
        if manual:
            return True
        component = self.owner.reference_component()
        if component is None:
            if show_error:
                QMessageBox.warning(self, "Reference ligand required",
                                    "This box-sizing method requires a selected reference ligand.")
            return False
        try:
            if method == "ligand_envelope_padding":
                _, size = envelope_box(component.atoms, self.padding.value())
            else:
                _, size = radius_of_gyration_box(component.atoms)
        except Exception as exc:
            if show_error:
                QMessageBox.warning(self, "Invalid docking box", str(exc))
            return False
        for spin, value in zip(self.size, size):
            spin.setValue(value)
        self.owner.box_size = size
        return True

    def validatePage(self) -> bool:
        method = str(self.method.currentData())
        try:
            if method == "ligand_envelope_padding":
                if not self.refresh_dimensions(show_error=True): return False
                size = tuple(spin.value() for spin in self.size)
            elif method == "radius_of_gyration":
                if not self.refresh_dimensions(show_error=True): return False
                size = tuple(spin.value() for spin in self.size)
            else: size = tuple(spin.value() for spin in self.size)
        except Exception as exc:
            QMessageBox.warning(self, "Invalid docking box", str(exc)); return False
        self.owner.box_size = size; self.owner.box_method = method
        self.owner.box_parameters = {"padding_angstrom": self.padding.value()} if method == "ligand_envelope_padding" else {}
        for spin, value in zip(self.size, size): spin.setValue(value)
        return True


class FinishPage(QWidget):
    def __init__(self) -> None:
        super().__init__(); self.setObjectName("wizardStep")
        self.activate = QCheckBox("Enable this receptor for docking"); self.activate.setChecked(True)
        self.preserve = QCheckBox("Preserve hydrogens already present")
        text = QLabel("Meeko will run in a temporary directory. Files are published only after preparation and validation succeed. Docking will not start automatically.")
        text.setWordWrap(True); layout = _page_surface(self); _page_heading(layout, "Prepare, validate, and save")
        layout.addWidget(text); layout.addWidget(self.preserve); layout.addWidget(self.activate); layout.addStretch()


class ReceptorPreparationWizard(QDialog):
    def __init__(self, parent=None) -> None:
        super().__init__(parent); self.setWindowTitle("Receptor Preparation"); self.resize(920, 540)
        self.setMinimumSize(820, 480)
        self.setObjectName("receptorPreparationDialog")
        self.setStyleSheet("""
            QDialog#receptorPreparationDialog, QWidget#wizardStep,
            QFrame#wizardPageSurface, QStackedWidget { background: #242424; color: #f3f4f6; }
            QFrame#dialogHeader, QFrame#dialogFooter { background: #2b2b2b; }
            QFrame#dialogHeader { border-bottom: 1px solid #3f3f46; }
            QFrame#dialogFooter { border-top: 1px solid #3f3f46; }
            QLabel#wizardPageHeading { color: #dbeafe; font-size: 17px; font-weight: 600; }
            QLabel#wizardPageSubtitle, QLabel[secondary="true"] { color: #a1a1aa; }
            QLabel { color: #f3f4f6; }
            QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
                background: #303030; color: #f3f4f6; border: 1px solid #666;
                border-radius: 4px; min-height: 30px; padding: 2px 7px;
            }
            QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
                border: 1px solid #60a5fa;
            }
            QLineEdit[readOnly="true"] { color: #d4d4d8; }
            QComboBox QAbstractItemView {
                background: #303030; color: #f3f4f6; selection-background-color: #2563eb;
            }
            QPushButton {
                background: #3b3b3b; color: #f3f4f6; border: 1px solid #71717a;
                border-radius: 4px; padding: 6px 14px; min-height: 28px; min-width: 78px;
            }
            QPushButton:hover:enabled { background: #4b4b4b; border-color: #93c5fd; }
            QPushButton#primaryButton { background: #2563eb; border-color: #60a5fa; }
            QPushButton:disabled { background: #303030; color: #71717a; border-color: #4b4b4b; }
            QRadioButton, QCheckBox { color: #f3f4f6; spacing: 7px; }
            QTableWidget {
                background: #2b2b2b; alternate-background-color: #323232; color: #f3f4f6;
                gridline-color: #555; border: 1px solid #666;
            }
            QHeaderView::section {
                background: #3a3a3a; color: #f3f4f6; padding: 6px; border: 0;
            }
        """)
        palette = self.palette()
        palette.setColor(QPalette.ColorRole.Window, QColor("#242424"))
        palette.setColor(QPalette.ColorRole.Base, QColor("#242424"))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor("#303030"))
        palette.setColor(QPalette.ColorRole.Text, QColor("#f3f4f6"))
        palette.setColor(QPalette.ColorRole.WindowText, QColor("#f3f4f6"))
        palette.setColor(QPalette.ColorRole.Button, QColor("#242424"))
        palette.setColor(QPalette.ColorRole.ButtonText, QColor("#f3f4f6"))
        palette.setColor(QPalette.ColorRole.Light, QColor("#454545"))
        palette.setColor(QPalette.ColorRole.Midlight, QColor("#393939"))
        palette.setColor(QPalette.ColorRole.Mid, QColor("#303030"))
        palette.setColor(QPalette.ColorRole.Dark, QColor("#171717"))
        self.setPalette(palette)
        self.analysis: StructureAnalysis | None = None; self.included_chains: tuple[str, ...] = ()
        self.box_center = (0.0, 0.0, 0.0); self.center_method = "manual"; self.center_parameters = {}
        self.box_size = (22.0, 22.0, 22.0); self.box_method = "manual"; self.box_parameters = {}
        self.source_page = SourcePage(self); self.components_page = ComponentsPage(self); self.integrity_page = IntegrityPage(self); self.binding_page = BindingSitePage(self)
        self.box_page = BoxPage(self); self.finish_page = FinishPage()
        self.pages = (self.source_page, self.components_page, self.integrity_page, self.binding_page, self.box_page, self.finish_page)
        self.stack = QStackedWidget()
        for page in self.pages:
            page.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
            page.setAutoFillBackground(True)
            page.setPalette(palette)
            self.stack.addWidget(page)

        header = QFrame(); header.setObjectName("dialogHeader")
        header_layout = QHBoxLayout(header); header_layout.setContentsMargins(18, 12, 18, 12)
        title = QLabel("Receptor Preparation"); title.setStyleSheet("font-size: 15px; font-weight: 600;")
        self.step_label = QLabel(); self.step_label.setProperty("secondary", True)
        header_layout.addWidget(title); header_layout.addStretch(); header_layout.addWidget(self.step_label)

        footer = QFrame(); footer.setObjectName("dialogFooter")
        footer_layout = QHBoxLayout(footer); footer_layout.setContentsMargins(18, 10, 18, 10)
        self.back_button = QPushButton("Back"); self.back_button.clicked.connect(self._back)
        self.next_button = QPushButton("Next"); self.next_button.setObjectName("primaryButton")
        self.next_button.setDefault(True); self.next_button.clicked.connect(self._advance)
        cancel_button = QPushButton("Cancel"); cancel_button.clicked.connect(self.reject)
        footer_layout.addStretch(); footer_layout.addWidget(self.back_button); footer_layout.addWidget(self.next_button); footer_layout.addWidget(cancel_button)

        layout = QVBoxLayout(self); layout.setContentsMargins(0, 0, 0, 0); layout.setSpacing(0)
        layout.addWidget(header); layout.addWidget(self.stack, 1); layout.addWidget(footer)
        self._update_navigation()

    def _update_navigation(self) -> None:
        index = self.stack.currentIndex(); total = len(self.pages)
        self.step_label.setText(f"Step {index + 1} of {total}")
        self.back_button.setEnabled(index > 0)
        self.next_button.setText("Create receptor" if index == total - 1 else "Next")

    def _back(self) -> None:
        index = self.stack.currentIndex()
        if index > 0:
            self.stack.setCurrentIndex(index - 1); self._update_navigation()

    def _advance(self) -> None:
        index = self.stack.currentIndex(); page = self.pages[index]
        validator = getattr(page, "validatePage", None)
        if callable(validator) and not validator():
            return
        if index == len(self.pages) - 1:
            self.accept(); return
        self.stack.setCurrentIndex(index + 1)
        if self.stack.currentWidget() is self.box_page:
            self.box_page.refresh_dimensions()
        self._update_navigation()

    def reference_component(self):
        if not self.analysis: return None
        for component, selector in zip(self.analysis.components, self.components_page.selectors):
            if selector.currentData() == ComponentRole.REFERENCE_LIGAND.value: return component
        return None

    def plan(self) -> ReceptorPreparationPlan:
        assert self.analysis
        roles = [ComponentRole(selector.currentData()) for selector in self.components_page.selectors]
        reference = next((component.key for component, role in zip(self.analysis.components, roles) if role == ComponentRole.REFERENCE_LIGAND), None)
        removed = tuple(component.key for component, role in zip(self.analysis.components, roles) if role == ComponentRole.REMOVE)
        retained = tuple(component.key for component, role in zip(self.analysis.components, roles) if role in {ComponentRole.RETAINED_COFACTOR, ComponentRole.RETAINED_ION, ComponentRole.RECEPTOR_COMPONENT})
        slug = re.sub(r"[^a-z0-9]+", "_", self.source_page.name.text().strip().lower()).strip("_") or "receptor"
        return ReceptorPreparationPlan(f"{slug}_{uuid.uuid4().hex[:8]}", self.source_page.name.text().strip(),
            Path(self.source_page.source.text()).resolve(), int(self.source_page.model.currentData()), self.included_chains,
            reference, removed, retained, self.box_center, self.box_size, self.box_method, self.box_parameters,
            center_method=self.center_method, center_parameters=self.center_parameters,
            altloc_choices=self.integrity_page.selected_altlocs(),
            excluded_receptor_residues=self.integrity_page.excluded_residues(),
            preserve_hydrogens=self.finish_page.preserve.isChecked(),
            chemistry_template_path=Path(self.components_page.template.text()).resolve() if self.components_page.template.text().strip() else None)


class ReceptorPreparationWorker(QObject):
    succeeded = pyqtSignal(str)
    failed = pyqtSignal(str)
    finished = pyqtSignal()

    def __init__(self, project_root: Path, plan: ReceptorPreparationPlan) -> None:
        super().__init__(); self.project_root = project_root; self.plan = plan

    @pyqtSlot()
    def run(self) -> None:
        try:
            self.succeeded.emit(str(prepare_receptor(self.project_root, self.plan)))
        except Exception as exc:
            self.failed.emit(str(exc))
        finally:
            self.finished.emit()
