"""Optional WebEngine frontend. Importing this module needs only ordinary Qt."""
from __future__ import annotations

import base64
from dataclasses import asdict
import json
from pathlib import Path
import re

from PyQt6.QtCore import QObject, QTimer, QUrl, Qt, pyqtSignal, pyqtSlot
from PyQt6.QtGui import QImage
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QFileDialog, QGridLayout, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QVBoxLayout, QWidget,
)

from ..services.structure_viewer import ViewerStructure
from ..services.interactions import TYPES


class ViewerBridge(QObject):
    command = pyqtSignal(str)
    ready_received = pyqtSignal()
    rendered_received = pyqtSignal(int)
    error_received = pyqtSignal(int, str)
    png_received = pyqtSignal(int, str)

    @pyqtSlot()
    def ready(self):
        self.ready_received.emit()

    @pyqtSlot(int)
    def rendered(self, token):
        self.rendered_received.emit(token)

    @pyqtSlot(int, str)
    def error(self, token, message):
        self.error_received.emit(token, message)

    @pyqtSlot(int, str)
    def png(self, token, data):
        self.png_received.emit(token, data)


def _web_view(parent: QWidget, bridge: ViewerBridge):
    # Keep WebEngine optional, including when the Results dialog is imported.
    from PyQt6.QtWebChannel import QWebChannel
    from PyQt6.QtWebEngineCore import QWebEnginePage, QWebEngineProfile, QWebEngineUrlRequestInterceptor
    from PyQt6.QtWebEngineWidgets import QWebEngineView

    assets = Path(__file__).parent / "assets" / "molecular_viewer"
    page_url = QUrl.fromLocalFile(str((assets / "viewer.html").resolve()))
    allowed_files = {str((assets / name).resolve()).casefold()
                     for name in ("viewer.html", "viewer.js", "3Dmol-2.5.5.min.js")}

    class LocalOnly(QWebEngineUrlRequestInterceptor):
        def interceptRequest(self, info):
            url = info.requestUrl()
            allowed = ((url.scheme() == "file" and str(Path(url.toLocalFile()).resolve()).casefold() in allowed_files)
                       or (url.scheme() == "qrc" and url.path() == "/qtwebchannel/qwebchannel.js"))
            info.block(not allowed)

    class LocalPage(QWebEnginePage):
        def acceptNavigationRequest(self, url, navigation_type, is_main_frame):
            return url == page_url

        def createWindow(self, window_type):
            return None

    view = QWebEngineView(parent)
    profile = QWebEngineProfile(view)  # Off-the-record: no disk cookies or cache.
    interceptor = LocalOnly(profile)
    profile.setUrlRequestInterceptor(interceptor)
    page = LocalPage(profile, view)
    # QWebEnginePage releases its engine state before QObject deletes children.
    # Keep the profile alive through that teardown, including parent destruction
    # at application exit; sibling ownership under the view deletes it too early.
    profile.setParent(page)
    view.setPage(page)
    channel = QWebChannel(page)
    channel.registerObject("bridge", bridge)
    page.setWebChannel(channel)
    view._viewer_objects = (profile, interceptor, page, channel)
    view.setContextMenuPolicy(Qt.ContextMenuPolicy.NoContextMenu)
    view.load(page_url)
    return view


class MolecularViewer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bridge = ViewerBridge(self)
        self._ready = False
        self._available = False
        self._token = 0
        self._payload = {"kind": "clear", "token": 0}
        self._structure = None
        self._export_path = None
        self._rendered_token = None
        self.title = QLabel("Select a completed docking result.")
        self.title.setWordWrap(True)
        self.status = QLabel("Starting offline 3D viewer…")
        self.status.setWordWrap(True)
        self.pose = QComboBox()
        self.pose.setEnabled(False)
        self.representation = QComboBox()
        self.representation.addItems(["Cartoon", "Sticks"])
        self.show_receptor = QCheckBox("Show receptor")
        self.show_receptor.setChecked(True)
        self.focus = QPushButton("Focus ligand")
        self.reset = QPushButton("Reset view")
        self.export = QPushButton("Export PNG…")
        controls = QHBoxLayout()
        for control in (self.representation, self.show_receptor, self.focus, self.reset, self.export):
            controls.addWidget(control)
        self.interaction_toggles = {}
        interaction_controls = QGridLayout()
        for index, (kind, (label, color)) in enumerate(TYPES.items()):
            toggle = QCheckBox(label)
            toggle.setStyleSheet(f"QCheckBox {{ color: {color}; }}")
            toggle.setChecked(kind != "hydrophobic")
            toggle.toggled.connect(self._style)
            self.interaction_toggles[kind] = toggle
            interaction_controls.addWidget(toggle, index // 3, index % 3)
        self.interaction_labels = QCheckBox("Contact labels")
        self.interaction_labels.setChecked(True)
        self.interaction_labels.toggled.connect(self._style)
        interaction_controls.addWidget(self.interaction_labels, 1, 2)
        self.interaction_info = QLabel()
        self.interaction_info.setWordWrap(True)
        self.interaction_info.setTextFormat(Qt.TextFormat.PlainText)
        self.criteria = QPushButton("Interaction criteria…")
        self.criteria.clicked.connect(self._show_criteria)
        layout = QVBoxLayout(self)
        layout.addWidget(self.title)
        layout.addWidget(self.pose)
        layout.addLayout(controls)
        layout.addLayout(interaction_controls)
        layout.addWidget(self.interaction_info)
        layout.addWidget(self.criteria)
        layout.addWidget(self.status)
        self._set_controls(False)
        self.bridge.ready_received.connect(self._on_ready)
        self.bridge.rendered_received.connect(self._on_rendered)
        self.bridge.error_received.connect(self._on_error)
        self.bridge.png_received.connect(self._save_png)
        self.representation.currentTextChanged.connect(self._style)
        self.show_receptor.toggled.connect(self._style)
        self.focus.clicked.connect(lambda: self._send("focus"))
        self.reset.clicked.connect(lambda: self._send("reset"))
        self.export.clicked.connect(self._export_png)
        self.timer = QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self._timeout)
        try:
            self.web = _web_view(self, self.bridge)
            layout.addWidget(self.web, 1)
            self.web.loadFinished.connect(self._page_loaded)
            self.web.page().renderProcessTerminated.connect(lambda *_: self._unavailable("The viewer process stopped. Reopen Results to retry."))
            self._available = True
            self.timer.start(30000)
        except (ImportError, RuntimeError, OSError) as exc:
            self.web = None
            self._unavailable(f"3D viewer unavailable: {exc}. Follow the optional 3D viewer installation in README, then restart.")

    @property
    def available(self):
        return self._available

    def _set_controls(self, enabled):
        for control in (self.representation, self.show_receptor, self.focus, self.reset, self.export,
                        self.interaction_labels, *self.interaction_toggles.values()):
            control.setEnabled(enabled)

    def _unavailable(self, message):
        self._available = False
        self._ready = False
        self.timer.stop()
        self._set_controls(False)
        self.status.setText(message)

    def _page_loaded(self, ok):
        if not ok:
            self._unavailable("The bundled viewer could not load. Check the application installation and reopen Results.")

    def _timeout(self):
        if not self._ready:
            self._unavailable("3D viewer did not start. Check WebEngine/WebGL support and reopen Results.")
        else:
            self._export_path = None
            self._on_error(self._token, "The viewer operation timed out. Select the result again to retry.")

    def _on_ready(self):
        if not self._available:
            return
        self.timer.stop()
        self._ready = True
        self.bridge.command.emit(json.dumps(self._payload))

    def _send(self, kind, **values):
        if self._ready:
            self.bridge.command.emit(json.dumps({"kind": kind, "token": self._token, **values}))

    def clear(self, message):
        self._token += 1
        self._structure = None
        self.interaction_info.clear()
        for kind, (label, _) in TYPES.items():
            self.interaction_toggles[kind].setText(label)
        self._export_path = None
        self._rendered_token = None
        self._payload = {"kind": "clear", "token": self._token}
        self.title.setText(message)
        self._set_controls(False)
        if self._available:
            self.status.setText(message)
        if self._ready:
            self.bridge.command.emit(json.dumps(self._payload))

    def set_structure(self, structure: ViewerStructure):
        self._token += 1
        self._structure = structure
        self._rendered_token = None
        self._export_path = None
        self.title.setText(f"{structure.receptor_name} · {structure.parent_id}\n"
                           f"{structure.state_id} · Mode {structure.mode_index} · {structure.affinity:g} kcal/mol")
        self._payload = {"kind": "load", "token": self._token, **asdict(structure),
                         **self._style_values()}
        for kind, (label, _) in TYPES.items():
            count = sum(c["kind"] == kind for c in structure.interactions)
            self.interaction_toggles[kind].setText(f"{label} ({count})")
        self.interaction_info.setText(
            structure.interaction_notes if structure.interaction_notes.startswith("Interaction analysis unavailable")
            else f"{len(structure.interactions)} predicted contacts. Dashed lines show geometry; H-bond candidates may lack donor H atoms.")
        self.interaction_info.setToolTip(structure.interaction_notes)
        self._set_controls(False)
        self.status.setText("Rendering receptor and ligand…")
        if self._ready:
            self.bridge.command.emit(json.dumps(self._payload))
        self.timer.start(30000)

    def _on_rendered(self, token):
        if token == self._token and self._structure is not None:
            self.timer.stop()
            self._rendered_token = token
            self._set_controls(True)
            self.status.setText("Drag to rotate · Scroll to zoom")

    def _on_error(self, token, message):
        if token == -1:
            self._unavailable(f"3D viewer unavailable: {message}. Check WebGL support and reopen Results.")
        elif token == self._token:
            self.timer.stop()
            self._export_path = None
            self._set_controls(self._rendered_token == token)
            self.status.setText(f"Viewer error: {message}")

    def _style(self, *_):
        self._send("style", **self._style_values())

    def _style_values(self):
        return dict(representation=self.representation.currentText(), show_receptor=self.show_receptor.isChecked(),
                    interaction_types=[kind for kind, toggle in self.interaction_toggles.items() if toggle.isChecked()],
                    interaction_labels=self.interaction_labels.isChecked())

    def _show_criteria(self):
        QMessageBox.information(self, "Predicted interaction criteria",
            "H-bonds: donor–acceptor ≤3.5 Å. With explicit donor H: H–acceptor ≤2.6 Å and D–H–A ≥120°. "
            "Without donor H: all heavy-neighbor–D–A angles ≥90°; labeled candidate.\n\n"
            "π–π: ring centers ≤5.5 Å, plane angle ≤30° (parallel) or ≥60° (T-shaped), minimum ring-plane offset ≤2 Å.\n"
            "Cation–π: positive group–ring center ≤6 Å, in-plane offset ≤2 Å.\n"
            "Salt bridges: opposite charge centers ≤4 Å.\n"
            "Hydrophobic: nonpolar C···C ≤4 Å, shortest contact per residue.\n\n"
            "Ligand donors/acceptors use RDKit features and charges use the SDF formal charges. "
            "Protein roles use standard amino-acid templates, assuming ASP/GLU negative and LYS/ARG positive. "
            "Histidine H-bond roles require explicit hydrogens. Waters, metals, cofactors and modified residues are excluded. "
            "No protonation, coordinates or docking scores are changed. These annotations are geometric predictions, "
            "not validated bonds or interaction energies. Full criteria and sources are in README.md.")

    def _export_png(self):
        if self._structure is None or self._rendered_token != self._token:
            return
        structure = self._structure
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{structure.receptor_name}_{structure.parent_id}_mode_{structure.mode_index}")
        path, _ = QFileDialog.getSaveFileName(self, "Export molecular view", safe + ".png", "PNG image (*.png)")
        if not path:
            return
        if not path.lower().endswith(".png"):
            path += ".png"
        self._export_path = path
        self._set_controls(False)
        self.status.setText("Rendering PNG…")
        self.timer.start(30000)
        self._send("export")

    def _save_png(self, token, data):
        if token != self._token or self._export_path is None:
            return
        self.timer.stop()
        path, self._export_path = self._export_path, None
        try:
            prefix = "data:image/png;base64,"
            if not data.startswith(prefix):
                raise ValueError("The viewer did not return a PNG image.")
            contents = base64.b64decode(data[len(prefix):], validate=True)
            if not contents.startswith(b"\x89PNG\r\n\x1a\n") or QImage.fromData(contents, "PNG").isNull():
                raise ValueError("The viewer returned an invalid PNG image.")
            # Atomic write preserves any existing image if writing fails.
            from PyQt6.QtCore import QIODevice, QSaveFile
            output = QSaveFile(path)
            if not output.open(QIODevice.OpenModeFlag.WriteOnly):
                raise OSError(output.errorString())
            if output.write(contents) != len(contents):
                output.cancelWriting()
                raise OSError(output.errorString())
            if not output.commit():
                raise OSError(output.errorString())
            self.status.setText(f"Saved PNG: {path}")
        except (ValueError, OSError) as exc:
            self.status.setText(f"PNG export failed: {exc}")
            QMessageBox.warning(self, "PNG export failed", str(exc))
        finally:
            self._set_controls(True)
