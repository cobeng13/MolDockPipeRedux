from __future__ import annotations

import sys
import codecs
from pathlib import Path

from PyQt6.QtCore import QProcess, QProcessEnvironment
from PyQt6.QtWidgets import (
    QApplication, QDialog, QFileDialog, QHBoxLayout, QLabel, QMessageBox,
    QPushButton, QTextEdit, QVBoxLayout,
)


class TroubleshooterDialog(QDialog):
    def __init__(self, root: Path, parent=None) -> None:
        super().__init__(parent)
        self.root = root
        self.setWindowTitle("MolDockPipe Troubleshooter")
        self.resize(850, 620)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Test RDKit, MolScrub, Meeko, Vina CPU and Vina GPU in the app's environment.\n"
            "Uses temporary test molecules. GPU checks run a small docking job; allow a few minutes.\n"
            "Vina GPU is optional for CPU users. Reports include local paths; review before sharing."
        ))
        self.report = QTextEdit()
        self.report.setReadOnly(True)
        layout.addWidget(self.report)
        buttons = QHBoxLayout()
        self.run_button = QPushButton("Run checks")
        self.run_button.clicked.connect(self.run_checks)
        buttons.addWidget(self.run_button)
        copy = QPushButton("Copy report")
        copy.clicked.connect(lambda: QApplication.clipboard().setText(self.report.toPlainText()))
        buttons.addWidget(copy)
        save = QPushButton("Save report…")
        save.clicked.connect(self.save_report)
        buttons.addWidget(save)
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        buttons.addWidget(close)
        layout.addLayout(buttons)
        self.process = QProcess(self)
        self.decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
        self.process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        environment = QProcessEnvironment.systemEnvironment()
        environment.insert("PYTHONIOENCODING", "utf-8")
        self.process.setProcessEnvironment(environment)
        self.process.readyReadStandardOutput.connect(self.read_output)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)

    def run_checks(self) -> None:
        self.decoder.reset()
        self.run_button.setEnabled(False)
        self.report.setPlainText(f"Running checks using {sys.executable}\nPlease wait…\n\n")
        self.process.start(sys.executable, ["-u", "-m", "moldockpipe", "doctor", "--project", str(self.root)])

    def read_output(self) -> None:
        cursor = self.report.textCursor()
        cursor.movePosition(cursor.MoveOperation.End)
        cursor.insertText(self.decoder.decode(bytes(self.process.readAllStandardOutput())))
        self.report.setTextCursor(cursor)

    def finished(self, code: int, status: QProcess.ExitStatus) -> None:
        self.read_output()
        self.run_button.setEnabled(True)
        if status == QProcess.ExitStatus.CrashExit:
            self.report.append("The diagnostic process crashed. See output above.")
        elif code:
            self.report.append("One or more checks failed. See details and suggested actions above.")

    def process_error(self, error: QProcess.ProcessError) -> None:
        self.report.append(self.process.errorString())
        if error == QProcess.ProcessError.FailedToStart:
            self.run_button.setEnabled(True)

    def save_report(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save diagnostic report", "moldockpipe-diagnostics.txt", "Text (*.txt)")
        if path:
            try:
                Path(path).write_text(self.report.toPlainText(), encoding="utf-8")
            except OSError as exc:
                QMessageBox.warning(self, "Could not save report", str(exc))

    def done(self, result: int) -> None:
        # Keep the dialog alive while its bounded subprocesses finish, including on Escape.
        if self.process.state() != QProcess.ProcessState.NotRunning:
            self.report.append("Checks are still running. Close this window after they finish.")
            return
        super().done(result)

    def closeEvent(self, event) -> None:
        if self.process.state() != QProcess.ProcessState.NotRunning:
            event.ignore()
            return
        super().closeEvent(event)
