from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QPainter
from PyQt6.QtWidgets import QWidget


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
