from __future__ import annotations

from PyQt6.QtWidgets import QDialog, QDialogButtonBox, QLabel, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget

from ..project import ProjectRepository


class DockingResultsDialog(QDialog):
    def __init__(self, repository: ProjectRepository, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.repository = repository
        self.setWindowTitle("Multidock Results")
        self.resize(900, 560)
        tabs = QTabWidget()
        profiles = repository.get_receptor_profiles(include_archived=True)
        for profile in profiles:
            tabs.addTab(self._profile_table(str(profile["id"])), str(profile.get("name", profile["id"])))
        if not profiles:
            tabs.addTab(QLabel("No receptor profiles configured."), "Results")
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close); buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self); layout.addWidget(tabs); layout.addWidget(buttons)

    def _profile_table(self, profile_id: str) -> QTableWidget:
        headers = ["Ligand", "Best state", "Affinity", "Mode", "Status", "Reason"]
        table = QTableWidget(0, len(headers)); table.setHorizontalHeaderLabels(headers)
        table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        with self.repository.connection() as conn:
            parents = conn.execute("SELECT parent_id FROM parent_ligands WHERE active=1 ORDER BY parent_id").fetchall()
            best = {row["parent_id"]: row for row in conn.execute("""WITH ranked AS (
                    SELECT s.parent_id, s.state_id, p.affinity, p.mode_index,
                           ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity, s.state_id, p.mode_index) rank
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    JOIN docking_poses p ON p.run_id=d.run_id
                    WHERE d.receptor_profile_id=? AND d.is_current=1 AND d.status='completed' AND s.active=1)
                SELECT * FROM ranked WHERE rank=1""", (profile_id,)).fetchall()}
            latest = {row["parent_id"]: row for row in conn.execute("""WITH ranked AS (
                    SELECT s.parent_id, d.status, COALESCE(d.reason, '') reason,
                           ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY d.started_at DESC, d.run_id DESC) rank
                    FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                    WHERE d.receptor_profile_id=? AND d.is_current=1 AND s.active=1)
                SELECT * FROM ranked WHERE rank=1""", (profile_id,)).fetchall()}
        parents = sorted(parents, key=lambda parent: (
            parent["parent_id"] not in best,
            best[parent["parent_id"]]["affinity"] if parent["parent_id"] in best else 0,
            parent["parent_id"],
        ))
        for parent in parents:
            parent_id = parent["parent_id"]
            score = best.get(parent_id); run = latest.get(parent_id)
            row = table.rowCount(); table.insertRow(row)
            values = [parent_id, score["state_id"] if score else "", score["affinity"] if score else "",
                      score["mode_index"] if score else "", "completed" if score else run["status"] if run else "pending",
                      "" if score else run["reason"] if run else ""]
            for column, value in enumerate(values):
                table.setItem(row, column, QTableWidgetItem(str(value)))
        return table
