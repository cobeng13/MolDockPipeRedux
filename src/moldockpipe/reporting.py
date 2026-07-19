from __future__ import annotations

import html
import json
import statistics
import uuid
from importlib.metadata import PackageNotFoundError, version
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from .fingerprints import file_sha256
from .project import ProjectRepository, utc_now


REPORT_SCHEMA = "moldockpipe-report-v1"


def _application_version() -> str:
    try:
        return version("moldockpipe")
    except PackageNotFoundError:
        return "development"


def _decode(value: object, fallback: Any) -> Any:
    if value in (None, ""):
        return fallback
    try:
        return json.loads(str(value))
    except (TypeError, ValueError, json.JSONDecodeError):
        return fallback


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    position = (len(values) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    weight = position - lower
    return values[lower] * (1 - weight) + values[upper] * weight


def _score_statistics(values: Iterable[float]) -> dict[str, object]:
    scores = sorted(float(value) for value in values)
    if not scores:
        return {"count": 0}
    return {
        "count": len(scores),
        "best_minimum": round(scores[0], 3),
        "worst_maximum": round(scores[-1], 3),
        "mean": round(statistics.fmean(scores), 3),
        "median": round(statistics.median(scores), 3),
        "q1": round(float(_percentile(scores, 0.25)), 3),
        "q3": round(float(_percentile(scores, 0.75)), 3),
        "threshold_counts": {
            "at_or_below_-7": sum(score <= -7 for score in scores),
            "at_or_below_-8": sum(score <= -8 for score in scores),
            "at_or_below_-9": sum(score <= -9 for score in scores),
        },
    }


class ReportDataBuilder:
    """Build a bounded report summary from structured project provenance."""

    def __init__(self, repository: ProjectRepository) -> None:
        self.repository = repository

    def build(self, *, report_id: str, title: str,
              source_manifest: dict[str, object]) -> dict[str, object]:
        settings = self.repository.get_settings()
        profiles = self.repository.get_receptor_profiles(include_archived=True)
        with self.repository.connection() as conn:
            project = conn.execute("SELECT * FROM projects ORDER BY created_at LIMIT 1").fetchone()
            tools = [dict(row) for row in conn.execute("SELECT * FROM tool_installations ORDER BY tool_name").fetchall()]
            environment_row = conn.execute("SELECT value_json FROM project_settings WHERE key='environment'").fetchone()
            flow = {
                "active_parent_compounds": conn.execute("SELECT COUNT(*) FROM parent_ligands WHERE active=1").fetchone()[0],
                "screening_completed": conn.execute("SELECT COUNT(*) FROM screening_results WHERE active=1 AND status='completed'").fetchone()[0],
                "screening_failed": conn.execute("SELECT COUNT(*) FROM screening_results WHERE active=1 AND status='failed'").fetchone()[0],
                "active_molecular_states": conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1").fetchone()[0],
                "pdbqt_ready_states": conn.execute("SELECT COUNT(*) FROM molecular_states WHERE active=1 AND status='pdbqt_ready'").fetchone()[0],
            }
            campaign_count = conn.execute("SELECT COUNT(*) FROM workflow_runs WHERE workflow_type='docking'").fetchone()[0]
            campaign_rows = conn.execute("""SELECT * FROM workflow_runs WHERE workflow_type='docking'
                ORDER BY started_at DESC LIMIT 25""").fetchall()
            campaigns = []
            for row in campaign_rows:
                stages = [dict(stage) for stage in conn.execute(
                    "SELECT * FROM workflow_run_stages WHERE workflow_run_id=? ORDER BY started_at",
                    (row["workflow_run_id"],)).fetchall()]
                for stage in stages:
                    stage["summary"] = _decode(stage.pop("summary_json"), {})
                campaign = dict(row)
                campaign["settings"] = _decode(campaign.pop("settings_json"), {})
                campaign["environment"] = _decode(campaign.pop("environment_json"), {})
                campaign["stages"] = stages
                campaigns.append(campaign)

        receptors = [self._receptor(profile) for profile in profiles]
        return {
            "schema": REPORT_SCHEMA,
            "report_id": report_id,
            "title": title,
            "generated_at": utc_now(),
            "scope": {
                "receptor_profile_ids": [str(profile.get("id")) for profile in profiles],
                "docking_runs": "current runs for active molecular states",
                "redocking_runs": "latest run per receptor profile",
                "docking_campaign_history_limit": 25,
            },
            "source_manifest": source_manifest,
            "project": {
                "project_id": project["project_id"] if project else None,
                "name": settings.get("name", self.repository.root.name),
                "created_at": project["created_at"] if project else None,
                "schema_version": settings.get("schema_version"),
                "application": {"name": "MolDockPipe Redux", "version": _application_version()},
                "environment": _decode(environment_row[0], {}) if environment_row else {},
                "tools": tools,
            },
            "pipeline_flow": flow,
            "docking_campaigns": {"total_count": campaign_count, "recent": campaigns},
            "receptors": receptors,
            "limitations": [
                "Vina affinities are ranking estimates and are not experimental binding free energies.",
                "Redocking remains ARTIFACTS_READY until an external RMSD result is recorded.",
                "Flexible-receptor, covalent, and metal-specific docking are outside the current workflow.",
                "Inferred ligand chemistry and deliberately excluded receptor residues should be reviewed.",
            ],
        }

    def _receptor(self, profile: dict[str, Any]) -> dict[str, object]:
        profile_id = str(profile.get("id"))
        with self.repository.connection() as conn:
            prep_row = conn.execute("""SELECT * FROM receptor_preparation_runs
                WHERE receptor_profile_id=? ORDER BY started_at DESC LIMIT 1""", (profile_id,)).fetchone()
            redock_row = conn.execute("""SELECT * FROM redocking_runs WHERE receptor_profile_id=?
                ORDER BY COALESCE(finished_at,interrupted_at,started_at) DESC LIMIT 1""", (profile_id,)).fetchone()
            docking_counts = conn.execute("""SELECT COUNT(*) total_runs,
                SUM(CASE WHEN d.status='completed' THEN 1 ELSE 0 END) completed_runs,
                SUM(CASE WHEN d.status='failed' THEN 1 ELSE 0 END) failed_runs,
                SUM(CASE WHEN d.status='interrupted' THEN 1 ELSE 0 END) interrupted_runs,
                COUNT(DISTINCT CASE WHEN d.status='completed' THEN s.parent_id END) docked_parent_compounds,
                COUNT(DISTINCT CASE WHEN d.status='completed' THEN d.state_id END) docked_states
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                WHERE d.receptor_profile_id=? AND d.is_current=1 AND s.active=1""", (profile_id,)).fetchone()
            pose_count = conn.execute("""SELECT COUNT(*) FROM docking_poses p JOIN docking_runs d ON d.run_id=p.run_id
                JOIN molecular_states s ON s.state_id=d.state_id WHERE d.receptor_profile_id=?
                AND d.is_current=1 AND d.status='completed' AND s.active=1""", (profile_id,)).fetchone()[0]
            best_rows = conn.execute("""WITH ranked AS (
                SELECT s.parent_id,s.state_id,d.run_id,p.mode_index,p.affinity,
                    ROW_NUMBER() OVER (PARTITION BY s.parent_id ORDER BY p.affinity ASC,d.ended_at DESC) rank
                FROM docking_poses p JOIN docking_runs d ON d.run_id=p.run_id
                JOIN molecular_states s ON s.state_id=d.state_id JOIN parent_ligands l ON l.parent_id=s.parent_id
                WHERE d.receptor_profile_id=? AND d.is_current=1 AND d.status='completed'
                    AND s.active=1 AND l.active=1)
                SELECT parent_id,state_id,run_id,mode_index,affinity FROM ranked WHERE rank=1
                ORDER BY affinity ASC,parent_id""", (profile_id,)).fetchall()
            failure_rows = conn.execute("""SELECT COALESCE(d.reason,'Unspecified failure') reason,COUNT(*) count
                FROM docking_runs d JOIN molecular_states s ON s.state_id=d.state_id
                WHERE d.receptor_profile_id=? AND d.is_current=1 AND d.status!='completed' AND s.active=1
                GROUP BY COALESCE(d.reason,'Unspecified failure') ORDER BY count DESC,reason""", (profile_id,)).fetchall()

            preparation = None
            if prep_row:
                preparation = dict(prep_row)
                for source, target, fallback in (("inventory_json", "inventory", {}), ("decisions_json", "decisions", {}),
                                                  ("command_json", "command", {}), ("artifacts_json", "artifacts", {}),
                                                  ("warnings_json", "warnings", [])):
                    preparation[target] = _decode(preparation.pop(source), fallback)

            redocking = None
            if redock_row:
                redocking = dict(redock_row)
                redocking["settings"] = _decode(redocking.pop("settings_json"), {})
                redocking["fingerprints"] = _decode(redocking.pop("fingerprints_json"), {})
                redocking["stages"] = []
                for stage in conn.execute("SELECT * FROM redocking_stages WHERE run_id=? ORDER BY started_at",
                                          (redock_row["run_id"],)).fetchall():
                    item = dict(stage); item["artifacts"] = _decode(item.pop("artifacts_json"), {})
                    redocking["stages"].append(item)
                poses = conn.execute("SELECT pose_rank,affinity FROM redocking_poses WHERE run_id=? ORDER BY pose_rank",
                                     (redock_row["run_id"],)).fetchall()
                redocking["pose_count"] = len(poses)
                redocking["affinities"] = [row["affinity"] for row in poses]
                redocking["dockrmsd"] = {"status": "Not calculated", "pose_rank": None, "rmsd_angstrom": None,
                                          "tool": None, "notes": None}
                transformation = self.repository.root / "inputs" / "receptors" / profile_id / "redocking" / str(redock_row["run_id"]) / "transformation.json"
                redocking["transformation"] = _decode(transformation.read_text(encoding="utf-8"), {}) if transformation.is_file() else {}

        best = [dict(row) for row in best_rows]
        parameters = {key: profile.get(key) for key in (
            "center_x", "center_y", "center_z", "size_x", "size_y", "size_z",
            "exhaustiveness", "num_modes", "energy_range", "seed", "cpu_count")}
        return {
            "profile": {key: value for key, value in profile.items() if key != "reference_ligand"},
            "reference_ligand": profile.get("reference_ligand"),
            "preparation": preparation,
            "redocking": redocking,
            "docking": {
                "counts": {**dict(docking_counts), "pose_count": pose_count},
                "parameters": parameters,
                "best_affinity_per_parent_statistics": _score_statistics(row["affinity"] for row in best),
                "top_compounds": best[:25],
                "failure_reasons": [dict(row) for row in failure_rows],
            },
        }


def _write_source_manifest(repository: ProjectRepository, destination: Path) -> dict[str, object]:
    queries = {
        "workflow_runs": "SELECT workflow_run_id FROM workflow_runs ORDER BY workflow_run_id",
        "receptor_preparation_runs": "SELECT preparation_run_id FROM receptor_preparation_runs ORDER BY preparation_run_id",
        "redocking_runs": "SELECT run_id FROM redocking_runs ORDER BY run_id",
        "current_docking_runs": "SELECT run_id FROM docking_runs WHERE is_current=1 ORDER BY run_id",
    }
    counts: dict[str, int] = {}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with repository.connection() as conn, destination.open("w", encoding="utf-8") as handle:
        handle.write('{"schema":"moldockpipe-source-runs-v1","runs":{')
        for query_index, (name, query) in enumerate(queries.items()):
            if query_index:
                handle.write(",")
            handle.write(json.dumps(name) + ":[")
            count = 0
            cursor = conn.execute(query)
            while True:
                rows = cursor.fetchmany(1000)
                if not rows:
                    break
                for row in rows:
                    if count:
                        handle.write(",")
                    handle.write(json.dumps(str(row[0])))
                    count += 1
            handle.write("]")
            counts[name] = count
        handle.write("}}")
    return {"path": destination.relative_to(repository.root).as_posix(),
            "sha256": file_sha256(destination), "counts": counts}


def _cell(value: object) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return f"{value:.3f}"
    if isinstance(value, (dict, list, tuple)):
        return html.escape(json.dumps(value, sort_keys=True))
    return html.escape(str(value))


def _table(headers: list[str], rows: Iterable[Iterable[object]]) -> str:
    body = "".join("<tr>" + "".join(f"<td>{_cell(value)}</td>" for value in row) + "</tr>" for row in rows)
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    return f"<div class='table-wrap'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def render_html(data: dict[str, object]) -> str:
    project = data["project"]
    flow = data["pipeline_flow"]
    receptor_sections = []
    for receptor in data["receptors"]:
        profile = receptor["profile"]
        prep = receptor["preparation"]
        redocking = receptor["redocking"]
        docking = receptor["docking"]
        preparation_html = "<p class='warning'>No receptor-preparation provenance is available for this imported profile.</p>"
        if prep:
            inventory = prep["inventory"]; decisions = prep["decisions"]
            counts = inventory.get("counts", {})
            preparation_html = (
                _table(["Source", "Model", "Chains", "Protein residues", "Waters", "Components"], [[
                    prep.get("source_path"), inventory.get("selected_model"), inventory.get("included_chains"),
                    counts.get("protein_residues"), counts.get("waters"), counts.get("nonpolymer_components")]])
                + _table(["Center method", "Center inputs", "Center", "Box method", "Box inputs", "Box size"], [[
                    decisions.get("center_method"), decisions.get("center_parameters"), decisions.get("box_center"),
                    decisions.get("box_method"), decisions.get("box_parameters"), decisions.get("box_size")]])
                + _table(["Decision", "Recorded value"], [
                    ["Reference ligand", decisions.get("reference_ligand")],
                    ["Components removed", decisions.get("removed_residues")],
                    ["Components retained", decisions.get("retained_components")],
                    ["Altloc choices", decisions.get("altloc_choices")],
                    ["Incomplete receptor residues excluded", decisions.get("excluded_receptor_residues")],
                    ["Preserve existing hydrogens", decisions.get("preserve_hydrogens")],
                ])
                + _table(["Component", "Category", "Suggested role", "Heavy atoms", "Reason"], [
                    [item.get("identity"), item.get("category"), item.get("suggested_role"),
                     item.get("heavy_atom_count"), item.get("classification_reason")]
                    for item in inventory.get("components", [])])
                + f"<p><strong>Status:</strong> {_cell(prep.get('status'))} &nbsp; <strong>Warnings:</strong> {_cell(prep.get('warnings'))}</p>"
            )
            incomplete = inventory.get("incomplete_protein_residues", [])
            if incomplete:
                preparation_html += _table(["Incomplete residue", "Missing atoms", "Altlocs", "Recommended"], [[
                    item.get("identity"), item.get("missing_atoms"), item.get("alternate_locations"), item.get("recommended_altloc")]
                    for item in incomplete])
        redocking_html = "<p>Not run.</p>"
        if redocking:
            redocking_html = _table(["Run", "Status", "Poses", "Settings", "DockRMSD"], [[
                redocking.get("run_id"), redocking.get("status"), redocking.get("pose_count"),
                redocking.get("settings"), redocking.get("dockrmsd", {}).get("status")]])
            transformation = redocking.get("transformation", {})
            reference = transformation.get("reference_ligand", {})
            if reference:
                redocking_html += _table(["Reference", "Source PDB", "Source SDF", "Chemistry source", "Heavy atoms", "Formal charge"], [[
                    reference.get("identity"), reference.get("source_pdb"), reference.get("source_sdf"),
                    reference.get("chemistry_source"), reference.get("heavy_atom_count"), reference.get("formal_charge")]])
            operations = transformation.get("operations", [])
            if operations:
                redocking_html += _table(["Stage", "Operation / implementation", "Coordinate or mapping policy"], [[
                    item.get("stage"), item.get("operation") or item.get("implementation"),
                    item.get("heavy_atom_coordinate_policy") or item.get("atom_mapping") or item.get("hydrogen_policy")]
                    for item in operations])
        stats = docking["best_affinity_per_parent_statistics"]
        docking_html = (
            _table(["Current runs", "Completed", "Failed", "Docked parents", "Docked states", "Poses"], [[
                docking["counts"].get("total_runs") or 0, docking["counts"].get("completed_runs") or 0,
                docking["counts"].get("failed_runs") or 0, docking["counts"].get("docked_parent_compounds") or 0,
                docking["counts"].get("docked_states") or 0, docking["counts"].get("pose_count") or 0]])
            + _table(["Best", "Median", "Mean", "Worst", "Q1", "Q3", "N"], [[
                stats.get("best_minimum"), stats.get("median"), stats.get("mean"), stats.get("worst_maximum"),
                stats.get("q1"), stats.get("q3"), stats.get("count")]])
            + _table(["Docking parameters", "Value"], [[key, value] for key, value in docking["parameters"].items()])
            + _table(["Affinity threshold", "Parent compounds"], [[key, value]
                for key, value in stats.get("threshold_counts", {}).items()])
            + _table(["Rank", "Parent", "State", "Affinity (kcal/mol)", "Mode"], [
                [index, item.get("parent_id"), item.get("state_id"), item.get("affinity"), item.get("mode_index")]
                for index, item in enumerate(docking["top_compounds"], 1)])
        )
        if docking["failure_reasons"]:
            docking_html += _table(["Failure reason", "Count"], [[item.get("reason"), item.get("count")]
                                                                   for item in docking["failure_reasons"]])
        receptor_sections.append(f"""
        <section><h2>{html.escape(str(profile.get('name', profile.get('id'))))}</h2>
        <p class='muted'>Profile ID: {_cell(profile.get('id'))} · Receptor: {_cell(profile.get('receptor'))}</p>
        <h3>Receptor preparation</h3>{preparation_html}
        <h3>Redocking</h3>{redocking_html}
        <h3>Docking</h3>{docking_html}</section>""")
    limitations = "".join(f"<li>{html.escape(str(item))}</li>" for item in data["limitations"])
    campaign_rows = []
    for campaign in data["docking_campaigns"]["recent"]:
        stage = next((item for item in campaign.get("stages", []) if item.get("stage_name") == "vina"), {})
        summary = stage.get("summary", {})
        campaign_rows.append([campaign.get("workflow_run_id"), campaign.get("status"), campaign.get("started_at"),
                              summary.get("total"), summary.get("succeeded"), summary.get("failed"), summary.get("reused")])
    return f"""<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(str(data['title']))}</title>
<style>
body{{font-family:Segoe UI,Arial,sans-serif;max-width:1180px;margin:32px auto;padding:0 24px;color:#202124;line-height:1.45}}
h1{{margin-bottom:4px}} h2{{border-bottom:2px solid #2563eb;padding-bottom:6px;margin-top:42px}} h3{{margin-top:26px}}
.muted{{color:#64748b}} .warning{{background:#fff7ed;border-left:4px solid #f97316;padding:10px 14px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:18px 0}}
.card{{border:1px solid #dbe1e8;border-radius:6px;padding:12px;background:#f8fafc}} .card strong{{display:block;font-size:1.35rem}}
.table-wrap{{overflow-x:auto;margin:10px 0 18px}} table{{border-collapse:collapse;width:100%;font-size:.91rem}}
th,td{{border:1px solid #dbe1e8;padding:7px 9px;text-align:left;vertical-align:top}} th{{background:#eff6ff}}
code{{font-family:Consolas,monospace}} @media print{{body{{max-width:none;margin:0}} section{{break-inside:avoid-page}}}}
</style></head><body>
<h1>{html.escape(str(data['title']))}</h1><p class='muted'>Generated {_cell(data['generated_at'])} · Report {_cell(data['report_id'])}</p>
<h2>Project summary</h2><p><strong>{_cell(project.get('name'))}</strong> · MolDockPipe {_cell(project.get('application', {}).get('version'))} · Schema {_cell(project.get('schema_version'))}</p>
<div class='cards'>{''.join(f"<div class='card'><strong>{_cell(value)}</strong>{html.escape(key.replace('_',' ').title())}</div>" for key,value in flow.items())}</div>
<h3>Recorded tools</h3>{_table(['Tool','Version','Location','SHA-256'], [[item.get('tool_name'),item.get('version'),item.get('location'),item.get('sha256')] for item in project.get('tools',[])])}
<h3>Recent docking campaigns</h3>{_table(['Campaign','Status','Started','Requested','Successful','Failed','Reused'], campaign_rows)}
{''.join(receptor_sections)}
<section><h2>Limitations</h2><ul>{limitations}</ul>
<p class='muted'>Exact included run identifiers are recorded in {_cell(data['source_manifest'].get('path'))}.</p></section>
</body></html>"""


def generate_project_report(repository: ProjectRepository, *, title: str | None = None,
                            output_root: Path | None = None) -> Path:
    report_id = uuid.uuid4().hex
    generated = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    folder = output_root or (repository.root / "exports" / "reports" / f"{generated}_{report_id[:8]}")
    folder.mkdir(parents=True, exist_ok=False)
    report_title = title or f"{repository.get_settings().get('name', repository.root.name)} — MolDockPipe Report"
    manifest_path = folder / "source_run_manifest.json"
    source_manifest = _write_source_manifest(repository, manifest_path)
    data = ReportDataBuilder(repository).build(report_id=report_id, title=report_title,
                                               source_manifest=source_manifest)
    data_path = folder / "report_data.json"
    data_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    html_path = folder / "report.html"
    html_path.write_text(render_html(data), encoding="utf-8")
    for path, artifact_type in ((manifest_path, "report_source_manifest"),
                                (data_path, "report_data_json"), (html_path, "project_report_html")):
        repository.add_artifact(path, artifact_type, "report")
    repository.record_report_snapshot(
        report_id=report_id,
        title=report_title,
        scope=dict(data["scope"]),
        source_runs=source_manifest,
        output_path=html_path,
        output_sha256=file_sha256(html_path),
    )
    return html_path
