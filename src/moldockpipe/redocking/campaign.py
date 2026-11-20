"""Plan and run a bounded redocking validation campaign."""

from __future__ import annotations

import math
from collections.abc import Callable

from ..execution_lock import project_execution_lock
from ..fingerprints import file_sha256
from ..project import ProjectRepository
from ..receptors.ad4zn import grid_definition, protocol_for
from .models import BOX_KEYS, SEARCH_SETTING_KEYS, RedockingSettings
from .runner import RedockingRunner

MAX_TRIALS = 21
TARGET_RMSD = 2.0
FINISHED_TRIALS = {"hit", "miss", "rmsd_missing"}


def docking_box(profile: dict) -> dict[str, float]:
    box = {key: float(profile[key]) for key in BOX_KEYS}
    if not all(math.isfinite(value) for value in box.values()) or min(box[key] for key in BOX_KEYS[3:]) <= 0:
        raise ValueError("The docking box needs finite coordinates and positive dimensions.")
    return box


def campaign_boxes(profile: dict, limit: int) -> list[dict[str, float]]:
    if not 1 <= limit <= MAX_TRIALS:
        raise ValueError(f"Trial cap must be between 1 and {MAX_TRIALS}.")
    baseline = docking_box(profile)
    zinc = protocol_for(profile) == "ad4zn"
    if zinc:
        grid_definition(profile)
    shifts = [(None, 0)] + [(key, delta) for key in BOX_KEYS[:3] for delta in (-2, 2)]
    boxes = []
    seen = set()
    for size_delta in (0, -2, 2):
        for center_key, shift in shifts:
            box = baseline.copy()
            if center_key is not None:
                box[center_key] += shift
            for key in BOX_KEYS[3:]:
                box[key] += size_delta
            try:
                docking_box(box)
                if zinc:
                    grid = grid_definition({**profile, **box})
                    identity = (tuple(grid["center"]), tuple(grid["npts"]))
                else:
                    identity = tuple(box.values())
            except ValueError:
                continue
            if identity in seen:
                continue
            seen.add(identity)
            boxes.append(box)
            if len(boxes) == limit:
                return boxes
    return boxes


def campaign_snapshot(repository: ProjectRepository, profile: dict) -> dict:
    reference = profile["reference_ligand"]
    paths = {"receptor": profile["receptor"], **{key: reference[key] for key in ("sdf", "mol2", "mapping")}}
    protocol = protocol_for(profile)
    return {
        "box": docking_box(profile),
        "protocol": protocol,
        "ad4zn": profile.get("ad4zn") if protocol == "ad4zn" else None,
        "reference_identity": reference["identity"],
        "inputs": {key: file_sha256(repository.root / str(path)) for key, path in paths.items()},
    }


def run_campaign(repository: ProjectRepository, item: dict, settings: RedockingSettings,
                 progress: Callable, cancel_event) -> dict:
    queue_id = item["queue_id"]
    with project_execution_lock(repository.root):
        trials = repository.redocking_campaign_trials(queue_id)
        for trial in trials:
            if trial["outcome"] == "hit":
                return {"campaign_id": queue_id, "run_id": trial["run_id"]}
            if trial["outcome"] in FINISHED_TRIALS:
                continue
            if cancel_event.is_set() or repository.redocking_cancel_requested(queue_id):
                raise InterruptedError("Redocking campaign stopped.")
            profile = next((p for p in repository.get_receptor_profiles()
                            if str(p["id"]) == item["receptor_profile_id"]), None)
            if profile is None or campaign_snapshot(repository, profile) != item["settings"]["campaign"]["snapshot"]:
                raise ValueError("The receptor, reference ligand, protocol, or saved box changed. Queue a new campaign.")
            trial_profile = {**profile, **trial["box"]}
            progress({"event": "campaign_trial", "trial": trial["trial_number"], "total": len(trials)})
            try:
                result = RedockingRunner(repository, trial_profile, settings, progress=progress,
                                         cancel_event=cancel_event, campaign_id=queue_id).run(trial["run_id"])
                rmsd = result["dockrmsd"]
                if not rmsd or rmsd["status"] == "Unavailable":
                    detail = rmsd["rank_1"].get("error") if rmsd else None
                    raise RuntimeError(f"DockRMSD is unavailable. {detail or 'Check the redocking RMSD log, then retry.'}")
                rank_1 = rmsd["rank_1"]
                value = rank_1["rmsd_angstrom"]
                if value is None or not math.isfinite(value) or value < 0:
                    outcome, value = "rmsd_missing", None
                else:
                    outcome = "hit" if value < TARGET_RMSD else "miss"
                repository.record_redocking_trial(trial["run_id"], outcome, rank_1["affinity"], value,
                                                  rank_1.get("error"))
            except InterruptedError:
                raise
            except Exception as exc:
                repository.record_redocking_trial(trial["run_id"], "error", error=str(exc))
                raise
            if outcome == "hit":
                return {**result, "campaign_id": queue_id}
        return {"campaign_id": queue_id}


def apply_winning_box(repository: ProjectRepository, queue_id: str) -> dict:
    with project_execution_lock(repository.root):
        item = next((row for row in repository.list_redocking_queue() if row["queue_id"] == queue_id), None)
        winner = next((trial for trial in repository.redocking_campaign_trials(queue_id)
                       if trial["outcome"] == "hit"), None)
        if not item or not winner or item["status"] != "completed":
            raise ValueError("This campaign has no completed winning trial.")
        profiles = repository.get_receptor_profiles(include_archived=True)
        profile = next((p for p in profiles if str(p["id"]) == item["receptor_profile_id"] and not p.get("archived")), None)
        if profile is None:
            raise ValueError("The receptor profile is missing or archived.")
        current = campaign_snapshot(repository, profile)
        starting = item["settings"]["campaign"]["snapshot"]
        if current != starting and current != {**starting, "box": winner["box"]}:
            raise ValueError("The receptor, reference ligand, protocol, or saved box changed. Queue a new campaign.")
        search_settings = {key: item["settings"][key] for key in SEARCH_SETTING_KEYS}
        updated = {**profile, **winner["box"], **search_settings}
        if updated != profile:
            repository.update_receptor_docking_defaults(profiles, updated)
        return updated
