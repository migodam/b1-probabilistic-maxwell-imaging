"""Read-only inventory of paid Research G caches; never import a solver.

Every returned state refers to an existing NPZ member. Repeated independent
validation states remain separate records, but carry an exact array-equality
link so they cannot be mistaken for additional material worlds. The software
phase-field demonstration is explicitly excluded from physical candidates.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import time

import numpy as np


def _reference(path, key):
    return {"path": str(path.resolve()), "array_key": key}


def _complex(value):
    return [float(np.real(value)), float(np.imag(value))]


def _change_metrics(delta):
    """Array-only DC/low-frequency diagnostics, not a prior acceptance test."""
    value = np.asarray(delta).reshape(24, 24)
    energy = float(np.vdot(value, value).real)
    yy, xx = np.meshgrid(np.arange(24) + .5, np.arange(24) + .5, indexing="ij")
    modes = [(a, b) for a in range(3) for b in range(3) if a + b <= 2]
    basis = []
    for a, b in modes:
        v = (np.cos(np.pi*a*yy/24)*np.cos(np.pi*b*xx/24)).reshape(-1)
        basis.append(v/np.linalg.norm(v))
    coeff = np.asarray(basis) @ value.reshape(-1)
    return {
        "mean_delta_complex": _complex(np.mean(value)),
        "rms_delta_complex": float(np.sqrt(energy/576)),
        "dc_energy_fraction": float(576*abs(np.mean(value))**2/energy) if energy > 0 else 0.,
        "dct_total_degree_le2_energy_fraction": float(np.vdot(coeff, coeff).real/energy) if energy > 0 else 0.,
        "changed_cells_abs_gt_1e_10": int(np.count_nonzero(abs(value) > 1e-10)),
        "interpretation": "POSTHOC_MATERIAL_CHANGE_DIAGNOSTIC_NOT_PRIOR_ACCEPTANCE",
    }


def inventory(source_stage: Path):
    """Return only whitelisted exposed-development G states and provenance."""
    stage = Path(source_stage).resolve()
    repo = stage.parents[1]
    summary_path = stage / "results/E1_LEGACY/RESULTS.json"
    summary = json.loads(summary_path.read_text())
    allowed = {"p250069_weak_a1_snr20", "p250073_weak_a1_snr20", "p250077_weak_a1_snr20"}
    if set(summary["cases"]) != allowed:
        raise ValueError("G_EXPOSED_CASE_WHITELIST_MISMATCH")
    commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    bases = {r["case_id"]: r for r in summary["bases"]}
    candidate_index = json.loads((stage/"CANDIDATE_MATERIALS/index.json").read_text())
    archived_ids = {r["candidate_id"] for r in candidate_index["candidates"]}
    locations = {b["candidate_id"]: (stage/b["array_path"], b["array_prefix"]) for b in bases.values()}
    episodes = [e for method in summary["methods"] for e in method["episodes"]]
    locations.update({e["endpoint_id"]: (stage/e["array_path"], "validation") for e in episodes})
    states, trajectories, mismatches = [], [], []

    def read_state(path, prefix, **metadata):
        # Never inspect truth, geometry, or clean-response arrays. These member
        # names are fixed by G's truth-free online state serializer.
        with np.load(path, allow_pickle=False) as saved:
            keys = {name: prefix+"__"+name for name in ("chi", "prediction", "current", "field")}
            arrays = {name: np.asarray(saved[key]).copy() for name, key in keys.items()}
            observation = np.asarray(saved["__observation"]).copy()
            background = np.asarray(saved["__background"]).copy()
        for name, shape in (("chi", (576,)), ("prediction", (32, 8)), ("current", (576, 8)), ("field", (576, 8))):
            if arrays[name].shape != shape or arrays[name].dtype != np.complex128:
                raise ValueError("G_CACHED_LAYOUT_INVALID:"+str(path)+":"+name)
        sigma = float(bases[metadata["case_id"]]["sigma"])
        residual = float(np.linalg.norm(arrays["prediction"]-observation)*np.sqrt(2.)/sigma)
        original = metadata.pop("saved_residual", None)
        if original is not None and abs(residual-original) > 1e-8:
            mismatches.append(dict(path=str(path), prefix=prefix, saved=original, replay=residual))
        row = dict(**metadata, array_path=str(path.resolve()), array_prefix=prefix,
                   arrays={name: _reference(path, key) for name, key in keys.items()},
                   observation=_reference(path, "__observation"), background=_reference(path, "__background"),
                   sigma=sigma, residual_recomputed=residual, ell_recomputed=.5*residual**2,
                   noise_ball_compatible=residual <= float(summary["threshold"]),
                   shape={name: list(value.shape) for name, value in arrays.items()},
                   dtype="complex128", prior_status="PRIOR_UNKNOWN", original_constraint="BOX_ONLY",
                   P0_EXACT_evidence=False, posterior_mass=None)
        states.append(row)
        return row, arrays, background

    for case, b in sorted(bases.items()):
        read_state(stage/b["array_path"], b["array_prefix"], state_id="G:"+b["candidate_id"],
                   candidate_id=b["candidate_id"], case_id=case, role="G_BASE", parent_id=None,
                   in_success_archive=True, saved_residual=b["validation"]["full_residual"],
                   source_summary=str(summary_path), historic_cost_states=1, historic_cost_rhs=8)
    for e in episodes:
        path = stage/e["array_path"]
        parent_id = e["action"]["parent_id"]
        parent_path, parent_prefix = locations[parent_id]
        with np.load(parent_path, allow_pickle=False) as p:
            parent = p[parent_prefix+"__chi"].copy()
        step_arrays = {}
        ids = []
        for prefix in ("raw", "correction1", "correction2", "validation"):
            step = next((r for r in e["steps"] if r["array_prefix"] == prefix), None)
            saved_residual = e["full_residual"] if prefix == "validation" else step["full_residual"]
            sid = "G:"+e["endpoint_id"]+":"+prefix
            row, arrays, background = read_state(
                path, prefix, state_id=sid, candidate_id=e["endpoint_id"], case_id=e["case_id"],
                role={"raw":"G_RAW", "correction1":"G_CORRECTION1", "correction2":"G_CORRECTION2", "validation":"G_VALIDATED_ENDPOINT"}[prefix],
                round_id=e["round_id"], method=e["method"], trajectory=e["trajectory"],
                parent_id=parent_id, parent_material=_reference(parent_path, parent_prefix+"__chi"),
                raw_proposal=_reference(path, "__raw_proposal"),
                in_success_archive=prefix == "validation" and e["endpoint_id"] in archived_ids,
                saved_residual=saved_residual, source_trajectory=str((stage/"results/E1_LEGACY/trajectories"/(e["endpoint_id"]+".json")).resolve()),
                accepted_correction_step=step.get("nonlinear_step_accepted") if step else None,
                independent_fresh_validation=prefix == "validation", historic_cost_states=1,
                historic_cost_rhs=40 if prefix in ("raw", "correction1") else 8,
                array_equal_to_step=None)
            if prefix == "validation":
                equal = [k for k, v in step_arrays.items() if np.array_equal(v["chi"], arrays["chi"]) and np.array_equal(v["prediction"], arrays["prediction"])]
                row["array_equal_to_step"] = equal
            else:
                step_arrays[prefix] = arrays
            row["change_from_parent"] = _change_metrics(arrays["chi"]-parent)
            row["change_from_raw"] = _change_metrics(arrays["chi"]-step_arrays["raw"]["chi"])
            ids.append(sid)
        raw = step_arrays["raw"]["chi"]
        correction = arrays["chi"]-raw
        edit = raw-parent
        # A broad, reproducible support diagnostic. The G actions were not
        # protected cores; the one-percent set is descriptive, not a new gate.
        extent = abs(edit) > max(float(np.max(abs(edit)))*.01, 1e-12)
        trajectories.append(dict(endpoint_id=e["endpoint_id"], case_id=e["case_id"],
            method=e["method"], round_id=e["round_id"], parent_id=parent_id, operator=e["action"]["operator"],
            source_trajectory=states[-1]["source_trajectory"], states=ids,
            action=e["action"], raw_compatible=states[-4]["noise_ball_compatible"],
            corrected_compatible=states[-1]["noise_ball_compatible"], raw_ell=states[-4]["ell_recomputed"],
            endpoint_ell=states[-1]["ell_recomputed"], corrector_cost=e["cost_delta"],
            correction_retention=e.get("retention"), original_offer_status=e.get("offer_status"),
            raw_edit_support_cells_one_percent=int(extent.sum()),
            correction_rms_inside_raw_extent=float(np.sqrt(np.mean(abs(correction[extent])**2))) if extent.any() else None,
            correction_rms_outside_raw_extent=float(np.sqrt(np.mean(abs(correction[~extent])**2))) if (~extent).any() else None,
            correction_global_diagnostics=_change_metrics(correction)))
    return dict(schema="prior_constrained_G_cache_inventory_v1", source_repository=str(repo),
                source_commit=commit, source_stage=str(stage), source_summary=str(summary_path),
                exposure="HISTORICALLY_EXPOSED_DEVELOPMENT", target=summary["target"],
                threshold=summary["threshold"], packing="source-major Re32 then Im32; sqrt(2)/sigma",
                gain=1, new_full_states=0, new_Maxwell_RHS=0, source_history_mutated=False,
                cost_claim="Existing 291 states/6936 RHS are provenance only; not recharged as new physics",
                states=states, trajectories=trajectories,
                counts=dict(states=len(states), trajectories=len(trajectories), roles=dict(Counter(r["role"] for r in states)),
                            archived_success_entries=len(archived_ids),
                            raw_noise_compatible=sum(r["raw_compatible"] for r in trajectories),
                            endpoint_noise_compatible=sum(r["corrected_compatible"] for r in trajectories),
                            raw_valid_to_endpoint_invalid=sum(r["raw_compatible"] and not r["corrected_compatible"] for r in trajectories)),
                residual_replay_mismatches=mismatches,
                excluded_material_only_arrays=[dict(path=str(stage/"results/CONTINUOUS_REPRESENTATION.npz"),
                    keys=["z", "fields"], reason="SYNTHETIC_SOFTWARE_REPRESENTATION_ONLY; no Maxwell prediction")],
                prior_note="All original online legal checks are box constraints; no fixed increment-phase or P0 generator retraction. No PHASE/P0_EXACT acceptance inferred.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-stage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    result = inventory(args.source_stage)
    result["postprocess_wall_seconds"] = time.perf_counter()-started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(counts=result["counts"], replay_mismatches=len(result["residual_replay_mismatches"]),
                          new_full_states=0, new_Maxwell_RHS=0, output=str(args.output))))


if __name__ == "__main__":
    main()
