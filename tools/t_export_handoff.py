"""Export paid Research T endpoints and replay saved arrays; ZERO new solves.

Only explicit qualified bases, explicit observation files, known background
physics arrays and genuine mechanism endpoint files are opened. Offline truth
audit archives, historical truth and G's workspace are not read.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from research_t.handoff import validate_bundle, whitened_residual, write_bundle
from research_t.structure import check_protected, material_distances

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_STAGE = ROOT / "research/research_t_20261015"
BACKGROUND_KEYS = ("chi_b", "eb", "yb", "Db", "G", "S", "incident", "volume", "grid",
                   "points", "source_points", "receiver_points", "descriptor")


def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def relative(path, base):
    return os.path.relpath(Path(path), Path(base))


def arrays(path, names):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key].copy() for key in names}


def relnorm(a, b):
    return float(np.linalg.norm(a) / max(float(np.linalg.norm(b)), 1e-300))


def cost_record(source, *, stage, bundle, wall, shared=None, selection_overhead=None):
    roles = source["RHS_by_role"]
    tangent = sum(v for key, v in roles.items() if "jvp" in key or "tangent" in key)
    adjoint = sum(v for key, v in roles.items() if "adjoint" in key)
    total = int(source["complex_full_solver_RHS"])
    return {"new_full_state_attempts": int(source["new_full_state_attempts"]),
            "factorizations": int(source["action_counts"].get("LU_factorizations", 0)),
            "forward_rhs": total - tangent - adjoint, "tangent_rhs": tangent,
            "adjoint_rhs": adjoint, "failed_attempts": int(source["failed_actions"]),
            "setup_seconds": None, "wall_seconds": wall,
            "setup_seconds_status": "not separately isolated; original shared setup receipt and production ledger preserved",
            "failed_attempts_semantics": "failed physical ACTION spans observed in the delimited production cost receipt",
            "cost_scope": "dedicated generation arm/trial; shared work and base selection overhead separately disclosed",
            "production_receipt": source,
            "production_ledger_ref": relative(stage / "COST_LEDGER.jsonl", bundle),
            "production_ledger_row_range_zero_based_half_open": source["ledger_rows"],
            "shared_setup_or_parent_replay_receipt": shared,
            "base_selection_other_arms": selection_overhead,
            "new_physics_during_export": {"full_state_attempts": 0, "complex_full_solver_RHS": 0},
            "historical_A30_A_B_P2_costs": "not reset, transferred, or included as new Research T costs"}


def saved_validation(state, obs, sigma, summary, background):
    residual = float(np.linalg.norm(whitened_residual(state["prediction"], obs, sigma=sigma)))
    x, current, field = state["chi"], state["current"], state["field"]
    constitutive = relnorm(current - x[:, None] * field, current)
    # Paid matrices only. These are algebraic applications, not inverse solves.
    field_from_current = background["incident"] + background["G"] @ current
    prediction_from_current = background["S"] @ current
    field_defect = relnorm(field - field_from_current, field)
    prediction_defect = relnorm(state["prediction"] - prediction_from_current, state["prediction"])
    legal = bool(np.all((x.real >= -.5-1e-12) & (x.real <= 4+1e-12)
                        & (x.imag >= -1e-12) & (x.imag <= .5+1e-12)))
    backward = float(state["backward_residual"])
    passed = legal and max(constitutive, field_defect, prediction_defect, backward) < 1e-10
    if not np.isclose(residual, summary["residual_norm"], rtol=1e-12, atol=1e-12):
        raise ValueError("saved endpoint and production summary residual differ")
    if not passed or residual > summary["tau"]:
        raise ValueError("selected endpoint is not actually qualified")
    return {"status": "FULL_VALIDATED_NUMERICAL", "matching_full_endpoint": True,
            "checks_passed": True, "observed_residual_norm": residual,
            "tau": float(summary["tau"]), "material_valid": legal,
            "constitutive_defect": constitutive, "constitutive_status": "CHECKED_FROM_SAVED_FULL_STATE",
            "backward_residual": backward, "output_error_upper_bound": None,
            "bound_provenance": None,
            "saved_field_equation_relative_defect": field_defect,
            "saved_prediction_equation_relative_defect": prediction_defect,
            "new_algebraic_G_columns": int(current.shape[1]),
            "new_algebraic_S_columns": int(current.shape[1]),
            "new_full_solver_rhs_for_replay": 0,
            "numerical_scope": "same finite 24x24 full system; no output enclosure or mesh-independent certificate"}


def physics_record(background, *, background_id, a30_revision, stage, bundle):
    desc = json.loads(str(background["descriptor"]))
    if desc["geometry"] != {"aperture": "full", "grid": 24, "n_rx": 32, "n_tx": 8,
                            "radius_rx": .8, "radius_tx": .75, "side": 1.0}:
        raise ValueError("frozen geometry metadata mismatch")
    return {"backend": "a30.physics.FrozenMaxwellAdapter -> unchanged vendored a24.physics.Physics2D/a2 VIE",
            "backend_revision": a30_revision, "dimension": 2, "field_model": "scalar TM electric field",
            "harmonic_sign": "exp(-i omega t)", "material_convention": "chi=epsilon_r-1; j=chi*field",
            "dispersion_law": "single frozen frequency only; no independently refitted frequency materials",
            "frequency_hz": float(desc["frequency_hz"]), "k_rad_per_m": 2*np.pi,
            "units": {"geometry": "metre", "chi": "dimensionless"},
            "grid": {"shape": [24, 24], "side_m": 1, "cell_mass_m2": float(background["volume"])},
            "green_convention": "g=i*Hankel1(0,k*r)/4; G includes k^2*cell area off diagonal",
            "quadrature": "off-diagonal cell midpoint collocation; original non-cell-integrated A2 default",
            "self_term": "k^2 times equal-area disk integral: .5j*pi*a/k*Hankel1(1,k*a)-1/k^2, a=h/sqrt(pi)",
            "state_normalization": "current=sqrt(cell_mass)*j; field=sqrt(cell_mass)*E_total",
            "response": "scattered channels S@current; no direct incident term in y",
            "source_coordinates_in_order": background["source_points"].tolist(),
            "receiver_coordinates_in_order": background["receiver_points"].tolist(),
            "background_id": background_id,
            "known_background_assets_ref": relative(stage/"T_PHYSICS_HANDOFF/backend_assets/background"/(background_id+".npz"), bundle),
            "material_shared_across_all_eight_sources": True}


def export_one(*, kind, case, candidate_id, state_path, state_summary,
               dedicated_cost, wall, shared, overhead, source_revision,
               a30_revision, stage, destination, background, lineage,
               core=None, guard=None, minimum_feature=None,
               parent_candidate_id=None):
    bundle = destination / "candidates" / candidate_id
    expected_state = arrays(state_path, ("chi", "current", "field", "prediction", "backward_residual"))
    data = arrays(stage/"qualification"/case/"observation.npz", ("observation", "sigma"))
    sigma = float(data["sigma"])
    validation = saved_validation(expected_state, data["observation"], sigma, state_summary, background)
    if bundle.exists():
        # Idempotent resume, never overwrite a prior published endpoint.
        existing = np.load(bundle/"material.npz", allow_pickle=False)
        if not np.array_equal(existing["chi"].reshape(-1), expected_state["chi"].reshape(-1)):
            raise ValueError(f"existing material differs; refusing overwrite {candidate_id}")
        result = validate_bundle(bundle)
        if not result["valid_bundle_semantics"]:
            raise ValueError(result)
        return {"candidate_id": candidate_id, "case": case, "kind": kind,
                "bundle_ref": relative(bundle, destination), "reused_exact_saved_bundle": True,
                "residual_norm": validation["observed_residual_norm"], "validation": result}
    if core is not None:
        protection = check_protected(expected_state["chi"], core, guard)
        if not protection["preserved"]:
            raise ValueError("genuine endpoint lost its protected core/guard")
    background_id = "weak" if "_weak_" in case else "ring"
    cost = cost_record(dedicated_cost, stage=stage, bundle=bundle, wall=wall,
                       shared=shared, selection_overhead=overhead)
    lineage.setdefault("generation_code_revision", "c752aa6fedd4d7f12d93062857013b721412135f")
    lineage.update(source_state_ref=relative(state_path, bundle),
                   original_observation_ref=relative(stage/"qualification"/case/"observation.npz", bundle),
                   export_source_revision=source_revision, truth_used_online=False,
                   verified_path_ref=None,
                   path_status="actual finite corrected endpoints retained; no feasible tube bound between endpoints",
                   producer_ledger_ref=relative(stage/"COST_LEDGER.jsonl", bundle))
    write_bundle(bundle, candidate_id=candidate_id,
                 parent_id=parent_candidate_id,
                 source_revision=source_revision, chi=expected_state["chi"],
                 observation=data["observation"], prediction=expected_state["prediction"],
                 coordinates=background["points"], mass=float(background["volume"]),
                 sigma=sigma, tau=state_summary["tau"],
                 physics=physics_record(background, background_id=background_id,
                                        a30_revision=a30_revision, stage=stage, bundle=bundle),
                 observation_id=case, noise_id="A25-original-noise-realization:"+case,
                 legal_domain={"real_chi": [-.5, 4], "imag_chi": [0, .5],
                               "shared_real_gain": [1, 1], "truth_or_topology_label_input": False},
                 validation=validation, cost=cost, lineage=lineage, core=core, guard=guard,
                 optional_state={"current": expected_state["current"], "field": expected_state["field"]},
                 minimum_feature_size=minimum_feature, minimum_feature_units="metre" if minimum_feature else None)
    result = validate_bundle(bundle)
    if not result["valid_bundle_semantics"]:
        raise ValueError(result)
    return {"candidate_id": candidate_id, "case": case, "kind": kind,
            "bundle_ref": relative(bundle, destination), "reused_exact_saved_bundle": False,
            "residual_norm": validation["observed_residual_norm"], "validation": result}


def archive_and_pairwise(records, destination):
    """Keep every real endpoint; report a deterministic material-space archive.

    Same-target candidates must meet BOTH frozen novelty criteria to create a
    new significant structural entry. This is a greedy archive rule, not an
    equivalence relation or a connected-component theorem.
    """
    pairs, materials, predictions = [], {}, {}
    for item in records:
        folder = destination/item["bundle_ref"]
        materials[item["candidate_id"]] = arrays(folder/"material.npz", ("chi",))["chi"]
        predictions[item["candidate_id"]] = arrays(folder/"prediction.npz", ("channels",))["channels"]
    for i, left in enumerate(records):
        for right in records[i+1:]:
            if left["case"] != right["case"]:
                continue
            a, b = left["candidate_id"], right["candidate_id"]
            d = material_distances(materials[a], materials[b], mass=1/576)
            data = arrays(destination/left["bundle_ref"]/"observation.npz", ("sigma_complex_rms",))
            pred_distance = float(np.linalg.norm(whitened_residual(predictions[a], predictions[b], sigma=float(data["sigma_complex_rms"]))))
            pairs.append({"case": left["case"], "candidate_a": a, "candidate_b": b,
                          "normalized_material_l2": d["normalized_material_l2"],
                          "phase_changed_cells": d["phase_changed_cells"],
                          "cavity_symmetric_difference_cells": d["enclosed_cavity_symmetric_difference_cells"],
                          "material_boundary_hausdorff_cells": d["material_boundary_hausdorff_cells"],
                          "pairwise_whitened_prediction_distance": pred_distance,
                          "distinct_under_frozen_structural_archive_rule": bool(d["normalized_material_l2"]>=.03 and d["phase_changed_cells"]>=6)})
    if pairs:
        with (destination/"PAIRWISE_DISTANCES.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(pairs[0]))
            writer.writeheader(); writer.writerows(pairs)
    lookup = {frozenset((p["candidate_a"],p["candidate_b"])):p for p in pairs}
    priority = {"qualified_minchange_base": 0, "qualified_base": 1, "protected_compensated_endpoint": 2}
    ordered = sorted(records, key=lambda r:(r["case"], priority[r["kind"]], r["candidate_id"]))
    retained, assignments = [], []
    for item in ordered:
        related = [r for r in retained if r["case"] == item["case"]]
        duplicate = next((r for r in related if not lookup[frozenset((r["candidate_id"],item["candidate_id"]))]["distinct_under_frozen_structural_archive_rule"]),None)
        if duplicate is None:
            retained.append(item)
        assignments.append({"candidate_id": item["candidate_id"], "case": item["case"],
                            "retained_as_structurally_distinct": duplicate is None,
                            "archive_representative": item["candidate_id"] if duplicate is None else duplicate["candidate_id"]})
    result = {"rule": "deterministic greedy archive, minchange default base first, then historical base, then chronological candidate id; novelty requires normalized_material_l2>=.03 AND phase_changed_cells>=6 against every retained same-target material",
              "not_an_equivalence_relation": True, "not_a_complete_world_or_connected_component_count": True,
              "raw_valid_endpoint_bundles": len(records), "retained_structural_representatives": len(retained),
              "per_case_retained_counts": {case:sum(r["case"]==case for r in retained) for case in sorted({r["case"] for r in records})},
              "assignments": assignments,
              "all_original_bundles_preserved": True, "pairwise_distance_ref": "PAIRWISE_DISTANCES.csv"}
    json_write(destination/"MATERIAL_ARCHIVE.json",result)
    return result


def preserve_actual_correction_constraints(records, destination):
    """Carry the actual source masks/boxes, including V2's extra phase lock."""
    for item in records:
        if item["kind"] != "protected_compensated_endpoint":
            continue
        folder = destination/item["bundle_ref"]
        lineage = json.loads((folder/"lineage.json").read_text())
        source = (folder/lineage["source_state_ref"]).resolve().parent/"edit_geometry.npz"
        with np.load(source, allow_pickle=False) as data:
            geometry = {key:data[key].copy() for key in data.files}
        output = folder/"correction_constraints.npz"
        if output.exists():
            with np.load(output, allow_pickle=False) as saved:
                if set(saved.files)!=set(geometry) or any(not np.array_equal(saved[k],v) for k,v in geometry.items()):
                    raise ValueError("existing correction geometry differs; not overwritten")
        else:
            np.savez_compressed(output,**geometry)
        chi = arrays(folder/"material.npz",("chi",))["chi"].reshape(-1)
        packed = np.r_[chi.real,chi.imag]
        violation = max(0.,float(np.max(geometry["lower"]-packed)),float(np.max(packed-geometry["upper"])))
        if violation>1e-8:
            raise ValueError("saved endpoint violates its actual correction bounds")
        lineage["original_correction_constraints_ref"] = relative(source,folder)
        lineage["fixed_legal_constraints_ref"] = "correction_constraints.npz"
        lineage["saved_correction_bounds_max_violation"] = violation
        json_write(folder/"lineage.json",lineage)
        manifest = json.loads((folder/"manifest.json").read_text())
        manifest["lineage"]["correction_constraints_ref"] = "correction_constraints.npz"
        json_write(folder/"manifest.json",manifest)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", type=Path, default=DEFAULT_STAGE)
    parser.add_argument("--a30-root", type=Path,
                        default=ROOT.parents[1]/"A30/a30-maxwell-material-manifold")
    parser.add_argument("--mode", choices=("bases", "all"), default="bases")
    args = parser.parse_args()
    begin = time.perf_counter()
    stage, a30 = args.stage.resolve(), args.a30_root.resolve()
    destination = stage/"T_PHYSICS_HANDOFF"
    source_revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    a30_revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=a30, text=True).strip()
    qualification = json.loads((stage/"qualification/SUMMARY.json").read_text())
    backgrounds = {}
    for background_id in ("weak", "ring"):
        source = a30/"inputs/a29_release/inputs/a25/background"/(background_id+".npz")
        backgrounds[background_id] = arrays(source, BACKGROUND_KEYS)
        path = destination/"backend_assets/background"/(background_id+".npz")
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(path, **backgrounds[background_id])
    records = []
    qualification_by_case = {row["case"]:row for row in qualification}
    for row in qualification:
        if not row["qualified"]:
            continue
        case = row["case"]
        selected = next(arm for arm in row["arms"] if arm["arm"] == row["selected_arm"])
        overhead = [{"arm": arm["arm"], "cost": arm["cost"], "wall_seconds": arm["wall_seconds"]}
                    for arm in row["arms"] if arm["arm"] != row["selected_arm"]]
        records.append(export_one(kind="qualified_base", case=case, candidate_id=case+"__base",
            state_path=stage/"qualification"/case/"qualified_base.npz",
            state_summary=row["selected_validation"], dedicated_cost=selected["cost"],
            wall=selected["wall_seconds"], shared=row["shared_setup"], overhead=overhead,
            source_revision=source_revision, a30_revision=a30_revision, stage=stage,
            destination=destination, background=backgrounds["weak" if "_weak_" in case else "ring"],
            lineage={"proposal_kind": "qualified_base_bounded_material_GN", "selected_arm": row["selected_arm"],
                     "base_role": "historical_lowest_residual_selection; default G base changed after failure analysis",
                     "selection_rule": row["selection"], "accepted_step_history": selected["path"],
                     "corrector": selected["diagnostics"]["config"],
                     "qualification_summary_ref": relative(stage/"qualification/SUMMARY.json", destination/"candidates"/(case+"__base"))}))
    minchange_file = stage/"qualification/MINCHANGE_SUMMARY.json"
    if minchange_file.exists():
        for row in json.loads(minchange_file.read_text()):
            case = row["case"]
            original = qualification_by_case[case]
            if row["selected_arm"] == "original_n3":
                selected = {"cost": {"new_full_state_attempts": 0, "complex_full_solver_RHS": 0,
                                     "RHS_by_role": {}, "action_counts": {}, "failed_actions": 0,
                                     "ledger_rows": [original["shared_setup"]["ledger_rows"][1]]*2},
                            "wall_seconds": None, "path": [], "diagnostics": {"config": None}}
            else:
                selected = next(a for a in original["arms"] if a["arm"]==row["selected_arm"])
            overhead = [{"arm":a["arm"],"cost":a["cost"],"wall_seconds":a["wall_seconds"]}
                        for a in original["arms"] if a["arm"]!=row["selected_arm"]]
            records.append(export_one(kind="qualified_minchange_base", case=case,
                candidate_id=case+"__minchange_base",
                state_path=stage/"qualification"/case/"qualified_minchange_base.npz",
                state_summary=row["validation"], dedicated_cost=selected["cost"],
                wall=selected["wall_seconds"],shared=original["shared_setup"],overhead=overhead,
                source_revision=source_revision,a30_revision=a30_revision,stage=stage,
                destination=destination,background=backgrounds["weak" if "_weak_" in case else "ring"],
                lineage={"proposal_kind":"qualified_minimum_change_base_selection",
                         "base_role":"DEFAULT_G_BASE", "selected_arm":row["selected_arm"],
                         "selection_rule":row["selection"],
                         "policy_revision_status":"changed after this exposed-development failure analysis; not an independent blind test",
                         "material_L2_change_from_original_n3":row["material_L2_change"],
                         "selection_new_physics_attempts":0,
                         "generation_code_revision":"a42845f",
                         "accepted_step_history":selected["path"],
                         "corrector":selected["diagnostics"]["config"]}))
    if args.mode == "all":
      for version, generation_revision in [(1,"c752aa6fedd4d7f12d93062857013b721412135f"),(2,"439773be7e663cca2cde3e33484e30fe0d849529"),(3,"a42845f")]:
        mechanism = stage/f"mechanism_v{version}/SUMMARY.json"
        rows = json.loads(mechanism.read_text()) if mechanism.exists() else []
        for row in rows:
            if not row["genuine_new_candidate"]:
                continue
            case, idx = row["case"], int(row["edit_id"])
            folder = stage/f"mechanism_v{version}"/case/f"edit_{idx:02}"
            state_path = folder/"new_valid_candidate.npz"
            if not state_path.exists():
                raise ValueError("completed genuine candidate summary lacks actual state")
            geometry = arrays(folder/"edit_geometry.npz", ("core", "guard"))
            base_name = "qualified_minchange_base.npz" if version==3 else "qualified_base.npz"
            base = arrays(stage/"qualification"/case/base_name, ("chi",))["chi"]
            endpoint = arrays(state_path, ("chi",))["chi"]
            distance = material_distances(base, endpoint, mass=1/576)
            if not distance["significant_phase_edit"] or distance["normalized_material_l2"] < .03:
                raise ValueError("genuine label contradicts saved material novelty")
            core = geometry["core"]
            core_new_air_cells = int(((abs(base.reshape(24,24))>.1) & core).sum())
            records.append(export_one(kind="protected_compensated_endpoint", case=case,
                candidate_id=case+f"__v{version}_edit_{idx:02}", state_path=state_path,
                state_summary=row["corrected_validation"], dedicated_cost=row["cost"],
                wall=row["wall_seconds"], shared=row["shared_base_cost"], overhead=None,
                source_revision=source_revision, a30_revision=a30_revision, stage=stage,
                destination=destination, background=backgrounds["weak" if "_weak_" in case else "ring"],
                parent_candidate_id=case+("__minchange_base" if version==3 else "__base"),
                core=core, guard=geometry["guard"], minimum_feature=row["edit_metadata"]["min_feature_radius_physical"],
                lineage={"proposal_kind": row["edit_metadata"]["kind"], "edit": row["edit_metadata"],
                         "mechanism_iteration":version,"generation_code_revision":generation_revision,
                         "extra_phase_protected_cells":row.get("extra_phase_protected_cells",0),
                         "fixed_legal_constraints_ref":relative(folder/"edit_geometry.npz",destination/"candidates"/(case+f"__v{version}_edit_{idx:02}")),
                         "compensation": {"raw_validation": row["raw_validation"],
                                          "compensated_validation": row["corrected_validation"],
                                          "finite_edited_local_score": row["compensation_score"],
                                          "gain_changed": False, "material_parameters_changed": True},
                         "material_distance_from_base": distance,
                         "core_nonair_to_air_cells": core_new_air_cells,
                         "core_cells_already_air_in_base": int(core.sum())-core_new_air_cells,
                         "genuine_new_candidate_as_recorded": True,
                         "invalid_to_valid_as_recorded": bool(row["invalid_to_valid"]),
                         "new_family_scope": "significant endpoint relative to qualified base; not a new connected feasible-component certificate",
                         "accepted_step_history": row["corrector"]["path"],
                         "corrector": row["corrector"]["diagnostics"]["config"]}))
    for item in records:
        item["base_role"] = ("DEFAULT_G_MINIMUM_CHANGE" if item["kind"]=="qualified_minchange_base"
                             else "HISTORICAL_LOWEST_RESIDUAL_SELECTION" if item["kind"]=="qualified_base"
                             else "COMPENSATED_ENDPOINT_FROM_DECLARED_ITERATION")
    preserve_actual_correction_constraints(records,destination)
    archive = archive_and_pairwise(records,destination)
    for phase in ["mechanism_v3","mechanism_v1"]:
        source = stage/phase/"SUMMARY.json"
        if source.exists():
            rows = json.loads(source.read_text())
            if rows:
                json_write(destination/"DEFAULT_CORRECTOR_CONFIG.json",rows[0]["corrector"]["diagnostics"]["config"])
                break
    json_write(destination/"BACKEND_ASSETS_PROVENANCE.json",{
        "backend_revision":a30_revision,
        "source_repository":"https://github.com/migodam/a30-maxwell-observable-material-manifold",
        "source_arrays":"inputs/a29_release/inputs/a25/background/{weak,ring}.npz",
        "copied_keys":list(BACKGROUND_KEYS),
        "not_copied":["truth", "object anatomy", "full teacher J", "step", "oracle", "locked data"],
        "known_background_prior_generation_cost":None,
        "historical_cost_status":"original upstream saved-physics setup; not re-counted as new T calls and not assumed to have been free",
        "Research_T_import_and_audit_cost":"shared setup receipts and the original production ledger, not zero-cost new physics",
        "new_physics_during_asset_copy":0,
        "consumer_reuse":"use original frozen A30 adapter; do not rebuild a different solver"})
    index = {"schema_version": "research-t-g/1.0", "producer": "Research T",
             "source_revision_read": source_revision, "backend_revision_read": a30_revision,
             "generation_code_revision": "c752aa6fedd4d7f12d93062857013b721412135f",
             "exported_at_utc": datetime.now(timezone.utc).isoformat(), "mode": args.mode,
             "qualified_base_count": sum(r["kind"]=="qualified_base" for r in records),
             "qualified_minchange_default_base_count": sum(r["kind"]=="qualified_minchange_base" for r in records),
             "default_bases_for_G":{r["case"]:r["candidate_id"] for r in records if r["kind"]=="qualified_minchange_base"},
             "compensated_endpoint_count": sum(r["kind"]=="protected_compensated_endpoint" for r in records),
             "raw_real_bundle_count":len(records),
             "deduplicated_structural_archive_count":archive["retained_structural_representatives"],
             "per_case_archive_counts":archive["per_case_retained_counts"],
             "archive_ref":"MATERIAL_ARCHIVE.json",
             "candidate_only": True, "probability_weight": None,
             "unexplored_structure_status": "UNQUANTIFIED; no complete coverage claim",
             "outside_model_support_mass": None, "new_physics_attempts": 0,
             "new_full_solver_RHS": 0,
             "producer_ledger_ref": "../COST_LEDGER.jsonl",
             "all_qualification_attempts_source":"../qualification/SUMMARY.json",
             "all_mechanism_attempts_sources":[f"../mechanism_v{v}/SUMMARY.json" for v in [1,2,3] if (stage/f"mechanism_v{v}/SUMMARY.json").exists()],
             "consumer_budget_rule": "G must provide a separate actual authorization and book; no use or reset of T allowance",
             "candidate_inventory": records}
    json_write(destination/"CANDIDATE_INDEX.json", index)
    receipt = {"utc": datetime.now(timezone.utc).isoformat(), "mode": args.mode,
               "source_revision_read": source_revision, "new_full_state_attempts": 0,
               "complex_full_solver_RHS": 0, "solver_imports": 0, "truth_data_reads": 0,
               "G_workspace_reads_or_writes": 0,
               "cached_matrix_columns_replayed_G_and_S_each": 8*len(records),
               "export_execution_wall_seconds": time.perf_counter()-begin,
               "qualified_bases": index["qualified_base_count"],
               "qualified_minchange_bases":index["qualified_minchange_default_base_count"],
               "compensated_endpoints": index["compensated_endpoint_count"],
               "raw_real_bundle_count":len(records),
               "archive_representatives":archive["retained_structural_representatives"],
               "per_case_archive_counts":archive["per_case_retained_counts"],
               "all_bundle_semantics_pass": all(r["validation"]["valid_bundle_semantics"] for r in records)}
    with (destination/"EXPORT_RECEIPTS.jsonl").open("a") as stream:
        stream.write(json.dumps(receipt, allow_nan=False)+"\n")
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
