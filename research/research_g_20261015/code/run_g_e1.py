"""Budgeted G search with the frozen historical A30 physics/corrector.

No solver is imported or constructed on import. Only ``--execute`` runs paid
physics. This online runner never imports a scorer or opens truth-labelled
arrays. Search density and discovery counts have no posterior interpretation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from g_archive import Archive, Candidate, material_distance, structure_features, structural_distinct
from g_budget import GBudget, BudgetStop, dump


REPO = Path(__file__).resolve().parents[3]
PHASE = Path(__file__).resolve().parents[1]
HISTORICAL = REPO / "research/research_a_20261015/results/E2"
A30 = REPO.parents[1] / "A30/a30-maxwell-material-manifold"
A_COMMIT = "f1b674eb4893e39126701ed48f3d3dcef3098c3c"
A30_COMMIT = "717d3d26a290248bd2b9aeb55f3d03694b407b84"
CASES = ("p250069_weak_a1_snr20", "p250073_weak_a1_snr20", "p250077_weak_a1_snr20")
METHODS = ("multistart", "static_scale", "QD", "bandit")
ROUNDS = ("first", "quality")
OPERATORS = ("global", "local_remove", "local_add", "distributed")
THRESHOLD = 23.785449004057032
SEED = 261015561


def project(field):
    """The inherited global-start material box, never a null-space claim."""
    value = np.asarray(field, np.complex128).reshape(576)
    return np.clip(value.real, -.5, 4.) + 1j*np.clip(value.imag, 0., .5)


def legal(field):
    x = np.asarray(field)
    return bool(x.shape == (576,) and x.dtype == np.complex128 and np.isfinite(x).all()
                and np.all((x.real >= -.5) & (x.real <= 4.) & (x.imag >= 0.) & (x.imag <= .5)))


def pack(data):
    """A30 source-major Re32 / Im32 packing; gain fixed to one."""
    y = np.asarray(data)
    if y.shape != (32, 8):
        raise ValueError("G_A30_32_RECEIVERS_8_SOURCES_REQUIRED")
    return np.concatenate((y.real.T, y.imag.T), axis=1).reshape(-1)


def residual(state, observation, white):
    return float(np.linalg.norm(pack(state.data-observation)*white))


def relative(first, second):
    return float(np.linalg.norm(np.asarray(first)-second)/max(np.linalg.norm(second), 1e-300))


def bounds_fraction(field):
    x = np.asarray(field)
    return float(np.mean((x.real <= -.5+1e-10) | (x.real >= 4.-1e-10)
                         | (x.imag <= 1e-10) | (x.imag >= .5-1e-10)))


def state_diagnostics(state, observation, white):
    data_norm = residual(state, observation, white)
    constitutive = relative(state.current, state.chi[:, None]*state.field)
    good_solver = bool(np.isfinite(state.residual) and state.residual <= 1e-10
                       and np.isfinite(constitutive) and constitutive <= 1e-10)
    return dict(full_residual=data_norm, solver_backward_residual=float(state.residual),
                constitutive_relative_defect=constitutive, material_legal=legal(state.chi),
                solver_validated=good_solver,
                compatible=bool(data_norm <= THRESHOLD and good_solver and legal(state.chi)))


def state_arrays(arrays, prefix, state):
    for name, attribute in (("chi", "chi"), ("prediction", "data"), ("current", "current"), ("field", "field")):
        arrays[prefix+"__"+name] = np.array(getattr(state, attribute), copy=True)


def assert_frozen_sources(a30_root):
    # Git commit identity/source equality are provenance, not digest checking.
    actual = subprocess.check_output(["git", "-C", str(a30_root), "rev-parse", "HEAD"], text=True).strip()
    if actual != A30_COMMIT:
        raise ValueError("G_A30_FROZEN_COMMIT_MISMATCH")
    files = ("src/research_a/e2.py",)
    for path in files:
        frozen = subprocess.check_output(["git", "-C", str(REPO), "show", A_COMMIT+":"+path])
        if frozen != (REPO/path).read_bytes():
            raise ValueError("G_A_CORRECTOR_SOURCE_DIFFERS_FROM_FROZEN_COMMIT:"+path)
    changed = subprocess.check_output(["git", "-C", str(a30_root), "diff", "HEAD", "--", "src/a30/physics.py", "src/a30/representations.py"], text=True)
    if changed:
        raise ValueError("G_A30_WORKTREE_PHYSICS_MODIFIED")


def load_whitelisted_inputs(a30_root):
    """Only four explicit A array keys per case and A30 sigma are accessed."""
    release = a30_root / "inputs/a29_release"
    metadata = {r["case_id"]: r for r in csv.DictReader((release/"results/a29/COVERAGE_CASES.csv").open())}
    base_log = {r["case_id"]: r for r in json.loads((HISTORICAL/"SEARCH_HISTORY.json").read_text()) if r["arm"] == "base"}
    result = {}
    with np.load(HISTORICAL/"MATERIAL_AND_PREDICTION_ARRAYS.npz", allow_pickle=False) as archive:
        for case in CASES:
            keys = [case+"__"+suffix for suffix in ("base", "base_prediction", "observation", "background")]
            if any(key not in archive.files for key in keys):
                raise ValueError("G_REQUIRED_HISTORICAL_BASE_ARRAY_MISSING:"+case)
            # Listing ZIP member names does not deserialize any other array.
            base, prediction, observation, background = [archive[key].copy() for key in keys]
            with np.load(release/"results/a29/pools"/(case+".npz"), allow_pickle=False) as pool:
                sigma = float(pool["sigma"])
            if not sigma > 0 or not np.isfinite(sigma):
                raise ValueError("G_INVALID_FROZEN_NOISE")
            white = np.sqrt(2.)/sigma
            if not legal(base) or prediction.shape != (32, 8) or observation.shape != (32, 8):
                raise ValueError("G_HISTORICAL_BASE_LAYOUT_OR_MATERIAL_INVALID")
            if not all(x.dtype == np.complex128 for x in (base, prediction, observation, background)):
                raise ValueError("G_COMPLEX128_REQUIRED")
            saved_norm = float(np.linalg.norm(pack(prediction-observation)*white))
            if not base_log[case]["compatible"] or saved_norm > THRESHOLD or abs(saved_norm-base_log[case]["residual"]) > 1e-8:
                raise ValueError("G_HISTORICAL_BASE_NOT_QUALIFIED:"+case)
            result[case] = dict(base=base, prediction=prediction, observation=observation,
                                background=background.reshape(576), sigma=sigma, white=white,
                                background_id=metadata[case]["background_id"], saved_residual=saved_norm,
                                read_array_keys=keys, sigma_source=str(release/"results/a29/pools"/(case+".npz")))
    return result


def choose_archive_parent(archive, uses, base):
    candidates = list(archive.active_archive)
    if not candidates:
        return base
    def key(candidate):
        others = [old for old in candidates if old.candidate_id != candidate.candidate_id]
        novelty = min((archive.pair(candidate, old)["dM"]+archive.pair(candidate, old)["dS"] for old in others), default=0.)
        return (uses.get(candidate.candidate_id, 0), -novelty, candidate.candidate_id)
    return min(candidates, key=key)


def propose(method, trajectory, case_index, archive, base, background, basis, rng, parent_uses, alpha, beta):
    if trajectory == 0 or method == "multistart":
        operator = "global"
    elif method == "static_scale":
        operator = "local_remove" if trajectory == 1 else "local_add"
    elif method == "QD":
        operator = OPERATORS[int(rng.integers(len(OPERATORS)))]
    elif method == "bandit":
        operator = OPERATORS[int(np.argmax(rng.beta(alpha, beta)))]
    else:
        raise ValueError("G_UNKNOWN_METHOD")
    parent = base if method in ("multistart", "static_scale") or operator == "global" else choose_archive_parent(archive, parent_uses, base)
    parent_uses[parent.candidate_id] = parent_uses.get(parent.candidate_id, 0)+1
    x = np.asarray(parent.field).reshape(576)
    details = dict(operator=operator, parent_id=parent.candidate_id,
                   forced_defensive_global=trajectory == 0, topology_label=None,
                   proposal_interpretation="SEARCH_PROPOSAL_NOT_POSTERIOR")
    if operator in ("global", "distributed"):
        coefficients = rng.standard_normal(basis.shape[1])
        update = basis@coefficients
        update *= .8/max(float(np.max(np.abs(update))), 1e-300)
        unprojected = x+update
        proposal = project(unprojected)
        details.update(coefficients=coefficients.tolist(), peak_unprojected_update=.8,
                       operator_note="base-centered correlated global restart, NOT uniform full-field prior" if operator == "global" else "archive-parent multiscale Haar material perturbation",
                       projection_displacement=float(np.linalg.norm(proposal-unprojected)))
    else:
        contrast = np.abs(x-background).reshape(24, 24)
        sites = np.argwhere(contrast > .18)
        if operator == "local_remove" and len(sites):
            center = sites[int(rng.integers(len(sites)))].astype(float)+rng.uniform(-.5, .5, 2)
            site_scope = "current_material_contrast_only"
        else:
            center = rng.uniform(0., 23., 2)
            site_scope = "entire_grid_no_truth"
        center = np.clip(center, 0., 23.)
        if method == "static_scale":
            scale = (1.5, 3., 6.)[(case_index+trajectory-1) % 3]
            amplitude = (.45, .7, .95)[(case_index+trajectory-1) % 3]
        else:
            scale = float(np.exp(rng.uniform(np.log(1.5), np.log(6.))))
            amplitude = float(rng.uniform(.4, .95))
        yy, xx = np.mgrid[:24, :24]
        mask = np.exp(-.5*((yy-center[0])**2+(xx-center[1])**2)/scale**2).reshape(576)
        target = background if operator == "local_remove" else background+complex(1.2, .12)
        unprojected = x+amplitude*mask*(target-x)
        proposal = project(unprojected)
        details.update(center_pixels=center.tolist(), gaussian_scale_pixels=scale,
                       amplitude=amplitude, signed_amplitude=-amplitude if operator == "local_remove" else amplitude,
                       site_scope=site_scope, known_background_used=True, nominal_added_contrast=[1.2, .12],
                       projection_displacement=float(np.linalg.norm(proposal-unprojected)))
    details["raw_bound_fraction"] = bounds_fraction(proposal)
    return proposal, details, parent


def archive_payload(archive, locations):
    return dict(archive_summary=archive.summary(),
                success_ledger=[dict(candidate_id=c.candidate_id, **locations[c.candidate_id],
                                     full_residual=c.full_residual, is_base=c.is_base,
                                     validation_receipt=c.validation_receipt) for c in archive.success_ledger],
                active_archive_ids=[c.candidate_id for c in archive.active_archive],
                material_representative_ids=[c.candidate_id for c in archive.material_representatives],
                structural_representative_ids=[c.candidate_id for c in archive.structural_representatives],
                frontier_ids=[c.candidate_id for c in archive.frontier],
                failures=list(archive.failure_ledger),
                top1_id=archive.top1().candidate_id if archive.top1() else None,
                topk_ids={str(k): [c.candidate_id for c in archive.topk(k)] for k in (2, 4, 8)},
                archive_events=[asdict(e) for e in archive.event_ledger])


def trajectory_run(model, common, basis, base, archive, method, trajectory, case_index,
                   round_id, case, book, output, locations, rng, parent_uses, alpha, beta,
                   gn_proposal):
    name = f"{round_id}__{case}__{method}__t{trajectory}"
    path = output/"arrays"/(name+".npz")
    relative_path = str(path.relative_to(PHASE))
    arrays = {"__observation": common["observation"], "__background": common["background"]}
    before = book.totals(); started = time.perf_counter()
    row = dict(round_id=round_id, case_id=case, method=method, trajectory=trajectory,
               endpoint_id=name, array_path=relative_path, array_prefix="validation", raw_array_prefix="raw",
               status="FAILED", steps=[], corrector_qp_receipts=[], posterior_weight=None,
               represented_mass=None, unmodeled_mass=None)
    parent = None
    try:
        raw, action, parent = propose(method, trajectory, case_index, archive, base,
                                      common["background"], basis, rng, parent_uses, alpha, beta)
        row["action"] = action
        arrays["__raw_proposal"] = raw.copy()
        current = None; trust = .35
        for step in range(3):
            prefix = ("raw", "correction1", "correction2")[step]
            state = model.full_state(raw, purpose=name+":"+prefix)
            diagnostics = state_diagnostics(state, common["observation"], common["white"])
            state_arrays(arrays, prefix, state)
            improved = current is None or diagnostics["full_residual"] < residual(current, common["observation"], common["white"])
            if improved:
                current = state
            else:
                trust *= .3
            row["steps"].append(dict(step=step, array_prefix=prefix, nonlinear_step_accepted=improved,
                                     trust_for_next_step=trust, independent_validation=False, **diagnostics))
            if step < 2:
                ledger_start = len(book.rows())
                raw = gn_proposal(model, current, basis, common["observation"], common["white"], trust=trust)
                receipt = [item for item in book.rows()[ledger_start:] if item.get("event") == "MATERIAL_QP"]
                displacement = raw-current.chi
                max_component = float(max(np.max(np.abs(displacement.real)), np.max(np.abs(displacement.imag))))
                row["corrector_qp_receipts"].append(dict(source_step=step, actual_source_receipts=receipt,
                                                        trust=trust, max_component_step=max_component,
                                                        trust_region_exit=bool(max_component > trust+1e-8)))
        # Independently paid fresh endpoint validation; no LU/cache replay.
        endpoint = model.full_state(current.chi.copy(), purpose=name+":fresh_independent_validation")
        state_arrays(arrays, "validation", endpoint)
        diagnostics = state_diagnostics(endpoint, common["observation"], common["white"])
        replay_error = relative(endpoint.data, current.data)
        diagnostics["independent_prediction_replay_error"] = replay_error
        validated = bool(diagnostics["compatible"] and replay_error <= 1e-10)
        receipt_id = "G-independent-full-validation:"+name
        book.append(dict(event="G_ENDPOINT_VALIDATION", candidate_id=name, receipt_id=receipt_id,
                         validation_passed=validated, **diagnostics))
        candidate = Candidate(name, endpoint.chi, endpoint.data, diagnostics["full_residual"],
                              validated, diagnostics["solver_validated"] and replay_error <= 1e-10,
                              receipt_id, is_base=False)
        old_success = list(archive.success_ledger)
        new_structure = bool(validated and all(archive.pair(candidate, old)["structural_distinct"] for old in old_success))
        new_material = bool(validated and all(archive.pair(candidate, old)["material_distinct"] for old in old_success))
        raw_feature = structure_features(arrays["__raw_proposal"], common["background"])
        endpoint_feature = structure_features(endpoint.chi, common["background"])
        parent_feature = structure_features(parent.field, common["background"])
        raw_structural_edit = structural_distinct(raw_feature, parent_feature)
        retained_structural_edit = structural_distinct(endpoint_feature, parent_feature)
        raw_prefix_distinct = all(structural_distinct(raw_feature, archive.features(old)) for old in old_success)
        edit = arrays["__raw_proposal"]-np.asarray(parent.field).reshape(576)
        corrected_edit = endpoint.chi-np.asarray(parent.field).reshape(576)
        denominator = float(np.vdot(edit, edit).real)
        retention = float(np.vdot(corrected_edit, edit).real/denominator) if denominator > 1e-24 else None
        offer = archive.offer(candidate)
        locations[name] = dict(array_path=relative_path, array_prefix="validation")
        if method == "bandit":
            index = OPERATORS.index(action["operator"])
            alpha[index] += int(new_structure); beta[index] += int(not new_structure)
        row.update(status="COMPLETE", valid=validated, full_residual=diagnostics["full_residual"],
                   validation=diagnostics, validation_receipt=receipt_id, offer_status=offer.status,
                   replaced_id=offer.replaced_id, new_structural_world=new_structure,
                   new_material_world=new_material, combined_new_world_reward=offer.new_world_reward,
                   retention=retention, raw_structural_edit=raw_structural_edit,
                   retained_structural_edit=retained_structural_edit,
                   correction_collapse_to_prefix=bool(validated and raw_prefix_distinct and not new_structure),
                   raw_material_delta=material_distance(arrays["__raw_proposal"], parent.field),
                   endpoint_material_delta=material_distance(endpoint.chi, parent.field),
                   correction_displacement=material_distance(endpoint.chi, arrays["__raw_proposal"]),
                   bounds_fraction=bounds_fraction(endpoint.chi),
                   bandit_alpha=alpha.tolist(), bandit_beta=beta.tolist(),
                   archive_after=archive.summary())
        if not validated:
            book.failure("G_ENDPOINT_INVALID", endpoint_id=name, case_id=case, method=method,
                         round_id=round_id, diagnostics=diagnostics)
    except BaseException as exc:
        row.update(error_type=type(exc).__name__, error=str(exc), valid=False)
        book.failure("G_TRAJECTORY_FAILURE", endpoint_id=name, error=repr(exc), steps_completed=len(row["steps"]))
        if method == "bandit" and "action" in row:
            # A failed physical trajectory has zero structural-discovery reward.
            beta[OPERATORS.index(row["action"]["operator"])] += 1
        if isinstance(exc, (BudgetStop, KeyboardInterrupt, SystemExit)):
            raise
        # No retries; failed paid attempts consume their actual budget.
    finally:
        # A30 transfer keys retain full LU states. No transfer is retained across
        # trajectories; within-trajectory legitimate same-state reuse is kept.
        model._transfers.clear()
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, **arrays)
        after = book.totals()
        row["cost_prefix"] = after
        row["cost_delta"] = {k: after[k]-before[k] for k in ("states", "rhs")}
        row["wall_seconds"] = time.perf_counter()-started
        dump(output/"trajectories"/(name+".json"), row)
    return row


def run(a30_root, job_name):
    protocol_path = PHASE/"E1_PROTOCOL_FROZEN.json"
    protocol = json.loads(protocol_path.read_text())
    if protocol.get("execution_authorized") is not True or protocol.get("backend_ready") is not True:
        raise ValueError("G_FROZEN_PHYSICAL_PROTOCOL_REQUIRED")
    output = PHASE/"results/E1_LEGACY"
    if (output/"RESULTS.json").exists():
        raise ValueError("G_COMPLETED_OR_FAILED_CAMPAIGN_ALREADY_EXISTS_NO_SILENT_RERUN")
    output.mkdir(parents=True, exist_ok=True)
    assert_frozen_sources(a30_root)
    common = load_whitelisted_inputs(a30_root)
    if str(REPO/"src") not in sys.path:
        sys.path.insert(0, str(REPO/"src"))
    from research_a.e2 import _modules, _bases, _gn_proposal
    book = GBudget(PHASE)
    summary = dict(status="FAILED", target="LEGACY_FROZEN_A30_24x24_scalar_TM",
                   a_commit=A_COMMIT, a30_commit=A30_COMMIT, cases=list(CASES), rounds=list(ROUNDS),
                   threshold=THRESHOLD, gain="FIXED_ONE", noise="sqrt(2)/original_sigma",
                   exposure="HISTORICALLY_EXPOSED_DEVELOPMENT", methods=[], bases=[],
                   T_guided="WAITING_FOR_T_NOT_IMPLEMENTED", posterior_weight=None,
                   represented_mass=None, unmodeled_mass=None,
                   primary_comparison="matched state attempts; actual RHS and wall reported separately",
                   max_expected_states=291, max_expected_RHS=6936, failures=[])
    models, bases, locations = {}, {}, {}
    try:
        with book.job(job_name, stage="E1"):
            physics, representations = _modules(a30_root)
            basis = np.asarray(_bases(representations)[1], dtype=np.complex128)
            if basis.shape != (576, 32):
                raise ValueError("G_INHERITED_HAAR_32REAL_CHART_REQUIRED")
            # Exactly three NEW full states establish independent legal bases.
            for case in CASES:
                data = common[case]
                model = physics.FrozenMaxwellAdapter.from_background_npz(
                    a30_root/"inputs/a29_release/inputs/a25/background"/(data["background_id"]+".npz"),
                    book, allowed_cache_paths=(), device="cpu", backend="scipy", purpose_prefix="Research_G_E1")
                models[case] = model
                if not np.array_equal(model.chi_b, data["background"]):
                    raise ValueError("G_KNOWN_BACKGROUND_IDENTITY_MISMATCH")
                state = model.full_state(data["base"], purpose=case+":independent_base_validation")
                diagnostics = state_diagnostics(state, data["observation"], data["white"])
                diagnostics["historical_prediction_relative_error"] = relative(state.data, data["prediction"])
                if not diagnostics["compatible"] or diagnostics["historical_prediction_relative_error"] > 1e-10:
                    raise ValueError("G_FRESH_BASE_REPLAY_OR_FEASIBILITY_FAILURE:"+case)
                candidate_id = "base__"+case
                receipt_id = "G-full-base-validation:"+case
                book.append(dict(event="G_BASE_VALIDATION", candidate_id=candidate_id, receipt_id=receipt_id, **diagnostics))
                candidate = Candidate(candidate_id, state.chi, state.data, diagnostics["full_residual"], True, True, receipt_id, is_base=True)
                bases[case] = candidate
                arrays = {"__observation": data["observation"], "__background": data["background"]}
                state_arrays(arrays, "base", state)
                path = output/"bases"/(case+".npz"); path.parent.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(path, **arrays)
                locations[candidate_id] = dict(array_path=str(path.relative_to(PHASE)), array_prefix="base")
                summary["bases"].append(dict(case_id=case, candidate_id=candidate_id, **locations[candidate_id],
                                              validation=diagnostics, validation_receipt=receipt_id,
                                              sigma=data["sigma"], historical_read_keys=data["read_array_keys"]))
                del state
            for round_index, round_id in enumerate(ROUNDS):
                for case_index, case in enumerate(CASES):
                    method_order = np.random.default_rng(SEED+100000+1000*round_index+case_index).permutation(METHODS).tolist()
                    book.append(dict(event="G_METHOD_ORDER", round_id=round_id, case_id=case, order=method_order))
                    for method in method_order:
                        archive = Archive(common[case]["background"], data_threshold=THRESHOLD, replacement_mode=round_id)
                        archive.offer(bases[case])
                        method_result = dict(round_id=round_id, case_id=case, method=method,
                                             replacement_mode=round_id, episodes=[], method_order=method_order)
                        summary["methods"].append(method_result)
                        parent_uses = {}; alpha = np.full(len(OPERATORS), .5); beta = np.full(len(OPERATORS), .5)
                        # No round term: paired raw RNG streams across archive repair.
                        rng = np.random.default_rng(SEED+1000*case_index+100*METHODS.index(method))
                        with book.arm(round_id, case, method, states=12, rhs=288):
                            for trajectory in range(3):
                                episode = trajectory_run(models[case], common[case], basis, bases[case], archive, method,
                                                         trajectory, case_index, round_id, case, book, output,
                                                         locations, rng, parent_uses, alpha, beta, _gn_proposal)
                                method_result["episodes"].append(episode)
                                method_result.update(archive_payload(archive, locations))
                                summary["resources"] = book.receipt()
                                dump(output/"PARTIAL_RESULTS.json", summary)
                        method_result["local_costs"] = book.totals(f"{round_id}/{case}/{method}")
                        print(json.dumps(dict(round=round_id, case=case, method=method,
                                              archive=archive.summary(), cumulative=book.totals())), flush=True)
            summary["status"] = "COMPLETE"
    except BaseException as exc:
        summary["failures"].append(dict(type=type(exc).__name__, message=str(exc)))
        raise
    finally:
        for model in models.values():
            model._transfers.clear()
        summary["resources"] = book.receipt()
        dump(output/"RESULTS.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--a30-root", type=Path, default=A30)
    parser.add_argument("--job", default="g-e1-legacy-search-v1")
    arguments = parser.parse_args()
    if not arguments.execute:
        raise SystemExit("No physics executed. Root review followed by --execute is required.")
    result = run(arguments.a30_root.resolve(), arguments.job)
    print(json.dumps(dict(status=result["status"], resources=result["resources"]), indent=2))
