"""Four real Research A handoffs under the frozen Research B E2 authorization.

No paid work happens on import or without --execute. The orchestrator must
review this runner before invoking it. Source A packages remain read-only.
Candidate-only acquisition never receives truth, future data or A probabilities.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from src.b1_joint.material import (BOUNDS, decode, geometry, event_ids,
                                   prior_logpdf, prior_sample)
from src.b1_joint.reference import Reference
from src.b1_joint.response import Response
from sensing import (Action, field_l1, normalized, resample_indices,
                     conditional_log_increment, expected_utility)
from active_experiment import candidate_scores
from e0_candidate_controls import profile
from ab_adapter import make_recommendation, updated_data_history, export_bundle
from cost import Cost, BudgetExceeded, ROOT, dump

FREQUENCIES = (1.2, 1.6, 2., 2.4, 2.8, 3.2)
OLD = 1.6
SEED_NOISE = 26101591
SEED_STARTS = 26101592
REFERENCE_ORDER = 13
REFERENCE_PARTICLES = 2048
REFERENCE_DECISIONS = 96
REFERENCE_OUTER = 128
MAX_NFEV = 80
ALPHA_COMPATIBILITY = .05  # Inherit A's chi2(.95, dimension) convention.
CANDIDATE_CATEGORY = "all_future_candidate_fields"
TRUTH_CATEGORY = "independent_truth_all6"
PROTECTED_CATEGORY = "protected_epsilon_corrector_all27_twofreq"
PRIOR_CATEGORY = "three_unlabelled_prior_multistarts_each4parents_twofreq"
REFERENCE_CATEGORY = "reference_selected_endpoint_checks"


def append(path, row):
    with Path(path).open("a") as stream:
        stream.write(json.dumps(row, default=lambda v: v.tolist() if hasattr(v, "tolist")
                                else str(v), allow_nan=False) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def key_for(x, k):
    return json.dumps([np.asarray(x, float).tolist(), float(k)], separators=(",", ":"))


def real_field(value):
    field = value[0] if isinstance(value, tuple) else value
    field = np.asarray(field)
    if np.iscomplexobj(field) or field.shape != (32,) or not np.isfinite(field).all():
        raise ValueError("Expected finite full-Mie Re16/Im16 field, shape(32,)")
    return field.astype(np.float64)


def epsilon_on_grid(x, edges):
    material = decode(np.asarray(x, float))
    mid = (edges[:-1] + edges[1:])/2
    return material["eps_layers"][np.searchsorted(material["radii"], mid)]


def load_bundle(parent):
    source = ROOT / "A_IMPORT" / (parent + "_B_menu_v3")
    package = json.loads((source / "MaterialHypothesisPackage.json").read_text())
    contract = json.loads((source / "contract.json").read_text())
    assets = {p.stem: np.load(p, allow_pickle=False) for p in source.glob("*.npy")}
    for p in source.glob("*.json"):
        if p.name not in ("MaterialHypothesisPackage.json", "contract.json",
                          "AcquisitionRecommendation.json", "updated_data_history.json"):
            assets[p.stem] = json.loads(p.read_text())
    x = np.asarray(assets["p2_latent_coordinates"], float)
    f0 = np.asarray(assets["unscaled_old_fields"], float)
    y = np.asarray(assets["observations"], float)
    cov = np.asarray(assets["covariance"], float)
    if (x.shape != (len(package["candidates"]), 4) or f0.shape != (len(x), 32)
            or y.shape != (32,) or not np.isfinite(f0).all()
            or not np.isfinite(prior_logpdf(x)).all()):
        raise ValueError("Malformed real A candidate/observation arrays")
    if (cov.shape != (32, 32) or not np.allclose(cov, np.diag(np.diag(cov)), atol=1e-15)
            or np.any(np.diag(cov) <= 0)):
        raise ValueError("This runner requires the explicitly diagonal 32-channel noise contract")
    if package["weights"]["kind"] != "none" or contract["allow_bayesian_risk"]:
        raise ValueError("Do not promote A candidates into posterior weights")
    if not contract.get("original_observation_exactly_preserved"):
        raise ValueError("A original observation identity is not verified")
    physics = contract["input_context"]["expected_physics"]
    if (physics["backend_commit"] != "7ed4faef9bc3dd214d4de689c2203cc383f7495f"
            or physics["solver_id"] != "B1-vector-Mie-tangential-boundary-matching"
            or physics["material_law_id"] != "passive-isotropic-relative-permittivity-mu1"
            or physics["time_convention"] != "exp(-i omega t)"):
        raise ValueError("Source A solver/material identity differs from the paid B adapter")
    actions = package["available_actions"]
    if (len(actions) != 6 or {float(a["frequency_kR"]) for a in actions} != set(FREQUENCIES)
            or not contract.get("B_menu_authorization_id")):
        raise ValueError("Verified independent B six-frequency menu extension required")
    edges = np.asarray(assets["radial_interval_edges"], float)
    for i, latent in enumerate(x):
        if not np.allclose(epsilon_on_grid(latent, edges), assets[f"material_{i}"], atol=1e-12, rtol=1e-12):
            raise ValueError("Latent no longer reconstructs A's original exact union-grid material")
    base = [i for i, c in enumerate(package["candidates"]) if c["base_solution"]]
    if len(base) != 1 or not np.array_equal(assets[f"material_{base[0]}"], assets["base_material"]):
        raise ValueError("A base-only control must retain the exact original base material")
    source_gains = np.asarray(assets["source_fitted_gains"])
    for i in range(len(x)):
        if not np.allclose(f0[i]*source_gains[i], assets[f"prediction_{i}"], atol=1e-13, rtol=1e-13):
            raise ValueError("A source gain was not correctly removed from the full field")
    return {"package": package, "contract": contract, "assets": assets,
            "source": source, "parent": parent, "base_index": base[0]}


class PaidFields:
    """Append-only local reservation cap on top of the global Cost callback.

    A callback attempt is reserved locally before even a cache lookup. Exact
    restored/seeded fields do not invoke the callback. Failed requests remain
    charged locally, including invalid output and retry attempts.
    """
    def __init__(self, cost, callback, output, frozen):
        self.cost, self.callback = cost, callback
        self.path = output / "PHYSICS_REQUESTS.jsonl"
        self.cap = min(240, int(frozen["new_attempt_subcap"]))
        self.caps = dict(frozen["allocation"])
        self.counts = {k: 0 for k in self.caps}
        self.attempts = 0
        self.cache, self.failed = {}, {}
        if self.path.exists():
            for line in self.path.read_text().splitlines():
                row = json.loads(line)
                if row["event"] == "BEGIN":
                    self.attempts += 1
                    self.counts[row["category"]] += 1
                    self.failed[row["key"]] = "Unfinished previous request"
                elif row["event"] == "END" and row["status"] == "PASS":
                    self.cache[row["key"]] = (np.asarray(row["field"]), row.get("info", {}))
                    self.failed.pop(row["key"], None)
                elif row["event"] == "END":
                    self.failed[row["key"]] = row.get("error", "Previous request failed")
        if self.attempts > self.cap or any(self.counts[k] > self.caps[k] for k in self.caps):
            raise BudgetExceeded("Existing E3 local ledger already exceeds its frozen subcap")

    def seed_original(self, latent, field, provenance):
        key = key_for(latent, OLD)
        field = real_field(field)
        if key in self.cache and not np.allclose(self.cache[key][0], field, atol=1e-12, rtol=1e-12):
            raise ValueError("Saved original A field conflicts with a restored full field")
        self.cache.setdefault(key, (field, {"source": provenance, "original_A_paid_field": True}))

    def call(self, latent, frequency, category, label):
        self.cost.check_time()
        key = key_for(latent, frequency)
        if key in self.cache:
            return self.cache[key][0].copy(), self.cache[key][1], "REUSED_PAID_FULL_FIELD"
        if key in self.failed:
            return None, {}, "PREVIOUS_FAILURE_PRESERVED: " + self.failed[key]
        if (category not in self.caps or self.attempts >= self.cap
                or self.counts[category] >= self.caps[category]):
            raise BudgetExceeded("E3 local total/category reservation cap reached before callback")
        self.attempts += 1
        self.counts[category] += 1
        header = dict(key=key, category=category, request=self.attempts,
                      latent=np.asarray(latent, float).tolist(), frequency=float(frequency), label=label)
        append(self.path, dict(event="BEGIN", **header))
        started = time.perf_counter()
        try:
            value = self.callback(np.asarray(latent, float), float(frequency), kind="E3_IMPORT_" + label)
            field = real_field(value)
            info = value[3] if isinstance(value, tuple) and len(value) > 3 else {}
        except BaseException as error:
            reason = type(error).__name__ + ": " + str(error)
            append(self.path, dict(event="END", status="FAIL", error=reason,
                                   wall_seconds=time.perf_counter()-started, **header))
            self.failed[key] = reason
            if isinstance(error, (BudgetExceeded, KeyboardInterrupt, SystemExit)):
                raise
            return None, {}, reason
        append(self.path, dict(event="END", status="PASS", field=field, info=info,
                               wall_seconds=time.perf_counter()-started, **header))
        self.cache[key] = (field.copy(), info)
        return field, info, "PAID_FULL_CALLBACK_PASS"

    def receipt(self):
        return dict(callback_reservations=self.attempts, callback_reservation_cap=self.cap,
                    by_category=self.counts, category_caps=self.caps,
                    failed_or_unfinished_keys=len(self.failed),
                    original_A1_6_fields_reused_without_callback=True,
                    total_global_ledger=str(self.cost.physics), local_ledger=str(self.path),
                    cache_hits_are_not_new_Maxwell_attempts=True)


class OfflineSimulator:
    """The only reader of A's confirmation truth and hidden persistent gain."""
    def __init__(self, parent_ids, output):
        source = REPO.parent / "b1-research-a/research/research_a_20261015/results"
        truth = json.loads((source / "offline_truth/CONFIRM_TRUTH.json").read_text())
        online = json.loads((source / "CONFIRM_OBSERVATIONS.json").read_text())
        self._truth = {r["id"]: r for r in truth if r["id"] in parent_ids}
        self._online = {r["id"]: r for r in online if r["id"] in parent_ids}
        if set(self._truth) != set(parent_ids) or set(self._online) != set(parent_ids):
            raise ValueError("Missing real A confirmation parent identity")
        self._fields, self._future, self._errors = {}, {}, []
        self.output = output

    def verify_original(self, parent, observation):
        if not np.array_equal(observation, np.asarray(self._online[parent]["y"])):
            raise ValueError("Imported A noisy observation differs from its original online record")

    def prepare(self, paid, sigma):
        for parent in sorted(self._truth):
            item = self._truth[parent]
            index = int(parent.split("_")[-1])
            for j, k in enumerate(FREQUENCIES):
                field, info, status = paid.call(item["latent"], k, TRUTH_CATEGORY, parent+"_TRUTH")
                if field is None:
                    self._errors.append(dict(parent=parent, frequency=k, status=status))
                    continue
                self._fields[(parent, k)] = field
                rng = np.random.default_rng(np.random.SeedSequence([SEED_NOISE, index, j]))
                self._future[(parent, k)] = float(item["gain"])*field + sigma[j]*rng.normal(size=32)
                if k == OLD:
                    original = np.asarray(item["measurement"], float)
                    discrepancy = float(np.linalg.norm(field-original)/max(np.linalg.norm(original), 1e-300))
                    if discrepancy > 1e-9:
                        raise ValueError("Real A truth response/solver identity does not match B simulator")
        rows = [(parent, k, field) for (parent, k), field in self._fields.items()]
        np.savez_compressed(self.output / "OFFLINE_SIMULATOR_FIELDS.npz",
                            parent=np.array([r[0] for r in rows]),
                            frequency=np.array([r[1] for r in rows]),
                            full_field=np.array([r[2] for r in rows]).reshape(-1, 32))
        dump(self.output / "OFFLINE_SIMULATOR_STATUS.json", dict(
            parents=sorted(self._truth), errors=self._errors, seed=SEED_NOISE,
            persistent_gain=True, original_measurements_replaced=False,
            source="A own confirmation population, not B E1 fresh parents"))

    def measure(self, parent, frequency):
        key = (parent, float(frequency))
        if key not in self._future:
            raise RuntimeError("No validated simulator response for selected action")
        return self._future[key].copy()

    def score(self, parent, selected, catalogue):
        truth = np.asarray(self._truth[parent]["latent"], float)
        best = float(field_l1(catalogue, truth).min())
        actual = float(field_l1(selected, truth)[0, 0])
        return dict(selected_MAE=actual, best_of_K_MAE=best, Coverage_at_K_005=best <= .05,
                    selected_event_id=int(event_ids(selected)),
                    selected_structure_correct=bool(event_ids(selected) == event_ids(truth)),
                    K=len(catalogue), material_metric="R1 exact componentwise complex-epsilon volume L1, s_ref1")


def candidate_menu(bundle, full_fields, sigma):
    actions = []
    for action in bundle["package"]["available_actions"]:
        k = float(action["frequency_kR"])
        j = FREQUENCIES.index(k)
        if np.isfinite(full_fields[:, j]).all():
            actions.append(Action(action["id"], full_fields[:, j], np.full(32, sigma[j])))
    return actions


def action_frequency(bundle, action_id):
    return float(next(a["frequency_kR"] for a in bundle["package"]["available_actions"]
                      if a["id"] == action_id))


def select_profile(catalogue, joint_fields, y, sigma):
    residual, gains = profile(joint_fields, y, sigma)
    index = int(np.argmin(residual))
    compatible = residual <= residual[index]+9.
    return index, dict(selected_profile_chi2=float(residual[index]), selected_gain=float(gains[index]),
                       compatible_candidate_count=int(compatible.sum()),
                       compatible_event_ids=np.unique(event_ids(catalogue[compatible])).tolist(),
                       candidate_weights_kind="none", posterior_probability=None)


def reference_acquisition(response, reference, bundle, sigma, parent_index, available_ids):
    """Independent B1 continuous prior; never a probability model of the A list."""
    y = bundle["assets"]["observations"]
    sd = np.sqrt(np.diag(bundle["assets"]["covariance"]))
    solved = reference.solve(y, sd, order=REFERENCE_ORDER)
    ids = resample_indices(solved["mass"], REFERENCE_PARTICLES, 26101594+parent_index)
    did = resample_indices(solved["mass"], REFERENCE_DECISIONS, 26101598+parent_index)
    x = solved["bank"]["x"][ids]
    decisions = solved["bank"]["x"][did]
    f = response.predict(x, FREQUENCIES).reshape(-1, 6, 32)
    w = np.full(len(x), 1/len(x))
    utilities = {}
    for a in bundle["package"]["available_actions"]:
        if a["id"] not in available_ids:
            continue
        j = FREQUENCIES.index(float(a["frequency_kR"]))
        action = Action(a["id"], f[:, j], np.full(32, sigma[j]))
        utility = expected_utility(x, w, f[:, 1], y, sd, action, decisions,
                                   seed=26101604+parent_index, outer=REFERENCE_OUTER)
        utilities[a["id"]] = {k: v for k, v in utility.items() if k != "risk_draws"}
    chosen = max(utilities, key=lambda aid: utilities[aid]["risk_reduction"])
    return dict(action=chosen, utilities=utilities, x=x, w=w, F=f, decisions=decisions,
                initial_event_mass=solved["event_mass"],
                target="independent B1 coefficient-uniform continuous prior and E1 six-frequency response proxy; not A candidate posterior")


def optimize_endpoint(response, start, frequencies, y, sd, start_gain, protected):
    """Ordinary bounded stochastic-start least squares; no truth arguments."""
    start = np.asarray(start, float)
    lower = np.r_[BOUNDS[2:, 0] if protected else BOUNDS[:, 0], .7]
    upper = np.r_[BOUNDS[2:, 1] if protected else BOUNDS[:, 1], 1.3]
    initial = np.r_[start[2:] if protected else start, start_gain]
    def latent(t):
        x = start.copy() if protected else np.asarray(t[:4]).copy()
        if protected:
            x[2:] = t[:2]
        return x
    def residual(t):
        fields = response.predict(latent(t), frequencies)[0]
        return (t[-1]*fields-y)/sd
    t0 = time.perf_counter()
    try:
        fit = least_squares(residual, initial, bounds=(lower, upper), max_nfev=MAX_NFEV,
                             xtol=1e-10, ftol=1e-10, gtol=1e-10)
        x = latent(fit.x)
        if not np.isfinite(fit.fun).all() or not np.isfinite(prior_logpdf(x)):
            raise FloatingPointError("Nonfinite or out-of-support optimizer endpoint")
        return dict(status="PROXY_ENDPOINT", latent=x, shared_proxy_gain=float(fit.x[-1]),
                    optimizer_success=bool(fit.success), optimizer_status=int(fit.status),
                    optimizer_message=str(fit.message), nfev=int(fit.nfev),
                    proxy_joint_chi2=float(fit.fun@fit.fun),
                    proxy_fields=response.predict(x, frequencies)[0].reshape(2, 32),
                    protected_geometry=protected, initial=start,
                    geometry_preserved=bool(np.array_equal(x[:2], start[:2])) if protected else None,
                    wall_seconds=time.perf_counter()-t0)
    except Exception as error:
        return dict(status="PROXY_OPTIMIZER_FAILURE", error=repr(error), initial=start,
                    optimizer_success=False, protected_geometry=protected,
                    wall_seconds=time.perf_counter()-t0)


def validate_endpoint(paid, fit, frequencies, y, sd, category, label, existing, log_path):
    row = {**fit, "label": label, "installed": False,
           "full_validation_frequencies": list(frequencies)}
    if fit["status"] != "PROXY_ENDPOINT":
        append(log_path, row)
        return row
    fields, statuses, infos = [], [], []
    try:
        for k in frequencies:
            field, info, status = paid.call(fit["latent"], k, category, label)
            statuses.append(status)
            infos.append(info)
            if field is None:
                row.update(status="FULL_VALIDATION_FAILURE", full_statuses=statuses)
                append(log_path, row)
                return row
            fields.append(field)
    except BudgetExceeded:
        row.update(status="BUDGET_EXHAUSTED_PARTIAL_VALIDATION", full_statuses=statuses,
                   validated_block_count=len(fields))
        append(log_path, row)
        raise
    joined = np.concatenate(fields)
    residual, gain = profile(joined[None], y, sd)
    threshold = float(chi2.ppf(1-ALPHA_COMPATIBILITY, len(y)))
    duplicate = bool(field_l1(fit["latent"], existing).min() <= 1e-12)
    compatible = bool(residual[0] <= threshold)
    row.update(status="FULL_VALIDATED_COMPATIBLE" if compatible else "FULL_VALIDATED_DATA_INCOMPATIBLE",
               full_fields=np.asarray(fields), full_statuses=statuses, full_solver_info=infos,
               full_joint_profile_chi2=float(residual[0]), full_shared_gain=float(gain[0]),
               nominal_A_style_joint_threshold=threshold, heuristic_compatibility=True,
               duplicate_material=duplicate, installed=bool(compatible and not duplicate),
               proxy_full_whitened_difference=float(np.linalg.norm(
                   (joined-np.asarray(fit["proxy_fields"]).ravel())/sd)))
    append(log_path, row)
    return row


def export_refresh(output, name, bundle, catalogue, fields, endpoint_rows, original_count=None):
    destination = output / name
    destination.mkdir(parents=True, exist_ok=True)
    original_edges = bundle["assets"]["radial_interval_edges"]
    edges = np.unique(np.r_[original_edges, *[decode(x)["radii"] for x in catalogue]])
    np.savez_compressed(destination / "FULL_VALIDATED_CATALOGUE.npz", latent=catalogue,
                        joint_unscaled_fields=fields, radial_interval_edges=edges,
                        radial_volume_weights=np.diff(edges**3),
                        complex_material_fields=np.array([epsilon_on_grid(x, edges) for x in catalogue]))
    dump(destination / "CATALOGUE_STATUS.json", dict(
        original_A_material_arrays=str(bundle["source"]), original_source_arrays_unchanged=True,
        K_original=(len(bundle["assets"]["p2_latent_coordinates"]) if original_count is None else original_count),
        K_full_validated=len(catalogue),
        installed_endpoints=sum(r.get("installed", False) for r in endpoint_rows),
        failed_optimizers=sum(r["status"] == "PROXY_OPTIMIZER_FAILURE" for r in endpoint_rows),
        nonconverged_optimizers=sum(not r.get("optimizer_success", False) for r in endpoint_rows),
        rejected_full_endpoints=sum(r["status"] == "FULL_VALIDATED_DATA_INCOMPATIBLE" for r in endpoint_rows),
        all_trials_path=str(output / "ENDPOINT_TRIALS.jsonl"), candidate_weights_kind="none"))


def write_exchange(output, method, bundle, action_id, new_y, new_sd, reference_policy=False):
    destination = output / "AcquisitionRecommendation" / method
    if destination.exists():
        return str(destination)
    criterion = ("independent_B1_continuous_reference_material_risk" if reference_policy
                 else "full_A_set_robust_action_then_base_capacity_control" if "same_action" in method
                 else "robust_finite_A_candidate_composite" if "robust" in method
                 else "fixed_frequency_control")
    recommendation = make_recommendation(bundle, [action_id], criterion=criterion,
                                           status="SIMULATION_DESIGN_ONLY")
    if reference_policy:
        recommendation["criterion_probability_status"] = "external_B1_continuous_proxy_reference; A candidates remain unweighted"
        recommendation["limitations"].append("Action requires external B1 prior/reference assets, not candidate-only deployment")
        recommendation["limitations"].append("New-history reference audit UNRESOLVED: approximate diagnostic, not a reliable upper bound")
    if "same_action" in method:
        recommendation["limitations"].append("Matched-action base capacity control: action chosen by the full A candidate set, not by the base alone")
    covariance = np.zeros((64, 64))
    covariance[:32, :32] = bundle["assets"]["covariance"]
    covariance[32:, 32:] = np.diag(new_sd**2)
    aid = bundle["package"]["package_id"] + ":E3:" + action_id
    record = dict(target_id=bundle["package"]["target_id"], solver_id=bundle["package"]["model"]["solver_id"],
                  units=bundle["contract"]["input_context"]["observation_units"], packing="all_real_then_all_imag",
                  observation_id=aid+":observation", noise_realization_id=aid+":noise-seed26101591",
                  provenance_id=str(output / "SIMULATED_NEW_OBSERVATIONS.npz")+"#"+action_id,
                  action_id=action_id, nuisance_group_by_parameter=bundle["contract"]["input_context"]["nuisance"]["group_by_parameter"],
                  source="simulation", noise_evidence="Independent parent/action SeedSequence[26101591,parent,frequency_index]; same persistent world gain",
                  independent_new_noise=True, values=new_y, joint_covariance=covariance)
    history = updated_data_history(bundle, record, updated_history_id=aid+":joint-history")
    if "corrector" in method:
        history["candidate_refresh_status"] = "B_INDEPENDENT_FULL_VALIDATED_REFRESH_IN_SEPARATE_CATALOGUE"
        history["refreshed_catalogue_uri"] = str(output / method / "FULL_VALIDATED_CATALOGUE.npz")
    export_bundle({k: bundle[k] for k in ("package", "contract", "assets")}, destination,
                  recommendation=recommendation, new_history=history)
    return str(destination)


def run(cost, output, resume=False):
    frozen = json.loads((ROOT / "configs/E3_IMPORT_FROZEN.json").read_text())
    if frozen["new_attempt_subcap"] > 240 or frozen["packages"] != [f"confirm_{i:02d}" for i in range(4)]:
        raise ValueError("E3 fixed first-four parent/subcap contract changed")
    if (output / "RESULTS.json").exists():
        raise FileExistsError("Completed real A handoff results must not be overwritten")
    if output.exists() and any(output.iterdir()) and not resume:
        raise FileExistsError("Partial output exists: review it and use --resume with a new Cost job identity")
    output.mkdir(parents=True, exist_ok=True)
    bundles = [load_bundle(parent) for parent in frozen["packages"]]
    if any(b["contract"]["B_menu_authorization_id"] != cost.overlay["authorization_id"] for b in bundles):
        raise ValueError("Real A B-menu snapshots do not match this independent authorization")
    if sum(len(b["assets"]["p2_latent_coordinates"]) for b in bundles) != frozen["K_total"]:
        raise ValueError("Real A first-four K total differs from frozen27")
    with np.load(ROOT / "results/E1_NOISE_ANCHORS.npz", allow_pickle=False) as data:
        if tuple(data["frequencies"]) != FREQUENCIES:
            raise ValueError("E1 noise/menu frequency identity mismatch")
        sigma = data["sigma"].copy()
    for bundle in bundles:
        for action in bundle["package"]["available_actions"]:
            frequency = float(action["frequency_kR"])
            covariance_name = Path(action["noise"]["covariance"]["uri"]).stem
            declared = bundle["assets"][covariance_name]
            expected = np.eye(32)*sigma[FREQUENCIES.index(frequency)]**2
            if not np.allclose(declared, expected, rtol=1e-12, atol=1e-15):
                raise ValueError("B menu covariance differs from the fixed-reference E1 action noise")
    details = dict(frozen=frozen, source_snapshots=[str(b["source"]) for b in bundles],
                   execution_stage="E2", reference_order=REFERENCE_ORDER,
                   reference_particles=REFERENCE_PARTICLES, reference_decisions=REFERENCE_DECISIONS,
                   reference_predictive_draws=REFERENCE_OUTER, optimization_max_nfev=MAX_NFEV,
                   future_noise_seed=SEED_NOISE, prior_multistart_seed=SEED_STARTS,
                   endpoint_acceptance="full fields on both history blocks; legal material; source-A-style chi2(.95,real_dimension) fit envelope; no truth-based acceptance",
                   reference_endpoint_check_scope="first3 parents: selected new-frequency field only; not installed into A catalogue",
                   new_history_reference_audit="UNRESOLVED; reference comparator/action diagnostic cannot be a reliable upper bound",
                   source_cost_scope="A global ledger snapshots repeated per package must not be summed",
                   reporting_scope="four actual A parents, descriptive only; not B's32-parent population")
    previous = output / "EXECUTION_DETAILS.json"
    if previous.exists() and json.loads(previous.read_text()) != details:
        raise ValueError("Resume cannot change frozen E3 implementation details")
    dump(previous, details)
    # The orchestrator's reviewed execution is the only entry that constructs
    # the paid backend. No solver is instantiated by the default CLI preview.
    from prepare_e1 import forward
    paid = PaidFields(cost, forward(cost), output, frozen)
    for bundle in bundles:
        for i, (x, f0) in enumerate(zip(bundle["assets"]["p2_latent_coordinates"],
                                       bundle["assets"]["unscaled_old_fields"])):
            paid.seed_original(x, f0, str(bundle["source"] / f"prediction_{i}.npy"))
    simulator = OfflineSimulator(frozen["packages"], output)
    for bundle in bundles:
        simulator.verify_original(bundle["parent"], bundle["assets"]["observations"])
    response = Response(ROOT / "results/E1_RESPONSE_BANK.npz")
    reference = Reference(response)
    records, parent_summaries = [], []
    simulator.prepare(paid, sigma)
    for parent_index, bundle in enumerate(bundles):
        cost.check_time()
        parent = bundle["parent"]
        directory = output / parent
        directory.mkdir(parents=True, exist_ok=True)
        x = bundle["assets"]["p2_latent_coordinates"].copy()
        y0 = bundle["assets"]["observations"]
        s0 = np.sqrt(np.diag(bundle["assets"]["covariance"]))
        full = np.full((len(x), 6, 32), np.nan)
        full[:, 1] = bundle["assets"]["unscaled_old_fields"]
        failed_predictions = []
        for j, k in enumerate(FREQUENCIES):
            if k == OLD:
                continue
            for i, latent in enumerate(x):
                field, _, status = paid.call(latent, k, CANDIDATE_CATEGORY, parent+"_CANDIDATE")
                if field is None:
                    failed_predictions.append(dict(candidate=i, frequency=k, status=status))
                else:
                    full[i, j] = field
        np.savez_compressed(directory / "A_FULL_ACTION_FIELDS.npz", latent=x, frequencies=FREQUENCIES,
                            full_fields=full, original_unscaled_fields=full[:, 1],
                            original_source_gains=bundle["assets"]["source_fitted_gains"])
        actions = candidate_menu(bundle, full, sigma)
        if not actions:
            raise RuntimeError("No fully validated candidate action remains")
        scores = {a.id: candidate_scores(x, full[:, 1], y0, s0, a.fields, a.sigma) for a in actions}
        robust = max(scores, key=lambda aid: scores[aid]["robust"])
        fixed = next((a.id for a in actions if action_frequency(bundle, a.id) == 3.2), None)
        external = reference_acquisition(response, reference, bundle, sigma, parent_index,
                                         {a.id for a in actions})
        selections = {"A_robust_action_then_profile": robust, "A_fixed3.2_then_profile": fixed,
                      "A_reference_action_then_profile": external["action"]}
        dump(directory / "ACTION_RANKINGS.json", dict(scores=scores, selected=selections,
             failed_candidate_fields=failed_predictions, candidate_weights_kind="none",
             external_reference={k: v for k, v in external.items() if k not in ("x", "w", "F", "decisions")}))
        initial_index, initial_fit = select_profile(x, full[:, 1], y0, s0)
        for method, index, catalogue in (("A_base_only_no_measurement", bundle["base_index"], x[[bundle["base_index"]]]),
                                         ("A_fixed_candidate_profile_no_measurement", initial_index, x)):
            record = dict(parent=parent, method=method, measurement_blocks=0, acquisition_input="A candidates only",
                          full_validated=True, generator_refresh=False, fixed_best_of_K_unchanged=True,
                          **simulator.score(parent, x[index], catalogue), **initial_fit)
            if method == "A_base_only_no_measurement":
                _, base_fit = select_profile(catalogue, full[[index], 1], y0, s0)
                record.update(base_fit)
            records.append(record)
        dump(output / "PARTIAL_RESULTS.json", dict(records=records, parents=parent_summaries,
                                                   local_cost=paid.receipt(), partial=True))
        observed = {}
        for action_id in dict.fromkeys(selections.values()):
            if action_id is None:
                continue
            k = action_frequency(bundle, action_id)
            new_y = simulator.measure(parent, k)
            new_sd = np.full(32, sigma[FREQUENCIES.index(k)])
            observed[action_id] = (new_y, new_sd)
        np.savez_compressed(directory / "SIMULATED_NEW_OBSERVATIONS.npz",
                            action_ids=np.array(list(observed)),
                            frequency=np.array([action_frequency(bundle, a) for a in observed]),
                            new_observations=np.array([observed[a][0] for a in observed]),
                            new_sigma=np.array([observed[a][1] for a in observed]))
        for method, action_id in selections.items():
            if action_id is None:
                records.append(dict(parent=parent, method=method, status="ACTION_NOT_AVAILABLE_PRESERVE_A",
                                    measurement_blocks=0, full_validated=True, generator_refresh=False,
                                    fixed_best_of_K_unchanged=True, candidate_weights_kind="none",
                                    **simulator.score(parent, x[initial_index], x)))
                continue
            k = action_frequency(bundle, action_id)
            new_y, ns = observed[action_id]
            joint_fields = np.column_stack((full[:, 1], full[:, FREQUENCIES.index(k)]))
            joint_y, joint_sd = np.r_[y0, new_y], np.r_[s0, ns]
            chosen, fit = select_profile(x, joint_fields, joint_y, joint_sd)
            handoff = write_exchange(directory, method, bundle, action_id, new_y, ns,
                                     reference_policy="reference" in method)
            records.append(dict(parent=parent, method=method, selected_action=action_id, frequency=k,
                                measurement_blocks=1, acquisition_input="external B1 continuous proxy reference" if "reference" in method else "A candidates or fixed menu",
                                full_validated=True, generator_refresh=False, fixed_best_of_K_unchanged=True,
                                acquisition_reference_status="UNRESOLVED_AFTER_NEW_HISTORY" if "reference" in method else "NOT_REQUIRED_CANDIDATE_ONLY_OR_FIXED",
                                handoff_dir=handoff, **simulator.score(parent, x[chosen], x), **fit))
        # Reference decision is a separate approximate Bayesian comparator,
        # never assigned as posterior weights on A's unweighted candidates.
        action_id = external["action"]
        k = action_frequency(bundle, action_id)
        new_y, ns = observed[action_id]
        increment = conditional_log_increment(external["F"][:, 1], y0, s0,
                                               external["F"][:, FREQUENCIES.index(k)], new_y, ns)[0]
        weights = normalized(np.log(external["w"])+increment)
        decision_losses = field_l1(external["decisions"], external["x"])
        chosen = int(np.argmin(decision_losses@weights))
        reference_selected = external["decisions"][chosen]
        audit = dict(status="NOT_IN_FIRST3_CHECKS", installed_into_A=False, full_history_validated=False)
        if parent_index < 3:
            fcheck, _, status = paid.call(reference_selected, k, REFERENCE_CATEGORY, parent+"_REFERENCE_SELECTED_CHECK")
            audit = dict(status=status, frequency=k, installed_into_A=False, full_history_validated=False)
            if fcheck is not None:
                fp = response.predict(reference_selected, (k,))[0]
                audit["new_frequency_proxy_full_whitened"] = float(np.linalg.norm((fcheck-fp)/ns))
        records.append(dict(parent=parent, method="declared_continuous_P2_reference_comparator", selected_action=action_id,
                            frequency=k, measurement_blocks=1, acquisition_input=external["target"],
                            full_validated=False, reference_audit=audit, installed_into_A=False,
                            reference_quality="finite independent proxy comparator, not a certified upper bound",
                            new_history_reference_audit="UNRESOLVED",
                            posterior_particle_ESS=float(1/(weights@weights)),
                            generator_refresh=False, candidate_weights_kind="A:none; separate independent B1 finite reference weights",
                            **simulator.score(parent, reference_selected, external["decisions"])))
        # Correctors operate on the robust policy's actual new data only.
        k = action_frequency(bundle, robust)
        new_y, ns = observed[robust]
        frequencies, joint_y, joint_sd = (OLD, k), np.r_[y0, new_y], np.r_[s0, ns]
        original_joint = np.column_stack((full[:, 1], full[:, FREQUENCIES.index(k)]))
        protected_pool, prior_pool = x.copy(), x.copy()
        protected_fields, prior_fields = original_joint.copy(), original_joint.copy()
        protected_trials, prior_trials = [], []
        log_path = directory / "ENDPOINT_TRIALS.jsonl"
        for i, latent in enumerate(x):
            cost.check_time()
            _, g = profile(original_joint[[i]], joint_y, joint_sd)
            fit = optimize_endpoint(response, latent, frequencies, joint_y, joint_sd, float(g[0]), True)
            trial = validate_endpoint(paid, fit, frequencies, joint_y, joint_sd, PROTECTED_CATEGORY,
                                      parent+f"_PROTECTED_{i:02d}", protected_pool, log_path)
            protected_trials.append(trial)
            if trial["installed"]:
                protected_pool = np.vstack((protected_pool, trial["latent"]))
                protected_fields = np.vstack((protected_fields, np.asarray(trial["full_fields"]).ravel()))
        # No new physical work: use exactly the already-validated protected
        # endpoint of the original A base. Acquisition still used the full set.
        base_index = bundle["base_index"]
        base_pool = x[[base_index]].copy()
        base_fields = original_joint[[base_index]].copy()
        base_trial = protected_trials[base_index]
        base_added = False
        if (base_trial["status"] == "FULL_VALIDATED_COMPATIBLE"
                and field_l1(base_trial["latent"], base_pool).min() > 1e-12):
            base_pool = np.vstack((base_pool, base_trial["latent"]))
            base_fields = np.vstack((base_fields, np.asarray(base_trial["full_fields"]).ravel()))
            base_added = True
        base_method = "A_base_only_same_action_eps_corrector"
        base_chosen, base_fit = select_profile(base_pool, base_fields, joint_y, joint_sd)
        base_score = simulator.score(parent, base_pool[base_chosen], base_pool)
        base_old_best = simulator.score(parent, x[base_index], x[[base_index]])["best_of_K_MAE"]
        base_handoff = write_exchange(directory, base_method, bundle, robust, new_y, ns)
        base_export_trial = {**base_trial, "installed": base_added}
        export_refresh(directory, base_method, bundle, base_pool, base_fields, [base_export_trial], original_count=1)
        records.append(dict(parent=parent, method=base_method, selected_action=robust, frequency=k,
            measurement_blocks=1, acquisition_input="full A set robust action; matched-action base capacity control",
            matched_action_control=True, pure_base_acquisition_claimed=False,
            acquisition_candidate_K=len(x), K_original=1, full_validated=True,
            candidate_refresh_attempted=True, generator_refresh=base_added,
            new_full_validated_candidates=int(base_added), fixed_candidate_best_of_K=base_old_best,
            fixed_best_of_K_unchanged=True,
            refreshed_best_of_K_difference=base_score["best_of_K_MAE"]-base_old_best,
            refresh_best_of_K_change_claimed=bool(base_score["best_of_K_MAE"] < base_old_best-1e-12),
            additional_callback_requests=0, endpoint_reused_from="protected base_index trial",
            handoff_dir=base_handoff, **base_score, **base_fit))
        starts = prior_sample(np.random.default_rng(SEED_STARTS+parent_index*101), 3)
        dump(directory / "PRIOR_MULTISTART_INITIALS.json", dict(seed=SEED_STARTS+parent_index*101,
              starts=starts, topology_labels_used=False, prior="original independent coefficient-uniform product"))
        for i, start in enumerate(starts):
            cost.check_time()
            fit = optimize_endpoint(response, start, frequencies, joint_y, joint_sd, 1., False)
            trial = validate_endpoint(paid, fit, frequencies, joint_y, joint_sd, PRIOR_CATEGORY,
                                      parent+f"_PRIOR_START_{i}", prior_pool, log_path)
            prior_trials.append(trial)
            if trial["installed"]:
                prior_pool = np.vstack((prior_pool, trial["latent"]))
                prior_fields = np.vstack((prior_fields, np.asarray(trial["full_fields"]).ravel()))
        for method, pool, fields, trials in (
            ("A_robust_then_protected_eps_corrector", protected_pool, protected_fields, protected_trials),
            ("A_robust_then_independent_prior_multistart_corrector", prior_pool, prior_fields, prior_trials)):
            chosen, fit = select_profile(pool, fields, joint_y, joint_sd)
            handoff = write_exchange(directory, method, bundle, robust, new_y, ns)
            export_refresh(directory, method, bundle, pool, fields, trials)
            fixed_best = simulator.score(parent, x[0], x)["best_of_K_MAE"]
            refreshed_score = simulator.score(parent, pool[chosen], pool)
            records.append(dict(parent=parent, method=method, selected_action=robust, frequency=k,
                measurement_blocks=1, acquisition_input="A finite candidates only",
                full_validated=True, candidate_refresh_attempted=True,
                generator_refresh=bool(len(pool)>len(x)), K_original=len(x),
                new_full_validated_candidates=int(len(pool)-len(x)),
                fixed_candidate_best_of_K=fixed_best,
                fixed_best_of_K_unchanged=True,
                refreshed_best_of_K_difference=refreshed_score["best_of_K_MAE"]-fixed_best,
                refresh_best_of_K_change_claimed=bool(refreshed_score["best_of_K_MAE"] < fixed_best-1e-12),
                optimizer_failures=sum(not r.get("optimizer_success", False) for r in trials),
                optimizer_exceptions=sum(r["status"] == "PROXY_OPTIMIZER_FAILURE" for r in trials),
                data_incompatible_endpoints=sum(r["status"] == "FULL_VALIDATED_DATA_INCOMPATIBLE" for r in trials),
                full_validation_failures=sum(r["status"] == "FULL_VALIDATION_FAILURE" for r in trials),
                handoff_dir=handoff, **refreshed_score, **fit))
        parent_summaries.append(dict(parent=parent, source_snapshot=str(bundle["source"]), K_original=len(x),
                                     selected_actions=selections, endpoint_trials=len(protected_trials)+len(prior_trials)))
        dump(output / "PARTIAL_RESULTS.json", dict(records=records, parents=parent_summaries,
                                                   local_cost=paid.receipt(), partial=True))
        print("E3", parent, "complete; globalattempts", cost.count(), "localrequests", paid.attempts, flush=True)
    summary = []
    for method in frozen["primary_methods"]:
        group = [r for r in records if r["method"] == method]
        summary.append(dict(method=method, parents=len(group),
                            selected_MAE=float(np.mean([r["selected_MAE"] for r in group])),
                            best_of_K_MAE=float(np.mean([r["best_of_K_MAE"] for r in group])),
                            Coverage_at_K_005=float(np.mean([r["Coverage_at_K_005"] for r in group])),
                            mean_K=float(np.mean([r["K"] for r in group])),
                            scope="four real A parents; descriptive, no large-population claim"))
    result = dict(status="REAL_A_CANDIDATE_HANDOFF_EXECUTED", frozen=frozen, records=records,
                  parent_summaries=parent_summaries, summary=summary, local_cost=paid.receipt(),
                  global_cost=cost.summary(), A_candidate_probabilities_assigned=False,
                  global_cost_scope="snapshot before this job's final compute receipt; FINAL_COST.json and COMPUTE_LEDGER include completed job/export",
                  source_A_modified=False, hardware_measurements=0, GPU=0,
                  method_cost_scope="One measured experiment portfolio shares all-menu candidate prediction and simulator setup; local categories retain cold costs. Do not sum counterfactual rows or claim isolated per-policy speed from these shared costs.",
                  limitations=["Four A parents differ from B E1 population; never pool them as a32-parent replication",
                     "B reference is an independent declared continuous proxy target, not posterior weights on A candidates",
                     "Raw source union-grid fields and original noisy observations retained in each handoff",
                     "Only explicitly full-validated compatible corrector endpoints change refreshed Best-of-K",
                     "Reference selected endpoints have at most one new-frequency audit; not installed as full-history validated candidates",
                     "New-history reference numerical audit is UNRESOLVED; neither reference action nor comparator is a reliable upper bound",
                     "All profile gains are one legal scalar per material world across history; no residual softmax probabilities",
                     "Local reservations include failed callbacks and do not borrow or reset historical balances"])
    dump(output / "RESULTS.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true", help="Orchestrator-reviewed paid execution under the existing frozen authorization")
    parser.add_argument("--job", default="e3-real-a-handoff-v1")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", type=Path, default=ROOT / "results/E3")
    args = parser.parse_args()
    output = args.output.resolve()
    allowed = (ROOT / "results/E3").resolve()
    if not output.is_relative_to(allowed):
        raise ValueError("E3 outputs must stay inside Research B results/E3")
    if not args.execute:
        print(json.dumps(dict(status="IMPLEMENTED_NOT_EXECUTED", frozen=str(ROOT / "configs/E3_IMPORT_FROZEN.json"),
                              source_snapshots="confirm_00..03_B_menu_v3", subcap=240,
                              paid_command="../../.venv_nn/bin/python research/research_b/code/e3_real_handoff.py --execute --job e3-real-a-handoff-v1",
                              new_Maxwell_attempts=0), indent=2))
        return
    cost = Cost()
    with cost.job(args.job, "E2"):
        run(cost, output, resume=args.resume)
    dump(output / "FINAL_COST.json", cost.summary())


if __name__ == "__main__":
    main()
