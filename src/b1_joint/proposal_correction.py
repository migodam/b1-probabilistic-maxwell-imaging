"""Matched full-Mie importance correction of frozen C1 material proposals.

Nothing runs on import, no solver/ledger is constructed, and no file is written.
The caller supplies the paid MaterialForward and meters the entire call.  Every
new response goes through forward.predict, with cache=False for equal budgets.
The C1 p2 directory must already be on sys.path; its portable NumPy density is
used for the exact proposal law.  The returned dictionaries can be saved by the
parent.  This is finite-sample numerical-target importance sampling, not a
posterior calibration, no-missed-mode, or continuum guarantee.
"""
from __future__ import annotations

import importlib
import math
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp
from scipy.stats import truncnorm

from .material import BOUNDS, WIDTH, event_ids, geometry, prior_logpdf, prior_sample
from .likelihood import conditional_gain, log_likelihood


METHODS = ("full_distill", "diag_distill", "gaussian_distill")
MODEL_SEED = 26101321
PROPOSAL_TARGET = "joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1"
FULL_TARGET = "full-vector-Mie-matching-n14:C:bounded-shared-gain"
DEFAULT_SAMPLE_SEED = 261013304
DEFAULT_BOOTSTRAP_SEED = 261013305
PRIOR_FRACTION = .1
SIGNIFICANT_MASS = .05
MASS_AWARE_DISCOVERY = .01
FALSE_EXCLUSION_MASS = .005


def _api():
    """Lazy C1 imports contain no physical response call."""
    density = importlib.import_module("c1p2.density")
    load = importlib.import_module("c1p2.train").load_model
    return density, load


def measurement_condition(y, sigma, gain_sd=.15):
    """The exact 34-coordinate C1 condition used for this one-frequency target."""
    y = np.asarray(y, dtype=float)
    sigma = np.broadcast_to(np.asarray(sigma, dtype=float), y.shape)
    if y.shape != (32,) or not np.isfinite(y).all() or not np.isfinite(sigma).all():
        raise ValueError("one-frequency C observation needs 32 finite Re/Im channels")
    if np.any(sigma <= 0) or not np.all(sigma == sigma[0]):
        raise ValueError("frozen C1 condition supports one positive shared real-channel sigma")
    if gain_sd != .15:
        raise ValueError("these frozen models were trained with gain_sd=0.15")
    return np.r_[y, np.log(sigma[0]), gain_sd][None, :]


def normalized_parameters(parameters):
    """Normalize the global categorical law, never an evaluation-node density.

    FP32 log-softmax exports can sum to 1 plus roundoff in FP64.  C1's NumPy
    categorical sampler already divides by that sum.  Subtracting this known
    global log normalizer also in the density makes the two laws identical.
    The export discrepancy is retained, and a non-roundoff discrepancy fails.
    """
    p = {key: np.asarray(value, dtype=float).copy() for key, value in parameters.items()}
    if not np.array_equal(p["lower"], BOUNDS[:, 0]) or not np.array_equal(p["upper"], BOUNDS[:, 1]):
        raise ValueError("C1 and B1 material support differ")
    means, factors, logw = p["means"], p["cholesky"], p["log_weights"]
    if means.ndim != 3 or means.shape[0] != 1 or means.shape[2] != 4:
        raise ValueError("exactly one four-coordinate material condition required")
    if factors.shape != (*means.shape[:2], 4, 4) or logw.shape != means.shape[:2]:
        raise ValueError("invalid exported mixture dimensions")
    if not all(np.isfinite(a).all() for a in (means, factors, logw)):
        raise FloatingPointError("nonfinite exported density parameters")
    if np.any(np.diagonal(factors, axis1=-2, axis2=-1) <= 0):
        raise ValueError("positive Cholesky diagonal required")
    if np.any(np.triu(factors, k=1) != 0):
        raise ValueError("exported Cholesky is not triangular")
    normalization = logsumexp(logw, axis=-1, keepdims=True)
    if np.max(np.abs(normalization)) > 1e-5:
        raise ValueError("mixture normalization differs beyond exported softmax roundoff")
    p["log_weights"] = logw - normalization
    p["weights"] = np.exp(p["log_weights"])
    return p, float(normalization.item())


def load_conditioned_proposal(model_path, condition, expected_target=PROPOSAL_TARGET):
    density, load = _api()
    model, receipt = load(Path(model_path), expected_target=expected_target)
    if receipt.get("seed") != MODEL_SEED or receipt.get("method") not in METHODS:
        raise ValueError("matched correction requires one of the three seed 26101321 frozen proposals")
    parameters, delta = normalized_parameters(model.export_parameters(condition))
    return parameters, dict(method=receipt["method"], model_seed=receipt["seed"],
                            selected_checkpoint=receipt.get("selected_checkpoint"),
                            training_updates=receipt.get("updates"),
                            proposal_target_identity=expected_target,
                            categorical_log_normalizer_removed=delta,
                            normalization_scope="explicit global categorical roundoff, not a material-grid or posterior normalization"), density


def sample_defensive(parameters, density, rng, n=64):
    """IID q=.9 bounded-MDN+.1 proper-uniform draws with exact dx logq."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError("at least two iid proposal draws required")
    mdn_branch = rng.random(n) >= PRIOR_FRACTION
    latent = np.empty((n, 4), dtype=float)
    transformed = np.empty_like(latent)
    component = np.full(n, -1, dtype=np.int64)
    rounded = np.zeros(n, dtype=bool)
    draw = density.sample_parameters(parameters, int(mdn_branch.sum()), rng)
    latent[mdn_branch] = draw["latent"][0]
    transformed[mdn_branch] = draw["transformed"][0]
    component[mdn_branch] = draw["component"][0]
    rounded[mdn_branch] = draw["rounded_endpoint"][0]
    latent[~mdn_branch] = prior_sample(rng, int((~mdn_branch).sum()))
    if np.any(rounded) or np.any((latent <= BOUNDS[:, 0]) | (latent >= BOUNDS[:, 1])):
        raise FloatingPointError("rounded/outside-support proposal draw; no clipping or replacement allowed")
    transformed[~mdn_branch] = density.to_transformed(latent[~mdn_branch])
    # Retained transformed coordinates avoid an inverse-logistic roundoff.
    mdn_logq = np.asarray(density.log_prob_batch(transformed, parameters,
                                               transformed=True, material_density=True))
    lp0 = np.asarray(prior_logpdf(latent))
    logq = np.logaddexp(np.log(.9) + mdn_logq, np.log(.1) + lp0)
    if logq.shape != (n,) or not np.isfinite(logq).all() or not np.isfinite(lp0).all():
        raise FloatingPointError("invalid exact defensive-mixture density")
    return dict(latent=latent, transformed=transformed, mdn_branch=mdn_branch,
                mdn_component=component, rounded_endpoint=rounded,
                log_q_mdn=mdn_logq, log_q=logq, log_p0=lp0,
                event_ids=event_ids(latent), geometry=geometry(latent))


def _raw_summary(draw):
    ids, n = draw["event_ids"], len(draw["event_ids"])
    count = np.bincount(ids, minlength=12)
    p = count / n
    mdn_n = int(draw["mdn_branch"].sum())
    mdn_count = np.bincount(ids[draw["mdn_branch"]], minlength=12)
    return dict(defensive_proposal_event_count=count, defensive_proposal_event_mass=p,
                defensive_proposal_event_SE=np.sqrt(p * (1 - p) / (n - 1)),
                mdn_branch_sample_count=mdn_n, prior_branch_sample_count=n - mdn_n,
                mdn_branch_event_count=mdn_count,
                mdn_branch_event_mass=mdn_count / mdn_n if mdn_n else None,
                status="uncorrected proposal probabilities, not physical posterior probabilities")


def _reference_case(reference, observation_id):
    r = reference["cases"][observation_id] if "cases" in reference else reference
    p = np.asarray(r["event_mass"], float)
    if p.shape != (12,) or not np.isfinite(p).all() or np.any(p < 0) or abs(p.sum() - 1) > 1e-10:
        raise ValueError("a normalized twelve-event reference is required")
    return r, p


def _posterior_statistics(logu, ids, rng, bootstrap_replicates):
    """SNIS ratio statistics and iid-row bootstrap, with no weight clipping."""
    n = len(logu)
    if not np.isfinite(logu).all():
        raise FloatingPointError("nonfinite full-likelihood importance weight")
    log_total = logsumexp(logu)
    weights = np.exp(logu - log_total)
    p = np.bincount(ids, weights=weights, minlength=12)
    indicator = ids[:, None] == np.arange(12)
    se = np.sqrt(n / (n - 1) * np.sum(weights[:, None] ** 2 * (indicator - p) ** 2, axis=0))
    if (not np.isfinite(p).all() or np.any(p < 0) or abs(p.sum() - 1) > 1e-10
            or not np.isfinite(se).all()):
        raise FloatingPointError("invalid importance posterior mass; no clipping or renormalization repair allowed")
    if isinstance(bootstrap_replicates, bool) or bootstrap_replicates < 2:
        raise ValueError("at least two bootstrap replicates required")
    rows = rng.integers(n, size=(bootstrap_replicates, n))
    b_logu = logu[rows]
    b_logtotal = logsumexp(b_logu, axis=1)
    b_w = np.exp(b_logu - b_logtotal[:, None])
    b_p = np.sum(b_w[:, :, None] * (ids[rows, None] == np.arange(12)), axis=1)
    # The normalization is across all twelve declared events, never Top-K.
    return dict(event_mass=p, event_mass_delta_method_SE=se,
                event_mass_bootstrap_SE=np.std(b_p, axis=0, ddof=1),
                event_mass_interval95=np.quantile(b_p, [.025, .975], axis=0).T,
                ESS=float(1 / np.sum(weights ** 2)), maximum_normalized_weight=float(weights.max()),
                log_evidence=float(log_total - np.log(n)),
                log_evidence_bootstrap_interval95=np.quantile(b_logtotal - np.log(n), [.025, .975]),
                evidence_relative_delta_method_SE=float(np.sqrt(np.sum((n * weights - 1) ** 2) / (n * (n - 1)))),
                bootstrap_replicates=bootstrap_replicates,
                underflowed_normalized_weight_count=int((weights == 0).sum())), weights, b_p


def _cost(ledger, before, requested, attempted, call_wall):
    rows = ledger.rows()
    begin = [row for row in rows if row["event"] == "BEGIN" and row["id"] > before]
    end = [row for row in rows if row["event"] == "END" and row["id"] > before]
    return dict(requested_material_frequency_evaluations=requested,
                forward_attempts=attempted, new_material_frequency_evaluations=len(begin),
                physical_ledger_ids=[row["id"] for row in begin],
                failed_evaluations=sum(row["status"] == "FAIL" for row in end),
                unfinished_evaluations=len(begin) - len(end),
                full_forward_wall_seconds=call_wall,
                kernel_wall_seconds=sum(row["wall_seconds"] for row in end),
                boundary_RHS=sum(row.get("boundary_RHS", 0) for row in end),
                cache_enabled=False, stage=ledger.stage,
                accounting="all predict calls use provided reserve ledger; setup, failures, repeats and RHS retained")


def correct_conditioned_proposal(forward, parameters, density, y, sigma, reference,
                                 *, method, sample_seed, bootstrap_seed,
                                 observation_id="large_cavity", n=64,
                                 bootstrap_replicates=2000):
    """One proposal correction; no I/O except the provided paid predict calls.

    Returns (summary, arrays).  A failed full call returns an unresolved summary
    and the raw proposal/partial predictions; failed draws are never dropped,
    given zero weight, retried, or replaced.  The parent records that outcome.
    """
    condition = measurement_condition(y, sigma)
    y = condition[0, :32]
    sigma = np.broadcast_to(np.asarray(sigma, float), (32,))
    ledger = forward.backend.ledger
    if ledger.stage != "reserve":
        raise ValueError("caller must select the existing reserve stage before correction")
    before = ledger.count()
    rng = np.random.default_rng(sample_seed)
    draw = sample_defensive(parameters, density, rng, n)
    raw = _raw_summary(draw)
    F_full = np.full((n, 32), np.nan)
    failures, attempted = [], 0
    call_wall = 0.
    for i, x in enumerate(draw["latent"]):
        start = time.perf_counter()
        attempted += 1
        try:
            result = forward.predict(x, 1.6, kind="C1_PROPOSAL_FULL_IS_" + method.upper(),
                                     nmax=14, implementation="matching", cache=False)
            F_full[i] = np.asarray(result[0])
            if not np.isfinite(F_full[i]).all():
                raise FloatingPointError("nonfinite full forward prediction")
        except Exception as exc:
            failures.append(dict(sample_index=i, type=type(exc).__name__, message=str(exc)))
            break
        finally:
            call_wall += time.perf_counter() - start
    summary = dict(method=method, observation_id=observation_id, target_identity=FULL_TARGET,
                   proposal_law="0.9 normalized bounded C1 MDN dx density +0.1 proper uniform latent prior",
                   n_requested=n, sample_seed=sample_seed, bootstrap_seed=bootstrap_seed,
                   prior_bounds=BOUNDS, gain_prior=dict(type="proper truncated Gaussian", mean=1.,
                                                      sd=.15, bounds=[.7, 1.3], shared_across_channels=True),
                   noise_convention="independent real Re/Im packed channels with given sigma^2",
                   raw_proposal=raw,
                   physical_cost=_cost(ledger, before, n, attempted, call_wall),
                   failures=failures)
    arrays = {**draw, "F_full": F_full, "y": y, "sigma": sigma}
    if failures:
        summary.update(status="UNRESOLVED_FULL_FORWARD_FAILURE", event_mass=None,
                       probability_status="No posterior estimate: failed/unevaluated proposal draws retained without replacement")
        return summary, arrays
    lf = np.asarray(log_likelihood(F_full, y, sigma, gain_sd=.15, gain_bounds=(.7, 1.3)))
    logu = lf + draw["log_p0"] - draw["log_q"]
    stats, weights, b_p = _posterior_statistics(logu, draw["event_ids"],
                                               np.random.default_rng(bootstrap_seed), bootstrap_replicates)
    r, reference_p = _reference_case(reference, observation_id)
    p = stats["event_mass"]
    significant = reference_p >= SIGNIFICANT_MASS
    observed = raw["defensive_proposal_event_count"] > 0
    gain_mean, gain_variance = conditional_gain(F_full, y, sigma, .15)
    gain_sd = np.sqrt(gain_variance)
    gain_a, gain_b = (.7 - gain_mean) / gain_sd, (1.3 - gain_mean) / gain_sd
    tmean, tvariance = truncnorm.stats(gain_a, gain_b, loc=gain_mean, scale=gain_sd, moments="mv")
    if not np.isfinite(tmean).all() or not np.isfinite(tvariance).all() or np.any(tvariance < 0):
        raise FloatingPointError("invalid truncated conditional gain moments")
    summary.update(status="FULL_LIKELIHOOD_CORRECTED_FINITE_SAMPLE", **stats,
                   comparison_to_full_Mie_reference=dict(
                       reference_event_mass=reference_p,
                       reference_event_mass_MC_SE=r.get("event_mass_MC_se"),
                       reference_interval95=r.get("event_mass_interval95"),
                       reference_numerical_status=r.get("numerical_status"),
                       corrected_event_MAE=float(np.mean(np.abs(p - reference_p))),
                       corrected_event_max_abs_error=float(np.max(np.abs(p - reference_p))),
                       raw_defensive_event_MAE=float(np.mean(np.abs(raw["defensive_proposal_event_mass"] - reference_p))),
                       raw_MDN_branch_event_MAE=float(np.mean(np.abs(raw["mdn_branch_event_mass"] - reference_p))) if raw["mdn_branch_event_mass"] is not None else None,
                       event_mass_difference=p - reference_p,
                       significant_event_ids=np.flatnonzero(significant),
                       sample_discovery_recall=float(observed[significant].mean()) if significant.any() else None,
                       mass_aware_recall=float((p[significant] >= MASS_AWARE_DISCOVERY).mean()) if significant.any() else None,
                       missed_significant_event_ids=np.flatnonzero(significant & ~observed),
                       false_exclusion_event_ids=np.flatnonzero(significant & (p <= FALSE_EXCLUSION_MASS)),
                       false_confidence=bool(p.max() >= .9 and significant.sum() >= 2),
                       log_evidence_difference=float(stats["log_evidence"] - r["log_evidence"])),
                   error_scope="iid importance ratio/bootstrap error conditional on frozen proposal and numerical full-Mie target; reference numerical/MC uncertainty separate",
                   finite_sample_caveats=["SNIS event masses are biased at finite n and consistent under integrability, not exact finite-sample probabilities",
                                          "ESS and defensive support do not certify discovery of rare/high-likelihood modes",
                                          "unvisited events yield zero empirical mass and potentially zero bootstrap intervals; this does not certify exclusion",
                                          "importance correction targets matching-n14 full-Mie likelihood, not continuum Maxwell or model adequacy",
                                          "no Top-K renormalization; all twelve disjoint modeled events are retained"],
                   outside_model_support_mass=None,
                   unexplored_mass_status="Unvisited modeled regions can retain unknown missed posterior mass; finite proposals provide no bound",
                   observation_exposure="preregistered development large_cavity/C, reused diagnostic; not fresh independent validation")
    arrays.update(log_likelihood_full=lf, log_importance_weight=logu, weights=weights,
                  bootstrap_event_mass=b_p, gain_mean_before_truncation=gain_mean,
                  gain_variance_before_truncation=gain_variance,
                  gain_conditional_mean=tmean, gain_conditional_variance=tvariance,
                  event_mass=p, event_mass_interval95=stats["event_mass_interval95"])
    return summary, arrays


def run_proposal_correction(forward, model_dir, y, sigma, reference, *, n=64,
                            sample_seed=DEFAULT_SAMPLE_SEED, bootstrap_seed=DEFAULT_BOOTSTRAP_SEED,
                            bootstrap_replicates=2000, observation_id="large_cavity"):
    """Three matched frozen proposals; caller meters/saves returned outputs.

    ``reference`` may be the loaded FULL_MIE_CORRECTED_REFERENCE.json or its
    one-case dictionary. Arrays are returned as {method: array_dictionary}.
    On an explicit full-call failure, later methods are NOT_RUN and no failed
    proposal is replaced. No density model is trained here.
    """
    condition = measurement_condition(y, sigma)
    _reference_case(reference, observation_id)
    ledger = forward.backend.ledger
    budget = ledger.summary()
    if ledger.stage != "reserve" or budget["remaining"] < n * len(METHODS) or 800 - budget["by_stage"]["reserve"] < n * len(METHODS):
        raise ValueError("matched three-proposal full correction does not fit existing reserve/total budget")
    start, cpu = time.perf_counter(), time.process_time()
    before = ledger.count()
    methods, arrays = {}, {}
    for i, method in enumerate(METHODS):
        load_start = time.perf_counter()
        path = Path(model_dir) / f"{method}-{MODEL_SEED}.pt"
        parameters, metadata, density = load_conditioned_proposal(path, condition)
        if metadata["method"] != method:
            raise ValueError("model filename and saved method differ")
        load_wall = time.perf_counter() - load_start
        row, raw = correct_conditioned_proposal(forward, parameters, density, y, sigma, reference,
                                                method=method, n=n, sample_seed=sample_seed + i,
                                                bootstrap_seed=bootstrap_seed + i,
                                                bootstrap_replicates=bootstrap_replicates,
                                                observation_id=observation_id)
        row.update(model_metadata=metadata, model_path=str(path), model_load_inference_wall_seconds=load_wall)
        methods[method], arrays[method] = row, raw
        if row["failures"]:
            for remaining in METHODS[i + 1:]:
                methods[remaining] = dict(status="NOT_RUN", reason="earlier full-forward failure retained", n_requested=n)
            break
    return dict(status="COMPLETE" if all(methods[m].get("status") == "FULL_LIKELIHOOD_CORRECTED_FINITE_SAMPLE" for m in METHODS) else "UNRESOLVED_PARTIAL_FAILURE",
                observation_id=observation_id, target_identity=FULL_TARGET,
                proposal_training_target_identity=PROPOSAL_TARGET, methods=methods,
                matched_budget_per_method=n, requested_total_physics=n * len(METHODS),
                actual_total_physics=ledger.count() - before,
                wall_seconds=time.perf_counter() - start, cpu_seconds=time.process_time() - cpu,
                GPU_training=0,
                setup_cost_disclosure="frozen C1 training/teacher and paid B1 bank are existing amortized costs; this summary records new correction cost only",
                physical_cost_comparison_scope="three C1 models matched at64 full calls; generic SMC has its separately disclosed actual cost"), arrays


def time_proposals(model_dir, y, sigma, *, repetitions=5, n=64, seed=261013306):
    """Rotated CPU timing; no physical calls, writes, retraining, or cache flush.

    Cold means each model is loaded and conditioned anew, not a flushed OS disk
    cache. Warm reuses the conditioned export, but regenerates all draws and
    exact defensive-mixture log probabilities. Imports are separately timed.
    The caller must meter this helper as a zero-new-physics controlled job.
    """
    if repetitions != 5:
        raise ValueError("matched timing protocol uses five rotations")
    condition = measurement_condition(y, sigma)
    import_start = time.perf_counter()
    _api()
    import_wall = time.perf_counter() - import_start
    rows = []
    for repetition in range(repetitions):
        order = METHODS[repetition % 3:] + METHODS[:repetition % 3]
        for position, method in enumerate(order):
            path = Path(model_dir) / f"{method}-{MODEL_SEED}.pt"
            start, cpu = time.perf_counter(), time.process_time()
            parameters, metadata, density = load_conditioned_proposal(path, condition)
            load_infer = time.perf_counter() - start
            sample_defensive(parameters, density, np.random.default_rng(seed + repetition * 3 + METHODS.index(method)), n)
            cold = time.perf_counter() - start
            cold_cpu = time.process_time() - cpu
            warm_start, warm_cpu_start = time.perf_counter(), time.process_time()
            sample_defensive(parameters, density, np.random.default_rng(seed + 100 + repetition * 3 + METHODS.index(method)), n)
            rows.append(dict(repetition=repetition, position=position, method=method,
                             order=list(order), n=n, model_seed=MODEL_SEED,
                             cold_load_inference_seconds=load_infer,
                             cold_load_infer_sample_logq_seconds=cold, cold_cpu_seconds=cold_cpu,
                             warm_sample_logq_seconds=time.perf_counter() - warm_start,
                             warm_cpu_seconds=time.process_time() - warm_cpu_start,
                             categorical_log_normalizer_removed=metadata["categorical_log_normalizer_removed"]))
    medians = {method: {key: float(np.median([row[key] for row in rows if row["method"] == method]))
                        for key in ("cold_load_inference_seconds", "cold_load_infer_sample_logq_seconds", "warm_sample_logq_seconds")}
               for method in METHODS}
    return dict(observation_id="large_cavity", methods=list(METHODS), model_seed=MODEL_SEED,
                repetitions=repetitions, sample_count=n, rows=rows, medians=medians,
                initial_C1_api_import_seconds=import_wall, new_material_frequency_evaluations=0,
                cold_definition="model reload + condition export + new64 defensive draws/exact logq; OS cache not flushed",
                warm_definition="conditioned density reused, fresh64 defensive draws and exact logq recomputed",
                timing_scope="proposal acquisition only; no likelihood correction or Maxwell in these timings")
