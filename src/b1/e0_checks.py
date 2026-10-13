"""Solver-free probability controls for the B1 E0 audit.

``run_checks`` reads supplied algebra and saved, already exposed arrays. It
never writes a source artifact, constructs a physical solver or ledger, or
evaluates a new material response. Numerical checks here are identities and
counterexamples; none is a Maxwell certificate or a new scientific gate.
"""
from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import math
from pathlib import Path
import warnings

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import logsumexp

from .inference import infer, log_likelihood, quadrature_nodes, top_k


SEED = 20261013
_ROUND_OFF = 512 * np.finfo(float).eps
_BOUNDS = np.array([[2.2, 2.8], [.04, .16]])


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {("pass" if k == "passed" else str(k)): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _warning_records(records):
    return [dict(category=w.category.__name__, message=str(w.message),
                 filename=w.filename, lineno=w.lineno) for w in records]


def _close(a, b):
    """Floating arithmetic consistency, not a posterior-accuracy tolerance."""
    return bool(np.allclose(a, b, rtol=_ROUND_OFF, atol=_ROUND_OFF))


def _capture(fn):
    stdout, stderr = io.StringIO(), io.StringIO()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                result = fn()
                error = None
            except Exception as exc:
                result = None
                error = dict(type=type(exc).__name__, message=str(exc))
    return result, _warning_records(caught), stdout.getvalue(), stderr.getvalue(), error


def _tree_differences(actual, saved, path=""):
    """Record platform-dependent replay differences without adding a gate."""
    rows = []
    if isinstance(actual, dict) and isinstance(saved, dict):
        for key in sorted(set(actual) | set(saved)):
            name = f"{path}.{key}" if path else str(key)
            if key not in actual or key not in saved:
                rows.append(dict(path=name, kind="missing_key"))
            else:
                rows.extend(_tree_differences(actual[key], saved[key], name))
    elif isinstance(actual, (list, tuple)) and isinstance(saved, (list, tuple)):
        if len(actual) != len(saved):
            rows.append(dict(path=path, kind="length", actual=len(actual), saved=len(saved)))
        for i, (a, b) in enumerate(zip(actual, saved)):
            rows.extend(_tree_differences(a, b, f"{path}[{i}]"))
    elif actual != saved:
        if (isinstance(actual, (int, float, np.number)) and
                isinstance(saved, (int, float, np.number)) and
                not isinstance(actual, (bool, np.bool_)) and
                not isinstance(saved, (bool, np.bool_))):
            rows.append(dict(path=path, kind="numeric", actual=float(actual),
                             saved=float(saved), absolute_difference=abs(float(actual-saved))))
        else:
            rows.append(dict(path=path, kind="value", actual=actual, saved=saved))
    return rows


def _supplied_theory(root):
    source = root / "inputs/b1_next_2026_10_13/verify_theory.py"
    saved_path = source.with_name("checks.json")
    spec = importlib.util.spec_from_file_location("b1_e0_supplied_verify_theory", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load supplied algebra file: {source}")
    module = importlib.util.module_from_spec(spec)
    # The reviewed file's __main__ block writes output; importing and calling
    # main() leaves that block inactive and only returns its algebra results.
    spec.loader.exec_module(module)
    result, recorded, stdout, stderr, error = _capture(module.main)
    saved = json.loads(saved_path.read_text(encoding="utf-8"))
    return dict(source=str(source.relative_to(root)), saved_source=str(saved_path.relative_to(root)),
                result=result, saved_checks_snapshot=saved, warnings=recorded,
                stdout=stdout, stderr=stderr, error=error,
                differences_from_saved_snapshot=[] if result is None else _tree_differences(result, saved),
                passed=bool(error is None and result["all_asserted_checks_pass"]),
                evidence_type="supplied_algebra_and_saved_number_arithmetic",
                interpretation="Supplied finite-dimensional identities and counterexamples; no physics evaluation.")


def _gain_dense_cholesky():
    rng = np.random.default_rng(SEED)
    rows = []
    for case in range(40):
        m = 32
        f, y = rng.normal(size=m), rng.normal(size=m)
        sigma, gain_sd = 10**rng.uniform(-1.5, .2), rng.uniform(0., .4)
        r = y-f
        covariance = sigma**2*np.eye(m) + gain_sd**2*np.outer(f, f)
        sign, dense_ld = np.linalg.slogdet(covariance)
        dense_q = r @ np.linalg.solve(covariance, r)
        cholesky = np.linalg.cholesky(covariance)
        whitened = np.linalg.solve(cholesky, r)
        cholesky_q = whitened @ whitened
        cholesky_ld = 2*np.log(np.diag(cholesky)).sum()
        dense = -.5*(dense_q+dense_ld+m*np.log(2*np.pi))
        chol = -.5*(cholesky_q+cholesky_ld+m*np.log(2*np.pi))
        fast, fast_q = log_likelihood(f, y, sigma, gain_sd)
        rows.append(dict(case=case, dimension=m, sigma=float(sigma), gain_sd=float(gain_sd),
                         covariance_determinant_sign=float(sign), finite=bool(np.all(np.isfinite(
                             [dense, chol, fast, dense_q, cholesky_q, fast_q, dense_ld, cholesky_ld]))),
                         dense_loglik=float(dense), cholesky_loglik=float(chol), rank_one_loglik=float(fast),
                         dense_logdet=float(dense_ld), cholesky_logdet=float(cholesky_ld),
                         fast_vs_dense_absolute_error=float(abs(fast-dense)),
                         fast_vs_cholesky_absolute_error=float(abs(fast-chol)),
                         dense_vs_cholesky_absolute_error=float(abs(dense-chol))))
    tolerance = 1e-7  # Retained from the supplied check, not an E0 scientific gate.
    errors = [max(r["fast_vs_dense_absolute_error"], r["fast_vs_cholesky_absolute_error"],
                  r["dense_vs_cholesky_absolute_error"]) for r in rows]
    return dict(evidence_type="analytic_random_matrices", seed=SEED, cases=rows, case_count=40,
                max_absolute_loglik_error=max(errors), supplied_identity_tolerance=tolerance,
                passed=bool(all(r["finite"] and r["covariance_determinant_sign"] > 0 for r in rows)
                          and max(errors) < tolerance),
                interpretation="The same 40 gain cases are replayed with dense and Cholesky cross-checks; determinant retained.")


def _missing_mode():
    width, height, n = 1e-4, 1e4, 1024
    samples = (np.arange(n)+.5)/n
    hit = samples < width
    likelihood = np.where(hit, height, 1.)
    weights = likelihood/likelihood.sum()
    exact = width*height/(width*height+1-width)
    estimate = weights[hit].sum()
    ess = 1/np.sum(weights**2)
    return dict(evidence_type="analytic_toy_counterexample", prior="uniform on [0,1]",
                likelihood="10000 on [0,0.0001); 1 elsewhere",
                proposal="uniform prior; deterministic midpoint realization with no rare-event hit",
                sample_count=n, hits=int(hit.sum()), true_event_mass=float(exact),
                estimated_event_mass=float(estimate), absolute_event_error=float(abs(exact-estimate)),
                ESS=float(ess), ESS_fraction=float(ess/n),
                iid_no_hit_probability=float((1-width)**n),
                iid_probability_scope="Probability for hypothetical independent uniform draws, not this deterministic realization.",
                passed=bool(not hit.any() and ess == n and estimate == 0 and exact > .49),
                interpretation="Maximum ESS within visited support coexists with a missed event carrying about half the true posterior mass.")


def _toy_bank(order=5):
    config = dict(material_prior=dict(eps_real=_BOUNDS[0].tolist(), eps_imag=_BOUNDS[1].tolist()))
    theta, weights = quadrature_nodes(order, config)
    return dict(H=np.repeat(np.arange(3), len(theta)), theta=np.tile(theta, (3, 1)),
                weights=np.tile(weights, 3), bounds=_BOUNDS.copy(), F=np.zeros((3*len(theta), 4)))


def _prior_only():
    bank = _toy_bank()
    prior = np.array([.52, .33, .15])
    result = infer(bank, np.ones(4), 1., prior)
    mean = _BOUNDS.mean(axis=1)
    covariance = np.diag(np.diff(_BOUNDS, axis=1).ravel()**2/12)
    quantiles = _BOUNDS[:, :1] + np.diff(_BOUNDS, axis=1)*np.array([.05, .5, .95])
    posterior_node_prior = prior[bank["H"]]*bank["weights"]
    summaries_ok = all(_close(row["mean"], mean) and _close(row["covariance"], covariance) and
                       _close(row["quantiles_05_50_95"], quantiles)
                       for row in result["conditional_material"])
    return dict(evidence_type="analytic_toy", likelihood="constant across topology and material",
                topology_prior=prior.tolist(), topology_posterior=result["topology_probabilities"],
                normalized_prior_weight_error=float(np.max(abs(result["weights"]-posterior_node_prior))),
                conditional_prior_mean=mean.tolist(), conditional_prior_covariance=covariance.tolist(),
                conditional_prior_quantiles=quantiles.tolist(), conditional_summaries=result["conditional_material"],
                passed=bool(_close(result["topology_probabilities"], prior) and
                          _close(result["weights"], posterior_node_prior) and summaries_ok),
                interpretation="An uninformative synthetic likelihood restores the declared proper prior.")


def _coordinate_transform():
    # u=r^3 is a coordinate change on r in [1,2], away from the endpoint
    # singularity. The transformed prior includes dr/du, so its measure is
    # still radius-uniform. Uniform-u is an explicitly different prior.
    lower, upper, cut = 1., 2., 1.5
    likelihood = lambda r: 1+r*r
    joint_u = lambda u: (1+u**(2/3))/(3*u**(2/3)*(upper-lower))
    z_exact = 1+(upper**3-lower**3)/(3*(upper-lower))
    event_exact = ((cut-lower)+(cut**3-lower**3)/3)/((upper-lower)*z_exact)
    mean_exact = ((upper**2-lower**2)/2+(upper**4-lower**4)/4)/((upper-lower)*z_exact)
    z_r, z_r_error = quad(likelihood, lower, upper, epsabs=1e-12, epsrel=1e-12)
    z_u, z_u_error = quad(joint_u, lower**3, upper**3, epsabs=1e-12, epsrel=1e-12)
    event_r = quad(likelihood, lower, cut, epsabs=1e-12, epsrel=1e-12)[0]/z_r
    event_u = quad(joint_u, lower**3, cut**3, epsabs=1e-12, epsrel=1e-12)[0]/z_u
    mean_u = quad(lambda u: u**(1/3)*joint_u(u), lower**3, upper**3,
                  epsabs=1e-12, epsrel=1e-12)[0]/z_u
    changed_joint = lambda u: (1+u**(2/3))/(upper**3-lower**3)
    changed_z = quad(changed_joint, lower**3, upper**3, epsabs=1e-12, epsrel=1e-12)[0]
    changed_event = quad(changed_joint, lower**3, cut**3, epsabs=1e-12, epsrel=1e-12)[0]/changed_z
    return dict(evidence_type="analytic_toy", original_prior="density 1 on r in [1,2]",
                coordinate_transform="u=r^3", transformed_prior="1/(3*u^(2/3)) on u in [1,8]",
                likelihood="L(r)=1+r^2", event="r<=1.5, equivalently u<=3.375",
                analytic_evidence=z_exact, original_coordinate_evidence=z_r, transformed_coordinate_evidence=z_u,
                integration_error_estimates=[z_r_error, z_u_error],
                analytic_event_probability=event_exact, original_event_probability=event_r,
                transformed_event_probability=event_u, analytic_radius_mean=mean_exact,
                transformed_radius_mean=mean_u,
                negative_control=dict(prior="uniform u on [1,8], a different measure",
                                      evidence=changed_z, event_probability=changed_event,
                                      absolute_event_change=abs(changed_event-event_exact)),
                passed=bool(_close([z_r, z_u, event_r, event_u, mean_u],
                                 [z_exact, z_exact, event_exact, event_exact, mean_exact]) and
                          not _close(changed_event, event_exact)),
                interpretation="Keeping the Jacobian preserves one probability measure; replacing it by a uniform transformed-coordinate prior changes the measure.")


def _radius_volume_priors():
    alpha, cut = 4., .6
    rows = []
    for power, name in [(0, "radius_uniform"), (2, "volume_uniform_3D")]:
        normalization = power+1
        z = 1+normalization*alpha/(power+3)
        event = ((1-cut**(power+1))+normalization*alpha*(1-cut**(power+3))/(power+3))/z
        first = normalization*(1/(power+2)+alpha/(power+4))/z
        second = normalization*(1/(power+3)+alpha/(power+5))/z
        integrand = lambda r: normalization*r**power*(1+alpha*r*r)
        quad_z, quad_error = quad(integrand, 0., 1., epsabs=1e-12, epsrel=1e-12)
        quad_event = quad(integrand, cut, 1., epsabs=1e-12, epsrel=1e-12)[0]/quad_z
        rows.append(dict(prior=name, prior_density=f"{normalization}*r^{power} on [0,1]",
                         prior_event_mass=1-cut**(power+1), evidence=z, log_evidence=math.log(z),
                         posterior_event_probability=event, posterior_radius_mean=first,
                         posterior_radius_variance=second-first*first,
                         quadrature_evidence=quad_z, quadrature_event_probability=quad_event,
                         quadrature_error_estimate=quad_error,
                         analytic_quadrature_agreement=_close([quad_z, quad_event], [z, event])))
    return dict(evidence_type="analytic_toy_prior_sensitivity", support="radius r in [0,1]",
                likelihood="L(r)=1+4*r^2 (positive analytic toy likelihood, no physical response)",
                event="r>=0.6", volume="V proportional to r^3; uniform V induces density 3*r^2",
                priors=rows, posterior_event_change=rows[1]["posterior_event_probability"]-rows[0]["posterior_event_probability"],
                passed=bool(all(r["analytic_quadrature_agreement"] for r in rows) and
                          not _close(rows[0]["posterior_event_probability"], rows[1]["posterior_event_probability"])),
                interpretation="This compares genuinely different radius and volume priors, not two coordinate representations of one prior.")


def _event_partition():
    bank = _toy_bank()
    prior = np.array([.5, .25, .25])
    result = infer(bank, np.zeros(4), 1., prior)
    report = top_k(result, bank, 2)
    masks = np.array([bank["H"] == h for h in range(3)])
    probe = np.array([-.01, 0., .25, .5, .75, 1., 1.01])
    labels = np.where((probe < 0) | (probe > 1), -1,
                      np.where(probe < .5, 0, np.where(probe < .75, 1, 2)))
    displayed = sum(row["posterior_mass"] for row in report["items"])
    remainder = report["remaining_modeled_mass"]
    boundary_ok = np.array_equal(labels, [-1, 0, 0, 1, 2, 2, -1])
    semantics_ok = (report["unmodeled_geometry_mass"] is None and
                    report["outside_adopted_prior_support_mass"] == 0. and
                    report["exact_image_probability"] == 0.)
    return dict(evidence_type="analytic_toy_event_semantics", events=["H=0", "H=1", "H=2"],
                disjoint_exhaustive_on_modeled_support=bool(np.all(masks.sum(axis=0) == 1)),
                boundary_rule_toy=dict(support="[0,1]", events=["[0,0.5)", "[0.5,0.75)", "[0.75,1]"],
                                       probes=probe.tolist(), labels=labels.tolist(), outside_label=-1),
                tie_rule="descending mass, stable original event order", top_k=report,
                displayed_mass=displayed, remaining_modeled_mass=remainder,
                total_modeled_mass=displayed+remainder,
                passed=bool(np.all(masks.sum(axis=0) == 1) and boundary_ok and semantics_ok and
                          [r["H"] for r in report["items"]] == [0, 1] and
                          _close([displayed, remainder, displayed+remainder], [.75, .25, 1.])),
                interpretation="Displayed events retain their unconditional modeled masses. The omitted modeled event has quantified mass; a broader unspecified physical model has unquantified mass. Support exclusion is a prior choice.")


def _tensor_prior_weights(theta, h_values, bounds):
    """Recover normalized GL weights for saved tensor grids, without physics."""
    weights = np.empty(len(h_values), dtype=float)
    rows = []
    for h in np.unique(h_values):
        ids = np.flatnonzero(h_values == h)
        x, ix = np.unique(theta[ids, 0], return_inverse=True)
        y, iy = np.unique(theta[ids, 1], return_inverse=True)
        codes = ix*len(y)+iy
        tensor_complete = len(ids) == len(x)*len(y) and len(np.unique(codes)) == len(ids)
        zx, wx = leggauss(len(x))
        zy, wy = leggauss(len(y))
        expected_x = bounds[0, 0]+(zx+1)/2*(bounds[0, 1]-bounds[0, 0])
        expected_y = bounds[1, 0]+(zy+1)/2*(bounds[1, 1]-bounds[1, 0])
        gauss_grid = _close(x, expected_x) and _close(y, expected_y)
        if not tensor_complete or not gauss_grid:
            raise ValueError(f"H={h}: saved theta is not the declared tensor Gauss-Legendre grid")
        weights[ids] = (wx[ix]/2)*(wy[iy]/2)
        rows.append(dict(H=int(h), orders=[len(x), len(y)], nodes=len(ids),
                         tensor_complete=tensor_complete, coordinate_grid_verified=gauss_grid,
                         normalized_prior_weight_sum=float(weights[ids].sum())))
    return weights, rows


def _posterior_from_loglik(loglik, h_values, weights, prior):
    evidence = np.array([logsumexp(loglik[h_values == h]+np.log(weights[h_values == h]))
                         for h in range(len(prior))])
    total = logsumexp(np.log(prior)+evidence)
    joint = np.exp(loglik+np.log(weights)+np.log(prior[h_values])-total)
    probabilities = np.exp(np.log(prior)+evidence-total)
    return probabilities, joint, evidence, float(total)


def _saved_gain_replay(root):
    config = json.loads((root/"configs/protocol.json").read_text(encoding="utf-8"))
    bounds = np.array([config["material_prior"]["eps_real"], config["material_prior"]["eps_imag"]])
    rows = []
    for case in ["gain_marginalized", "refined_gain_marginalized"]:
        source = root/f"results/cases/{case}/posterior.npz"
        probability_source = source.with_name("probabilities.json")
        saved_report = json.loads(probability_source.read_text(encoding="utf-8"))
        with np.load(source, allow_pickle=False) as archive:
            arrays = {key: archive[key] for key in archive.files}
        if "gain_sd" not in arrays:
            raise ValueError(f"Corrected saved case lacks gain_sd: {source.relative_to(root)}")
        f, y, h_values, theta = arrays["forward"], arrays["data"], arrays["H"], arrays["theta"]
        sigma, gain_sd, prior = float(arrays["noise_real_sd"]), float(arrays["gain_sd"]), arrays["topology_prior"]
        if sigma <= 0 or gain_sd < 0 or not np.all(np.isfinite([sigma, gain_sd])):
            raise ValueError(f"Invalid saved likelihood metadata: {case}")
        weights, grid = _tensor_prior_weights(theta, h_values, bounds)
        loglik, quadratic = log_likelihood(f, y, sigma, gain_sd)
        p, joint, evidence, total = _posterior_from_loglik(loglik, h_values, weights, prior)
        candidate_logdet = f.shape[1]*np.log(sigma*sigma)+np.log1p(gain_sd**2*np.sum(f*f, axis=1)/sigma**2)
        without_logdet = -.5*(quadratic+f.shape[1]*np.log(sigma*sigma)+f.shape[1]*np.log(2*np.pi))
        wrong_p = _posterior_from_loglik(without_logdet, h_values, weights, prior)[0]
        frozen = saved_report["result"]
        expected_gain = saved_report["likelihood"]["gain_sd"]
        metadata_match = (gain_sd == expected_gain and sigma == saved_report["likelihood"]["real_noise_sd"])
        probability_agreement = _close(p, frozen["topology_probabilities"])
        evidence_agreement = _close(evidence, frozen["log_evidence"]) and _close(total, frozen["log_total_evidence"])
        weight_agreement = _close(joint, arrays["posterior_weights"])
        history_source = root/f"results/history/{case}_posterior_INITIAL.npz"
        historical = dict(source=str(history_source.relative_to(root)), available=history_source.exists())
        if history_source.exists():
            with np.load(history_source, allow_pickle=False) as old:
                historical["missing_gain_sd"] = "gain_sd" not in old.files
                same_target_arrays = all(np.array_equal(old[key], arrays[key])
                                         for key in ["data", "noise_real_sd", "topology_prior", "forward", "theta", "H"])
                historical["same_saved_data_fields_and_prior"] = same_target_arrays
                # Quantify the old default's target change using only the
                # already saved arrays; do not import the CLI/forward stack.
                legacy_gain = float(old["gain_sd"]) if "gain_sd" in old.files else 0.
                old_ll, _ = log_likelihood(old["forward"], old["data"], float(old["noise_real_sd"]), legacy_gain)
                old_weights, _ = _tensor_prior_weights(old["theta"], old["H"], bounds)
                old_p = _posterior_from_loglik(old_ll, old["H"], old_weights, old["topology_prior"])[0]
                historical.update(unsafe_legacy_default_gain_sd=legacy_gain,
                                  probabilities_under_legacy_default=old_p.tolist(),
                                  maximum_probability_change_from_corrected_target=float(np.max(abs(old_p-p))),
                                  interpretation="Missing metadata made default fixed-gain replay a different likelihood target.")
        rows.append(dict(case=case, source=str(source.relative_to(root)), probability_source=str(probability_source.relative_to(root)),
                         gain_sd=gain_sd, noise_real_sd=sigma, metadata_matches_saved_probability_target=metadata_match,
                         prior_weight_reconstruction="normalized tensor Gauss-Legendre weights verified against declared prior bounds",
                         grid=grid, topology_probabilities=p.tolist(), saved_topology_probabilities=frozen["topology_probabilities"],
                         maximum_probability_replay_error=float(np.max(abs(p-np.asarray(frozen["topology_probabilities"])))),
                         maximum_posterior_weight_replay_error=float(np.max(abs(joint-arrays["posterior_weights"]))),
                         maximum_log_evidence_replay_error=float(np.max(abs(evidence-np.asarray(frozen["log_evidence"])))),
                         log_determinant_retained=True, candidate_logdet_range=[float(candidate_logdet.min()), float(candidate_logdet.max())],
                         omitted_candidate_logdet_negative_control=dict(topology_probabilities=wrong_p.tolist(),
                             maximum_probability_change=float(np.max(abs(wrong_p-p))),
                             interpretation="Omitting the candidate-dependent determinant defines a wrong likelihood target, not a reference improvement."),
                         historical_metadata=historical,
                         passed=bool(metadata_match and probability_agreement and weight_agreement and evidence_agreement and
                                   np.all(np.isfinite(loglik)) and _close(weights[h_values == 0].sum(), 1.))))
    cli_source = root/"results/cli_gain_regression_corrected/probabilities.json"
    cli = json.loads(cli_source.read_text(encoding="utf-8"))
    cli_p = cli["result"]["topology_probabilities"]
    refined_p = rows[1]["topology_probabilities"]
    cli_match = _close(cli_p, refined_p) and cli["measurement_target"]["gain_sd"] == rows[1]["gain_sd"]
    return dict(evidence_type="saved_array_replay", cases=rows,
                corrected_cli_replay=dict(source=str(cli_source.relative_to(root)), gain_sd=cli["measurement_target"]["gain_sd"],
                                          topology_probabilities=cli_p, matches_saved_refined_target=bool(cli_match)),
                passed=bool(all(row["passed"] for row in rows) and cli_match),
                interpretation="Replays normalization on exposed saved response arrays. Agreement checks reproducibility of that saved target, not its Maxwell accuracy or full posterior convergence.")


def _distance_to_interval(mean, lower, upper):
    return max(lower-mean, 0., mean-upper)


def _summary_coherence(root):
    config = json.loads((root/"configs/protocol.json").read_text(encoding="utf-8"))
    bounds = np.array([config["material_prior"]["eps_real"], config["material_prior"]["eps_imag"]])
    sources = ["results/PRIMARY_SUMMARY.json", "results/refinement/SUMMARY.json"]
    rows, coverage, structural_errors = [], [], []
    for source in sources:
        document = json.loads((root/source).read_text(encoding="utf-8"))
        source_cases = []
        for case, item in sorted(document["cases"].items()):
            summaries = item["result"]["conditional_material"]
            source_cases.append(case)
            h_seen = set()
            for summary in summaries:
                h = int(summary["H"])
                h_seen.add(h)
                means, cov, quantiles = map(np.asarray, (summary["mean"], summary["covariance"], summary["quantiles_05_50_95"]))
                if means.shape != (2,) or cov.shape != (2, 2) or quantiles.shape != (2, 3):
                    structural_errors.append(dict(source=source, case=case, H=h, error="Unexpected material summary dimensions"))
                    continue
                for j, coordinate in enumerate(["Re(epsilon)", "Im(epsilon)"]):
                    mean, variance = float(means[j]), float(cov[j, j])
                    q05, q50, q95 = map(float, quantiles[j])
                    lower, upper = bounds[j]
                    finite = bool(np.all(np.isfinite([mean, variance, q05, q50, q95])))
                    ordered = q05 <= q50 <= q95
                    tail_bound = .05*(max(mean-q05, 0.)**2+max(q95-mean, 0.)**2)
                    median_bound = (tail_bound+.45*_distance_to_interval(mean, q05, q50)**2+
                                    .45*_distance_to_interval(mean, q50, q95)**2)
                    slack = 1e-12*(upper-lower)**2
                    coordinate_slack = 1e-12*(upper-lower)
                    violations = []
                    if not finite:
                        violations.append("nonfinite_summary")
                    if not ordered:
                        violations.append("quantiles_not_ordered")
                    if variance < -slack:
                        violations.append("negative_variance")
                    if mean < lower-coordinate_slack or mean > upper+coordinate_slack:
                        violations.append("mean_outside_adopted_prior_support")
                    if q05 < lower-coordinate_slack or q95 > upper+coordinate_slack:
                        violations.append("quantile_outside_adopted_prior_support")
                    if finite and ordered and variance+slack < median_bound:
                        violations.append("variance_below_necessary_quantile_lower_bound")
                    support_upper = float((upper-mean)*(mean-lower))
                    if lower <= mean <= upper and variance > support_upper+slack:
                        violations.append("variance_above_bounded_support_upper_bound")
                    rows.append(dict(source=source, case=case, H=h, coordinate=coordinate, mean=mean,
                                     reported_variance=variance, quantiles_05_50_95=[q05, q50, q95],
                                     necessary_two_tail_variance_lower_bound=tail_bound,
                                     necessary_median_aware_variance_lower_bound=median_bound,
                                     bounded_support_variance_upper_bound=support_upper,
                                     two_tail_lower_bound_divided_by_reported_variance=None if variance <= 0 else tail_bound/variance,
                                     median_aware_lower_bound_divided_by_reported_variance=None if variance <= 0 else median_bound/variance,
                                     floating_roundoff_slack=slack, violations=violations,
                                     status="INCOMPATIBLE_NECESSARY_CONDITION" if violations else "NO_NECESSARY_CONDITION_VIOLATION_DETECTED"))
            if h_seen != {0, 1, 2} or len(summaries) != 3:
                structural_errors.append(dict(source=source, case=case, error="Expected one summary for each H=0,1,2"))
        coverage.append(dict(source=source, cases=source_cases, case_count=len(source_cases),
                             diagnostic_coordinate_rows=sum(row["source"] == source for row in rows)))
    native = {(r["case"], r["H"], r["coordinate"]): r for r in rows if r["source"] == sources[0]}
    refined = {(r["case"], r["H"], r["coordinate"]): r for r in rows if r["source"] == sources[1]}
    changes = []
    for key in sorted(native.keys() & refined.keys()):
        a, b = native[key], refined[key]
        ratio = None if a["reported_variance"] <= 0 else b["reported_variance"]/a["reported_variance"]
        changes.append(dict(case=key[0], H=key[1], coordinate=key[2],
                            original_variance=a["reported_variance"], refined_variance=b["reported_variance"],
                            refined_to_original_variance_ratio=ratio,
                            refined_to_original_sd_ratio=None if ratio is None or ratio < 0 else math.sqrt(ratio),
                            original_status=a["status"], refined_status=b["status"]))
    complete = (not structural_errors and native.keys() == refined.keys() and
                all(c["case_count"] == 6 and c["diagnostic_coordinate_rows"] == 36 for c in coverage))
    return dict(evidence_type="saved_summary_necessary_condition_diagnostic", coverage=coverage, rows=rows,
                structural_errors=structural_errors, complete_case_H_coordinate_coverage=complete,
                incompatible_row_count=sum(bool(row["violations"]) for row in rows),
                changes_between_saved_summaries=changes,
                necessary_condition="For exact quantiles of one law: Var(X)>=.05*(mu-q05)_+^2+.05*(q95-mu)_+^2. Median bins give additional nonnegative lower-bound terms.",
                quantile_convention="generalized-inverse CDF quantiles; the tests apply if the reported values are interpreted as exact quantiles",
                interpretation="Violations identify incompatible approximate moment/CDF summaries. No detected violation is not certification. Refinement differences are descriptive; neither saved summary is an exact reference. This is not a likelihood-normalization failure.")


def run_checks(root) -> dict:
    """Execute necessary E0 probability controls and return JSON-safe results.

    The caller owns timing, receipts and output files. This function is read
    only and has no physical-call allowance.
    """
    root = Path(root).resolve()
    checks, warnings_out, errors = {}, [], []
    jobs = [
        ("supplied_verify_theory", lambda: _supplied_theory(root)),
        ("gain_40_dense_cholesky", _gain_dense_cholesky),
        ("high_ESS_missing_mode_counterexample", _missing_mode),
        ("prior_only_posterior", _prior_only),
        ("same_measure_coordinate_transform", _coordinate_transform),
        ("radius_uniform_vs_volume_uniform", _radius_volume_priors),
        ("deterministic_event_partition_and_remainder", _event_partition),
        ("gain_metadata_saved_case_replay", lambda: _saved_gain_replay(root)),
        ("saved_complete_case_moment_quantile_coherence", lambda: _summary_coherence(root)),
    ]
    for name, job in jobs:
        result, caught, stdout, stderr, error = _capture(job)
        warnings_out.extend(dict(check=name, **warning) for warning in caught)
        if error is not None:
            checks[name] = dict(execution_status="ERROR", error=error, stdout=stdout, stderr=stderr)
            errors.append(dict(check=name, **error))
            continue
        result["execution_status"] = "EXECUTED"
        if stdout:
            result["captured_stdout"] = stdout
        if stderr:
            result["captured_stderr"] = stderr
        checks[name] = result
        if name == "supplied_verify_theory":
            warnings_out.extend(dict(check=name, **warning) for warning in result["warnings"])
            if result["error"] is not None:
                errors.append(dict(check=name, **result["error"]))
    diagnostic = checks.get("saved_complete_case_moment_quantile_coherence", {})
    all_asserted = (not errors and all(check.get("passed", True) for check in checks.values()) and
                    diagnostic.get("complete_case_H_coordinate_coverage", False))
    result = dict(date="2026-10-13", seed=SEED, new_physical_evaluations=0, neural_training=False,
                  source_artifacts_modified=False, old_locked_truth_read=False, whole_old_suite_executed=False,
                  scope="finite-dimensional algebra, analytic toy probabilities, and exposed saved-array/summary replay only",
                  checks=checks, warnings=warnings_out, errors=errors, all_asserted_checks_pass=bool(all_asserted),
                  status_meaning="Passing controls verifies identities, counterexample detection and saved-target replay; it supplies no Maxwell approximation, mixing, calibration, or scientific acceptance certificate.")
    result = _jsonable(result)
    json.dumps(result, allow_nan=False)
    return result
