"""Necessary solver-free probability controls for B1 x C1 Joint P2.

Only shared material algebra, gain likelihood, analytic probability examples,
and the exposed B1 E0 audit are read.  There is no physical forward or resource
ledger interface here.  Run under the parent's shared controlled-job meter.
JSON is returned on stdout; --report writes only the owned E0_AUDIT.md.
"""
from __future__ import annotations

import argparse
import importlib
import json
import math
from pathlib import Path
import sys
import time
import warnings

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm, truncnorm


ROOT = Path(__file__).resolve().parents[1]
SEED = 2026101302
CENTRE_SCALE, SHELL_SCALE = .45, .18
ROOT_FACTOR = math.sqrt(1 - 1 / math.sqrt(3))
RMAX, HMAX = CENTRE_SCALE * ROOT_FACTOR, SHELL_SCALE * ROOT_FACTOR
CENTRE_CUTS = [0., .10, .20, RMAX]
SHELL_CUTS = [0., .06, HMAX]


def _native(value):
    if isinstance(value, np.ndarray):
        return _native(value.tolist())
    if isinstance(value, np.generic):
        return _native(value.item())
    if isinstance(value, dict):
        return {str(k): _native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _a_of_radius(r, scale):
    return (1 - (np.asarray(r) / scale) ** 2) ** -2


def _radius(z, scale):
    z = np.asarray(z, float)
    return scale * np.sqrt(np.maximum(0., 1 - 1 / np.sqrt(np.maximum(-z, 1.))))


def _positive_radius_density(r, scale):
    """The original dz/4 measure, excluding its radius-zero atom."""
    return r / scale ** 2 * (1 - (r / scale) ** 2) ** -3


def _bin_masses(cuts, scale):
    return np.r_[.5, np.diff(_a_of_radius(np.asarray(cuts), scale)) / 4]


def _split_prior_nodes(cuts, scale, order=3):
    # Split at each material-event boundary.  Integrating a constant likelihood
    # in coefficient space is exact on each segment, without a radius-bin prior.
    zcuts = np.unique(np.r_[-_a_of_radius(cuts, scale), 1.])
    v, w = leggauss(order)
    x, weight = [], []
    for lo, hi in zip(zcuts[:-1], zcuts[1:]):
        x.extend((lo + (v + 1) * (hi - lo) / 2).tolist())
        weight.extend((w * (hi - lo) / 8).tolist())
    return np.asarray(x), np.asarray(weight)


def _shared_material(material):
    expected_bounds = np.array([[-3., 1.], [-3., 1.], [2.2, 2.8], [.04, .16]])
    expected = np.outer(_bin_masses(CENTRE_CUTS, CENTRE_SCALE),
                        _bin_masses(SHELL_CUTS, SHELL_SCALE)).ravel()
    c, cw = _split_prior_nodes(CENTRE_CUTS, CENTRE_SCALE)
    s, sw = _split_prior_nodes(SHELL_CUTS, SHELL_SCALE)
    xx, yy = np.meshgrid(c, s, indexing="ij")
    theta = np.column_stack([xx.ravel(), yy.ravel(),
                             np.full(xx.size, 2.5), np.full(xx.size, .1)])
    weights = (cw[:, None] * sw[None, :]).ravel()
    geometry = np.asarray(material.geometry(theta))
    ids = np.asarray(material.event_ids(theta), int)
    independent_geometry = np.column_stack([_radius(theta[:, 0], CENTRE_SCALE),
                                            _radius(theta[:, 1], SHELL_SCALE)])
    independent_ids = (np.where(independent_geometry[:, 0] == 0., 0,
                               np.searchsorted(CENTRE_CUTS[1:-1],
                                               independent_geometry[:, 0], side="left") + 1) * 3
                       + np.where(independent_geometry[:, 1] == 0., 0,
                                  np.searchsorted(SHELL_CUTS[1:-1],
                                                  independent_geometry[:, 1], side="left") + 1))
    actual = np.asarray(material.prior_event_masses(), float)
    integrated = np.bincount(ids, weights=weights, minlength=12)
    probes = np.array([[1, 1, 2.5, .1], [-1, -1, 2.5, .1],
                       [-3, 1, 2.5, .1], [1, -3, 2.5, .1], [-3, -3, 2.5, .1]])
    probe_geometry = np.asarray(material.geometry(probes))
    probe_ids = np.asarray(material.event_ids(probes), int)
    matches = np.array_equal(ids, independent_ids)
    err = float(np.max(np.abs(actual - expected)))
    qerr = float(np.max(np.abs(integrated - expected)))
    return dict(passed=bool(np.array_equal(np.asarray(material.BOUNDS), expected_bounds)
                            and geometry.shape == independent_geometry.shape
                            and np.allclose(geometry, independent_geometry, rtol=1e-13, atol=1e-14)
                            and matches and len(actual) == 12 and err < 1e-12
                            and qerr < 1e-12 and np.all((ids >= 0) & (ids < 12))
                            and np.array_equal(probe_ids, [0, 0, 9, 2, 11])),
                prior="Independent coefficient-uniform and complex-epsilon-uniform proper prior",
                event_order="centre_bin * 3 + annular_bin",
                event_prior_mass=actual, independent_analytic_prior_mass=expected,
                split_coefficient_quadrature_prior_mass=integrated,
                maximum_analytic_error=err, maximum_quadrature_error=qerr,
                radius_zero_atoms=[.5, .5], solid_prior_mass=.25,
                event_partition_disjoint_and_exhaustive=bool(matches),
                independent_decoder_max_error=float(np.max(np.abs(geometry - independent_geometry))),
                boundary_probes=dict(latent=probes, geometry=probe_geometry, event_ids=probe_ids),
                latent_fiber_semantics="Inactive coefficients retain their integrated prior mass; no fixed topology-label generator.")


def _prior_only(material, likelihood):
    c, cw = _split_prior_nodes(CENTRE_CUTS, CENTRE_SCALE)
    s, sw = _split_prior_nodes(SHELL_CUTS, SHELL_SCALE)
    xx, yy = np.meshgrid(c, s, indexing="ij")
    theta = np.column_stack([xx.ravel(), yy.ravel(), np.full(xx.size, 2.5),
                             np.full(xx.size, .1)])
    weights = (cw[:, None] * sw[None, :]).ravel()
    F = np.zeros((len(theta), 16))
    y, sigma = np.linspace(-.4, .5, 16), .8
    ll = np.asarray(likelihood.log_likelihood(F, y, sigma, gain_sd=.15,
                                             gain_bounds=(.7, 1.3)))
    posterior = weights * np.exp(ll - np.max(ll))
    posterior /= posterior.sum()
    event_mass = np.bincount(material.event_ids(theta), weights=posterior, minlength=12)
    expected = np.asarray(material.prior_event_masses())
    gaussian_noise = -.5 * (y @ y / sigma ** 2 + len(y) * math.log(2 * math.pi * sigma ** 2))
    return dict(passed=bool(np.max(np.abs(event_mass - expected)) < 1e-12
                            and np.max(np.abs(ll - gaussian_noise)) < 1e-12),
                likelihood="Synthetic F=0 for every latent: a constant observation likelihood, no Maxwell evaluation",
                posterior_event_mass=event_mass, prior_event_mass=expected,
                maximum_event_error=float(np.max(np.abs(event_mass - expected))),
                loglikelihood_minus_known_constant=float(np.max(np.abs(ll - gaussian_noise))),
                geometry_zero_atom_preserved=bool(abs(event_mass[0] - .25) < 1e-12))


def _joint_loglik(g, F, y, sigma, sd, bounds):
    prior_normalizer = ndtr((bounds[1] - 1) / sd) - ndtr((bounds[0] - 1) / sd)
    residual = y - g * F
    return (-.5 * (residual @ residual / sigma ** 2 + len(F) * math.log(2 * math.pi * sigma ** 2))
            -.5 * ((g - 1) / sd) ** 2 - math.log(sd * math.sqrt(2 * math.pi))
            - math.log(prior_normalizer))


def _gain_reference(F, y, sigma, sd, bounds):
    F, y = np.asarray(F, float), np.asarray(y, float)
    precision = 1 / sd ** 2 + F @ F / sigma ** 2
    variance = 1 / precision
    mean = variance * (1 / sd ** 2 + F @ y / sigma ** 2)
    mode = float(np.clip(mean, *bounds))
    maximum = _joint_loglik(mode, F, y, sigma, sd, bounds)
    integrand = lambda g: math.exp(_joint_loglik(g, F, y, sigma, sd, bounds) - maximum)
    points = [mean] if bounds[0] < mean < bounds[1] else None
    integral, error = quad(integrand, *bounds, points=points, epsabs=1e-12,
                           epsrel=1e-12, limit=150)
    log_quad = maximum + math.log(integral)
    covariance = sigma ** 2 * np.eye(len(F)) + sd ** 2 * np.outer(F, F)
    cholesky = np.linalg.cholesky(covariance)
    whitened = np.linalg.solve(cholesky, y - F)
    logdet = 2 * np.log(np.diag(cholesky)).sum()
    normal_loglik = -.5 * (whitened @ whitened + logdet + len(F) * math.log(2 * math.pi))
    posterior_sd = math.sqrt(variance)
    prior_z = ndtr((bounds[1] - 1) / sd) - ndtr((bounds[0] - 1) / sd)
    posterior_z = ndtr((bounds[1] - mean) / posterior_sd) - ndtr((bounds[0] - mean) / posterior_sd)
    log_cholesky = normal_loglik + math.log(posterior_z) - math.log(prior_z)
    return dict(log_quad=log_quad, log_cholesky=log_cholesky,
                untruncated_conditional_mean=mean, untruncated_conditional_variance=variance,
                quadrature_relative_error_estimate=error / integral,
                logdet=float(logdet), prior_normalizer=float(prior_z),
                posterior_truncation_probability=float(posterior_z))


def _truncated_gain(likelihood):
    rng = np.random.default_rng(SEED)
    bounds, sd = (.7, 1.3), .15
    rows, F_rows, y_rows = [], [], []
    for i, physical_gain in enumerate([.7, .85, 1., 1.18, 1.3]):
        F = rng.normal(size=16)
        sigma = .28 + .03 * i
        y = physical_gain * F + sigma * rng.normal(size=16)
        reference = _gain_reference(F, y, sigma, sd, bounds)
        value = float(np.asarray(likelihood.log_likelihood(F, y, sigma, gain_sd=sd,
                                                          gain_bounds=bounds)))
        conditional_mean, conditional_variance = likelihood.conditional_gain(F, y, sigma, sd)
        conditional_mean, conditional_variance = float(np.asarray(conditional_mean)), float(np.asarray(conditional_variance))
        rows.append(dict(case=i, declared_synthetic_gain=physical_gain, sigma=sigma,
                         shared_loglik=value, **reference,
                         shared_untruncated_conditional_mean=conditional_mean,
                         shared_untruncated_conditional_variance=conditional_variance,
                         quad_absolute_error=abs(value - reference["log_quad"]),
                         cholesky_absolute_error=abs(value - reference["log_cholesky"])))
        F_rows.append(F)
        y_rows.append(y)
    # The same data and shared scalar gain are compared with the deliberately
    # different two-independent-gain likelihood.  This is a semantic negative control.
    F, y, sigma = F_rows[2], y_rows[2], rows[2]["sigma"]
    shared = float(likelihood.log_likelihood(F, y, sigma, gain_sd=sd, gain_bounds=bounds))
    independent = float(likelihood.log_likelihood(F[:8], y[:8], sigma, gain_sd=sd, gain_bounds=bounds)
                        + likelihood.log_likelihood(F[8:], y[8:], sigma, gain_sd=sd, gain_bounds=bounds))
    known = float(likelihood.log_likelihood(F, y, sigma, gain_sd=0., gain_bounds=bounds))
    residual = y - F
    known_reference = -.5 * (residual @ residual / sigma ** 2 + len(F) * math.log(2 * math.pi * sigma ** 2))
    errors = [max(row["quad_absolute_error"], row["cholesky_absolute_error"]) for row in rows]
    conditional_errors = [max(abs(row["shared_untruncated_conditional_mean"] - row["untruncated_conditional_mean"]),
                              abs(row["shared_untruncated_conditional_variance"] - row["untruncated_conditional_variance"]))
                          for row in rows]
    return dict(passed=bool(max(errors) < 1e-9 and max(conditional_errors) < 1e-12
                            and abs(known - known_reference) < 1e-10
                            and abs(shared - independent) > 1e-4),
                noise_convention="Each real packed channel has independent N(0,sigma^2) noise; no complex-variance factor silently substituted",
                prior="Gain N(1,0.15^2) truncated to [0.7,1.3], normalized",
                reference="Scaled independent one-dimensional quadrature and dense Cholesky rank-one marginal + posterior/prior truncation ratio",
                cases=rows, maximum_absolute_loglik_error=max(errors),
                maximum_conditional_gain_error=max(conditional_errors),
                known_gain_absolute_error=abs(known - known_reference),
                two_frequency_shared_gain=dict(shared_gain_loglik=shared,
                                              independent_gains_loglik=independent,
                                              difference=shared - independent,
                                              interpretation="Stacking both frequency channels uses one gain prior. Summing separately marginalized frequencies defines another target."))


def _same_measure():
    scale, rmax, cut = CENTRE_SCALE, RMAX, .20
    L = lambda r: 1 + 3 * (r / scale) ** 2
    z_joint = lambda z: L(float(_radius(z, scale))) / 4
    r_joint = lambda r: L(r) * _positive_radius_density(r, scale)
    atom = .5 * L(0.)
    z_evidence = atom + quad(z_joint, -3., -1., epsabs=1e-12, epsrel=1e-12)[0]
    r_evidence, r_error = quad(r_joint, 0., rmax, epsabs=1e-12, epsrel=1e-12)
    r_evidence += atom
    z_event = (atom + quad(z_joint, -float(_a_of_radius(cut, scale)), -1.,
                           epsabs=1e-12, epsrel=1e-12)[0]) / z_evidence
    r_event = (atom + quad(r_joint, 0., cut, epsabs=1e-12, epsrel=1e-12)[0]) / r_evidence
    z_mean = quad(lambda z: float(_radius(z, scale)) * z_joint(z), -3., -1.,
                  epsabs=1e-12, epsrel=1e-12)[0] / z_evidence
    r_mean = quad(lambda r: r * r_joint(r), 0., rmax,
                  epsabs=1e-12, epsrel=1e-12)[0] / r_evidence
    # Preserve the atom/support while changing only the active-radius measure.
    changed_joint = lambda r: .5 * L(r) / rmax
    changed_evidence = atom + quad(changed_joint, 0., rmax)[0]
    changed_event = (atom + quad(changed_joint, 0., cut)[0]) / changed_evidence
    return dict(passed=bool(np.allclose([z_evidence, z_event, z_mean],
                                        [r_evidence, r_event, r_mean], rtol=1e-11, atol=1e-12)
                            and abs(changed_event - z_event) > 1e-3),
                original_measure="z_c uniform[-3,1] with a 0.5 probability atom at radius=0",
                transformed_continuous_density="r/(0.45^2) * [1-(r/0.45)^2]^-3 on (0,rmax] plus atom 0.5",
                synthetic_likelihood="1+3*(r/0.45)^2", event="central radius <=0.20",
                coefficient_coordinate=dict(evidence=z_evidence, event_mass=z_event, mean=z_mean),
                radius_coordinate=dict(evidence=r_evidence, event_mass=r_event, mean=r_mean,
                                       quadrature_error_estimate=r_error),
                different_prior_negative_control=dict(prior="0.5 atom at 0 plus 0.5 uniform radius on (0,rmax]",
                                                      evidence=changed_evidence, event_mass=changed_event,
                                                      event_mass_difference=changed_event - z_event),
                interpretation="A coordinate transform retains the Jacobian and atom. Radius-uniform is a different prior, not the same target.")


def _continuous_coherence():
    # Known normalized mixture: disconnected narrow and broad components, on
    # an explicit common continuous target, not node/bin hybrid summaries.
    probabilities, means, sds = np.array([.55, .45]), np.array([.22, .78]), np.array([.008, .045])
    aa, bb = -means / sds, (1 - means) / sds
    normalizers = ndtr(bb) - ndtr(aa)
    density = lambda x: float(np.sum(probabilities * norm.pdf((x - means) / sds) / sds / normalizers))
    cdf = lambda x: float(np.sum(probabilities * (ndtr((x - means) / sds) - ndtr(aa)) / normalizers))
    component_mean, component_variance = truncnorm.stats(aa, bb, loc=means, scale=sds, moments="mv")
    mean = float(probabilities @ component_mean)
    variance = float(probabilities @ (component_variance + component_mean ** 2) - mean ** 2)
    quantiles = np.array([brentq(lambda x: cdf(x) - q, 0., 1., xtol=1e-13) for q in [.05, .5, .95]])
    grid = np.linspace(0., 1., 65)
    subdivisions = sorted(set([0., 1.] + [float(np.clip(m + k * s, 0., 1.))
                                           for m, s in zip(means, sds) for k in [-8, -3, 0, 3, 8]]))
    integral = lambda f, cut=1.: sum(quad(f, lo, min(hi, cut), epsabs=1e-12, epsrel=1e-12)[0]
                                   for lo, hi in zip(subdivisions[:-1], subdivisions[1:]) if lo < cut)
    numeric_cdf = np.array([integral(density, cut) for cut in grid])
    cdf_error = float(np.max(np.abs(numeric_cdf - np.array([cdf(x) for x in grid]))))
    numeric_mean = integral(lambda x: x * density(x))
    numeric_variance = integral(lambda x: (x - numeric_mean) ** 2 * density(x))
    lower_bound = .05 * (max(mean - quantiles[0], 0.) ** 2 + max(quantiles[2] - mean, 0.) ** 2)
    event_mass = cdf(.5)
    rho = .995
    covariance = np.array([[1., rho], [rho, 1.]])
    # The whitening transformation includes |det C| and keeps one normalized
    # correlated Gaussian measure. Its covariance/marginal checks are exact algebra.
    chol = np.linalg.cholesky(covariance)
    transformed_covariance = chol @ chol.T
    return dict(passed=bool(abs(integral(density) - 1) < 1e-11 and cdf_error < 1e-10
                            and abs(numeric_mean - mean) < 1e-11
                            and abs(numeric_variance - variance) < 1e-11
                            and variance + 1e-12 >= lower_bound
                            and np.all(np.diff(numeric_cdf) >= -1e-12)
                            and abs(event_mass - .55) < 1e-9
                            and np.allclose(transformed_covariance, covariance)),
                evidence_scope="Analytic normalized toy distributions; no Maxwell or generator claim",
                target="0.55 truncated N(0.22,.008^2) + 0.45 truncated N(0.78,.045^2), support [0,1]",
                normalized_integral=integral(density), event="x<=0.5", event_mass=event_mass,
                cdf_grid=grid, cdf=numeric_cdf, max_cdf_error=cdf_error,
                mean=mean, variance=variance, SD=math.sqrt(variance), quantiles_05_50_95=quantiles,
                necessary_variance_lower_bound=lower_bound,
                independent_mean=numeric_mean, independent_variance=numeric_variance,
                correlated_ridge=dict(correlation=rho, covariance=covariance,
                                      whitening_log_jacobian=float(np.log(np.diag(chol)).sum()),
                                      interpretation="An elongated connected Gaussian is a ridge, not two distinct modes."))


def _historical_negative_controls(root):
    path = root / "research/b1_next/e0/results/CONTROLS.json"
    saved = json.loads(path.read_text(encoding="utf-8"))
    checks = saved["checks"]
    rare = checks["high_ESS_missing_mode_counterexample"]
    old = checks["supplied_verify_theory"]["result"]["checks"]["B1_saved_moment_quantile_coherence"]
    lower = .05 * (max(old["native_mean"] - old["native_q05"], 0.) ** 2
                   + max(old["native_q95"] - old["native_mean"], 0.) ** 2)
    mismatch = lower > old["reported_variance"] + .6 ** 2 * 1e-12
    rare_detected = (rare["ESS_fraction"] == 1. and rare["estimated_event_mass"] == 0.
                     and rare["true_event_mass"] > .49)
    return dict(passed=bool(mismatch and rare_detected), source=str(path.relative_to(root)),
                new_campaign=False, rare_mode=rare,
                historical_coherence=dict(case="distinguishable/H0/Re(epsilon)",
                                          mean=old["native_mean"], q05=old["native_q05"], q95=old["native_q95"],
                                          reported_variance=old["reported_variance"],
                                          independent_necessary_lower_bound=lower,
                                          violation_factor=lower / old["reported_variance"],
                                          detected=mismatch,
                                          interpretation="Incompatible numerical summaries; this is not an evidence-normalization diagnosis."),
                interpretation="Existing exposed E0 counterexamples are reused without repeating any physics or oracle campaign.")


def _mass_semantics(material):
    masses = np.asarray(material.prior_event_masses(), float)
    order = np.argsort(-masses, kind="stable")
    selected = order[:3]
    represented = float(masses[selected].sum())
    remainder = float(masses[order[3:]].sum())
    package = dict(hypotheses=[dict(event_id=int(i), posterior_mass=float(masses[i])) for i in selected],
                   represented_mass=represented, unshown_modeled_mass=remainder,
                   unexplored_mass_status="All twelve events analytically integrated for this prior-only example; no external-support guarantee",
                   outside_model_support_mass=None,
                   exact_continuous_image_probability=None)
    return dict(passed=bool(abs(represented + remainder - 1.) < 1e-12 and remainder > 0.
                            and represented < 1. and package["outside_model_support_mass"] is None),
                package=package,
                interpretation="Top-K masses are not renormalized. Modeled remainder is quantified; unspecified material worlds remain unquantified. A representative is not a point probability.")


def run_controls(root=ROOT):
    root = Path(root).resolve()
    if str(root / "src") not in sys.path:
        sys.path.insert(0, str(root / "src"))
    wall, cpu = time.perf_counter(), time.process_time()
    loaded_before = set(sys.modules)
    material = importlib.import_module("b1_joint.material")
    likelihood = importlib.import_module("b1_joint.likelihood")
    checks, errors, warnings_out = {}, [], []
    jobs = [("material_prior_decoder_event_partition", lambda: _shared_material(material)),
            ("proper_prior_only_posterior", lambda: _prior_only(material, likelihood)),
            ("truncated_gain_marginalization_shared_frequency", lambda: _truncated_gain(likelihood)),
            ("same_measure_coordinate_change_prior_sensitivity", _same_measure),
            ("continuous_CDF_moment_quantile_coherence", _continuous_coherence),
            ("reused_high_ESS_and_incoherent_summary_counterexamples", lambda: _historical_negative_controls(root)),
            ("Top_K_remainder_outside_model_semantics", lambda: _mass_semantics(material))]
    for name, fn in jobs:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                checks[name] = dict(execution_status="EXECUTED", **fn())
            except Exception as exc:
                checks[name] = dict(execution_status="ERROR", passed=False,
                                    error=dict(type=type(exc).__name__, message=str(exc)))
                errors.append(dict(check=name, type=type(exc).__name__, message=str(exc)))
        warnings_out.extend(dict(check=name, category=w.category.__name__, message=str(w.message)) for w in caught)
    unwanted = sorted(name for name in set(sys.modules) - loaded_before
                      if name.startswith(("b1.physics", "b1_joint.physics", "b1.ledger", "b1_joint.ledger")))
    return _native(dict(stage="B1-C1-Joint-P2-E0", date="2026-10-13", seed=SEED,
                        new_material_frequency_evaluations=0,
                        physical_forward_interface=False, resource_ledger_write_interface=False,
                        neural_training=False, locked_truth_access=False,
                        historical_results_modified=False,
                        checks=checks, errors=errors, warnings=warnings_out,
                        unintended_physics_or_ledger_imports=unwanted,
                        all_required_checks_pass=bool(not unwanted and not errors and all(row["passed"] for row in checks.values())),
                        wall_seconds=time.perf_counter() - wall, cpu_seconds=time.process_time() - cpu,
                        status_meaning="Finite probability/algebra controls only. No Maxwell approximation, posterior mixing, calibration, continuum, or model-adequacy certificate."))


def render_report(result):
    checks = result["checks"]
    status = "PASS" if result["all_required_checks_pass"] else "FAIL/需定位"
    rows = ["# B1 × C1 Joint P2 E0：零新增物理概率审计", "",
            f"日期：2026-10-13。实际运行状态：{status}。新增材料–频率评估 **0**，无物理 forward/ledger 写入接口。",
            f"本脚本实际 wall {result['wall_seconds']:.6f} 秒、CPU {result['cpu_seconds']:.6f} 秒；阶段共享 job meter 另记总受控费用。", "",
            "这些是共同 prior、decoder、likelihood 与概率摘要的必要核对，不能替代 full-Mie 求解器资格、连续参考收敛、采样 mixing 或现实概率校准。", "",
            "| 必要核对 | 实际结果 | 证据 |", "|---|---|---|"]
    for name, item in checks.items():
        if item["execution_status"] == "ERROR":
            evidence = f"{item['error']['type']}: {item['error']['message']}"
        elif name == "material_prior_decoder_event_partition":
            evidence = f"解析事件质量误差 {item['maximum_analytic_error']:.3g}；分段 coefficient 积分误差 {item['maximum_quadrature_error']:.3g}；solid prior mass=.25"
        elif name == "proper_prior_only_posterior":
            evidence = f"合成常数 likelihood 恢复 prior；事件误差 {item['maximum_event_error']:.3g}"
        elif name == "truncated_gain_marginalization_shared_frequency":
            evidence = f"独立一维积分/Cholesky loglik 最大误差 {item['maximum_absolute_loglik_error']:.3g}；保留 logdet、截断 prior 归一化和共享两频 gain"
        elif name == "same_measure_coordinate_change_prior_sensitivity":
            evidence = f"原测度事件质量 {item['coefficient_coordinate']['event_mass']:.8f}；换成 radius-uniform prior 后 {item['different_prior_negative_control']['event_mass']:.8f}"
        elif name == "continuous_CDF_moment_quantile_coherence":
            evidence = f"窄/宽双峰同密度 CDF 误差 {item['max_cdf_error']:.3g}；SD、quantile、event mass 一致；相关 ridge 不当作两个 modes"
        elif name == "reused_high_ESS_and_incoherent_summary_counterexamples":
            evidence = f"复用既有高 ESS 漏约半数质量负例；检测历史方差下界违反 {item['historical_coherence']['violation_factor']:.3f} 倍，未重跑历史物理"
        else:
            evidence = "保留 modeled remainder；模型外质量 null/unquantified；未将 Top-K 重新归一化"
        rows.append(f"| {name} | {'PASS' if item['passed'] else 'FAIL'} | {evidence} |")
    rows.extend(["", "共同材料 prior 由连续系数生成，中心和环形几何分别具有 0.5 的零半径/半宽原子质量；联合实心质量为 0.25。正半径的密度由 coefficient-uniform 测度与解析 root Jacobian 诱导，不能用均匀半径替代。inactive latent fibers 在积分中保留，拓扑只由材料派生。", "",
                 "噪声约定为 Re 后接 Im 的各个实通道独立方差 σ²。gain proper prior 为截断 N(1,.15²)、范围 [.7,1.3]；两频先合并通道再边缘化一个共享 gain。分别对两频边缘化意味着两个 gain，是另一个 target。补偿最小残差并未被当成 evidence。", "",
                 "原生 B1 distinguishable/H0/Re(ε) 的 node moment 与 bin quantile 不相容仍保留为历史诊断，不称为 evidence 归一化错误。高 ESS 负例说明 visited support 的良好权重并不保证全 posterior mode coverage。", "",
                 "本次未读取 locked truth，未训练 NN，未调用 Maxwell，未执行旧 campaign、额外 smoke 或 SHA-256 检查。共同物理参考、C1 密度归一化/训练与 fresh-panel 结论由后续实际实验另列。", ""])
    if result["warnings"] or result["errors"] or result["unintended_physics_or_ledger_imports"]:
        rows += ["## 保留的运行问题", "", "```json", json.dumps(dict(errors=result["errors"], warnings=result["warnings"], imports=result["unintended_physics_or_ledger_imports"]), ensure_ascii=False, indent=2), "```", ""]
    return "\n".join(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--report", action="store_true", help="Write the owned research/joint_p2/E0_AUDIT.md after the controlled run")
    args = parser.parse_args()
    result = run_controls(args.root)
    if args.report:
        path = args.root / "research/joint_p2/E0_AUDIT.md"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_report(result), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2))
    return 0 if result["all_required_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
