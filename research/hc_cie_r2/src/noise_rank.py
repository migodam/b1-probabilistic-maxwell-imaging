"""Propagate the archived receive noise into the projected reciprocal equations."""
import argparse
import csv
import json
from pathlib import Path
import time

import spectra as sp  # single-thread native BLAS
import numpy as np
from scipy.linalg import svd
from scipy.stats import beta, chi2

from spectra_inputs import load_observations


HERE = Path(__file__).resolve().parents[1]
DESIGN = dict(r13_cases=["r13_fresh_0000", "r13_fresh_0001"], r15_cases=[], r16_cases={},
              srap_objects=["s2d_010"], srap_layouts=["pw360", "pw180", "ls_line_same"])
SCALES = (0.01, 0.10, 1.00)
ALPHA = 0.01
CUT = 0.01


def skew_rows(m, n):
    d, s = m.shape
    a, b = np.triu_indices(s, 1)
    rows = np.zeros((len(a), d, s), dtype=complex)
    for index, (i, j) in enumerate(zip(a, b, strict=True)):
        rows[index, :, i] = m[:, j]
        rows[index, :, j] = -m[:, i]
    return rows.reshape(len(a), d*s) / np.sqrt(n)


def form_m(basis, incident):
    # A conservative complex-dot-product envelope for the stored basis and B.
    value = basis.T @ incident
    envelope = sp._gamma(incident.shape[0]) * np.linalg.norm(abs(basis).T @ abs(incident))
    return value, float(envelope)


def exterior_rows(v0):
    s, q = v0.shape
    a, b = np.triu_indices(s, 1)
    i, j = np.triu_indices(q, 1)
    return v0[a[None, :], i[:, None]] * v0[b[None, :], j[:, None]] - \
           v0[b[None, :], i[:, None]] * v0[a[None, :], j[:, None]]


def factor_rows(operator, input_error):
    u, s, vh = svd(operator, full_matrices=False, check_finite=False)
    # Add a conservative numerical scale floor to the assembly perturbation.
    floor = input_error + sp._gamma(max(operator.shape)) * np.linalg.norm(operator)
    return u, s, vh, float(floor)


def whitened_ball(obs, detector_dimension, output_rank):
    if output_rank == 0:
        return 0.0
    if obs.sigma is not None:
        return float(obs.sigma * np.sqrt(chi2.isf(ALPHA, 2*output_rank)/2))
    if output_rank == detector_dimension:
        return obs.delta
    return float(obs.delta * np.sqrt(beta.isf(ALPHA, output_rank, detector_dimension-output_rank)))


def one_mode_family_radius(obs, detector_dimension, pair_count):
    if obs.sigma is not None:
        return float(obs.sigma * np.sqrt(-np.log(ALPHA/pair_count)))
    return float(obs.delta * np.sqrt(beta.isf(ALPHA/pair_count, 1, detector_dimension-1)))


def reference_material(obs):
    if obs.cohort == "r13":
        file = HERE / f"results/real_weights/{obs.key}.npz"
        with np.load(file, allow_pickle=False) as archive:
            return archive["anchor_material"], f"results/real_weights/{obs.key}.npz:anchor_material"
    file = HERE / "results/spectra/feasible_followup.npz"
    with np.load(file, allow_pickle=False) as archive:
        return archive["base_material"], "results/spectra/feasible_followup.npz:base_material (same material on all SRAP layouts)"


def assay(obs, recv_left, receive, vh, count, label, current_norm, reference_residual):
    n, s = obs.b.shape
    dark = vh[count:].conj().T
    light = vh[:count].conj().T / receive[:count]
    m_l, m_l_error = form_m(dark, obs.b)
    m_q, m_q_error = form_m(light, obs.b)
    l = skew_rows(m_l, n)
    q = skew_rows(m_q, n)
    root2n = np.sqrt(2/n)
    ul, sig_l, _vl, l_floor = factor_rows(l, root2n*m_l_error)
    uq, sig_q, _vq, q_floor = factor_rows(q, root2n*m_q_error)
    p, ds = l.shape
    detector_dimension = obs.y.size
    noise_sigma = obs.sigma if obs.sigma is not None else obs.delta/np.sqrt(detector_dimension)
    projected_q = ul.conj().T @ q
    projection_noise = np.linalg.norm(projected_q, axis=1)
    mode_upper = projection_noise + q_floor
    mode_lower = np.maximum(sig_l - l_floor, 0)
    mode_covariance = noise_sigma**2 * (projected_q @ projected_q.conj().T)
    q_covariance = noise_sigma**2 * (q @ q.conj().T)
    std = np.sqrt(np.maximum(np.diag(mode_covariance).real, 0))
    norm_scale = np.outer(std, std)
    row_correlation = np.divide(mode_covariance, norm_scale, out=np.zeros_like(mode_covariance), where=norm_scale>0)
    off_diagonal = row_correlation - np.diag(np.diag(row_correlation))
    a, b = np.triu_indices(s, 1)
    modal_y = recv_left[:, :count].conj().T @ obs.y
    actual_j0 = light @ modal_y
    direct_reciprocity = actual_j0.T @ obs.b
    direct_reciprocity = (direct_reciprocity-direct_reciprocity.T)[a,b]/np.sqrt(n)
    propagated_y = q @ modal_y.ravel()
    rng = np.random.default_rng(20261128089)
    raw_probe = rng.standard_normal(obs.y.shape) + 1j*rng.standard_normal(obs.y.shape)
    raw_probe *= obs.delta/np.linalg.norm(raw_probe)
    modal_probe = recv_left[:, :count].conj().T @ raw_probe
    direct_probe = (light @ modal_probe).T @ obs.b
    direct_probe = (direct_probe-direct_probe.T)[a,b]/np.sqrt(n)
    propagated_probe = q @ modal_probe.ravel()

    # The source null can be structural (fewer retained rows than sources). Build it
    # from M_q itself, before using the noisy positive row spectrum.
    _um, m_values, vm = svd(m_q, full_matrices=True, check_finite=False)
    padded = np.pad(m_values, (0, max(s-len(m_values), 0)))
    analytic_q = np.sort(np.hypot(padded[a], padded[b]))[::-1] / np.sqrt(n)
    source_null_dimensions = max(s-count, 0)
    noise_zero_rows = source_null_dimensions*(source_null_dimensions-1)//2
    r0 = exterior_rows(vm[-source_null_dimensions:].conj().T) if noise_zero_rows else np.empty((0,p), complex)
    q0_leak = float(np.linalg.norm(r0 @ q)) if noise_zero_rows else 0.0
    l0 = r0 @ l
    if noise_zero_rows:
        _u0, sig_l0, v0, l0_floor = factor_rows(l0, l_floor)
        gap = max(m_values[-1]-m_q_error, 0.0)
        source_angle_bound = min(2.0, 2*m_q_error/max(gap, 1e-300))
        l0_floor += source_angle_bound*sig_l[0]
    else:
        sig_l0, v0, l0_floor = np.empty(0), np.empty((0,ds), complex), 0.0
    # Never whiten the structural zeros or the positive values buried in the
    # computed noise-map formation floor.
    noisy = (np.arange(p) < p-noise_zero_rows) & (sig_q > q_floor)
    q_min = float(sig_q[noisy][-1]) if np.any(noisy) else None
    white = uq[:, noisy].conj().T / sig_q[noisy,None]
    inflated = 1 + q_floor/q_min if q_min else 1.0
    wl = white @ l
    uncertainty_white = l_floor/q_min if q_min else 0.0
    whitened_99 = whitened_ball(obs, detector_dimension, int(np.count_nonzero(noisy)))
    one_mode_99 = one_mode_family_radius(obs, detector_dimension, p)
    report, arrays = [], dict(l_singular_values=sig_l, q_gain_singular_values=sig_q,
                              detector_covariance_of_reciprocity_rhs=q_covariance,
                              detector_covariance_on_l_singular_rows=mode_covariance,
                              l_mode_noise_correlation=row_correlation, mode_noise_rms=std, mode_noise_gain=projection_noise,
                              filtered_y_reciprocity_rhs=propagated_y, q_check_direct_noise_rhs=direct_probe,
                              q_check_propagated_noise_rhs=propagated_probe,
                              mode_noise_gain_upper=mode_upper, mode_linear_gain_lower=mode_lower,
                              noise_zero_signal_singular_values=sig_l0, stable_noisy_rows=noisy)
    for fraction in SCALES:
        radius = fraction * current_norm
        solved0 = radius*np.maximum(sig_l0-l0_floor, 0) > 2*obs.delta*(q0_leak+q_floor)
        condensed = wl - (wl @ v0[solved0].conj().T) @ v0[solved0]
        if len(condensed):
            whitened_singular = svd(condensed, compute_uv=False, check_finite=False)
            floor_white = uncertainty_white + sp._gamma(max(condensed.shape))*np.linalg.norm(condensed)
            if np.any(solved0):
                # The projection is formed from the nominal exact-noise-free rows.
                # Account for uncertainty in their retained right singular subspace.
                safe_gap = max(sig_l0[solved0][-1]-l0_floor, 1e-300)
                floor_white += np.linalg.norm(wl, ord=2)*min(2.0, 2*l0_floor/safe_gap)
            white_lower = np.maximum(whitened_singular-floor_white, 0)
        else:
            whitened_singular, white_lower, floor_white = np.empty(0), np.empty(0), 0.0
        arrays[f"{fraction:g}_joint_gain_singular_values"] = whitened_singular
        arrays[f"{fraction:g}_joint_gain_lower"] = white_lower
        counts = dict(
            directional_bounded_1delta=int(np.count_nonzero(radius*mode_lower > obs.delta*mode_upper)),
            directional_bounded_2delta=int(np.count_nonzero(radius*mode_lower > 2*obs.delta*mode_upper)),
            directional_family_99=int(np.count_nonzero(radius*mode_lower > one_mode_99*mode_upper)),
            joint_bounded_1delta=int(np.count_nonzero(radius*white_lower > obs.delta*inflated)),
            joint_bounded_2delta=int(np.count_nonzero(radius*white_lower > 2*obs.delta*inflated)),
            joint_99=int(np.count_nonzero(radius*white_lower > whitened_99*inflated)),
            qualified_structural_noise_zero=int(np.count_nonzero(solved0)),
        )
        for key in ("joint_bounded_1delta", "joint_bounded_2delta", "joint_99"):
            counts[key+"_with_structural"] = counts[key]+counts["qualified_structural_noise_zero"]
        report.append(dict(observation=obs.key, truncation=label, reference_scale_fraction=fraction,
                           reference_current_fro=current_norm, dark_perturbation_fro=radius,
                           reference_forward_residual_over_delta=reference_residual,
                           dark_coordinate_columns=ds, retained_receive_modes=count,
                           source_count=s, complex_reciprocity_rows=p, receive_noise_rms=noise_sigma,
                           saved_receive_ball_radius=obs.delta, family_bonferroni_99_preprop_radius=one_mode_99,
                           joint_whitened_ball_99=whitened_99, l_absolute_floor=l_floor,
                           l_sigma_max=float(sig_l[0]), q_absolute_floor=q_floor, q_sigma_max=float(sig_q[0]),
                           l_rows_above_formation_floor=int(np.count_nonzero(sig_l>l_floor)),
                           max_off_diagonal_l_noise_correlation=float(np.max(abs(off_diagonal))),
                           q_analytic_spectrum_max_difference=float(np.max(abs(sig_q-analytic_q))),
                           q_independent_noise_action_relative_error=float(np.linalg.norm(direct_probe-propagated_probe)/np.linalg.norm(direct_probe)),
                           q_independent_observed_y_relative_error=float(np.linalg.norm(direct_reciprocity-propagated_y)/np.linalg.norm(direct_reciprocity)),
                           filtered_j0_full_norm=float(np.linalg.norm(actual_j0)),
                           smallest_stable_q=q_min, noisy_stable_rank=int(noisy.sum()),
                           noisy_unresolved_rows=p-int(noisy.sum())-noise_zero_rows,
                           structurally_zero_noise_rows=noise_zero_rows,
                           direct_structural_q_leakage_bound=q0_leak + q_floor,
                           structural_l_absolute_floor=l0_floor, whitened_noise_inflation=inflated,
                           whitened_gain_absolute_floor=float(floor_white), **counts,
                           **{f"{key}_fraction": value/ds for key,value in counts.items()}))
    return report, arrays


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--workspace", type=Path, default=Path.cwd())
    cli.add_argument("--output", type=Path, default=HERE / "results/spectra")
    args = cli.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    prior = json.loads((HERE / "results/spectra/spectrum.json").read_text())
    start, cpu = time.perf_counter(), time.process_time()
    published, arrays, costs = [], {}, []
    for obs in load_observations(args.workspace.resolve(), DESIGN):
        clock, c0 = time.perf_counter(), time.process_time()
        material, reference = reference_material(obs)
        current = sp._current(obs.g, obs.b, material)
        reference_residual = float(np.linalg.norm(obs.s@current-obs.y)/obs.delta)
        current_norm = float(np.linalg.norm(current))
        u, receive, vh = svd(obs.s, full_matrices=True, check_finite=False)
        for label, cut in (("roundoff_cut", max(obs.s.shape)*sp.EPS), ("tsvd_1em2", CUT)):
            count = int(np.count_nonzero(receive > cut*receive[0]))
            match = next(r for r in prior if r["observation"]==obs.key and r["mode"]==label
                         and r["source_count"]==obs.y.shape[1])
            if match["locked"]!=count or match["gamma_rank"]!=len(obs.b)-1:
                raise ValueError(f"{obs.key}: prior constant-only projection does not match")
            rows, payload = assay(obs, u, receive, vh, count, label, current_norm, reference_residual)
            published.extend(rows)
            arrays.update({f"{obs.key}__{label}__{k}": v for k,v in payload.items()})
            print(obs.key,label,[(r["reference_scale_fraction"],r["directional_bounded_1delta"],
                                     r["joint_bounded_2delta_with_structural"]) for r in rows], flush=True)
        costs.append(dict(observation=obs.key, source=obs.provenance, physical_current_reference=reference,
                          wall_s=time.perf_counter()-clock, cpu_s=time.process_time()-c0,
                          physical_forward_solves=1, forward_rhs=obs.y.shape[1]))
    np.savez_compressed(output / "noise_rank_arrays.npz", **arrays)
    with (output / "noise_rank.csv").open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(published[0]))
        writer.writeheader()
        writer.writerows(published)
    (output / "noise_rank.json").write_text(json.dumps(dict(rows=published, costs=costs,
        wall_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu, alpha=ALPHA,
        scale_reference_fractions=SCALES, receive_relative_cut=CUT), ensure_ascii=False, indent=2)+"\n")


if __name__ == "__main__":
    main()
