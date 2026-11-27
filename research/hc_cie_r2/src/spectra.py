"""Real-grid HC-CIE spectrum and exact reciprocal affine elimination."""
import argparse
import csv
import json
import os
from pathlib import Path
import platform
import resource
import time

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.linalg import eigh, solve, svd
from scipy.stats import beta, chi2

from spectra_inputs import load_observations

EPS = np.finfo(float).eps


def _gamma(k):
    """Conservative operation count for a complex BLAS dot product."""
    scale = (8 * k + 8) * EPS
    return scale / (1 - scale)


def _wedge(x, y, a, b):
    return x[:, a] * y[:, b] - x[:, b] * y[:, a]


def _orthogonal_constant(n):
    """Householder direction mapping the last coordinate to the constant vector."""
    h = np.full(n, -1 / np.sqrt(n))
    h[-1] += 1
    return h / np.linalg.norm(h)


def _gamma_adjoint_norm(v, dark, gdark):
    """Norm of Gamma* v via a skew matrix, without a matrix of all minor outputs."""
    weighted = v.conj()[:, None]
    product = dark.T @ (weighted * gdark)
    return np.linalg.norm(product - product.T) / np.sqrt(2)


def gamma_spectrum(g, dark, seed):
    """Complete spectrum from the restricted Gram and a forward error bound."""
    n, d = dark.shape
    gd = g @ dark
    p = dark @ dark.conj().T
    k = gd @ gd.conj().T
    c = dark @ gd.conj().T
    cross = c * c.conj().T
    gram = p * k - cross
    gram = (gram + gram.conj().T) * 0.5

    # Entrywise bounds by Cauchy--Schwarz and the complex dot-product model.
    nd = np.linalg.norm(dark, axis=1)
    ng = np.linalg.norm(gd, axis=1)
    ng_error = _gamma(n) * np.linalg.norm(g, axis=1) * np.linalg.norm(dark)
    pe = _gamma(d) * np.outer(nd, nd)
    ke = _gamma(d) * np.outer(ng, ng) + np.outer(ng_error, ng) + np.outer(ng, ng_error)
    ke += np.outer(ng_error, ng_error)
    ce = _gamma(d) * np.outer(nd, ng) + np.outer(nd, ng_error)
    formation = np.abs(p) * ke + np.abs(k) * pe + pe * ke
    formation += np.abs(c) * ce.T + np.abs(c.T) * ce + ce * ce.T
    formation += _gamma(1) * (np.abs(p) * np.abs(k) + np.abs(cross) + np.abs(gram))
    gram_error = float(np.linalg.norm(formation))

    # Use the exactly known constant left null direction, not sqrt of its roundoff eigenvalue.
    h = _orthogonal_constant(n)
    gh = gram @ h
    restricted = gram - 2 * (np.outer(h, gh.conj()) + np.outer(gh, h))
    restricted += 4 * (h @ gh).real * np.outer(h, h)
    restricted = restricted[:-1, :-1]
    restricted = (restricted + restricted.conj().T) * 0.5
    values, basis = eigh(restricted, check_finite=False)
    basis_err = float(np.linalg.norm(basis.conj().T @ basis - np.eye(n - 1)))
    eigen_residual = float(np.linalg.norm(restricted @ basis - basis * values))
    restricted_norm = float(np.linalg.norm(restricted))
    comp_error = _gamma(n) * (12 * np.linalg.norm(gram) + 8 * restricted_norm * np.sqrt(n))
    eig_error = (restricted_norm + 2 * max(abs(values))) * basis_err
    eig_error += eigen_residual * np.sqrt(1 + basis_err)
    error = float(gram_error + comp_error + eig_error)
    positive = values > error
    lower_gram = np.maximum(values - error, 0)
    upper_gram = np.maximum(values + error, 0)
    # Preserve signed eigenvalues and intervals if any direction fails to resolve.
    spectrum = np.r_[np.where(positive, np.sqrt(np.maximum(values, 0)), np.nan)[::-1], 0.0]
    interval = np.r_[np.c_[np.sqrt(lower_gram), np.sqrt(upper_gram)][::-1], [[0.0, 0.0]]]
    image_basis = np.vstack([basis, np.zeros((1, n - 1))])
    image_basis -= 2 * h[:, None] * (h @ image_basis)
    reference = np.ones(n) / np.sqrt(n)
    constant_action = _gamma_adjoint_norm(reference, dark, gd)
    test_rng = np.random.default_rng(seed)
    z = test_rng.standard_normal((d, 2)) + 1j * test_rng.standard_normal((d, 2))
    az, gz = dark @ z, gd @ z
    wedge = az[:, 0] * gz[:, 1] - az[:, 1] * gz[:, 0]
    v = test_rng.standard_normal(n) + 1j * test_rng.standard_normal(n)
    v /= np.linalg.norm(v)
    tests = {}
    for name, direction in (("smallest", image_basis[:, 0]), ("largest", image_basis[:, -1]), ("random", v)):
        direct = _gamma_adjoint_norm(direction, dark, gd)
        gram_norm = np.sqrt(max(0.0, float(np.vdot(direction, gram @ direction).real)))
        tests[name] = dict(direct=float(direct), gram=float(gram_norm),
                           absolute_difference=float(abs(direct - gram_norm)))
    q_bound = float(np.sqrt(upper_gram[-1]) * np.sqrt((min(d, 2) - 1) / (2 * min(d, 2))))
    if np.linalg.norm(wedge) > q_bound * np.linalg.norm(z) ** 2:
        raise ValueError("Rank-two direct quadratic action exceeded its upper bound")
    return dict(
        spectrum=spectrum, intervals=interval, signed_eigenvalues=values[::-1], basis=image_basis,
        rank_lower=int(np.count_nonzero(positive)), radius=error, gram_error=gram_error,
        eig_error=float(eig_error), compression_error=float(comp_error), action=tests,
        constant_action=float(constant_action),
        wedge_sum_relative=float(abs(wedge.sum()) / np.linalg.norm(wedge)),
        cols=d * (d - 1) // 2,
    )


def _receive_noise(obs, receive_sv, measured, alpha):
    m, s = measured.shape
    if obs.sigma is not None:
        modal_noise_rms = obs.sigma
        cutoff = modal_noise_rms * np.sqrt(chi2.isf(alpha / m, 2 * s) / 2)
    else:
        modal_noise_rms = obs.delta / np.sqrt(m * s)
        cutoff = obs.delta * np.sqrt(beta.isf(alpha / m, s, (m - 1) * s))
    significance = np.linalg.norm(measured, axis=1) > cutoff
    failures = np.flatnonzero(~significance)
    prefix = int(failures[0]) if len(failures) else m
    return dict(cutoff=float(cutoff), prefix=prefix, noise_rms=float(modal_noise_rms),
                significant=significance, whitened_receive=receive_sv / modal_noise_rms)


def _reciprocal_model(dark, b):
    n, s = b.shape
    m = dark.T @ b
    left, values, vh = svd(m, full_matrices=False, check_finite=False)
    v = vh.conj().T
    a, z = np.triu_indices(s, 1)
    pair = np.hypot(values[a], values[z])
    sig_l = pair / np.sqrt(n)
    tol = max(m.shape) * EPS * values[0]
    source_rank = int(np.count_nonzero(values > tol))
    # The ambient input has ds coordinates. All additional singular values are exactly zero
    # even when M has full column rank. Sub-floor entries here are recorded separately.
    rank = int(np.count_nonzero(pair > tol))
    return dict(left=left, v=v, values=values, a=a, z=z, norm_pair=pair,
                sigma=np.sort(sig_l)[::-1], rank=rank, rank_m=source_rank,
                tol=float(tol), rows=s * (s - 1) // 2, cols=m.size)


def reciprocal_affine_part(model, rhs):
    """Minimum-norm solution of the constant-row skew equation, if compatible."""
    l, v, values = (model[k] for k in ("left", "v", "values"))
    a, b = (model[k] for k in ("a", "z"))
    pair_norm = model["norm_pair"]
    mask = pair_norm > model["tol"]
    t = (v.T @ rhs @ v)[a, b]
    dual = np.zeros_like(t)
    dual[mask] = t[mask] / pair_norm[mask] ** 2
    w = np.zeros((len(values), len(values)), dtype=complex)
    w[b, a] = values[b] * dual
    w[a, b] = -values[a] * dual
    return l.conj() @ w @ v.conj().T


def _current(g, b, material):
    lhs = -material[:, None] * g
    lhs.flat[:: len(material) + 1] += 1
    return solve(lhs, material[:, None] * b, check_finite=False, assume_a="gen", overwrite_a=True)


def _skew(current, b):
    cross = current.T @ b
    return cross - cross.T


def _contrast_partial(obs, dark, reciprocal, j0, b, y, g_off, gamma, config):
    a, c = reciprocal["a"], reciprocal["z"]
    n = b.shape[0]
    rhs = -_skew(j0, b)
    z_affine = reciprocal_affine_part(reciprocal, rhs)
    j_affine = j0 + dark @ z_affine
    e_affine = b + g_off @ j_affine
    wedge_affine = _wedge(j_affine, e_affine, a, c)
    full_illuminations = b.shape[1] == obs.b.shape[1]
    assert full_illuminations
    contrast = config["contrast_scales"]
    saved_scales = config.get("save_contrast_scales", (0.8, 1.0, 1.2))
    base_current = _current(obs.g, b, obs.base)
    output = dict(base_material=obs.base, base_current=base_current, affine_current=j_affine)
    rows = []
    for scale in contrast:
        material = np.clip(obs.lower + scale * (obs.base - obs.lower), obs.lower, obs.upper)
        current = base_current if scale == 1 else _current(obs.g, b, material)
        diff = current - base_current
        projected_dark = dark.conj().T @ (current - j0)
        removable = reciprocal_affine_part(reciprocal, _skew(dark @ projected_dark, b))
        null_component = dark @ (projected_dark - removable)
        alternative = j_affine + null_component
        g_null = g_off @ null_component
        linear = _wedge(null_component, e_affine, a, c) + _wedge(j_affine, g_null, a, c)
        quadratic = _wedge(null_component, g_null, a, c)
        direct = _wedge(alternative, b + g_off @ alternative, a, c)
        physical_q = _wedge(diff, g_off @ diff, a, c)
        q_norm = np.linalg.norm(quadratic)
        q_fraction = {}
        for fraction in config["partial_keep_fractions"]:
            n_keep = max(1, int(np.floor(fraction * (n - 1))))
            retained = gamma["basis"][:, :n_keep].conj().T
            q_keep = np.linalg.norm(retained @ quadratic)
            linear_keep = np.linalg.norm(retained @ linear)
            q_fraction[str(fraction)] = dict(
                retained_rows=n_keep, full_curvature_fraction=float(q_keep / max(q_norm, 1e-300)),
                linear_norm=float(linear_keep), quadratic_norm=float(q_keep),
                sigma_max_upper=float(gamma["intervals"][n - 1 - n_keep, 1]),
            )
        field = b + obs.g @ alternative
        denom = np.einsum("ij,ij->i", field.conj(), field).real
        readout = np.clip(np.einsum("ij,ij->i", field.conj(), alternative).real / denom,
                          obs.lower, obs.upper)
        readout_current = _current(obs.g, b, readout)
        rows.append(dict(
            scale=scale, material_min=float(material.min()), material_max=float(material.max()),
            material_change=float(np.linalg.norm(material - obs.base)),
            physical_current_change=float(np.linalg.norm(diff)),
            physical_port_change_over_delta=float(np.linalg.norm(obs.s @ diff) / obs.delta),
            physical_port_over_delta=float(np.linalg.norm(obs.s @ current - y) / obs.delta),
            physical_quadratic_norm=float(np.linalg.norm(physical_q)),
            alternative_full_port_over_delta=float(np.linalg.norm(obs.s @ alternative - y) / obs.delta),
            alternative_physical_readout_port_over_delta=float(np.linalg.norm(obs.s @ readout_current - y) / obs.delta),
            alternative_closure=float(np.linalg.norm(alternative - readout[:, None] * field) / np.linalg.norm(alternative)),
            projected_derivative_norm=float(np.linalg.norm(linear.sum(axis=0)) / np.sqrt(n)),
            projected_quadratic_norm=float(np.linalg.norm(quadratic.sum(axis=0)) / np.sqrt(n)),
            projected_affine_residual_norm=float(np.linalg.norm(direct.sum(axis=0)) / np.sqrt(n)),
            projected_affine_residual_relative=float(np.linalg.norm(direct.sum(axis=0)) / (np.sqrt(n) * max(1e-300, np.linalg.norm(direct)))),
            full_remainder_norm=float(np.linalg.norm(direct)),
            full_quadratic_norm=float(q_norm), full_linear_norm=float(np.linalg.norm(linear)),
            decomposition_error=float(np.linalg.norm(direct - wedge_affine - linear - quadratic) / np.linalg.norm(direct)),
            partial=q_fraction,
        ))
        if scale in saved_scales:
            tag = str(scale).replace(".", "p")
            output.update({f"material_{tag}": material, f"current_{tag}": current,
                           f"alternative_{tag}": alternative, f"quadratic_{tag}": quadratic})
    return rows, output


def _direct_reciprocal_spectrum(dark, b, model):
    """Build only the p by ds skew map as an independent check of its closed spectrum."""
    m = dark.T @ b
    d, s = m.shape
    rows = np.zeros((model["rows"], d, s), dtype=complex)
    for row, (a, c) in enumerate(zip(model["a"], model["z"], strict=True)):
        rows[row, :, a] = m[:, c]
        rows[row, :, c] = -m[:, a]
    return svd(rows.reshape(model["rows"], d * s) / np.sqrt(b.shape[0]),
               compute_uv=False, check_finite=False)


def _linear_action(g, dark, j0, b, model, seed):
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((dark.shape[1], b.shape[1]))
    z = z + 1j * rng.standard_normal(z.shape)
    h = dark @ z
    a, c = model["a"], model["z"]
    first = _wedge(h, b + g @ j0, a, c)
    second = _wedge(j0, g @ h, a, c)
    divisor = np.sqrt(len(b))
    direct = (first + second).sum(axis=0) / divisor
    skew = _skew(h, b)[a, c] / divisor
    return dict(
        relative_error=float(np.linalg.norm(direct - skew) / np.linalg.norm(skew)),
        cancellation_scaled_error=float(np.linalg.norm(direct - skew) /
                                        (np.linalg.norm(first) + np.linalg.norm(second))),
    )


def _cut_name(relative):
    return f"tsvd_{relative:.0e}".replace("-0", "m").replace("-", "m")


def process(obs, config, cache):
    clock, cpu = time.perf_counter(), time.process_time()
    g_off = (obs.g + obs.g.T) * 0.5
    g_change = float(np.linalg.norm(g_off - obs.g) / np.linalg.norm(obs.g))
    if g_change > 1e-13:
        raise ValueError(f"{obs.key}: Green fails reciprocity")
    np.fill_diagonal(g_off, 0)
    u, receive, vh = svd(obs.s, full_matrices=True, check_finite=False)
    measure = u.conj().T @ obs.y
    m, n = obs.s.shape
    full_sources = obs.y.shape[1]
    noise = _receive_noise(obs, receive, measure, config["noise_familywise_alpha"])
    raw_tol = max(m, n) * EPS * receive[0]
    full_rank = int(np.count_nonzero(receive > raw_tol))
    cuts = {"roundoff_cut": full_rank}
    cuts.update({_cut_name(t): int(np.count_nonzero(receive > t * receive[0]))
                 for t in config["relative_receive_cuts"]})
    cuts["noise_detected_prefix"] = min(full_rank, noise["prefix"])
    parts = dict(receive_singular_values=receive, receive_whitened=noise["whitened_receive"],
                 receive_modal_data_norm=np.linalg.norm(measure, axis=1),
                 incident_singular_values=svd(obs.b, compute_uv=False, check_finite=False))
    rows, contrast_stats, gamma_costs = [], [], []
    partial_arrays = {}
    data_rows = list(config["source_counts"])
    for mode, num_locked in cuts.items():
        dark = vh[num_locked:].conj().T
        light = vh[:num_locked].conj().T
        if num_locked:
            j0 = (light / receive[:num_locked]) @ measure[:num_locked]
            inverse_max = 1 / receive[num_locked - 1]
        else:
            j0 = np.zeros_like(obs.b)
            inverse_max = 0.0
        # Array equality ensures the shared spectra really are from the same operators.
        shared = next((entry for entry in cache if num_locked == entry["rank"]
                       and np.array_equal(g_off, entry["g"]) and np.array_equal(obs.s, entry["s"])), None)
        if shared is None:
            start, start_cpu = time.perf_counter(), time.process_time()
            gamma = gamma_spectrum(g_off, dark, config["action_seed"])
            gamma_elapsed = time.perf_counter() - start
            gamma_cpu = time.process_time() - start_cpu
            shared = dict(rank=num_locked, g=g_off, s=obs.s, gamma=gamma, source=f"{obs.key}:{mode}")
            cache.append(shared)
        else:
            gamma = shared["gamma"]
            gamma_elapsed, gamma_cpu = 0.0, 0.0
        if gamma["rank_lower"] != n - 1:
            raise ValueError(f"{obs.key}: unqualified constant-only projection at {mode}")
        locked_error = np.linalg.norm(u[:, :num_locked].conj().T @ (obs.s @ j0 - obs.y))
        effective_dark_leak = float(np.linalg.norm(u[:, :num_locked].conj().T @ (obs.s @ dark)))
        gamma_costs.append(dict(mode=mode, locked=num_locked, source=shared["source"],
                                wall_s=gamma_elapsed, cpu_s=gamma_cpu))
        parts[f"{mode}_gamma_singular_values"] = gamma["spectrum"]
        parts[f"{mode}_gamma_sigma_intervals"] = gamma["intervals"]
        parts[f"{mode}_gamma_signed_gram_eigenvalues"] = gamma["signed_eigenvalues"]
        for stop in [s for s in data_rows if s <= full_sources]:
            indices = np.rint(np.linspace(0, full_sources - 1, stop)).astype(int)
            b, y = obs.b[:, indices], obs.y[:, indices]
            rec = _reciprocal_model(dark, b)
            sig_l = rec["sigma"]
            parts[f"{mode}_{stop}_sources_reciprocity_l_singular_values"] = sig_l
            parts[f"{mode}_{stop}_sources_incident_dark_singular_values"] = rec["values"]
            direct_fit = -_skew(j0[:, indices], b)
            affine = reciprocal_affine_part(rec, direct_fit)
            j_affine = j0[:, indices] + dark @ affine
            residual = _skew(j_affine, b)
            details = None
            if stop == full_sources:
                independent = _direct_reciprocal_spectrum(dark, b, rec)
                parts[f"{mode}_direct_full_l_singular_values"] = independent
                details = _linear_action(g_off, dark, j0, b, rec, config["action_seed"])
                details["direct_reduced_spectrum_max_abs_difference"] = float(np.max(np.abs(independent - sig_l)))
                details["direct_reduced_spectrum_relative_error"] = float(np.linalg.norm(independent - sig_l) / np.linalg.norm(sig_l))
                details["L_full_ambient_spectral_entries"] = min(n * len(sig_l), rec["cols"])
                if mode == _cut_name(config["relative_receive_cuts"][0]):
                    contrast_stats, partial_arrays = _contrast_partial(obs, dark, rec, j0, b, y, g_off, gamma, config)
            rows.append(dict(
                observation=obs.key, physical_id=obs.physical_id, cohort=obs.cohort,
                mode=mode, spatial_n=n, receivers=m, locked=num_locked, dark=dark.shape[1],
                source_count=stop, source_indices=indices.tolist(), gamma_rank=gamma["rank_lower"],
                gamma_columns=gamma["cols"], gamma_right_nullity=gamma["cols"] - gamma["rank_lower"],
                gamma_sigma_max=float(gamma["spectrum"][0]), gamma_sigma_min_positive=float(gamma["spectrum"][-2]),
                gamma_sigma_min_lower=float(gamma["intervals"][-2, 0]),
                gamma_relative_rank={str(t): int(np.count_nonzero(gamma["spectrum"] > t * gamma["spectrum"][0]))
                                     for t in config["relative_rank_report"]},
                gamma_squared_error_radius=gamma["radius"],
                reciprocity_rank=rec["rank"], reciprocity_m_rank=rec["rank_m"],
                reciprocity_m_floor=rec["tol"], reciprocity_spectrum_rows=rec["rows"],
                reciprocity_columns=rec["cols"], reciprocity_right_nullity=rec["cols"] - rec["rank"],
                reciprocity_structural_ambient_spectrum_zeros=min(n * rec["rows"], rec["cols"]) - rec["rows"],
                reciprocity_alpha=0.0, reciprocity_sigma_max=float(sig_l[0]),
                reciprocity_sigma_min_reduced=float(sig_l[-1]),
                reciprocity_relative_rank={str(t): int(np.count_nonzero(sig_l > t * sig_l[0]))
                                           for t in config["relative_rank_report"]},
                incident_rank=int(np.count_nonzero(parts["incident_singular_values"] > max(obs.b.shape) * EPS *
                                                   parts["incident_singular_values"][0])) if stop == full_sources else None,
                incident_dark_fraction=float(np.linalg.norm(rec["values"]) / np.linalg.norm(b)),
                l_over_columns=rec["rank"] / rec["cols"],
                gamma_action=gamma["action"], gamma_constant_action=gamma["constant_action"],
                gamma_direct_wedge_sum_relative=gamma["wedge_sum_relative"],
                gamma_formation_error=gamma["gram_error"], gamma_eigen_error=gamma["eig_error"],
                gamma_compression_error=gamma["compression_error"],
                reduced_action=details, green_reciprocity_correction_relative=g_change,
                noise_amplification=float(obs.delta * inverse_max),
                noise_rms_amplification=float(noise["noise_rms"] * np.sqrt(stop) * np.linalg.norm(1 / receive[:num_locked])),
                locked_inverse_norm=float(inverse_max),
                receive_dropped_data_norm=float(np.linalg.norm(measure[num_locked:, indices])),
                receive_dropped_forward_norm=float(receive[num_locked]) if num_locked < m else 0.0,
                receive_solver_floor=float(raw_tol),
                receive_floor_rank=full_rank, receive_noise_detected_prefix=int(noise["prefix"]),
                receive_noise_significant_count=int(np.count_nonzero(noise["significant"])),
                receive_noise_detection_threshold=noise["cutoff"],
                receive_min_max_ratio=float(receive[-1] / receive[0]),
                locked_data_absolute_error=float(locked_error), locked_dark_norm=effective_dark_leak,
                j0_norm=float(np.linalg.norm(j0[:, indices])),
                relative_linear_consistency=float(np.linalg.norm(residual) /
                                                  max(np.linalg.norm(j0[:, indices]) * np.linalg.norm(b), 1e-300)),
                affine_norm=float(np.linalg.norm(affine)), delta=obs.delta,
            ))
    spent = dict(observation=obs.key, source=obs.provenance,
                 per_observation_wall_s=time.perf_counter() - clock,
                 per_observation_cpu_s=time.process_time() - cpu,
                 gamma_costs=gamma_costs, partial_forward_solves=2 * len(contrast_stats),
                 forward_rhs=2 * len(contrast_stats) * full_sources)
    return rows, parts, partial_arrays, contrast_stats, spent


def write_csv(path, rows):
    column_names = [k for k, val in rows[0].items() if not isinstance(val, (dict, list))]
    with path.open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=column_names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--workspace", type=Path, default=Path.cwd())
    cli.add_argument("--output", type=Path)
    cli.add_argument("--states-output", type=Path, help="Directory for the bulk finite-perturbation field arrays")
    args = cli.parse_args()
    workspace = args.workspace.resolve()
    here = Path(__file__).resolve().parents[1]
    output = args.output or here / "results/spectra"
    output.mkdir(parents=True, exist_ok=True)
    states_output = args.states_output or output / "states"
    states_output.mkdir(parents=True, exist_ok=True)
    config = json.loads((here / "configs/spectra.json").read_text())
    cache, all_rows, partial, all_spent, indexed = [], [], [], [], []
    start, start_cpu = time.perf_counter(), time.process_time()
    for observation in load_observations(workspace, config):
        rows, arrays, states, contrast, spent = process(observation, config, cache)
        np.savez_compressed(output / f"{observation.key}.npz", **arrays)
        np.savez_compressed(states_output / f"{observation.key}.npz", **states)
        indexed.append(dict(observation=observation.key, spectrum_file=f"{observation.key}.npz",
                            field_file=f"{observation.key}.npz", spectrum_keys=list(arrays),
                            field_keys=list(states)))
        all_rows.extend(rows)
        partial.extend(dict(observation=observation.key, **line) for line in contrast)
        all_spent.append(spent)
        sample = [row for row in rows if row["mode"] == "roundoff_cut"][-1]
        print(observation.key, f"Γ {sample['gamma_rank']} L {sample['reciprocity_rank']}/{sample['reciprocity_columns']}",
              f"σΓ_min/max {sample['gamma_sigma_min_positive']/sample['gamma_sigma_max']:.3g}",
              f"{spent['per_observation_wall_s']:.2f} s", flush=True)
    field_path = str(states_output.relative_to(workspace) if states_output.is_relative_to(workspace) else states_output)
    (output / "array_index.json").write_text(json.dumps(dict(spectra=indexed, field_output_option="--states-output",
                                                            local_fields_path=field_path), indent=2) + "\n")
    (output / "spectrum.json").write_text(json.dumps(all_rows, ensure_ascii=False, indent=2) + "\n")
    (output / "partial.json").write_text(json.dumps(partial, ensure_ascii=False, indent=2) + "\n")
    write_csv(output / "spectrum.csv", all_rows)
    write_csv(output / "partial.csv", partial)
    receipt = dict(platform=platform.platform(), numpy=np.__version__, blas_threads=1,
                   python=platform.python_version(), wall_s=time.perf_counter() - start,
                   cpu_s=time.process_time() - start_cpu,
                   peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   observation_count=len(all_spent), spectra_rows=len(all_rows),
                   physical_forward_solves=sum(s["partial_forward_solves"] for s in all_spent),
                   distinct_gamma_computations=len(cache), observation_costs=all_spent)
    (output / "receipt.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
