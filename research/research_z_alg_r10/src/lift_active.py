"""Exact response-row lifting and coherent-source benchmarks for the scalar TM VIE."""

import json
import os
from pathlib import Path
import sys
import time

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares
from scipy.special import j0

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent / "research_z_alg_r7" / "src"))
from physical_experiment import Model, ReducedKernel, Response, pack

BOUNDS = (0.5, 12.0)
SEEDS = {12: 47850, 24: 54601}
NOISE = (0.0, 2e-5)
RANK_RTOL = 1e-11


def design(k):
    """Distinct apertures and independent 3x4 and 4x6 active grids."""
    start = time.process_time()
    rng = np.random.default_rng(SEEDS[k])
    nx, ny, pitch, source_count, radius, disk = {
        12: (4, 3, .52, 28, 1.80, .065),
        24: (6, 4, .56, 40, 2.54, .066),
    }[k]
    x = pitch * (np.arange(nx) - (nx - 1) / 2)
    y = pitch * (np.arange(ny) - (ny - 1) / 2)
    xx, yy = np.meshgrid(x, y)
    active_xy = np.column_stack((xx.ravel(), yy.ravel())) + rng.uniform(-.047, .047, (k, 2))
    corners_x = x[:-1] + pitch / 2
    corners_y = y[:-1] + pitch / 2
    hx, hy = np.meshgrid(corners_x, corners_y)
    host_xy = np.column_stack((hx.ravel(), hy.ravel()))[::2]
    centers = np.concatenate((active_xy, host_xy))
    host = np.r_[np.zeros(k), rng.uniform(1.4, 3.9, len(host_xy))]
    theta = 2 * np.pi * (np.arange(source_count) + .18) / source_count
    theta += .042 * np.sin(3 * theta) if k == 24 else 0
    rx = np.array([[.37, 1.31] if k == 12 else [-.85, 1.79]])
    sources = radius * np.column_stack((np.cos(theta), np.sin(theta)))
    if k == 24:
        sources *= np.array([1.08, .86])
    model = Model(centers=centers, sources=sources, receivers=rx, radius=disk, k=2*np.pi)
    kernel = ReducedKernel(model, host, active=np.arange(k))
    gram = j0(model.k * np.linalg.norm(sources[:, None] - sources[None, :], axis=-1))
    g, u = np.linalg.eigh(gram)
    keep = g >= 1/16
    whitened = u[:, keep] / np.sqrt(g[keep])
    _, s, vh = np.linalg.svd(kernel.q @ whitened, full_matrices=False)
    if len(s) < k or s[-1] < RANK_RTOL * s[0]:
        raise ValueError(f"Insufficient physical incident rank for {k} cells: {s}")
    feeds = whitened @ vh.conj().T
    return model, kernel, gram, whitened, feeds, s, time.process_time()-start


def random_basis(rng, shape):
    a = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    q, r = np.linalg.qr(a)
    return q * (np.diag(r) / abs(np.diag(r)))


def make_materials(rng, k, sample_count):
    """Interior real contrasts on each separate active disk."""
    return rng.uniform(1.1, 10.9, (sample_count, k))


def full_data(model, kernel, material):
    full = kernel.m0.copy()
    full[kernel.active] = material
    return (model.forward(full, independent=True)["y"] - kernel.y0).ravel()


def affine_response(B, observed):
    """Particular row w0 and columns N with N.T @ B = 0 (no conjugation)."""
    u, s, vh = np.linalg.svd(B.T, full_matrices=True)
    rank = int(np.sum(s > RANK_RTOL * s[0])) if len(s) else 0
    V = vh.conj().T
    w0 = V[:, :rank] @ ((u[:, :rank].conj().T @ observed) / s[:rank])
    return w0, V[:, rank:], s, rank


def directional_feed(kernel, source_basis, direction):
    sensitivity = direction @ kernel.q @ source_basis
    return source_basis @ (sensitivity.conj() / np.linalg.norm(sensitivity))


def make_case(k):
    model, kernel, gram, whitened, modes, singular, setup_cpu = design(k)
    rng = np.random.default_rng(SEEDS[k] + 1)
    material = make_materials(rng, k, 14)
    data = np.array([full_data(model, kernel, z) for z in material])

    # Only this intentionally adversarial old array uses the two training endpoints.
    incident_reference_difference = (data[0] - data[1]) @ modes
    basis = null_space(incident_reference_difference[None, :])
    design_coordinates = basis @ random_basis(rng, (k-1, k-1))
    old_feeds = modes @ design_coordinates
    old_B = kernel.q @ old_feeds
    old_data = data @ old_feeds
    _, N, old_sigma, rank = affine_response(old_B, old_data[0])
    if rank != k - 1:
        raise ValueError("Designed array failed to leave one complex missing direction")
    p_new = directional_feed(kernel, modes, N[:, 0])
    generic_feeds = modes @ random_basis(rng, (k, k))
    export = dict(
        k=k, centers=model.centers, sources=model.sources, receivers=model.receivers,
        disk_radius=model.radius, wavenumber=model.k, host=kernel.m0, active=kernel.active,
        H=kernel.K, R=kernel.L.ravel(), q=kernel.q, y_host=kernel.y0.ravel(),
        gram=gram, whitened_source_basis=whitened, observable_modes=modes,
        p_old=old_feeds, B_old=old_B, p_new=p_new, b_new=kernel.q @ p_new,
        pair_materials=material[:2], heldout_materials=material[2:],
        pair_old_data=old_data[:2], heldout_old_data=old_data[2:],
        pair_new_data=data[:2] @ p_new, heldout_new_data=data[2:] @ p_new,
        independent_all_source_data=data,
        generic_feeds=generic_feeds, generic_B=kernel.q @ generic_feeds,
        generic_data=data @ generic_feeds, old_sigma=old_sigma, full_sigma=singular,
        design_cpu_s=setup_cpu,
    )
    np.savez_compressed(ROOT / "results" / f"lift_k{k}.npz", **export)
    print(f"k={k} prod-cond={singular[0]/singular[-1]:.3g} "
          f"old-cond={old_sigma[0]/old_sigma[-1]:.3g} "
          f"pair-old-gap={np.linalg.norm(old_data[0]-old_data[1]):.3g} "
          f"pair-new-gap={abs((data[0]-data[1])@p_new):.3g}", flush=True)
    return model, kernel, export, material, data, rng


def lifted_coefficients(w0, N, H, R):
    """Q(w0 + N t) = constant + M features(t) with t t* Hermitian."""
    a = N.T
    b = a @ H
    v0 = R + w0 @ H
    w_linear = a * v0.conj()
    v_linear = w0 * b.conj()
    diagonal = (a * b.conj()).imag
    columns = [(w_linear + v_linear).imag.T, (w_linear - v_linear).real.T, diagonal.T]
    for j in range(len(a)):
        for k in range(j + 1, len(a)):
            a_jk = a[j] * b[k].conj()
            a_kj = a[k] * b[j].conj()
            columns.extend([((a_jk + a_kj).imag)[:, None], ((a_jk - a_kj).real)[:, None]])
    return (w0 * v0.conj()).imag, np.column_stack(columns)


def features(t):
    mixed = [part for j in range(len(t)) for k in range(j + 1, len(t))
             for part in ((t[j] * t[k].conj()).real, (t[j] * t[k].conj()).imag)]
    return np.r_[t.real, t.imag, abs(t)**2, mixed]


def normalized_gap(s):
    return float(s[-1] / s[0]) if len(s) and s[0] else 0.0


def reality(w, H, R):
    return (w * (R + w @ H).conj()).imag


def reality_jacobian(w, N, H, R):
    a = N.T * (R + w @ H).conj()
    b = w * (N.T @ H).conj()
    return np.column_stack(((a + b).imag.T, (a - b).real.T))


def material_from_row(w, H, R):
    v = R + w @ H
    return w / v


def recover(B, observed, H, R, rng, *, rank2=False, nstarts=8):
    """Reconstruct physical rows. The nonlinear fallback is a candidate search."""
    begin = time.process_time()
    w0, N, s_B, rank_B = affine_response(B, observed)
    d = N.shape[1]
    scale = max(np.linalg.norm(w0) / np.sqrt(len(w0)), 1e-8)
    constant, raw_operator = lifted_coefficients(w0, N, H, R)
    scales = np.r_[np.full(2*d, scale), np.full(d*d, scale**2)]
    M = raw_operator * scales
    u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > RANK_RTOL * s[0])) if len(s) else 0
    used_rank = 2 if d == 1 and rank2 else rank
    lift = vh[:used_rank].T @ ((u[:, :used_rank].T @ -constant) / s[:used_rank])
    t_lift = (lift[:d] + 1j * lift[d:2*d]) * scale
    lift_error = np.linalg.norm(features(t_lift) - lift * scales) / max(np.linalg.norm(lift * scales), 1e-14)
    t_trials = [(t_lift, "linear lift")]
    orbit_roots = None
    nfev = 0
    if d == 0:
        method = "direct complex row"
    elif used_rank == len(scales):
        method = "full rank lift"
    elif d == 1 and used_rank == 2:
        method = "rank-2 line/paraboloid"
        q = vh[2]
        # Dimensionless coordinates x=Re(t)/scale, y=Im(t)/scale, h=|t|²/scale².
        coeff = [q[0]**2 + q[1]**2, 2 * np.dot(lift[:2], q[:2]) - q[2],
                 np.dot(lift[:2], lift[:2]) - lift[2]]
        if coeff[0] <= 1e-24:
            coeff = coeff[1:]
        roots = np.roots(coeff)
        orbit_roots = [float(r.real) for r in roots if abs(r.imag) < 1e-10]
        t_trials = [
            (scale * np.array([(lift[0]+rho*q[0]) + 1j*(lift[1]+rho*q[1])]), method)
            for rho in orbit_roots
        ]
    else:
        method = "rank-one candidate search"

    if d > 0 and method != "rank-2 line/paraboloid":
        # Solve the structured constraints in t, avoiding free quadratic features.
        # For a full-rank lift use its answer and one optional noisy-data refinement.
        guess = np.r_[t_lift.real, t_lift.imag] / scale
        guesses = [guess]
        if method == "rank-one candidate search":
            guesses = [guess, np.zeros(2*d)]
            guesses.extend(rng.normal(size=(nstarts-2, 2*d)))
        q_scale = max(np.linalg.norm(constant), np.linalg.norm(M), 1e-16)

        def row(v):
            return w0 + N @ ((v[:d] + 1j * v[d:]) * scale)

        for initial in guesses:
            fit = least_squares(lambda v: reality(row(v), H, R) / q_scale, initial,
                                jac=lambda v: reality_jacobian(row(v), N, H, R) * scale / q_scale,
                                ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=100)
            nfev += fit.nfev
            t_trials.append((scale * (fit.x[:d] + 1j * fit.x[d:]), "rank-one refinement"))

    raw_material = material_from_row(w0 + N @ t_lift, H, R)
    candidates = []
    for t, label in t_trials:
        w = w0 + N @ t
        z = material_from_row(w, H, R)
        if not np.all(np.isfinite(z)):
            continue
        material = np.clip(z.real, *BOUNDS)
        Hm = np.eye(len(z)) - material[:, None] * H
        v = np.linalg.solve(Hm.T, R)
        actual_row = material * v
        predicted = actual_row @ B
        candidate = dict(
            w=w, material=material, unprojected_material=z, label=label,
            residual=np.linalg.norm(predicted - observed), prediction=predicted,
            reality_gap=np.linalg.norm(z.imag) / max(np.linalg.norm(z), 1e-14),
            projection_gap=np.linalg.norm(actual_row - w) / max(np.linalg.norm(w), 1e-14),
            box_inside=bool(np.all((z.real >= BOUNDS[0]) & (z.real <= BOUNDS[1]))),
            material_row=actual_row, resolvent_sigma_min=np.linalg.svd(Hm, compute_uv=False)[-1],
        )
        candidates.append(candidate)

    unique = []
    for cand in sorted(candidates, key=lambda c: c["residual"]):
        if not any(np.linalg.norm(cand["material"] - seen["material"]) < 1e-5 for seen in unique):
            unique.append(cand)
    exact_fits = [c for c in unique if c["box_inside"] and c["reality_gap"] < 1e-7
                  and c["residual"] < 1e-8 * max(np.linalg.norm(observed), 1e-6)]
    t_check = scale * (np.arange(1, 1+d) + 1j * np.arange(d))
    formula_gap = np.linalg.norm(constant + raw_operator @ features(t_check)
                                 - reality(w0 + N @ t_check, H, R))
    best = unique[0] if unique else None
    fit_jac_s = np.linalg.svd(reality_jacobian(best["w"], N, H, R), compute_uv=False) if best else []
    return dict(best=best, candidates=unique, exact_candidates=exact_fits, method=method,
                d=d, B_rank=rank_B, B_singular_values=s_B, lift_size=len(scales), lift_rank=rank,
                lift_singular_values=s, assumed_rank=used_rank, lift_cone_discrepancy=lift_error,
                local_reality_normalized_gap=normalized_gap(fit_jac_s),
                lift_formula_discrepancy=formula_gap, nfev=nfev, line_roots=orbit_roots,
                cpu_s=time.process_time()-begin, N=N, particular=w0, raw_lift_material=raw_material)


def material_fit(kernel, feeds, data, initials):
    """Ordinary bounded nonlinear least squares with a reused R7 analytic Jacobian."""
    begin = time.process_time()
    f = Response(kernel, feeds)
    norm = max(np.linalg.norm(data), 1e-6)
    fit_rows = []
    for z0 in initials:
        ans = least_squares(lambda z: pack(f.predict(z) - data) / norm, z0,
                            jac=lambda z: f.predict(z, derivative=True) / norm,
                            bounds=BOUNDS, ftol=1e-12, xtol=1e-12, gtol=1e-12,
                            max_nfev=110)
        fit_rows.append(dict(material=ans.x, residual=np.linalg.norm(ans.fun) * norm,
                             nfev=ans.nfev, njev=ans.njev, status=ans.status))
    best = min(fit_rows, key=lambda row: row["residual"])
    return dict(best=best, starts=fit_rows, nlu=f.nlu, nrhs=f.nrhs,
                cpu_s=time.process_time()-begin)


def disc_noise(rng, count, sigma):
    return sigma * np.sqrt(rng.uniform(size=count)) * np.exp(2j*np.pi*rng.uniform(size=count))


def error(material, truth):
    return float(np.linalg.norm(material - truth) / np.sqrt(len(truth)))


def score(rec, truth):
    keys = ("d", "B_rank", "B_singular_values", "lift_size", "lift_rank", "lift_singular_values",
            "assumed_rank", "method", "lift_cone_discrepancy", "cpu_s", "nfev",
            "lift_formula_discrepancy", "local_reality_normalized_gap", "line_roots")
    out = {key: rec[key] for key in keys}
    raw = rec["raw_lift_material"]
    raw_rmse = error(raw.real, truth) if np.all(np.isfinite(raw)) else None
    out.update(raw_linear_lift_rmse=raw_rmse, incident_normalized_gap=normalized_gap(rec["B_singular_values"]),
               lift_normalized_gap=normalized_gap(rec["lift_singular_values"]),
               exact_candidate_count=len(rec["exact_candidates"]),
               exact_candidate_rmse=[error(c["material"], truth) for c in rec["exact_candidates"]],
               best_candidate_rmse=error(rec["best"]["material"], truth) if rec["best"] else None,
               best_candidate_data_residual=rec["best"]["residual"] if rec["best"] else None,
               best_candidate_reality_gap=rec["best"]["reality_gap"] if rec["best"] else None,
               candidates=[{"material": c["material"], "residual": c["residual"],
                            "reality_gap": c["reality_gap"], "inside_box": c["box_inside"],
                            "resolvent_sigma_min": c["resolvent_sigma_min"]} for c in rec["candidates"]])
    return out


def score_trf(fit, truth):
    return dict(rmse=error(fit["best"]["material"], truth), residual=fit["best"]["residual"],
                cpu_s=fit["cpu_s"], nlu=fit["nlu"], nrhs=fit["nrhs"], starts=[
                    dict(row, rmse=error(row["material"], truth)) for row in fit["starts"]
                ])


def feed_power(feeds, gram):
    return np.einsum("is,ij,js->s", feeds.conj(), gram, feeds).real


def active_feed(kernel, modes, rec):
    """Resolve the widest observed-data candidate pair or the weakest local direction."""
    candidates = rec["exact_candidates"]
    if rec["method"] == "rank-2 line/paraboloid" and len(candidates) < 2:
        candidates = [candidate for candidate in rec["candidates"] if candidate["box_inside"]]
    if len(candidates) >= 2:
        differences = [(a["w"] - b["w"]) for i, a in enumerate(candidates)
                       for b in candidates[i+1:]]
        direction = max(differences, key=lambda v: np.linalg.norm(v @ kernel.q @ modes))
        return directional_feed(kernel, modes, direction), "separate recovered material candidates"

    N = rec["N"]
    w = rec["best"]["w"] if rec["best"] else rec["particular"]
    vh = np.linalg.svd(reality_jacobian(w, N, kernel.K, kernel.L.ravel()), full_matrices=False)[2]
    direction = N @ (vh[-1, :N.shape[1]] + 1j * vh[-1, N.shape[1]:])
    return directional_feed(kernel, modes, direction), "sense weakest local reality direction"


def score_sequence(kernel, info, feeds, measured, truth, initial, rng, old_rec=None):
    B = kernel.q @ feeds
    rec = old_rec if old_rec is not None else recover(B, measured, info["H"], info["R"], rng)
    fit = material_fit(kernel, feeds, measured, initial)
    power = feed_power(feeds, info["gram"])
    return dict(lift=score(rec, truth), trf=score_trf(fit, truth),
                shots=feeds.shape[1], radiated_power_sum=np.sum(power),
                max_feed_l2_norm=np.linalg.norm(feeds, axis=0).max(),
                max_unit_power_error=np.max(abs(power-1)))


def designed_pair_experiment(kernel, info, material, all_data, rng, initials):
    B, H, R = info["B_old"], info["H"], info["R"]
    old_feeds, modes = info["p_old"], info["observable_modes"]
    k = len(R)
    rows = []
    for truth_id in (0, 1):
        truth = material[truth_id]
        complete_y = all_data[truth_id]
        clean_old = complete_y @ old_feeds
        noise_unit = disc_noise(rng, k, 1.0)
        random_feed = modes @ random_basis(rng, (k, 1))[:, 0]
        for sigma in NOISE:
            noisy_old = clean_old + sigma*noise_unit[:-1]
            # A rank-2 fit tests the designed two-branch hypothesis under nonzero noise.
            # The numerical rank of the untruncated operator is still recorded.
            old_rec = recover(B, noisy_old, H, R, rng, rank2=sigma > 0)
            selected, rule = active_feed(kernel, modes, old_rec)
            setup_start = time.process_time()
            _, N, _, _ = affine_response(B, noisy_old)
            completion = directional_feed(kernel, modes, N[:, 0])
            align = np.vdot(completion, selected)
            completion *= align / abs(align)
            completion_cpu = time.process_time()-setup_start
            experiments = {"selected": selected, "random": random_feed,
                           "repeat_old": old_feeds[:, 0], "optimal_completion": completion}
            protocols = {}
            endpoint_old_max_errors = np.max(abs(all_data[:2] @ old_feeds - noisy_old), axis=1)
            for name, feed in experiments.items():
                combined = np.column_stack((old_feeds, feed))
                y = np.r_[noisy_old, complete_y @ feed + sigma*noise_unit[-1]]
                # Repeated feed preserves the lone missing direction in the stress array.
                rec = recover(kernel.q @ combined, y, H, R, rng, rank2=sigma > 0)
                report = score_sequence(kernel, info, combined, y, truth, initials, rng, rec)
                if name == "optimal_completion":
                    report["feed_difference_vs_selected_after_phase_alignment"] = np.linalg.norm(
                        completion-selected)
                report["added_shot_true_pair_gap"] = abs((all_data[0]-all_data[1]) @ feed)
                endpoint_errors = np.maximum(endpoint_old_max_errors, abs(all_data[:2] @ feed - y[-1]))
                report["bounded_noise_feasible_prototype_count"] = int(np.sum(endpoint_errors <= sigma + 1e-12))
                if old_rec["candidates"]:
                    incident = kernel.q @ feed
                    best = min(old_rec["candidates"], key=lambda c: abs(c["material_row"] @ incident-y[-1]))
                    local = material_fit(kernel, combined, y, [best["material"]])
                    report.update(branch_selection_rmse=error(best["material"], truth),
                                  branch_plus_polish=score_trf(local, truth))
                protocols[name] = report
            rows.append(dict(truth_id=truth_id, sigma=sigma, selected_rule=rule,
                             completion_cpu_s=completion_cpu,
                             bounded_noise_feasible_prototype_count_old=int(np.sum(
                                 endpoint_old_max_errors <= sigma + 1e-12)),
                             old=score_sequence(
                                 kernel, info, old_feeds, noisy_old, truth, initials, rng, old_rec),
                             protocols=protocols))

    # Isolate discrimination from old-data uncertainty: a higher-noise single added shot.
    # The radius is chosen from the achieved test gap for the stress sweep only.
    recovered = recover(B, all_data[0] @ old_feeds, H, R, rng)
    start = time.process_time()
    best, rule = active_feed(kernel, modes, recovered)
    design_cpu = time.process_time()-start
    optimum_gap = abs((all_data[0] - all_data[1]) @ best)
    stress_rows = []
    for fraction in (0.10, 0.30, 0.48, 0.65):
        sigma = fraction * optimum_gap
        draws = 256
        noise = disc_noise(rng, draws, sigma)
        results = {name: [] for name in ("selected", "random", "repeat_old")}
        for i in range(draws):
            randomized = modes @ random_basis(rng, (k, 1))[:, 0]
            true_id = i % 2
            for name, feed in zip(results, (best, randomized, old_feeds[:, 0])):
                y_added = all_data[true_id] @ feed + noise[i]
                predictions = np.array([c["material_row"] @ (kernel.q @ feed)
                                        for c in recovered["exact_candidates"]])
                distance = abs(predictions - y_added)
                true_cand = min(range(len(predictions)), key=lambda x: error(
                    recovered["exact_candidates"][x]["material"], material[true_id]))
                inside = distance <= sigma + 1e-12
                difference = abs(predictions[0]-predictions[1])
                if difference < 1e-12:
                    winner = 0.5
                else:
                    winner = int(np.argmin(distance) == true_cand)
                results[name].append((winner, int(sum(inside)), int(inside[true_cand]), difference))
        stress_rows.append(dict(noise_fraction_of_optimal_gap=fraction, sigma=sigma, draws=draws,
                                protocols={name: dict(
                                    nearest_accuracy=float(np.mean(np.array(values)[:, 0])),
                                    unique_witness_rate=float(np.mean(np.array(values)[:, 1] == 1)),
                                    physical_truth_witness_containment=float(np.mean(np.array(values)[:, 2])),
                                    gap_median=float(np.median(np.array(values)[:, 3])),
                                    energy_to_guarantee_separation_at_median_gap=
                                    float((2*sigma / np.median(np.array(values)[:, 3]))**2)
                                    if np.median(np.array(values)[:, 3]) > 1e-12 else None,
                                ) for name, values in results.items()}))

    return dict(clean_candidate_design_cpu_s=design_cpu,
                unit_power_optimal_gap=optimum_gap,
                true_material_separation_rmse=error(material[0], material[1]),
                true_old_pair_gap=np.linalg.norm((all_data[0]-all_data[1]) @ old_feeds),
                selected_rule=rule, rows=rows, isolated_new_shot_stress=stress_rows)


def generic_experiment(kernel, info, material, all_data, rng, initials):
    k = len(info["R"])
    H, R = info["H"], info["R"]
    modes, feeds = info["observable_modes"], info["generic_feeds"]
    d_threshold = int(np.sqrt(k+1)-1)
    before_rank_shots = k-d_threshold-1
    heldout = []
    sweeps = []
    baseline_rng = np.random.default_rng(SEEDS[k] + 490)
    strong_initials = np.vstack((initials, [np.full(k, j) for j in (1, 2.5, 4, 5.5, 7, 8.5, 10, 11.5)],
                                 make_materials(baseline_rng, k, 8)))
    for i, true in enumerate(material[2:], 2):
        complete_y = all_data[i]
        noise_unit = disc_noise(rng, k, 1.0)
        common = feeds[:, :before_rank_shots]
        for sigma in NOISE:
            old_obs = complete_y @ common + sigma*noise_unit[:before_rank_shots]
            old_rec = recover(kernel.q @ common, old_obs, H, R, rng)
            start = time.process_time()
            selected, rule = active_feed(kernel, modes, old_rec)
            design_cpu = time.process_time()-start
            paths = {"old": common,
                     "selected": np.column_stack((common, selected)),
                     "random": feeds[:, :before_rank_shots+1],
                     "full": modes}
            readings = {}
            for name, matrix in paths.items():
                measured = complete_y @ matrix + sigma*noise_unit[:matrix.shape[1]]
                rec = old_rec if name == "old" else recover(kernel.q @ matrix, measured, H, R, rng)
                report = score_sequence(kernel, info, matrix, measured, true, initials, rng, rec)
                if rec["best"]:
                    report["response_seeded_trf"] = score_trf(
                        material_fit(kernel, matrix, measured, [rec["best"]["material"]]), true)
                if k == 24 and name in ("selected", "full"):
                    # Small-contrast physical least-squares initialization, then all 20 starts fit the same data.
                    born = (R[:, None] * (kernel.q @ matrix)).T
                    born_real = np.stack((born.real, born.imag), axis=1).reshape(-1, k)
                    born_start = np.linalg.lstsq(born_real, pack(measured), rcond=None)[0]
                    starts = [*strong_initials, np.clip(born_start, *BOUNDS)]
                    report["strong_20_start_trf"] = score_trf(
                        material_fit(kernel, matrix, measured, starts), true)
                readings[name] = report
            # The same designed old array also sees independently generated material.
            adv = info["p_old"]
            old_adversarial_data = complete_y @ adv + sigma*noise_unit[:k-1]
            rec = recover(kernel.q @ adv, old_adversarial_data, H, R, rng)
            heldout.append(dict(
                id=i-2, truth=true, sigma=sigma, selected_rule=rule,
                online_design_cpu_s=design_cpu, shots=readings,
                adversarial_array_on_unseen_material=score(rec, true),
            ))
        if i < 6:
            for count in range(before_rank_shots-2, k+1):
                matrix = feeds[:, :count]
                rec = recover(kernel.q @ matrix, complete_y @ matrix, H, R, rng)
                sweeps.append(dict(id=i-2, shots=count, lift=score(rec, true)))
    return dict(first_count_with_full_size_lift=k-d_threshold, heldout=heldout, rank_sweep=sweeps)


def aggregate(rows, sigma, method):
    methods = [row["shots"][method] for row in rows if row["sigma"] == sigma]
    rmse = [method["lift"]["best_candidate_rmse"] for method in methods]
    valid_rmse = [v for v in rmse if v is not None]
    raw_rmse = [m["lift"]["raw_linear_lift_rmse"] for m in methods]
    return dict(
        shots=methods[0]["shots"], lifts_with_estimate=len(valid_rmse),
        lift_rmse_median=np.median(valid_rmse) if valid_rmse else None,
        raw_lift_rmse_median=np.median(raw_rmse),
        trf_rmse_median=np.median([m["trf"]["rmse"] for m in methods]),
        lift_cpu_median_s=np.median([m["lift"]["cpu_s"] for m in methods]),
        trf_cpu_median_s=np.median([m["trf"]["cpu_s"] for m in methods]),
        trf_nlu_median=np.median([m["trf"]["nlu"] for m in methods]),
        incident_ratio_median=np.median([m["lift"]["incident_normalized_gap"] for m in methods]),
        lift_ratio_median=np.median([m["lift"]["lift_normalized_gap"] for m in methods]),
        exact_candidate_counts=[m["lift"]["exact_candidate_count"] for m in methods],
        median_radiated_power_sum=np.median([m["radiated_power_sum"] for m in methods]),
        response_seeded_trf_rmse_median=np.median([m["response_seeded_trf"]["rmse"] for m in methods]),
        response_seeded_trf_cpu_median=np.median([m["response_seeded_trf"]["cpu_s"] for m in methods]),
        strong_20_start_trf_rmse_median=(np.median([m["strong_20_start_trf"]["rmse"] for m in methods])
                                         if "strong_20_start_trf" in methods[0] else None),
        strong_20_start_trf_cpu_median=(np.median([m["strong_20_start_trf"]["cpu_s"] for m in methods])
                                        if "strong_20_start_trf" in methods[0] else None),
    )


def plots(summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    palette = {"selected": "#225F93", "random": "#E07B40", "repeat_old": "#666666",
               "old": "#787878", "full": "#18876B"}
    label = {"selected": "Selected", "random": "Random", "repeat_old": "Repeated",
             "old": "Old only", "full": "Full span"}
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.8), layout="constrained")
    for row, (k, task) in enumerate(summary.items()):
        rank_ax, stress_ax, noise_ax = axes[row]
        sweeps = task["generic"]["rank_sweep"]
        shots = sorted({case["shots"] for case in sweeps})
        lifts = [[case["lift"] for case in sweeps if case["shots"] == s] for s in shots]
        rank_ax.plot(shots, [cases[0]["lift_rank"] for cases in lifts], "o-", color="#225F93",
                     ms=4, label="Measured rank")
        rank_ax.plot(shots, [cases[0]["lift_size"] for cases in lifts], "s--", color="#E07B40",
                     ms=3.5, label="Feature count")
        rank_ax.axvline(task["generic"]["first_count_with_full_size_lift"], color="gray", lw=.8)
        rank_ax.set(title=f"{k} independent cells", xlabel="Generic old shots", ylabel="Rank or feature count")
        if row == 0:
            rank_ax.legend(frameon=False, fontsize=8)

        probes = task["designed"]["isolated_new_shot_stress"]
        ratio = [p["noise_fraction_of_optimal_gap"] for p in probes]
        for name in ("selected", "random", "repeat_old"):
            stress_ax.plot(ratio, [p["protocols"][name]["unique_witness_rate"] for p in probes],
                           "o-", color=palette[name], label=label[name], ms=4)
        stress_ax.set(ylim=(-.04, 1.07), xlabel="Noise bound / best achievable pair gap",
                      ylabel="Fraction excluding the other witness", title="Designed two-material stress")
        if row == 0:
            stress_ax.legend(frameon=False, fontsize=8)

        x = np.arange(4)
        names = ("old", "selected", "random", "full")
        report = task["aggregate"][str(NOISE[-1])]
        noise_ax.bar(x-.15, [report[name]["lift_rmse_median"] for name in names],
                     width=.30, color=[palette[name] for name in names], label="Response lift")
        noise_ax.bar(x+.15, [report[name]["trf_rmse_median"] for name in names], width=.30,
                     edgecolor=[palette[name] for name in names], color="white", hatch="///",
                     linewidth=1.2, label="3-start TRF")
        noise_ax.scatter(x, [report[name]["response_seeded_trf_rmse_median"] for name in names],
                         label="Lift then 1-start TRF", marker="d", c="#663889", s=24, zorder=4)
        if row == 1:
            noise_ax.scatter([1, 3], [report[name]["strong_20_start_trf_rmse_median"]
                                       for name in ("selected", "full")], marker="x",
                             label="20-start TRF", c="#1D1D1D", s=35, zorder=5)
        noise_ax.set_yscale("log")
        noise_ax.set(xticks=x, xticklabels=[f"{label[name]}\n({report[name]['shots']} shots)" for name in names],
                     ylabel="Median material RMSE", title="12 unseen materials, noisy data")
        noise_ax.set_ylim(top=.0065 if row == 0 else 13)
        noise_ax.tick_params(axis="x", labelsize=9)
        noise_ax.legend(frameon=False, fontsize=7, loc="upper right")
    for ext in ("png", "pdf"):
        fig.savefig(ROOT / "figures" / f"lift_active.{ext}", dpi=210)
    plt.close(fig)


def json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported result type {type(value)}")


def main():
    summary = {}
    for k in SEEDS:
        model, kernel, info, material, all_data, rng = make_case(k)
        initials = np.vstack((np.full(k, np.mean(BOUNDS)), make_materials(rng, k, 2)))
        designed = designed_pair_experiment(kernel, info, material, all_data, rng, initials)
        generic = generic_experiment(kernel, info, material, all_data, rng, initials)
        summary[k] = dict(
            design=dict(k=k, active_cells=k, known_host_cells=model.n-k, external_sources=model.ns,
                        receivers=model.nr, Q_floor=1/16,
                        retained_power_modes=info["whitened_source_basis"].shape[1],
                        max_power_whitening_error=np.linalg.norm(
                            info["observable_modes"].conj().T @ info["gram"] @ info["observable_modes"]-np.eye(k)),
                        design_cpu_s=info["design_cpu_s"], incidence_singular_values=info["full_sigma"],
                        independent_simulation=dict(model.counts)),
            designed=designed, generic=generic,
            aggregate={str(sigma): {key: aggregate(generic["heldout"], sigma, key)
                                    for key in ("old", "selected", "random", "full")}
                       for sigma in NOISE},
        )
        print(f"k={k}, two-candidate gap={designed['unit_power_optimal_gap']:.5g}", flush=True)
        print(json.dumps(summary[k]["aggregate"], default=json_value), flush=True)
        (ROOT / "results" / "lift_active.json").write_text(
            json.dumps(summary, indent=2, default=json_value, allow_nan=False) + "\n")
    plots(summary)


if __name__ == "__main__":
    main()
