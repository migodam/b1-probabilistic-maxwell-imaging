"""Compare one real-weight optical front end against a common filtered fiber."""
import argparse
import csv
import json
from pathlib import Path
import time

import spectra as sp  # set native BLAS threads before loading NumPy
import numpy as np
from scipy.linalg import lu_factor, lu_solve, svd
from scipy.optimize import minimize
from scipy.stats import norm as gaussian

from real_weights import imh, matrix_h, weak_weights


HERE = Path(__file__).resolve().parents[1]
CONFIG = json.loads((HERE / "configs/real_frontend.json").read_text())


def reals(z):
    return np.r_[z.real.ravel(), z.imag.ravel()]


def complex_coordinates(x, shape):
    part = len(x) // 2
    return (x[:part] + 1j * x[part:]).reshape(shape)


def hermitian_coordinates(x):
    a, b = np.triu_indices(x.shape[0], 1)
    return np.r_[x.diagonal().real, np.sqrt(2) * x[a, b].real, np.sqrt(2) * x[a, b].imag]


def weighted_coordinates(j, e, weight):
    return hermitian_coordinates(imh(j.conj().T @ (weight[:, None] * e)))


def weight_rows(t):
    s = t.shape[1]
    basis = np.zeros((s * s, s, s), complex)
    for i in range(s):
        basis[i, i, i] = 1.0
    a, b = np.triu_indices(s, 1)
    diag = np.arange(len(a)) + s
    basis[diag, a, b] = 1 / np.sqrt(2)
    basis[diag, b, a] = 1 / np.sqrt(2)
    basis[diag + len(a), a, b] = 1j / np.sqrt(2)
    basis[diag + len(a), b, a] = -1j / np.sqrt(2)
    vectors = (-1j * (t @ basis)).reshape(s*s, -1)
    return np.concatenate((vectors.real, vectors.imag), axis=1)


def reciprocity_rows(m, spatial):
    d, s = m.shape
    a, b = np.triu_indices(s, 1)
    coefficients = np.zeros((len(a), d, s), complex)
    for k, (i, j) in enumerate(zip(a, b, strict=True)):
        coefficients[k, :, i] = m[:, j]
        coefficients[k, :, j] = -m[:, i]
    c = coefficients.reshape(len(a), -1) / np.sqrt(spatial)
    return np.block([[c.real, -c.imag], [c.imag, c.real]])


def reciprocity_offset(j, b):
    c = j.T @ b
    a, z = np.triu_indices(c.shape[0], 1)
    return reals((c - c.T)[a, z]) / np.sqrt(j.shape[0])


def factor_rows(a, b_norm):
    u, s, vh = svd(a, full_matrices=False, check_finite=False)
    cutoff = max(CONFIG["affine_relative_row_cut"] * s[0], CONFIG["numerical_source_floor"] * b_norm)
    return u, s, vh, s > cutoff, cutoff


def minimum_norm_step(offset, mode_u, mode_s, mode_vh, keep):
    return -mode_vh[keep].T @ ((mode_u[:, keep].T @ offset) / mode_s[keep])


def best_quadratic_step(offset, lin, quad):
    if not np.any(lin):
        return 0.0
    lo, hi = CONFIG["quadratic_step_interval"]
    coefficients = [2 * (quad @ quad), 3 * (lin @ quad), lin @ lin + 2 * (offset @ quad), offset @ lin]
    roots = np.roots(coefficients)
    candidates = [lo, hi, 1.0, *(float(x.real) for x in roots if abs(x.imag) < 1e-8 and lo <= x.real <= hi)]
    return min(candidates, key=lambda a: np.linalg.norm(offset + a * lin + a * a * quad))


def array_readout(j, b, g, lower, upper):
    e = b + g @ j
    chi = np.sum((e.conj() * j).real, axis=1) / np.sum(abs(e)**2, axis=1)
    m = np.clip(chi, lower, upper)
    return m, float(np.linalg.norm(j - m[:, None] * e) / np.linalg.norm(j))


def corrector(g, b, s, y, material, lower, upper):
    wall, cpu = time.perf_counter(), time.process_time()
    y_norm_squared = np.linalg.norm(y)**2
    history = []
    best = None

    def physical(x):
        nonlocal best
        system = -x[:, None] * g
        system.flat[::len(x) + 1] += 1
        lu = lu_factor(system, overwrite_a=True, check_finite=False)
        j = lu_solve(lu, x[:, None] * b, check_finite=False)
        pred = s @ j
        diff = pred - y
        adjoint = lu_solve(lu, s.conj().T @ diff, trans=2, check_finite=False)
        loss = float(np.linalg.norm(diff))
        history.append(loss)
        if best is None or loss < best[0]:
            best = (loss, x.copy(), pred.copy())
        gradient = np.sum((adjoint.conj() * (b + g @ j)).real, axis=1) / y_norm_squared
        return 0.5 * loss**2 / y_norm_squared, gradient

    fit = minimize(physical, material, jac=True, method="L-BFGS-B", bounds=[(lower, upper)] * len(material),
                   options=dict(maxiter=CONFIG["corrector_iterations"], maxfun=CONFIG["corrector_maxfun"],
                                ftol=1e-13, gtol=1e-9))
    return best[1], best[2], dict(corrector_wall_s=time.perf_counter() - wall,
                                  corrector_cpu_s=time.process_time() - cpu,
                                  accepted=int(fit.nit), stop=str(fit.message), evaluations=len(history),
                                  forward_rhs=len(history) * b.shape[1], adjoint_rhs=len(history) * b.shape[1],
                                  residual_history=history)


def optical_frontend(g, b, s, y, sigma, initial, weight):
    u, singular, vh = svd(s, full_matrices=True, check_finite=False)
    r_full = np.count_nonzero(singular > max(s.shape)*sp.EPS*singular[0])
    r_strong = np.count_nonzero(singular > CONFIG["receive_relative_cut"]*singular[0])
    dark, light = vh[r_full:].conj().T, vh[:r_strong].conj().T / singular[:r_strong]
    shape, b_norm = (dark.shape[1], b.shape[1]), np.linalg.norm(b)
    gn, gl = g @ dark, g @ light
    h_dark = matrix_h(weight, dark, gn)
    bp = initial + light @ (u[:, :r_strong].conj().T @ (y - s @ initial))
    ar = reciprocity_rows(dark.T @ b, len(b))
    off_rec = reciprocity_offset(bp, b)
    ur, sr, vr, keep_r, cut_r = factor_rows(ar, b_norm)
    xr = minimum_norm_step(off_rec, ur, sr, vr, keep_r)
    jr = bp + dark @ complex_coordinates(xr, shape)
    er = b + g @ jr
    phi = weighted_coordinates(jr, er, weight)
    t_dark = dark.conj().T @ (weight[:, None] * er) - gn.conj().T @ (weight[:, None] * jr)
    t_light = light.conj().T @ (weight[:, None] * er) - gl.conj().T @ (weight[:, None] * jr)
    aw = weight_rows(t_dark)
    awr = aw - (aw @ vr[keep_r].T) @ vr[keep_r]
    uw, sw, vw, keep_w, cut_w = factor_rows(awr, b_norm)
    receive_jac = weight_rows(t_light)
    rsvd_light = reciprocity_rows(light.T @ b, len(b))
    map_back = -(ur[:, keep_r].T @ rsvd_light) / sr[keep_r, None]
    receive_jac += (aw @ vr[keep_r].T) @ map_back
    b_modes = uw.T @ phi
    std = sigma / np.sqrt(2) * np.linalg.norm(uw.T @ receive_jac, axis=1)
    zcrit = gaussian.isf(CONFIG["noise_family_alpha"] / (2 * len(std)))
    significant = abs(b_modes) > zcrit * std
    qualified = significant & keep_w
    rng = np.random.default_rng(CONFIG["noise_seed"])
    linear_noise, quadratic_noise = [], []
    for _ in range(CONFIG["noise_second_order_probes"]):
        sample = sigma / np.sqrt(2) * rng.standard_normal(receive_jac.shape[1])
        correction = dark @ complex_coordinates(vr[keep_r].T @ (map_back @ sample), shape)
        h = light @ complex_coordinates(sample, (r_strong, b.shape[1])) + correction
        linear_noise.append(float(np.linalg.norm(receive_jac @ sample)))
        quadratic_noise.append(float(np.linalg.norm(weighted_coordinates(h, g @ h, weight))))
    correction_diagnostics = {}
    currents = dict(initial=initial, bp=bp, reciprocity=jr)
    for label, select in (("weighted_exact", keep_w), ("weighted_qualified_exact", qualified)):
        dx = minimum_norm_step(phi, uw, sw, vw, select)
        dz = complex_coordinates(dx, shape)
        h = dark @ dz
        q = hermitian_coordinates(dz.conj().T @ (h_dark @ dz))
        l = aw @ dx
        coeff = best_quadratic_step(phi, l, q)
        if label == "weighted_exact":
            currents["weighted_linear"] = jr + h
        currents[label] = jr + coeff * h
        direct_q = weighted_coordinates(h, g @ h, weight)
        r0_plus = weighted_coordinates(jr + h, b + g @ (jr + h), weight)
        correction_diagnostics[label] = dict(step=float(coeff), current_correction=float(np.linalg.norm(h)),
                                             current_correction_over_recip_current=float(np.linalg.norm(h)/np.linalg.norm(jr)),
                                             quadratic=float(np.linalg.norm(q)), linear=float(np.linalg.norm(l)),
                                             quadratic_over_linear=float(np.linalg.norm(q)/np.linalg.norm(l)) if np.any(l) else None,
                                             remainder_direct_difference=float(np.linalg.norm(q - direct_q)),
                                             decomposition_error=float(np.linalg.norm(r0_plus - phi - l - q)),
                                             measured_weighted_norm_at_alpha1=float(np.linalg.norm(r0_plus)),
                                             measured_weighted_norm_at_given_alpha=float(np.linalg.norm(phi + coeff*l + coeff**2*q)))
    info = dict(receive_rank=int(r_full), strong_receive_rank=int(r_strong),
                recomputed_sigma=sigma, input_projector_for_noise_norm=float(1/singular[r_strong-1]),
                reciprocity_rows_retained=int(np.count_nonzero(keep_r)), reciprocity_cutoff=float(cut_r),
                absolute_reciprocity_svd_max=float(sr[0]),
                weight_stable_rows=int(np.count_nonzero(keep_w)), weight_significant_rows=int(np.count_nonzero(qualified)),
                weight_cutoff=float(cut_w), raw_largest_weight_linear=float(np.linalg.norm(t_dark, ord=2)),
                residualized_svd_max=float(sw[0]), weight_svd_relative_to_incident_fro=float(sw[0]/b_norm),
                phi_before=float(np.linalg.norm(phi)), bonferroni_z=float(zcrit),
                curvature_fro=float(np.linalg.norm(h_dark)),
                noise_linear_rms=float(np.linalg.norm(receive_jac)*sigma/np.sqrt(2)),
                sampled_noise_linear=linear_noise, sampled_noise_quadratic=quadratic_noise,
                corrections=correction_diagnostics)
    arrays = dict(reciprocity_singular_values=sr, weighted_residualized_singular_values=sw,
                  weighted_offset_modes=b_modes, weighted_noisy_mode_std=std,
                  weighted_stable_modes=keep_w, weighted_qualified_modes=qualified)
    return currents, info, arrays


def original_root(workspace):
    local = workspace / "research/research_max_struct_r13"
    if local.is_dir() and (local / "data_public/fresh48").is_dir():
        return local
    return workspace / "research/delegated/max_struct_r1/r13_worktree/research/research_max_struct_r13"


def _open_operators(data, case):
    with np.load(data / f"{case}.npz", allow_pickle=False) as saved:
        return saved["G"], saved["S"]


def run(workspace, output):
    study = original_root(workspace)
    data = study / "data_public/fresh48"
    pred_dir = output / "predictions"
    pred_dir.mkdir(exist_ok=True)
    start, c = time.perf_counter(), time.process_time()
    rows, tests = [], []
    geometry = {}
    for case in CONFIG["cases"]:
        with np.load(data / f"{case}.npz", allow_pickle=False) as saved:
            g, b, s, y, start_m = (saved[k] for k in ("G", "E0", "S", "Y", "m_init"))
            sigma, radius, lower, upper = (float(saved[k]) for k in ("sigma", "delta", "lower", "upper"))
            acquisition = "half" if s.shape[0] == 24 else "full"
        wall, cpu = time.perf_counter(), time.process_time()
        if acquisition not in geometry:
            reference = "r13_fresh_0000" if acquisition == "half" else "r13_fresh_0001"
            g_ref, s_ref = _open_operators(data, reference)
            packaged = workspace / f"research/hc_cie_r2/replay_inputs/real_frontend/{reference}.npz"
            stored_path = packaged if packaged.is_file() else HERE / f"results/real_weights/r13_{reference}.npz"
            with np.load(stored_path, allow_pickle=False) as stored:
                geometry[acquisition] = (g_ref, s_ref, stored["minimum_real_weight"], reference)
        g_ref, s_ref, w, weight_ref = geometry[acquisition]
        if not np.array_equal(g, g_ref) or not np.array_equal(s, s_ref):
            _, ss, vh = svd(s, full_matrices=True, check_finite=False)
            n = vh[np.count_nonzero(ss > max(s.shape) * sp.EPS * ss[0]):].conj().T
            w = weak_weights(g, n, g @ n)[0]
            weight_ref = "same-case full receiver dark subspace"
        j0 = sp._current(g, b, start_m)
        states, meta, arrays = optical_frontend(g, b, s, y, sigma, j0, w)
        tests.append(dict(observation=case, array=acquisition, source=str((data/f"{case}.npz").relative_to(workspace)),
                          weight_reference=weight_ref, initial_chi_min=float(np.min(start_m)),
                          initial_chi_max=float(np.max(start_m)), frontend_wall_s=time.perf_counter()-wall,
                          frontend_cpu_s=time.process_time()-cpu, initial_forward_rhs=b.shape[1], **meta))
        for method in CONFIG["arms"]:
            j = states[method]
            wall, cpu = time.perf_counter(), time.process_time()
            material, closure = array_readout(j, b, g, lower, upper)
            native_pred = s @ sp._current(g, b, material)
            pure_audit_wall_s, pure_audit_cpu_s = time.perf_counter()-wall, time.process_time()-cpu
            final_m, final_y, correction = corrector(g, b, s, y, material, lower, upper)
            tag = f"{case}__{method}"
            np.savez_compressed(pred_dir / f"{tag}.npz", raw_m=material, raw_y=native_pred,
                                corrected_m=final_m, corrected_y=final_y)
            rows.append(dict(observation=case, array=acquisition, arm=method, stored=tag,
                             algebraic_data_residual_over_delta=float(np.linalg.norm(s @ j - y)/radius),
                             algebraic_bounded_constitutive_closure=closure,
                             reciprocity_offset_norm=float(np.linalg.norm(reciprocity_offset(j, b))),
                             optical_offset_norm=float(np.linalg.norm(weighted_coordinates(j, b+g@j,w))),
                             raw_residual_over_delta=float(np.linalg.norm(native_pred-y)/radius),
                             corrected_residual_over_delta=float(np.linalg.norm(final_y-y)/radius),
                             raw_audit_wall_s=pure_audit_wall_s, raw_audit_cpu_s=pure_audit_cpu_s,
                             raw_audit_rhs=b.shape[1], **correction))
        np.savez_compressed(output / f"{case}__qualifications.npz", minimum_real_weight=w, **arrays)
        print(case, acquisition, "rows", meta["weight_stable_rows"], meta["weight_significant_rows"],
              "quad/lin", meta["corrections"]["weighted_exact"]["quadratic_over_linear"], flush=True)
    (output / "physical_receipt.json").write_text(json.dumps(dict(rows=rows, tests=tests,
       wall_s=time.perf_counter()-start, cpu_s=time.process_time()-c), ensure_ascii=False, indent=2)+"\n")


def score(workspace, output):
    study = original_root(workspace)
    paid = json.loads((output / "physical_receipt.json").read_text())
    assert len(paid["rows"]) == len(CONFIG["cases"]) * len(CONFIG["arms"])
    scores = []
    for row in paid["rows"]:
        case = row["observation"]
        with np.load(study / f"sealed_truth/fresh48/{case}.npz", allow_pickle=False) as saved:
            truth = saved["material_coarse"].ravel()
        true_support = truth > 1
        amplitude = np.max(truth) - 0.1
        with np.load(output / f"predictions/{row['stored']}.npz", allow_pickle=False) as saved:
            for view in ("raw", "corrected"):
                m = saved[f"{view}_m"]
                e = m - truth
                mask = m > 1
                scores.append(dict(observation=case, array=row["array"], arm=row["arm"], view=view,
                                   NRMSE=float(np.sqrt(np.mean(e**2))/amplitude),
                                   interior_NRMSE=float(np.sqrt(np.mean(e[true_support]**2))/amplitude),
                                   background_NRMSE=float(np.sqrt(np.mean(e[~true_support]**2))/amplitude),
                                   L2=float(np.linalg.norm(e)/np.linalg.norm(truth)),
                                   IoU=float((mask & true_support).sum()/(mask | true_support).sum()),
                                   port_residual_over_delta=row[f"{view}_residual_over_delta"]))
    (output / "scores.json").write_text(json.dumps(scores, ensure_ascii=False, indent=2)+"\n")
    with (output / "scores.csv").open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scores[0]))
        writer.writeheader()
        writer.writerows(scores)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("generate", "score", "all"))
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=HERE / "results/real_frontend")
    args = parser.parse_args()
    root, output = args.workspace.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.stage in ("generate", "all"):
        run(root, output)
    if args.stage in ("score", "all"):
        score(root, output)


if __name__ == "__main__":
    main()
