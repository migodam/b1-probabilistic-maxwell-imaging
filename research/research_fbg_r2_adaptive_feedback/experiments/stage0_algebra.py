#!/usr/bin/env python3
"""One deterministic finite-matrix audit of feedback gauge and A17 algebra.

This is complex linear algebra only: no Maxwell solve, optimizer, GPU, GT, or
precision sweep. Run with the local NumPy/SciPy environment if available.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "audit" / "STAGE0_MATH.json"


def relerr(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(1.0, np.linalg.norm(a), np.linalg.norm(b)))


def objective(m: np.ndarray, g: np.ndarray, e: np.ndarray, s: np.ndarray,
              y: np.ndarray, b: np.ndarray) -> tuple[float, dict[str, np.ndarray]]:
    n = len(m)
    x = np.diag(m)
    q = np.linalg.solve(np.eye(n) - g @ x, e)
    j = x @ q
    r = s @ j - y
    f = q + b @ j  # transformed TOTAL field, not a new incident source
    h = f.conj().T @ f
    k = np.linalg.solve(h, np.eye(e.shape[1]))
    psi = float(np.trace(r @ k @ r.conj().T).real)
    return psi, dict(x=x, q=q, j=j, r=r, f=f, h=h, k=k)


def run() -> dict:
    rng = np.random.default_rng(20261030)

    def cm(shape, scale=1.0):
        return scale * (rng.normal(size=shape) + 1j * rng.normal(size=shape)) / np.sqrt(2.0)

    n, ns, nr, d = 6, 2, 4, 3
    m = np.array([0.31, 0.47, 0.60, 0.72, 0.83, 1.02])
    g = cm((n, n), 0.12)
    e = cm((n, ns))
    s = cm((nr, n))
    b = cm((n, n), 0.11)
    y = cm((nr, ns), 0.3)  # fixed observation matrix, independent of m
    direction = np.array([0.41, -0.25, 0.19, 0.32, -0.27, 0.36])
    eps = 1e-6
    psi, z = objective(m, g, e, s, y, b)
    x, q, j, r, f, k = (z[key] for key in ("x", "q", "j", "r", "f", "k"))
    ident = np.eye(n)
    mf = ident + b @ x
    a = ident + x @ b
    xt = np.linalg.solve(a, x)
    gt = g + b
    qt = np.linalg.solve(ident - gt @ xt, e)
    jt = xt @ qt

    # Physical directional current and total-field derivatives, fixed B.
    l = ident - x @ g
    t = np.linalg.solve(l, np.diag(direction) @ q)
    qh = g @ t
    fh_direct = b @ np.diag(direction) @ q + mf @ qh
    fh_short = (g + b) @ t
    w = r @ k
    v = f @ k @ (r.conj().T @ r) @ k
    adj = np.linalg.solve(l.conj().T, s.conj().T @ w - (g + b).conj().T @ v)
    gradient = 2.0 * np.real(np.sum(adj.conj() * q, axis=1))
    fd = (objective(m + eps * direction, g, e, s, y, b)[0]
          - objective(m - eps * direction, g, e, s, y, b)[0]) / (2 * eps)
    analytic = float(np.dot(gradient, direction))

    # Source-span projector derivative; unlike the gauge current, this can move.
    pinv_f = k @ f.conj().T
    proj = f @ pinv_f
    dproj = (ident - proj) @ fh_short @ pinv_f
    dproj += pinv_f.conj().T @ fh_short.conj().T @ (ident - proj)
    fp = objective(m + eps * direction, g, e, s, y, b)[1]["f"]
    fm = objective(m - eps * direction, g, e, s, y, b)[1]["f"]
    pp = fp @ np.linalg.solve(fp.conj().T @ fp, fp.conj().T)
    pm = fm @ np.linalg.solve(fm.conj().T @ fm, fm.conj().T)
    dproj_fd = (pp - pm) / (2 * eps)

    # Unrestricted static B=-G removes all transformed-field material variation.
    cancel = -g
    psi_cancel, z_cancel = objective(m, g, e, s, y, cancel)
    cancel_gradient_field_norm = float(np.linalg.norm((g + cancel) @ t))
    cancel_fd = (objective(m + eps * direction, g, e, s, y, cancel)[0]
                 - objective(m - eps * direction, g, e, s, y, cancel)[0]) / (2 * eps)

    # Another pathology: B=tI makes the proposal objective vanish by scaling.
    psi_zero_b, _ = objective(m, g, e, s, y, np.zeros_like(g))
    psi_scaled, z_scaled = objective(m, g, e, s, y, 1e4 * ident)

    # Noncommutation and conditioning examples separate algebraic from safe.
    wrong_order = np.linalg.solve(ident + b @ x, x)
    near_b = np.zeros((n, n), dtype=complex)
    near_b[0, 0] = -(1.0 - 1e-6) / m[0]
    near_mf = ident + near_b @ x
    nilpotent = np.zeros((n, n), dtype=complex)
    nilpotent[0, 1] = 20.0
    nonnormal_b = nilpotent @ np.linalg.inv(x)
    nonnormal_mf = ident + nonnormal_b @ x

    # Invertible source mixing preserves psi and P, while raw Frobenius data
    # loss changes for this nonunitary C.
    mix = np.array([[1.3, 0.4j], [0.1, 0.7]], dtype=complex)
    fmix, rmix = f @ mix, r @ mix
    kmix = np.linalg.solve(fmix.conj().T @ fmix, np.eye(ns))
    psi_mix = float(np.trace(rmix @ kmix @ rmix.conj().T).real)
    proj_mix = fmix @ kmix @ fmix.conj().T
    # R4-style production scaling is not invariant when recomputed after an
    # arbitrary nonunitary source mix; a frozen original-acquisition factor is.
    m_ref = np.full(n, 0.55)
    q_ref = objective(m_ref, g, e, s, y, np.zeros_like(g))[1]["q"]
    factor = 0.5 * (np.linalg.norm(q_ref)**2 / ns) / np.linalg.norm(y)**2
    factor_mix = 0.5 * (np.linalg.norm(q_ref @ mix)**2 / ns) / np.linalg.norm(y @ mix)**2
    regularizer = 0.2
    normalized = factor * psi + regularizer
    normalized_mix_recomputed = factor_mix * psi_mix + regularizer
    normalized_mix_frozen = factor * psi_mix + regularizer

    # Actual A17 shared-anchor Schur/Galerkin exchange algebra. Its material
    # injection is a separate matrix, not the feedback-design B above.
    basis, _ = np.linalg.qr(cm((n, 4)))
    core, old, new = basis[:, :2], basis[:, 2:3], basis[:, 3:4]
    inj = cm((n, d))
    core_res = core @ np.linalg.solve(core.conj().T @ l @ core, core.conj().T)
    k0 = core_res @ inj

    def schur(c, injection):
        kbase = core_res @ injection
        cbar = (ident - core_res @ l) @ c
        yc = s @ cbar
        nc = c.conj().T @ (injection - l @ kbase)
        gamma = c.conj().T @ l @ cbar
        delta = yc @ np.linalg.solve(gamma, nc)
        delta_k = cbar @ np.linalg.solve(gamma, nc)
        u = np.column_stack((core, c))
        ku = u @ np.linalg.solve(u.conj().T @ l @ u, u.conj().T @ injection)
        return delta, delta_k, s @ ku, ku, gamma

    old_delta, old_delta_k, jo, ko, gamma_old = schur(old, inj)
    new_delta, new_delta_k, jn, kn, gamma_new = schur(new, inj)
    anchor = s @ k0
    exchange = new_delta - old_delta
    current_exchange = new_delta_k - old_delta_k
    # Finite A17 exchange: compute projector endpoints from the exact common-
    # anchor current difference and independently from direct Galerkin solves.
    # Use two of the three real material directions, so the tangent span does
    # not automatically equal the entire three-dimensional current trial space.
    def current_projector(a):
        a = a[:, :2]
        return a @ np.linalg.solve(a.conj().T @ a, a.conj().T)

    tangent_span_projector_delta_direct = current_projector(kn) - current_projector(ko)
    tangent_span_projector_delta_exchange = current_projector(ko + current_exchange) - current_projector(ko)

    # Requested A17-to-feedback bridge: one material direction at the SAME
    # exact physical Q anchor. Approximate only its current tangent through
    # old/new Galerkin spaces. The illumination projector P_Q is held at that
    # exact anchor, so its derivative is linear in the supplied current T.
    injection_h = direction[:, None] * q
    _, old_delta_h, _, t_old, _ = schur(old, injection_h)
    _, new_delta_h, _, t_new, _ = schur(new, injection_h)
    current_exchange_h = new_delta_h - old_delta_h
    qdagger = np.linalg.solve(q.conj().T @ q, q.conj().T)
    p_exact_q = q @ qdagger

    def illumination_projector_jvp(current_tangent):
        part = (ident - p_exact_q) @ g @ current_tangent @ qdagger
        return part + part.conj().T

    illumination_delta_direct = (illumination_projector_jvp(t_new)
                                 - illumination_projector_jvp(t_old))
    illumination_delta_exchange = illumination_projector_jvp(current_exchange_h)
    illumination_full_jvp = illumination_projector_jvp(t)
    # Real material coordinates: stack complex data rather than treating each
    # illumination as an independent material parameter.
    jo_real = np.vstack((jo.real, jo.imag))
    jn_real = np.vstack((jn.real, jn.imag))
    residual_real = rng.normal(size=2 * nr)
    lam = np.diag([0.5, 0.8, 1.1])
    so = -np.linalg.solve(jo_real.T @ jo_real + lam, jo_real.T @ residual_real)
    sn = -np.linalg.solve(jn_real.T @ jn_real + lam, jn_real.T @ residual_real)
    dd = jn_real - jo_real
    delta_info = jo_real.T @ dd + dd.T @ jo_real + dd.T @ dd
    delta_offset = dd.T @ residual_real
    step_formula = -np.linalg.solve(jn_real.T @ jn_real + lam,
                                    delta_offset + delta_info @ so)

    checks = {
        "gauge_current_relative_error": relerr(jt, j),
        "gauge_total_field_relative_error": relerr(qt, mf @ q),
        "gauge_material_order_identity_relative_error": relerr(xt, x @ np.linalg.inv(ident + b @ x)),
        "noncommuting_wrong_order_gap": relerr(xt, wrong_order),
        "feedback_directional_identity_relative_error": relerr(fh_direct, fh_short),
        "gradient_directional_absolute_error": abs(analytic - fd),
        "gradient_directional_relative_error": abs(analytic - fd) / max(1.0, abs(analytic), abs(fd)),
        "projector_directional_relative_error": relerr(dproj, dproj_fd),
        "projector_directional_norm": float(np.linalg.norm(dproj)),
        "source_mixing_psi_relative_error": abs(psi_mix - psi) / max(1.0, abs(psi)),
        "source_mixing_projector_relative_error": relerr(proj_mix, proj),
        "source_mixing_raw_loss_ratio": float(np.linalg.norm(rmix)**2 / np.linalg.norm(r)**2),
        "source_mixing_recomputed_R4_factor_ratio": float(factor_mix / factor),
        "source_mixing_recomputed_normalized_objective_relative_change": float(
            abs(normalized_mix_recomputed - normalized) / max(1.0, abs(normalized))),
        "source_mixing_frozen_factor_objective_relative_error": float(
            abs(normalized_mix_frozen - normalized) / max(1.0, abs(normalized))),
        "cancel_field_minus_incident_relative_error": relerr(z_cancel["f"], e),
        "cancel_field_directional_norm": cancel_gradient_field_norm,
        "cancel_objective_directional_fd": float(cancel_fd),
        "scaled_psi_ratio_to_original": float(psi_scaled / psi_zero_b),
        "scaled_physical_residual_relative_error": relerr(z_scaled["r"], r),
        "near_pole_Mf_condition": float(np.linalg.cond(near_mf)),
        "nonnormal_BX_spectral_radius": float(max(abs(np.linalg.eigvals(nonnormal_b @ x)))),
        "nonnormal_Mf_condition": float(np.linalg.cond(nonnormal_mf)),
        "a17_old_schur_relative_error": relerr(jo, anchor + old_delta),
        "a17_new_schur_relative_error": relerr(jn, anchor + new_delta),
        "a17_exchange_relative_error": relerr(jn - jo, exchange),
        "a17_step_difference_relative_error": relerr(sn - so, step_formula),
        "a17_old_schur_condition": float(np.linalg.cond(gamma_old)),
        "a17_new_schur_condition": float(np.linalg.cond(gamma_new)),
    }
    # Direct independent CPU check of the production helper's Torch algebra.
    # It invokes no CUDA lane or Maxwell workload; the comparison uses the
    # NumPy forward, finite difference, and adjoint above as independent oracle.
    import torch
    sys.path.insert(0, str(ROOT.parents[1]))
    from research.research_fbg_r2_adaptive_feedback.src.feedback_fit import (
        FeedbackOperator, feedback_value_gradient, fit_feedback)

    tc = lambda arr: torch.as_tensor(arr, dtype=torch.complex128)
    tr = torch.as_tensor(m, dtype=torch.float64)
    tg, te, ts, ty = map(tc, (g, e, s, y))
    tl = torch.eye(n, dtype=torch.complex128) - tr[:, None] * tg
    lu, pivots = torch.linalg.lu_factor(tl)
    tj = torch.linalg.lu_solve(lu, pivots, tr[:, None] * te)
    tq = te + tg @ tj
    torch_state = {"m": tr, "q": tq, "r": ts @ tj - ty, "lu": lu, "pivots": pivots}
    torch_shared = SimpleNamespace(GH=tg.mH, SH=ts.mH)
    torch_feedback = FeedbackOperator.dense(tc(b))
    tv, tgrad, tinfo = feedback_value_gradient(torch_shared, torch_state, torch_feedback, 1.0)
    checks["production_helper_value_relative_error"] = abs(float(tv) - psi) / max(1.0, psi)
    checks["production_helper_gradient_relative_error"] = relerr(tgrad.detach().numpy(), gradient)
    checks["production_helper_cholesky_info"] = int(tinfo)
    checks["a17_current_exchange_relative_error"] = relerr(kn - ko, current_exchange)
    checks["a17_tangent_span_projector_delta_bridge_relative_error"] = relerr(
        tangent_span_projector_delta_direct, tangent_span_projector_delta_exchange)
    checks["a17_tangent_span_projector_delta_norm"] = float(np.linalg.norm(
        tangent_span_projector_delta_direct))
    checks["a17_illumination_DP_Q_current_exchange_relative_error"] = relerr(
        t_new - t_old, current_exchange_h)
    checks["a17_illumination_DP_Q_bridge_relative_error"] = relerr(
        illumination_delta_direct, illumination_delta_exchange)
    checks["a17_illumination_DP_Q_delta_norm"] = float(np.linalg.norm(
        illumination_delta_direct))
    checks["a17_illumination_DP_Q_full_tail_relative"] = float(
        np.linalg.norm(illumination_full_jvp - illumination_projector_jvp(t_new))
        / max(np.linalg.norm(illumination_full_jvp), 1e-300))

    # One n=16 production fit audit: normal equations for complex diagonal and
    # disjoint 2x2 blocks, plus exact source-basis invariance of the whitened
    # ridge fit. All columns and directions are from one retained physical L.
    nf, side, nfit_dir = 16, 4, 4
    gf, ef = cm((nf, nf), 0.06), cm((nf, ns))
    mf = np.linspace(0.23, 0.97, nf)
    lf = np.eye(nf) - mf[:, None] * gf
    jf = np.linalg.solve(lf, mf[:, None] * ef)
    qf = ef + gf @ jf
    directions = rng.normal(size=(nfit_dir, nf))
    directions /= np.max(np.abs(directions), axis=1, keepdims=True)
    tf = np.stack([np.linalg.solve(lf, h[:, None] * qf) for h in directions])
    fit_diag = fit_feedback(gf, qf, tf, kind="complex", ridge=0.08)
    fit_block = fit_feedback(gf, qf, tf, kind="block", side=side, ridge=0.08)
    fit_diag_mix = fit_feedback(gf, qf @ mix, tf @ mix, kind="complex", ridge=0.08)
    fit_block_mix = fit_feedback(gf, qf @ mix, tf @ mix, kind="block", side=side, ridge=0.08)
    cf = np.linalg.cholesky(qf.conj().T @ qf)
    whiten = np.linalg.solve(cf.conj().T, np.eye(ns))
    fit_t = (tf @ whiten).transpose(1, 0, 2).reshape(nf, -1)
    fit_y = gf @ fit_t
    diag_num = np.sum(fit_y * fit_t.conj(), axis=1)
    diag_den = np.sum(np.abs(fit_t)**2, axis=1) + fit_diag.fit_lambda
    diag_normal_res = fit_diag.data * diag_den + diag_num
    block_res = []
    for block, idx in zip(fit_block.data, fit_block.indices):
        ti, yi = fit_t[idx], fit_y[idx]
        block_res.append(block @ (ti @ ti.conj().T) + fit_block.fit_lambda * block
                         + yi @ ti.conj().T)
    checks["production_fit_diag_normal_relative_residual"] = float(
        np.linalg.norm(diag_normal_res) / max(1.0, np.linalg.norm(diag_num)))
    checks["production_fit_block_normal_relative_residual"] = float(
        np.linalg.norm(block_res) / max(1.0, np.linalg.norm(fit_block.data)))
    checks["production_fit_diag_source_mixing_relative_error"] = relerr(
        fit_diag.data, fit_diag_mix.data)
    checks["production_fit_block_source_mixing_relative_error"] = relerr(
        fit_block.data, fit_block_mix.data)
    checks["production_fit_lambda_source_mixing_relative_error"] = max(
        abs(fit_diag.fit_lambda - fit_diag_mix.fit_lambda),
        abs(fit_block.fit_lambda - fit_block_mix.fit_lambda)) / max(1.0, fit_diag.fit_lambda,
                                                                   fit_block.fit_lambda)
    tolerances = {
        "algebra_relative": 1e-10,
        "directional_relative": 1e-5,
        "gradient_absolute": 1e-5,
    }
    passes = {
        "gauge_and_order": max(checks[k] for k in (
            "gauge_current_relative_error", "gauge_total_field_relative_error",
            "gauge_material_order_identity_relative_error")) < tolerances["algebra_relative"],
        "noncommutation_witness": checks["noncommuting_wrong_order_gap"] > 1e-3,
        "directional_identity": checks["feedback_directional_identity_relative_error"] < tolerances["algebra_relative"],
        "gradient": checks["gradient_directional_relative_error"] < tolerances["directional_relative"]
                    and checks["gradient_directional_absolute_error"] < tolerances["gradient_absolute"],
        "projector_derivative": checks["projector_directional_relative_error"] < tolerances["directional_relative"]
                                and checks["projector_directional_norm"] > 1e-4,
        "source_basis_invariance": max(checks["source_mixing_psi_relative_error"],
                                       checks["source_mixing_projector_relative_error"]) < tolerances["algebra_relative"],
        "raw_loss_basis_caveat": abs(checks["source_mixing_raw_loss_ratio"] - 1.0) > 1e-2,
        "production_normalization_caveat": abs(checks["source_mixing_recomputed_R4_factor_ratio"] - 1.0) > 1e-2
                                           and checks["source_mixing_frozen_factor_objective_relative_error"]
                                           < tolerances["algebra_relative"],
        "cancel_counterexample": checks["cancel_field_minus_incident_relative_error"] < tolerances["algebra_relative"]
                                 and checks["cancel_field_directional_norm"] < tolerances["algebra_relative"]
                                 and abs(checks["cancel_objective_directional_fd"]) > 1e-4,
        "scale_degeneracy": checks["scaled_psi_ratio_to_original"] < 1e-5
                            and checks["scaled_physical_residual_relative_error"] < tolerances["algebra_relative"],
        "a17_schur_exchange_step": max(checks[k] for k in (
            "a17_old_schur_relative_error", "a17_new_schur_relative_error",
            "a17_exchange_relative_error", "a17_step_difference_relative_error")) < tolerances["algebra_relative"],
        "a17_tangent_span_projector_delta_bridge": max(
            checks["a17_current_exchange_relative_error"],
            checks["a17_tangent_span_projector_delta_bridge_relative_error"])
            < tolerances["algebra_relative"]
            and checks["a17_tangent_span_projector_delta_norm"] > 1e-4,
        "a17_illumination_DP_Q_bridge": max(
            checks["a17_illumination_DP_Q_current_exchange_relative_error"],
            checks["a17_illumination_DP_Q_bridge_relative_error"])
            < tolerances["algebra_relative"]
            and checks["a17_illumination_DP_Q_delta_norm"] > 1e-4,
        "production_helper_CPU": max(checks["production_helper_value_relative_error"],
                                     checks["production_helper_gradient_relative_error"])
                                 < tolerances["algebra_relative"]
                                 and checks["production_helper_cholesky_info"] == 0,
        "production_fit_normal_equations": max(checks["production_fit_diag_normal_relative_residual"],
                                               checks["production_fit_block_normal_relative_residual"])
                                           < tolerances["algebra_relative"],
        "production_fit_source_mixing": max(checks["production_fit_diag_source_mixing_relative_error"],
                                            checks["production_fit_block_source_mixing_relative_error"],
                                            checks["production_fit_lambda_source_mixing_relative_error"])
                                        < tolerances["algebra_relative"],
    }
    result = {
        "status": "PASS" if all(passes.values()) else "FAIL",
        "scope": "one deterministic small complex matrix audit; no Maxwell physics or imaging inference",
        "seed": 20261030,
        "dimensions": {"core_state": n, "sources": ns, "receivers": nr,
                       "real_material_directions_A17": d, "production_fit_state": nf,
                       "production_fit_directions": nfit_dir, "production_fit_block_side": side},
        "finite_difference_step": eps,
        "checks": checks,
        "passes": passes,
        "tolerances": tolerances,
        "objective_values": {"feedback": psi, "B_zero_original": psi_zero_b,
                             "B_minus_G": psi_cancel, "B_large_positive": psi_scaled},
        "safety_boundary": "conditioning examples demonstrate risk, not a general admissible B set",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "passes": passes,
                      "audit": str(OUT)}, indent=2))
    return result


if __name__ == "__main__":
    run()
