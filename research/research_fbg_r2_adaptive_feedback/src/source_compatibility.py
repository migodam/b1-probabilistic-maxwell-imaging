"""Public-source row compatibility diagnostics for finite current proposals.

At a public anchor m, L=I-diag(m)G, J=L^-1 diag(m)E, Q=E+GJ and
T=L^-1 diag(h)Q. A finite current proposal w=tT is representable by one
shared *real* material increment u only if
    Lw=diag(u)(Q+Gw).
This finite condition is stronger than exact first-order tangency. It is a
mechanism diagnostic inspired by source-row compatibility, not an
implementation or performance claim about A17-X or a recovery algorithm.
"""
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve, solve_triangular


def analyze(G, E, h, *, anchor=0.1, step=0.5, real_pass_tolerance=1e-6):
    """Analyze one acquisition with one physical LU and two source RHS blocks.

    Sources are whitened by the exact Cholesky H=Q*Q at the anchor. The same
    invertible right mixing is applied to Q and the current tangent T, so the
    row compatibility and projector geometry do not depend on source basis.
    No Cholesky jitter or hidden rank regularization is used. Zero-excitation
    rows remain explicit in the receipt.
    """
    start = time.perf_counter()
    G, E, h = np.asarray(G), np.asarray(E), np.asarray(h, dtype=np.float64)
    n, p = E.shape
    if G.shape != (n, n) or h.shape != (n,):
        raise ValueError("operator and direction shapes are inconsistent")
    m = np.full(n, anchor, dtype=np.float64)
    L = np.eye(n, dtype=np.complex128) - m[:, None] * G
    fact = lu_factor(L, check_finite=False)
    J = lu_solve(fact, m[:, None] * E, check_finite=False)
    Q = E + G @ J
    T = lu_solve(fact, h[:, None] * Q, check_finite=False)
    costs = {"physical_LU": 1, "initial_source_RHS": p,
             "tangent_source_RHS": p, "source_Gram_Cholesky": 1,
             "full_Maxwell_candidate_checks": 0}

    field_gram = Q.conj().T @ Q
    eigenvalues = np.linalg.eigvalsh(field_gram)
    gram_cond = float(eigenvalues[-1] / eigenvalues[0]) if eigenvalues[0] > 0 else None
    chol = np.linalg.cholesky(field_gram)
    # Qw=Q chol(H)^-*; triangular solve avoids an explicit inverse.
    right = solve_triangular(chol.conj(), np.eye(p, dtype=np.complex128), lower=True, check_finite=False).T
    Qw, Tw = Q @ right, T @ right
    whitened_gram = Qw.conj().T @ Qw
    dQ = G @ Tw
    # Low-rank exact tangent dP=(I-P)dQ Q+ + adjoint. Since Qw*Qw≈I,
    # ||dP||_F²=2 tr(U*U Hwhite^-1), U=(I-P)dQ.
    projection_coeff = np.linalg.solve(whitened_gram, Qw.conj().T @ dQ)
    U = dQ - Qw @ projection_coeff
    projector_derivative_norm_squared = 2 * float(np.trace(np.linalg.solve(
        whitened_gram, U.conj().T @ U)).real)

    w = step * Tw
    V = Qw + G @ w
    Lw = L @ w
    v_power = np.sum(np.abs(V) ** 2, axis=1)
    q_power = np.sum(np.abs(Qw) ** 2, axis=1)
    gt_power = np.sum(np.abs(dQ) ** 2, axis=1)
    zero_q = np.flatnonzero(q_power == 0)
    zero_v = np.flatnonzero(v_power == 0)
    valid_v = v_power > 0
    u = np.zeros(n, dtype=np.float64)
    u[valid_v] = np.real(np.sum(V[valid_v].conj() * Lw[valid_v], axis=1)) / v_power[valid_v]
    real_defect = Lw - u[:, None] * V
    finite_compatibility_residual = float(np.linalg.norm(real_defect) / np.linalg.norm(Lw))
    complex_u = np.zeros(n, dtype=np.complex128)
    complex_u[valid_v] = np.sum(V[valid_v].conj() * Lw[valid_v], axis=1) / v_power[valid_v]
    complex_defect = Lw - complex_u[:, None] * V
    complex_compatibility_residual = float(np.linalg.norm(complex_defect) / np.linalg.norm(Lw))

    valid_q = q_power > 0
    gamma = np.zeros(n, dtype=np.complex128)
    gamma[valid_q] = np.sum(Qw[valid_q].conj() * dQ[valid_q], axis=1) / q_power[valid_q]
    row_propagation_defect = dQ.copy()
    row_propagation_defect[valid_q] -= gamma[valid_q, None] * Qw[valid_q]
    propagation_dispersion = float(np.linalg.norm(row_propagation_defect) / np.linalg.norm(dQ))
    rational_denominator = 1 + step * gamma[valid_q]
    near_zero_rational = np.abs(rational_denominator) <= 1e-12
    rational_u = np.full(n, np.nan + 1j * np.nan)
    usable_q_indices = np.flatnonzero(valid_q)[~near_zero_rational]
    rational_u[usable_q_indices] = step * h[usable_q_indices] / (
        1 + step * gamma[usable_q_indices])
    finite_rational = rational_u[np.isfinite(rational_u)]
    rational_imag_relative = (float(np.linalg.norm(finite_rational.imag) /
                                    np.linalg.norm(finite_rational))
                              if finite_rational.size and np.linalg.norm(finite_rational) else None)
    apparent_real_pass = bool(finite_compatibility_residual <= real_pass_tolerance
                              and len(zero_q) == 0 and len(zero_v) == 0)
    return {
        "sources": p, "cells": n, "anchor_material": anchor, "finite_step": step,
        "direction_min": float(h.min()), "direction_max": float(h.max()),
        "field_gram_condition": gram_cond,
        "source_whitened_gram_condition": float(np.linalg.cond(whitened_gram)),
        "projector_derivative_frobenius_norm_squared": projector_derivative_norm_squared,
        "projector_derivative_frobenius_norm": float(np.sqrt(projector_derivative_norm_squared)),
        "source_row_real_compatibility_residual": finite_compatibility_residual,
        "source_row_complex_compatibility_residual": complex_compatibility_residual,
        "complex_row_gamma_dispersion": propagation_dispersion,
        "gamma_imag_max_abs": float(np.max(np.abs(gamma[valid_q].imag))) if valid_q.any() else None,
        "gamma_imag_rms": float(np.sqrt(np.mean(gamma[valid_q].imag ** 2))) if valid_q.any() else None,
        "rational_increment_imag_relative": rational_imag_relative,
        "rational_near_zero_denominator_rows": int(near_zero_rational.sum()),
        "recovered_real_increment_min": float(u.min()),
        "recovered_real_increment_max": float(u.max()),
        "recovered_real_increment_mean": float(u.mean()),
        "zero_excitation_Q_row_count": int(len(zero_q)),
        "zero_excitation_Q_row_indices": zero_q.tolist(),
        "zero_finite_V_row_count": int(len(zero_v)),
        "zero_finite_V_row_indices": zero_v.tolist(),
        "apparent_real_pass": apparent_real_pass,
        "real_pass_tolerance": real_pass_tolerance,
        "costs": costs, "wall_seconds": time.perf_counter() - start,
        "interpretation_scope": "local public-anchor mechanism; no algorithm, A17 implementation, or recovery claim"
    }
