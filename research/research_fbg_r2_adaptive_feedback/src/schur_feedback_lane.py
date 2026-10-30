"""External exploratory Schur-current feedback constructor for one Q9 arm.

This uses the common-core scalar-TM Galerkin/Schur algebra associated with
A17-style current reduction, but it does not implement or certify A17's full
Jacobian method. It changes only the fixed diagonal feedback fitted at each
public anchor. The parent adaptive solver's true Maxwell forward, objective,
adjoint, physical corrector and budget logic remain unchanged.
"""
import time

import numpy as np
import torch

from .feedback_fit import fit_feedback
from .adaptive_solver import _spatial_directions


def _projector_tangent_norm_squared(Q, dQ):
    """||dP_Q||_F²=2 tr(U*U H^-1), U=(I-P_Q)dQ, at fixed Q."""
    H = Q.mH @ Q
    C = torch.linalg.cholesky(H)
    right = torch.linalg.solve_triangular(C.mH, torch.eye(H.shape[0], dtype=H.dtype, device=H.device), upper=True)
    Qw, dQw = Q @ right, dQ @ right
    Hw = Qw.mH @ Qw
    U = dQw - Qw @ torch.linalg.solve(Hw, Qw.mH @ dQw)
    return 2 * torch.trace(torch.linalg.solve(Hw, U.mH @ U)).real


def make_schur_feedback(lane, x, kind, cfg):
    """Return ``(fixed B, JSON details)`` for ``kind='schur_diag'``.

    At the shared public anchor, U=orth(Q), C=orth((I-UU*)[GQ,G²Q]). With
    L=I-diag(m)G and F=diag(h)Q for seven public spatial directions, let
    A=U*LU, Tret=U A^-1 U*F, Cbar=C-U A^-1 U*LC, Gamma=C*L Cbar.
    Then Taug=Tret+Cbar Gamma^-1 C*(F-L Tret). Small Galerkin systems are
    solved directly, without a numerical epsilon or full tangent solve. One
    sampled exact p-source tangent uses the anchor's retained Maxwell LU only
    for a diagnostic tail. The fit itself uses Taug and the same source-
    normalized .05 ridge and whole-box safety as the full tangent arm.
    """
    if kind != "schur_diag":
        raise ValueError("Schur constructor only handles schur_diag")
    s = lane.shared
    started = time.perf_counter()
    with torch.cuda.stream(lane.stream), torch.no_grad():
        z = lane._forward(x)
        m_host = lane._pull(z["m"])
        Q, L = z["q"], z["a"]
        p = s.ns
        hs = torch.as_tensor(_spatial_directions(s.side, m_host), dtype=torch.float64, device=s.device)
        ndirections = len(hs)
        F = (hs[:, :, None] * Q[None, :, :]).permute(1, 0, 2).reshape(s.n, -1)
        U, Ru = torch.linalg.qr(Q, mode="reduced")
        GQ = s.G @ Q
        G2Q = s.G @ GQ
        candidate = torch.cat((GQ, G2Q), dim=1)
        candidate = candidate - U @ (U.mH @ candidate)
        C, Rc = torch.linalg.qr(candidate, mode="reduced")
        LU, LC = L @ U, L @ C
        A = U.mH @ LU
        A_for_rhs = torch.linalg.solve(A, U.mH @ F)
        A_for_C = torch.linalg.solve(A, U.mH @ LC)
        Tret = U @ A_for_rhs
        Cbar = C - U @ A_for_C
        Gamma = C.mH @ (L @ Cbar)
        schur_rhs = C.mH @ (F - L @ Tret)
        Taug = Tret + Cbar @ torch.linalg.solve(Gamma, schur_rhs)
        tangents = Taug.reshape(s.n, ndirections, p).permute(1, 0, 2)
        # Only the first public direction receives an exact tangent solve.
        Texact = torch.linalg.lu_solve(z["lu"], z["pivots"], F[:, :p])
        Tapprox = Taug[:, :p]
        current_tail = torch.linalg.vector_norm(Tapprox - Texact) / torch.linalg.vector_norm(Texact)
        dQ_exact = s.G @ Texact
        dQ_error = s.G @ (Tapprox - Texact)
        projector_tail = torch.sqrt(_projector_tangent_norm_squared(Q, dQ_error) /
                                    _projector_tangent_norm_squared(Q, dQ_exact))
        feedback = fit_feedback(s.G, Q, tangents, kind="complex", side=s.side,
                                ridge=float(cfg.get("feedback_ridge", .05)))
        fit_lambda = float(feedback.fit_lambda)
        feedback.name = kind
        feedback, box_receipt = feedback.box_certify(lower_m=.02, upper_m=20., perturbation_cap=.8)
        reference = torch.full((s.n,), .1, dtype=torch.float64, device=s.device)
        reference_q = s.feedback_reference_q
        reference_v = feedback.field(reference_q, reference)
        feedback.reference_scale = (reference_v.abs().square().sum() /
                                    reference_q.abs().square().sum()).real
        safety = feedback.safety(z["m"], lower=.02, upper=20.)
        G_tangents = s.G @ tangents
        base_energy = G_tangents.abs().square().sum()
        feedback_energy = (G_tangents + feedback.apply(tangents)).abs().square().sum()
        metrics = torch.stack((current_tail.real, projector_tail.real,
                               (feedback_energy / base_energy).real,
                               (torch.diagonal(Ru).abs().min() / torch.diagonal(Ru).abs().max()).real,
                               (torch.diagonal(Rc).abs().min() / torch.diagonal(Rc).abs().max()).real))
        sampled_current_tail, sampled_projector_tail, energy, u_qr_ratio, c_qr_ratio = map(float, lane._pull(metrics))
        lane.counts.update({"schur_setups": 1, "schur_QR": 2,
                            "schur_small_solves": 3,
                            "schur_small_RHS": int((2 * p) + (ndirections * p) * 2),
                            "schur_ancillary_small_solves": 7,
                            "schur_ancillary_small_RHS": 7 * p,
                            "schur_G_products": 5,
                            "schur_source_Gram_Cholesky": 3,
                            "tangent_RHS": p, "sampled_exact_tangent_solves": 1,
                            "feedback_fits": 1, "feedback_reference": 1})
        details = {"kind": kind, "approximation": "common_core_U_plus_2p_Schur_current",
                   "common_core_rank": p, "candidate_rank": 2 * p,
                   "directions": ndirections, "sampled_exact_tangent_direction": 0,
                   "sampled_current_tail_relative": sampled_current_tail,
                   "sampled_projector_linear_map_tail_relative": sampled_projector_tail,
                   "U_QR_minmax_diagonal_ratio": u_qr_ratio,
                   "C_QR_minmax_diagonal_ratio": c_qr_ratio,
                   "tangent_energy_ratio": energy,
                   "fit_lambda": fit_lambda,
                   "safety": safety, "whole_box": box_receipt,
                   "fit_wall_s": time.perf_counter() - started,
                   "schur_costs": {"QR": 2, "small_solves": 3,
                                   "small_RHS": int((2 * p) + (ndirections * p) * 2),
                                   "ancillary_projector_and_fit_small_solves": 7,
                                   "ancillary_projector_and_fit_small_RHS": 7 * p,
                                   "sampled_exact_tangent_RHS": p,
                                   "source_Gram_Cholesky": 3},
                   "normalization": "fixed public material .1; same source-normalized ridge and prior as full diag arm",
                   "scientific_scope": "approximate current fit only; no A17 full-Jacobian certificate or basin claim"}
        feedback.details = details
        return feedback, details
