"""Frozen joint-GN profile algebra and coordinate-invariant soft release."""
from __future__ import annotations

import numpy as np
from scipy import linalg, optimize

try:
    from .current_normal import thin_basis
except ImportError:
    from current_normal import thin_basis


def solve_symmetric(H, b):
    try:
        return linalg.cho_solve(linalg.cho_factor(H, check_finite=False), b, check_finite=False)
    except linalg.LinAlgError:
        return linalg.lstsq(H, b, cond=1e-12, lapack_driver="gelsy", check_finite=False)[0]


def least_squares(A, b):
    if A.shape[1] == 0:
        return np.empty(0)
    if A.shape[0] == 0:
        return np.zeros(A.shape[1])
    return linalg.lstsq(A, b, cond=1e-11, lapack_driver="gelsy", check_finite=False)[0]


def box_quadratic_step(H, g, m=None, bounds=None):
    """Solve the actual constrained QP; never assign an unclipped score to clip."""
    p = -solve_symmetric(H, g)
    if m is None or bounds is None:
        return p
    lo, hi = bounds
    lb = np.broadcast_to(lo, np.shape(m)) - m
    ub = np.broadcast_to(hi, np.shape(m)) - m
    if np.all(p >= lb) and np.all(p <= ub):
        return p
    result = optimize.minimize(lambda h: 0.5 * h @ H @ h + g @ h,
                               np.clip(p, lb, ub), jac=lambda h: H @ h + g,
                               bounds=list(zip(lb, ub)), method="L-BFGS-B",
                               options={"maxiter": 180, "ftol": 1e-13, "gtol": 1e-9, "maxls": 30})
    return result.x


def profile_schur(b, B, A, Lambda=None, ell=None, lm=0.0, m=None, bounds=None, cost=None):
    b, B, A = np.asarray(b, float), np.asarray(B, float), np.asarray(A, float)
    n = B.shape[1]
    Lambda = np.zeros((n, n)) if Lambda is None else np.asarray(Lambda, float)
    ell = np.zeros(n) if ell is None else np.asarray(ell, float)
    Q, _ = thin_basis(A)
    residual_B = B - Q @ (Q.T @ B)
    residual_b = b - Q @ (Q.T @ b)
    H = B.T @ residual_B + Lambda + float(lm) * np.eye(n)
    H = 0.5 * (H + H.T)
    g = B.T @ residual_b + ell
    p = box_quadratic_step(H, g, m, bounds)
    c = -least_squares(A, b + B @ p)
    if cost:
        cost.add("rank_revealing_qr")
        cost.add("material_qp_solves")
        cost.add("small_schur_solves")
    return {"b": b, "B": B, "A": A, "Q_old": Q, "H": H, "g": g,
            "p": p, "c": c, "Lambda": Lambda, "ell": ell, "lm": lm,
            "m": m, "bounds": bounds, "real_dof": A.shape[1], "residual_rank": Q.shape[1]}


def release_candidate(profile, A_C, release="soft", eta_target=0.1, gamma=None, cost=None):
    """Release in residual-QR coordinates u=R c, penalty gamma*||u||^2/2.

    R retains original candidate coefficient coordinates, so C T + U K is a
    true reparameterization.  Nonzero basis scaling is never used as release.
    """
    A_C = np.asarray(A_C, float)
    Q0 = profile["Q_old"]
    Z = A_C - Q0 @ (Q0.T @ A_C)
    Q, R = thin_basis(Z, reference_scale=float(np.linalg.norm(A_C)))
    D, v = Q.T @ profile["B"], Q.T @ profile["b"]
    if len(v):
        DHinvDt = D @ solve_symmetric(profile["H"], D.T)
        eta0 = max(0.0, float(linalg.eigvalsh(0.5 * (DHinvDt + DHinvDt.T), check_finite=False)[-1]))
    else:
        eta0 = 0.0
    if release == "hard":
        gamma = 0.0
    elif release == "none":
        gamma = np.inf
    elif release == "soft":
        if not 0.0 < eta_target < 1.0:
            raise ValueError("eta_target must lie in (0,1)")
        gamma = max(0.0, eta0 / eta_target - 1.0) if gamma is None else float(gamma)
        if gamma < 0:
            raise ValueError("soft gamma must be nonnegative")
    else:
        raise ValueError(f"Unknown release: {release}")
    alpha = 1.0 / (1.0 + gamma)
    H = profile["H"] - alpha * D.T @ D
    g = profile["g"] - alpha * D.T @ v
    p = box_quadratic_step(H, g, profile["m"], profile["bounds"])
    u = -alpha * (v + D @ p)
    c_C = least_squares(R, u)
    c_U = -least_squares(profile["A"], profile["b"] + profile["B"] @ p + A_C @ c_C)
    dp = p - profile["p"]
    Hnorm = float(np.sqrt(max(0.0, dp @ profile["H"] @ dp)))
    bound = (alpha * np.sqrt(eta0) / (1.0 - alpha * eta0) * np.linalg.norm(v + D @ profile["p"])) if alpha * eta0 < 1.0 else None
    if cost:
        cost.add("rank_revealing_qr")
        cost.add("material_qp_solves")
        cost.add("small_schur_solves")
    return {"p": p, "delta_p": dp, "c_U": c_U, "c_C": c_C,
            "H": H, "g": g, "D": D, "v": v, "Q": Q, "R": R,
            "eta0": eta0, "eta": alpha * eta0, "gamma": gamma, "alpha": alpha,
            "real_dof": A_C.shape[1], "residual_rank": Q.shape[1],
            "duplicate": Q.shape[1] == 0, "step_Hnorm": Hnorm,
            "unconstrained_step_bound": bound, "A_C": A_C}


def replacement_steps(core_profile, A_old, A_new, **kwargs):
    """Same-capacity replacement, both increments relative to a retained core."""
    old = release_candidate(core_profile, A_old, **kwargs)
    new = release_candidate(core_profile, A_new, **kwargs)
    if old["real_dof"] != new["real_dof"]:
        raise ValueError("Replacement requires equal actual real current DOF")
    return {"old": old, "new": new, "delta_p": new["p"] - old["p"],
            "delta_H": new["H"] - old["H"], "delta_g": new["g"] - old["g"]}


def physical_task_gain(r, Jp, Jdp, p, dp, ell=None, Lambda=None):
    gain = -float((r + Jp) @ Jdp) - 0.5 * float(Jdp @ Jdp)
    if ell is not None:
        gain -= float(np.asarray(ell) @ dp)
    if Lambda is not None:
        gain -= float(p @ Lambda @ dp + 0.5 * dp @ Lambda @ dp)
    return gain
