"""Shared-real-material current geometry; no ground-truth inputs.

A mode tensor has shape (n, p, d) and d **real** coefficients.  A complex
spatial atom with independent coefficients for p sources has d = 2*p.
"""
from __future__ import annotations

import numpy as np
from scipy import linalg


def real_pack(z):
    z = np.asarray(z)
    return np.concatenate((z.real.ravel(), z.imag.ravel()))


def real_unpack(x, shape):
    x = np.asarray(x, dtype=float)
    k = int(np.prod(shape))
    return (x[:k] + 1j * x[k:]).reshape(shape)


def spatial_modes(U, p):
    """Lift spatial complex columns to all independent real source coordinates."""
    U = np.asarray(U, dtype=complex)
    if U.ndim == 1:
        U = U[:, None]
    n, k = U.shape
    V = np.zeros((n, p, 2 * k * p), dtype=complex)
    for a in range(k):
        for s in range(p):
            V[:, s, a * p + s] = U[:, a]
            V[:, s, k * p + a * p + s] = 1j * U[:, a]
    return V


def mode_matrix(V, pack=real_pack):
    V = np.asarray(V, dtype=complex)
    return np.column_stack([pack(V[:, :, k]) for k in range(V.shape[2])]) if V.shape[2] else np.empty((2 * V.shape[0] * V.shape[1], 0))


def current_action(V, c):
    return np.einsum("npd,d->np", V, np.asarray(c, dtype=float), optimize=True)


def thin_basis(A, rtol=1e-11, reference_scale=None):
    """Rank revealing QR, with R in the original coefficient coordinates."""
    A = np.asarray(A, dtype=float)
    if A.shape[1] == 0:
        return np.empty((A.shape[0], 0)), np.empty((0, 0))
    Q, R, _ = linalg.qr(A, mode="economic", pivoting=True, check_finite=False)
    diag = np.abs(np.diag(R))
    scale = max(float(diag.max(initial=0.0)), float(reference_scale or 0.0))
    rank = int(np.count_nonzero(diag > rtol * scale))
    Q = Q[:, :rank]
    return Q, Q.T @ A


def normal_projection(Q, V, node_tol=1e-14):
    """Pi_Q^{perp,R} V, for a current (n,p) or tensor (n,p,d).

    Nodes have no defined one-dimensional material direction; leave their
    residual unchanged and return their mask rather than adding an epsilon.
    """
    Q = np.asarray(Q, complex)
    V = np.asarray(V, complex)
    power = np.sum(np.abs(Q) ** 2, axis=1)
    nodes = power <= node_tol * max(1.0, float(power.max(initial=0.0)))
    denominator = np.where(nodes, 1.0, power)
    if V.ndim == 2:
        coeff = np.real(np.sum(Q.conj() * V, axis=1)) / denominator
        coeff[nodes] = 0.0
        out = V - coeff[:, None] * Q
    else:
        coeff = np.real(np.einsum("np,npd->nd", Q.conj(), V)) / denominator[:, None]
        coeff[nodes] = 0.0
        out = V - Q[:, :, None] * coeff[:, None, :]
    return out, nodes


def normal_diagnostic(G, m, Q, V, beta=0.0, gradient=None, cost=None):
    LV = V - np.asarray(m)[:, None] * (G @ V)
    if cost:
        cost.add("green_vector_actions", V.shape[1])
    N, nodes = normal_projection(Q, LV)
    total = float(np.vdot(LV, LV).real)
    normal = float(np.vdot(N, N).real)
    power = np.sum(np.abs(Q) ** 2, axis=1)
    bias = np.zeros(len(m))
    ok = ~nodes
    if beta:
        bias[ok] = beta * np.sum(np.abs(N[ok]) ** 2, axis=1) / (power[ok] * (1.0 + beta * np.asarray(m)[ok]))
    return {"normal_fraction": normal / total if total > 1e-28 else None,
            "normal_energy": normal, "state_action_energy": total,
            "current_energy": float(np.vdot(V, V).real),
            "no_op": total <= 1e-28, "node_rows": np.flatnonzero(nodes).tolist(),
            "cie_bias": bias, "bias_alignment": float(-np.asarray(gradient) @ bias) if gradient is not None else None}


def quotient_normal_signature(G, m, Q, retained, candidate, rtol=1e-11, cost=None):
    """Invariant normal trace after quotienting in ||L V||_F current metric."""
    def apply(V):
        return V - np.asarray(m)[:, None, None] * np.einsum("ij,jpd->ipd", G, V, optimize=True)
    AU, AC = mode_matrix(apply(retained)), mode_matrix(apply(candidate))
    if cost:
        cost.add("green_vector_actions", retained.shape[1] * retained.shape[2] + candidate.shape[1] * candidate.shape[2])
        cost.add("rank_revealing_qr", 2)
    QU, _ = thin_basis(AU, rtol)
    Z = AC - QU @ (QU.T @ AC)
    # A pure retained-space candidate leaves only cancellation roundoff. Rank
    # must be measured against the ORIGINAL candidate, not against that noise.
    W, _ = thin_basis(Z, rtol, reference_scale=float(np.linalg.norm(AC)))
    if W.shape[1] == 0:
        return {"normal_trace": None, "real_dof": 0, "duplicate": True, "node_rows": []}
    normals = []
    nodes = np.zeros(Q.shape[0], bool)
    for k in range(W.shape[1]):
        Nk, nodes = normal_projection(Q, real_unpack(W[:, k], Q.shape))
        normals.append(float(np.vdot(Nk, Nk).real))
    return {"normal_trace": float(np.mean(normals)), "real_dof": W.shape[1],
            "duplicate": False, "node_rows": np.flatnonzero(nodes).tolist()}
