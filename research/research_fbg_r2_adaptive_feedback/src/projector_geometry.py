"""Small, independent tangent/projector audit primitives.

These routines operate on explicit complex matrices. They are intended for
small scientific audits, not for the CUDA objective callback. All derivatives
are analytic and valid for locally constant tangent rank. They describe local
source/tangent geometry; no basin or global identifiability follows from them.
"""
import numpy as np


def current_tangent(G, E, m, direction):
    """Exact scalar-TM current tangent j_h=L^-1 diag(h)Q.

    Returns (j_h,Q,j,L), retaining the field and state needed to audit the
    feedback identity. ``direction`` is real and may be a single n-vector.
    """
    G, E, m, direction = map(np.asarray, (G, E, m, direction))
    L = np.eye(len(m), dtype=G.dtype) - m[:, None] * G
    j = np.linalg.solve(L, m[:, None] * E)
    Q = E + G @ j
    return np.linalg.solve(L, direction[:, None] * Q), Q, j, L


def transformed_field_tangent(G, B, E, m, direction):
    """Return V_h=(G+B)j_h and the expanded expression independently.

    Expanded: B diag(h)Q + (I+B diag(m))G j_h. Equality with the compact
    expression follows from L j_h=diag(h)Q. ``B`` is an explicit n×n matrix.
    """
    jh, Q, j, L = current_tangent(G, E, m, direction)
    compact = (G + B) @ jh
    expanded = B @ (direction[:, None] * Q) + (np.eye(len(m), dtype=G.dtype) + B * m[None, :]) @ (G @ jh)
    return compact, expanded


def orthonormal_tangent(T, rank_tol=None):
    """Thin orthonormal basis for the complex column span of T via SVD.

    A numerical rank threshold is explicit. A projector differential requires
    all columns to be independent; use ``projector_jvp`` for that case only.
    """
    U, singular, _ = np.linalg.svd(np.asarray(T), full_matrices=False)
    if rank_tol is None:
        rank_tol = max(T.shape) * np.finfo(singular.dtype).eps * (singular[0] if len(singular) else 0.)
    rank = int(np.count_nonzero(singular > rank_tol))
    return U[:, :rank], singular, rank


def tangent_projector(T, rank_tol=None):
    """Hermitian orthogonal projector onto span(T), with SVD rank metadata."""
    U, singular, rank = orthonormal_tangent(T, rank_tol)
    return U @ U.conj().T, {"rank": rank, "singular_values": singular}


def projector_jvp(T, dT):
    """Exact dP for full-column-rank T and a tangent perturbation dT.

    With P=T(T*T)^-1 T*, T+=(T*T)^-1 T*,
    dP=(I-P)dT T+ + [that term]*. This is a low-rank product, requiring only
    a k×k solve (k=T.shape[1]); it avoids an n×n inverse or finite differences.
    Raises on rank loss because the projector can be nondifferentiable there.
    """
    T, dT = np.asarray(T), np.asarray(dT)
    if T.shape != dT.shape or T.ndim != 2:
        raise ValueError("T and dT must have the same two-dimensional shape")
    if T.shape[1] > T.shape[0] or np.linalg.matrix_rank(T) != T.shape[1]:
        raise np.linalg.LinAlgError("projector JVP requires full column rank")
    gram = T.conj().T @ T
    tp = np.linalg.solve(gram, T.conj().T)
    P = T @ tp
    A = dT @ tp - P @ (dT @ tp)
    return A + A.conj().T


def projector_jvp_factors(T, dT):
    """Return U,V with dP=U V*+V U*, without forming an n×n matrix.

    U=(I-P)dT and V=T(T*T)^-1. This is convenient for low-rank trace and
    Galerkin audits when the observation dimension is larger than rank.
    """
    T, dT = np.asarray(T), np.asarray(dT)
    if T.shape != dT.shape or T.ndim != 2:
        raise ValueError("T and dT must have the same two-dimensional shape")
    if np.linalg.matrix_rank(T) != T.shape[1]:
        raise np.linalg.LinAlgError("projector JVP requires full column rank")
    V = np.linalg.solve(T.conj().T @ T, T.conj().T).conj().T
    U = dT - T @ (V.conj().T @ dT)
    return U, V


def galerkin_exchange(T_left, T_right, operator=None, rank_tol=None):
    """Compare two tangent subspaces and their exact Galerkin restrictions.

    Returns principal-angle cosines, tr(P_left P_right), and squared chordal
    distance ||P_left-P_right||_F². If ``operator`` is given, also returns the
    self/cross reduced matrices U_left* A U_left, U_right* A U_right,
    U_left* A U_right. This is an exchange of local tangent trial spaces,
    not a claim that the underlying nonlinear response is unchanged.
    """
    U0, s0, r0 = orthonormal_tangent(T_left, rank_tol)
    U1, s1, r1 = orthonormal_tangent(T_right, rank_tol)
    overlap = U0.conj().T @ U1
    cosines = np.linalg.svd(overlap, compute_uv=False)
    trace_overlap = float(np.linalg.norm(overlap, "fro") ** 2)
    result = {"left_rank": r0, "right_rank": r1,
              "left_singular_values": s0, "right_singular_values": s1,
              "principal_cosines": cosines, "projector_trace_overlap": trace_overlap,
              "projector_chordal_distance_squared": float(r0 + r1 - 2 * trace_overlap),
              "basis_overlap": overlap}
    if operator is not None:
        A = np.asarray(operator)
        result.update({"left_galerkin": U0.conj().T @ A @ U0,
                       "right_galerkin": U1.conj().T @ A @ U1,
                       "cross_galerkin": U0.conj().T @ A @ U1})
    return result
