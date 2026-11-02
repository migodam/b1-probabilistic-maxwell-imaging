"""Unmodified selected A9 kernels; source Gaussian/A9/code/a9_engine.py.
Copied additively for standalone study; no new Green backend.
"""
import numpy as np

def dyadic_green(dv, k):
    """Free-space dyadic Green with the A8 `_receiver` normalisation.

    Parameters
    ----------
    dv : (..., 3) float array
        Receiver minus source displacement(s); `dv = 0` is not allowed (the
        self term belongs to the polarizability).
    k : float
        Background wavenumber.

    Returns
    -------
    (..., 3, 3) complex128 with
        G_cd = exp(i k d)/(4 pi d) [ k^2 (d_cd - n_c n_d)
                                     + (i k/d - 1/d^2)(d_cd - 3 n_c n_d) ].
    """
    dv = np.asarray(dv, dtype=float)
    squeeze = dv.ndim == 1
    if squeeze:
        dv = dv[None, :]
    k = float(k)
    d = np.linalg.norm(dv, axis=-1)
    if np.any(d <= 0.0):
        raise ValueError("dyadic_green: zero separation (self term excluded)")
    inv = 1.0 / d
    n = dv * inv[..., None]
    nn = n[..., :, None] * n[..., None, :]
    eye = np.eye(3)
    pref = np.exp(1j * k * d) * inv / (4.0 * np.pi)
    out = pref[..., None, None] * (
        k * k * (eye - nn)
        + (1j * k * inv - inv * inv)[..., None, None] * (eye - 3.0 * nn)
    )
    return out[0] if squeeze else out

def build_goff(points, k, block=None, counters=None):
    """Dense off-diagonal free-space Green operator G_off (3N x 3N).

    Layout: index 3*i + c is cell i, cartesian component c (cell-major,
    component-minor) -- the same layout as A8 `_receiver` and as `_inc`.
    The diagonal 3x3 blocks are exactly zero.
    """
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    N = len(pts)
    n = 3 * N
    if block is None:
        # keep the transient (block, N, 3, 3) complex128 buffer below ~50 MB
        block = max(1, int(5.0e7 / (16.0 * 9.0 * max(N, 1))))
    G = np.empty((n, n), dtype=complex)
    for i0 in range(0, N, block):
        i1 = min(N, i0 + block)
        dv = pts[i0:i1, None, :] - pts[None, :, :]
        d = np.linalg.norm(dv, axis=-1)
        safe = d > 0.0
        inv = np.where(safe, 1.0 / np.where(safe, d, 1.0), 0.0)
        nhat = dv * inv[..., None]
        nn = nhat[..., :, None] * nhat[..., None, :]
        eye = np.eye(3)
        pref = np.where(safe, np.exp(1j * k * d) * inv / (4.0 * np.pi), 0.0)
        gb = pref[..., None, None] * (
            k * k * (eye - nn)
            + (1j * k * inv - inv * inv)[..., None, None] * (eye - 3.0 * nn)
        )
        gb = np.where(safe[..., None, None], gb, 0.0)
        # gb is (block, cell_i, c, cell_j, d); put (cell, component) contiguous
        # on each side so that row 3*i+c, column 3*j+d
        G[3 * i0 : 3 * i1, :] = gb.transpose(0, 2, 1, 3).reshape((i1 - i0) * 3, n)
    if counters is not None:
        counters.bump("kernel_build_blocks", int(np.ceil(N / block)))
        counters.bump("kernel_entries", n * n)
    return G
