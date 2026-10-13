"""Unmodified function extracts from local A13/REVISION_R2/src/mie_validation.py.
Provenance retained in SOURCE_MANIFEST.json; no upstream file modified.
Used for coefficient/angle/convention audits, not running its DDA campaign.
No SHA checks. Redistribution license not assumed.
"""

import numpy as np

import math

from scipy.special import spherical_jn,spherical_yn

BACKWARD_DEGENERACY_SIN=1e-10

def _riccati(n: int, z: complex):
    """Riccati-Bessel psi_n, psi'_n, xi_n, xi'_n in the e^{-i omega t} convention."""
    jn = spherical_jn(n, z)
    dj = spherical_jn(n, z, derivative=True)
    yn = spherical_yn(n, z)
    dy = spherical_yn(n, z, derivative=True)
    h = jn + 1j * yn
    dh = dj + 1j * dy
    return z * jn, jn + z * dj, z * h, h + z * dh

def mie_coefficients(nmax: int, m: complex, x: float):
    """Mie a_n, b_n coefficients (Bohren-Huffman, e^{-i omega t}).

    This is a local implementation used to check the existing
    ``sphere_reversal.mie_coefficients`` without importing it.
    """
    a = np.zeros(int(nmax), dtype=complex)
    b = np.zeros(int(nmax), dtype=complex)
    mx = m * x
    for n in range(1, int(nmax) + 1):
        psi_x, dpsi_x, xi_x, dxi_x = _riccati(n, x)
        psi_m, dpsi_m, _, _ = _riccati(n, mx)
        a[n - 1] = (
            (m * psi_m * dpsi_x - psi_x * dpsi_m)
            / (m * psi_m * dxi_x - xi_x * dpsi_m)
        )
        b[n - 1] = (
            (psi_m * dpsi_x - m * psi_x * dpsi_m)
            / (psi_m * dxi_x - m * xi_x * dpsi_m)
        )
    return a, b

def mie_s1_s2(mu: float, a: np.ndarray, b: np.ndarray):
    """Vector-Mie angular functions S1(mu), S2(mu) via pi/tau recurrences."""
    nmax = int(len(a))
    pi_prev = 0.0
    pi_cur = 1.0
    s1 = 0.0 + 0.0j
    s2 = 0.0 + 0.0j
    for n in range(1, nmax + 1):
        tau = n * mu * pi_cur - (n + 1) * pi_prev
        factor = (2.0 * n + 1.0) / (n * (n + 1.0))
        s1 += factor * (a[n - 1] * pi_cur + b[n - 1] * tau)
        s2 += factor * (a[n - 1] * tau + b[n - 1] * pi_cur)
        pi_next = ((2.0 * n + 1.0) * mu * pi_cur - (n + 1.0) * pi_prev) / n
        pi_prev, pi_cur = pi_cur, pi_next
    return s1, s2

def _unit_vector(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(3)
    norm = float(np.linalg.norm(v))
    if norm == 0.0:
        raise ValueError("zero vector has no direction")
    return v / norm

def mie_farfield_vector(
    u_inc,
    e_inc,
    n_obs,
    a: np.ndarray,
    b: np.ndarray,
    k: float,
    radius: float,
) -> np.ndarray:
    """Physical vector far field from a unit-incidence Mie sphere.

    Parameters
    ----------
    u_inc : (3,) real
        Unit propagation direction of the incident plane wave
        ``exp(+i k u_inc . x)``.
    e_inc : (3,) complex
        Incident electric polarization (linearly scaled; it is projected onto
        the plane transverse to ``u_inc``).
    n_obs : (3,) real
        Unit observation direction.
    a, b : (nmax,) complex
        Mie coefficients for the background wavenumber ``k``.
    k : float
        Background wavenumber.
    radius : float
        Receiver radius at which the physical far field is returned.

    Returns
    -------
    (3,) complex
        ``E_s`` at ``radius * n_obs`` including the ``exp(+i k radius)/radius``
        factor and the standard ``1/(-i k)`` amplitude convention.
    """
    u = _unit_vector(u_inc)
    n = _unit_vector(n_obs)
    e = np.asarray(e_inc, dtype=complex).reshape(3)
    e = e - np.dot(e, u) * u  # algebraic transverse projection

    c = float(np.dot(u, n))
    c = max(-1.0, min(1.0, c))
    s = math.sqrt(max(0.0, 1.0 - c * c))
    s1, s2 = mie_s1_s2(c, a, b)

    # 1/(-i k) = i/k.
    pref = 1j * np.exp(1j * k * radius) / (k * radius)

    if s <= BACKWARD_DEGENERACY_SIN:
        # Forward is non-degenerate after summing S1=S2=S0.  Exact backward
        # is a caustic of the scattering-plane basis; use the deterministic
        # limit whose transverse direction is the incident polarization.
        nmax = int(len(a))
        nvec = np.arange(1, nmax + 1, dtype=float)
        s0 = 0.5 * np.sum((2.0 * nvec + 1.0) * (a + b))
        if c >= 0.0:
            return pref * s0 * e
        return -pref * s2 * e

    # Scattering-plane basis.  e_perp is the same vector on the input and
    # output sides; e_in_par is transverse to u_inc, e_out_par to n_obs.
    e_perp = np.cross(u, n) / s
    e_in_par = (n - c * u) / s
    e_out_par = (c * n - u) / s

    comp_perp = np.dot(e, e_perp)  # algebraic dot product
    comp_par = np.dot(e, e_in_par)
    return pref * (s1 * comp_perp * e_perp + s2 * comp_par * e_out_par)

def mie_efficiencies(a: np.ndarray, b: np.ndarray, x: float):
    """Extinction/scattering/absorption efficiencies for the Mie coefficients."""
    n = np.arange(1, len(a) + 1, dtype=float)
    q_ext = float((2.0 / x**2) * np.sum((2.0 * n + 1.0) * np.real(a + b)))
    q_sca = float(
        (2.0 / x**2)
        * np.sum((2.0 * n + 1.0) * (np.abs(a) ** 2 + np.abs(b) ** 2))
    )
    return q_ext, q_sca, q_ext - q_sca
