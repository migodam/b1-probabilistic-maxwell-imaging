"""Scalar TM volume-integral model for the finite dielectric plate."""

from functools import lru_cache
from typing import NamedTuple

import numpy as np
from scipy.fft import fft2, ifft2, next_fast_len
from scipy.linalg import lu_factor, lu_solve
from scipy.special import hankel1, jv


C0 = 299_792_458.0
PLATE_SIDE = 0.024
PLATE_PERMITTIVITY = 2.25
BACKGROUND = PLATE_PERMITTIVITY - 1.0
BOUNDS = (-1.25, 1.5)
FREQUENCIES_GHZ = np.array(
    [2, 5.5, 6, 9, 10, 12, 12.5, 15.68, 16, 16.32, 19.5, 20, 23, 24, 26.5, 28, 30]
)


def ring(count, radius, shift=0.0):
    angle = (np.arange(count) + shift) * (2 * np.pi / count)
    return radius * np.column_stack((np.cos(angle), np.sin(angle)))


SOURCES = ring(8, 0.045)
RECEIVERS = ring(16, 0.050, 0.5)


class Host(NamedTuple):
    Eh: np.ndarray
    Sh: np.ndarray
    Gh: np.ndarray
    xy: np.ndarray


def grid(n):
    axis = ((np.arange(n) + 0.5) / n - 0.5) * PLATE_SIDE
    x, y = np.meshgrid(axis, axis, indexing="xy")
    return np.column_stack((x.ravel(), y.ravel()))


def distance(points, centers):
    return np.linalg.norm(points[:, None, :] - centers[None, :, :], axis=-1)


class Air:
    """Equal-area disk quadrature for square-cell pulse collocation."""

    def __init__(self, n, f_ghz):
        self.n = n
        self.xy = grid(n)
        self.width = PLATE_SIDE / n
        self.k = 2 * np.pi * f_ghz * 1e9 / C0
        ka = self.k * self.width / np.sqrt(np.pi)
        self.prefactor = 0.5j * np.pi * ka * jv(1, ka)
        self.self_term = 0.5j * np.pi * ka * hankel1(1, ka) - 1
        # The on-shell amplitude of H_0^(1) decays as 1/sqrt(k*r).
        # This fixes the outgoing scalar cylindrical-wave power norm.
        self.E0 = 0.25j * np.sqrt(f_ghz / 16.0) * hankel1(
            0, self.k * distance(self.xy, SOURCES)
        )
        self.S = self.external(distance(RECEIVERS, self.xy))
        offsets = np.arange(n) * self.width
        distances = np.hypot(offsets[:, None], offsets[None, :])
        distances[0, 0] = 1.0
        self.kernel = self.external(distances)
        self.kernel[0, 0] = self.self_term

    def external(self, r):
        return self.prefactor * hankel1(0, self.k * r)

    def dense(self):
        row, col = np.indices((self.n, self.n)).reshape(2, -1)
        dy = np.abs(row[:, None] - row[None, :])
        dx = np.abs(col[:, None] - col[None, :])
        return self.kernel[dy, dx]

    def convolution(self):
        """FFT apply for the same finite-grid Green matrix as dense()."""
        n = self.n
        padded = next_fast_len(2 * n - 1)
        signed = np.arange(1 - n, n)
        embedded = np.zeros((padded, padded), dtype=complex)
        embedded[np.ix_(signed, signed)] = self.kernel[np.ix_(abs(signed), abs(signed))]
        spectrum = fft2(embedded)

        def apply(vector):
            cells = np.asarray(vector).reshape(n, n)
            transformed = fft2(cells, s=(padded, padded))
            return ifft2(spectrum * transformed)[:n, :n].ravel()

        return apply


def _host(n, f_ghz, background):
    air = Air(n, f_ghz)
    green = air.dense()
    plate = np.eye(n * n, dtype=complex) - background * green
    factor = lu_factor(plate, check_finite=False)
    eh = lu_solve(factor, air.E0, check_finite=False)
    sh = lu_solve(factor, air.S.T, trans=1, check_finite=False).T
    gh = lu_solve(factor, green, check_finite=False)
    return Host(eh, sh, gh, air.xy)


@lru_cache(maxsize=32)
def host(n, f_ghz):
    """Eliminate the full known plate and retain all n*n unknown deviations.

    For delta = epsilon_r - 2.25, the perturbation current p obeys
    p = diag(delta) (Eh + Gh p), with differential receiver data Sh p.
    The coordinate array is indexed as (x, y) in meters, with x fastest.
    """
    return _host(n, f_ghz, BACKGROUND)


def debye_epsilon(f_ghz):
    """Relative permittivity for the passive Debye host with exp(-i*omega*t)."""
    return 2.10 + 0.35 / (1 - 1j * 2 * np.pi * np.asarray(f_ghz) * 1e-2)


def dielectric_factor(f_ghz):
    """Multiply real reference material q by this factor before a Debye solve."""
    return (debye_epsilon(f_ghz) - 1) / BACKGROUND


@lru_cache(maxsize=32)
def host_debye(n, f_ghz):
    """Same geometry and interface as host() for the known lossy Debye plate."""
    return _host(n, f_ghz, debye_epsilon(f_ghz) - 1)


def differential(model, delta):
    delta = np.asarray(delta).reshape(-1, 1)
    matrix = np.eye(delta.size, dtype=complex) - delta * model.Gh
    current = lu_solve(
        lu_factor(matrix, check_finite=False), delta * model.Eh, check_finite=False
    )
    return model.Sh @ current


def read_data(path=None):
    """Read the stored case/frequency/receiver/transmitter arrays."""
    from pathlib import Path

    if path is None:
        path = Path(__file__).resolve().parents[1] / "results" / "acquisition.npz"
    with np.load(path, allow_pickle=False) as source:
        return {name: source[name] for name in source.files}
