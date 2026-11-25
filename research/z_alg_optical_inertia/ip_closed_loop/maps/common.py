"""Frozen acquisition and the exact discrete convolution convention."""
from pathlib import Path
from time import perf_counter, process_time
import json

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.special import hankel1, jv

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
FREEZE = json.loads((HERE / "contract.json").read_text())


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return {"wall_s": perf_counter() - start[0], "cpu_s": process_time() - start[1]}


def save_json(path, entry):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp")
    temp.write_text(json.dumps(entry, indent=2, allow_nan=False, ensure_ascii=False) + "\n")
    temp.replace(path)


def pattern(side, shape):
    scale = side / 32
    offset = np.arange(-int(6 * scale), int(6 * scale) + 1)
    y, x = np.meshgrid(offset, offset, indexing="ij")
    xh, yh = x / scale, y / scale

    def d(xc, yc):
        return (xh - xc)**2 + (yh - yc)**2

    if shape == "solid":
        core = d(0, 0) <= 3.15**2
        extended = d(0, 0) <= 3.80**2
    elif shape == "ring":
        core = (d(0, 0) >= 1.90**2) & (d(0, 0) <= 3.95**2)
        extended = (d(0, 0) >= 1.60**2) & (d(0, 0) <= 4.31**2)
    elif shape == "double":
        z = np.minimum(d(-2.9, 0), d(2.9, 0))
        core, extended = z <= 2.18**2, z <= 2.55**2
    elif shape == "austria":
        ring = d(0, -0.9)
        satellites = np.minimum(d(-2.35, 3.65), d(2.35, 3.65))
        core = ((ring >= 1.47**2) & (ring <= 2.65**2)) | (satellites <= 1.34**2)
        extended = ((ring >= 1.2**2) & (ring <= 2.95**2)) | (satellites <= 1.63**2)
    else:
        raise ValueError(shape)
    assert np.all(~core | extended)
    return np.stack((x[extended], y[extended]), axis=1), np.where(core[extended], 1.0, 0.5)


def place(side, center, offsets, profile):
    xy = offsets + np.array(center)
    valid = np.all((xy >= 0) & (xy < side), axis=1)
    cells = xy[valid, 1] * side + xy[valid, 0]
    return cells, profile[valid]


def grid_positions(side):
    y, x = np.indices((side, side))
    return np.stack((x.ravel(), y.ravel()), axis=1)


def geometry32():
    settings = FREEZE["model"]
    aperture = FREEZE["aperture"]
    side = 32
    dx = settings["side_m"] / side
    xy = (grid_positions(side) + 0.5) * dx - settings["side_m"] / 2
    k = 2 * np.pi * settings["frequency_hz"] / settings["host_speed_m_per_s"]
    radius = dx / np.sqrt(np.pi)
    tx = np.deg2rad(aperture["tx_angles_degrees"])
    rx = np.deg2rad(aperture["rx_angles_degrees"])
    tx = np.stack((np.cos(tx), np.sin(tx)), axis=1)
    rx = aperture["sensor_radius_m"] * np.stack((np.cos(rx), np.sin(rx)), axis=1)
    B = aperture["transmitter_amplitude"] * np.exp(1j * k * (xy @ tx.T))
    constant = 0.5j * np.pi * k * radius * jv(1, k * radius)
    S = constant * hankel1(0, k * np.linalg.norm(rx[:, None] - xy[None, :], axis=2))

    offset = np.arange(2 * side)
    offset = np.where(offset < side, offset, offset - 2 * side)
    dy, dx_ = np.meshgrid(offset, offset, indexing="ij")
    distance = dx * np.hypot(dx_, dy)
    distance[0, 0] = 1.0
    green = constant * hankel1(0, k * distance)
    green[0, 0] = 0.5j * np.pi * k * radius * hankel1(1, k * radius) - 1
    # Lag ±side is never used; zero these rows to tighten the embedding.
    green[side, :] = 0
    green[:, side] = 0
    return {"side": side, "B": B, "S": S, "spectrum": fft2(green)}


class Geometry:
    def __init__(self, archive):
        with np.load(archive, allow_pickle=False) as arrays:
            self.side = int(arrays["side"])
            self.B = arrays["B"]
            self.S = arrays["S"]
            self.symbol = arrays["spectrum"]
        self.n = self.side**2
        self.xy = grid_positions(self.side)
        self.lag = ifft2(self.symbol)
        self.gram_lag = ifft2(np.abs(self.symbol)**2)

    def embed(self, u):
        view = u.reshape(self.side, self.side, -1)
        return fft2(view, s=self.symbol.shape, axes=(0, 1))

    def apply(self, u):
        field = ifft2(self.symbol[:, :, None] * self.embed(u), axes=(0, 1))
        return field[:self.side, :self.side].reshape(u.shape)

    def restrict(self, cells, lag=None):
        lag = self.lag if lag is None else lag
        xy = self.xy[cells]
        shift = (xy[:, None] - xy[None, :]) % lag.shape[0]
        return lag[shift[:, :, 1], shift[:, :, 0]]

    def forward(self, x):
        support = np.flatnonzero(x)
        if len(support) <= 700:
            g = self.restrict(support)
            j = np.linalg.solve(np.eye(len(g)) - x[support, None] * g,
                                x[support, None] * self.B[support])
            current = np.zeros_like(self.B)
            current[support] = j
            return self.S @ current, current, 0
        operator = LinearOperator((self.n, self.n),
                                  matvec=lambda u: u - x * self.apply(u),
                                  dtype=np.complex128)
        current = np.empty_like(self.B)
        max_info = 0
        for source, incident in enumerate(self.B.T):
            current[:, source], info = gmres(operator, x * incident, rtol=2e-12,
                                              maxiter=60, restart=60)
            max_info = max(max_info, abs(info))
        return self.S @ current, current, max_info


def geom_path(side):
    return HERE / "public" / f"geometry_{side}.npz"
