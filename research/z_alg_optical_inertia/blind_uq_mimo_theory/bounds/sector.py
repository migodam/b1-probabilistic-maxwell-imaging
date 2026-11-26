"""Native-grid Theorem 1 state/affine blocks and reduced dark-current proposals."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, null_space, qr

HERE = Path(__file__).resolve().parent
PUBLIC = HERE.parent / "ports" / "public"
UPPER = 4.5


def adj(value):
    return value.conj().T


def herm(value):
    return (value + adj(value)) / 2


@dataclass
class Green:
    spectrum: np.ndarray
    side: int
    left: np.ndarray
    right: np.ndarray

    @classmethod
    def read(cls, operator, side):
        data = np.load(operator)
        n = side * side
        left = data["low_rank_left"] if "low_rank_left" in data else np.empty((n, 0), complex)
        right = data["low_rank_right"] if "low_rank_right" in data else np.empty((0, n), complex)
        return cls(data["spectrum"], side, left, right)

    @property
    def n(self):
        return self.side * self.side

    @property
    def pad(self):
        return self.spectrum.shape[0]

    def apply(self, value, transpose=False):
        z = value.reshape(self.side, self.side, -1)
        spec = self.spectrum.conj() if transpose else self.spectrum
        result = ifft2(fft2(z, s=self.spectrum.shape, axes=(0, 1)) * spec[..., None], axes=(0, 1))
        bare = np.ascontiguousarray(result[:self.side, :self.side]).reshape(value.shape)
        if transpose:
            return bare + adj(self.right) @ (adj(self.left) @ value)
        return bare + self.left @ (self.right @ value)

    def matrix(self, rows, cols):
        kernel = ifft2(self.spectrum)
        dy = rows[:, None] // self.side - cols[None, :] // self.side
        dx = rows[:, None] % self.side - cols[None, :] % self.side
        return kernel[dy % self.pad, dx % self.pad] + self.left[rows] @ self.right[:, cols]

    def waves(self, count):
        """Approximate subspace only, ordered by high real Green symbol."""
        mode = np.argsort(-self.spectrum.real.ravel(), kind="stable")[:count]
        fy, fx = np.unravel_index(mode, self.spectrum.shape)
        y, x = np.indices((self.side, self.side))
        phase = (y.reshape(-1, 1) * fy + x.reshape(-1, 1) * fx) * (2j * np.pi / self.pad)
        return np.exp(phase) / self.side


@dataclass
class Profile:
    low: np.ndarray
    high: np.ndarray
    weight: np.ndarray
    optical: np.ndarray
    receive: np.ndarray

    @property
    def mix(self):
        return (self.weight * (self.low + self.high) + 1j * self.optical) / 2

    @property
    def prod(self):
        return self.weight * self.low * self.high


def blocks(g, s, b, y, p):
    mix, prod = p.mix[:, None], p.prod[:, None]
    receiver = p.receive[:, None]
    mg = mix * g
    k = herm(adj(s) @ (receiver * s) - mg - adj(mg))
    k.flat[::len(k)+1] += p.weight
    rhs = adj(s) @ (receiver * y) + mix * b
    h0 = adj(y) @ (receiver * y)
    if np.any(prod):
        k += adj(g) @ (prod * g)
        rhs -= adj(g) @ (prod * b)
        h0 += adj(b) @ (prod * b)
    return herm(k), rhs, herm(h0)


def tail_from_blocks(k, rhs, h0):
    eig, modes = eigh(k, driver="evd", check_finite=False)
    gap = float(np.min(np.abs(eig)))
    nu = int((eig < 0).sum())
    if nu >= h0.shape[0] or gap < 1e-12:
        return 0.0, nu, gap, []
    source = adj(modes) @ rhs
    sh = herm(h0 - adj(source) @ (source / eig[:, None]))
    sh_eigs = eigh(sh, eigvals_only=True, check_finite=False)[::-1]
    return float(np.maximum(sh_eigs[nu:], 0).sum()**0.5), nu, gap, sh_eigs.tolist()


def dark_trial(green, s, indices, rank=144):
    if not len(indices):
        return np.empty((green.n,0),complex),np.empty((green.n,0),complex)
    trial = green.waves(rank)
    inactive = np.ones(green.n, bool)
    inactive[indices] = False
    trial[inactive] = 0
    frame, r, _ = qr(trial[indices], mode="economic", pivoting=True, check_finite=False)
    keep = np.abs(r.diagonal()) > abs(r[0, 0]) * 1e-9
    frame = frame[:, keep]
    null = null_space(s[:, indices] @ frame, rcond=1e-13)
    trial = np.zeros((green.n, null.shape[1]), complex)
    trial[indices] = frame @ null
    if not trial.shape[1]:
        return trial,trial.copy()
    scattered=np.hstack([green.apply(trial[:,j:j+16]) for j in range(0,trial.shape[1],16)])
    return trial,scattered


def dark_compression(trial, scattered, s, p):
    mix = adj(trial) @ (p.mix[:, None] * scattered)
    measured = s @ trial
    return herm(adj(measured) @ (p.receive[:, None] * measured)
                + adj(trial) @ ((p.weight)[:, None] * trial)
                + adj(scattered) @ (p.prod[:, None] * scattered) - mix - adj(mix))


def surrogate_uniform(green, s, b, y, p):
    """Lower K envelope using *all* physical currents, diagonal Fourier symbol and exact loading.

    Common sector endpoint and common optical/sector weights required. If l>0, we
    discard the PSD G*diag(w l u)G from K only. L and C stay exact.
    """
    if any(np.ptp(v) for v in (p.low, p.high, p.weight, p.optical)):
        raise ValueError("Uniform envelope requires spatially uniform sector and optical terms")
    mix = p.mix[0]
    core = p.weight[0] - 2 * np.real(mix * green.spectrum)
    alpha = (p.weight[0]) / 2
    keep = np.flatnonzero(core.ravel() < alpha)
    fy, fx = np.unravel_index(keep, green.spectrum.shape)
    dy, dx = np.indices((green.side, green.side))
    phase = (dy.reshape(-1, 1) * fy + dx.reshape(-1, 1) * fx) * (2j * np.pi / green.pad)
    fourier = np.exp(phase) / green.pad
    remove = p.receive == 0
    s = s[~remove]
    y = y[~remove]
    receiver = p.receive[~remove]
    f = np.column_stack((fourier, green.left, adj(green.right), adj(s)))
    nkr = green.left.shape[1]
    d = np.diag(np.r_[core.ravel()[keep] - alpha, np.zeros(2*nkr), receiver]).astype(complex)
    sl = slice(len(keep), len(keep)+nkr)
    sr = slice(sl.stop, sl.stop+nkr)
    d[sl, sr] = -mix * np.eye(nkr)
    d[sr, sl] = -mix.conjugate() * np.eye(nkr)
    if nkr and not mix:
        valid = np.r_[np.arange(len(keep)), np.arange(sr.stop, len(d))]
        d = d[np.ix_(valid, valid)]
        f = f[:, valid]
    if not len(d):
        eig = np.array([])
        ds = eig
        solved = None
    else:
        ds = eigh(d, eigvals_only=True, check_finite=False)
        small = herm(np.linalg.inv(d) + adj(f) @ f / alpha)
        eig, modes = eigh(small, driver="evd", check_finite=False)
        solved = modes
    nu = int((eig > 0).sum() - (ds > 0).sum())
    rhs = adj(s) @ (receiver[:, None] * y) + mix * b
    h0 = adj(y) @ (receiver[:, None] * y)
    if p.prod[0]:
        rhs -= green.apply(p.prod[:, None] * b, transpose=True)
        h0 += adj(b) @ (p.prod[:, None] * b)
    if nu >= y.shape[1]:
        return 0.0, nu, [], dict(keep=int(len(keep)), proof_dimension=int(len(d)))
    projected = adj(f) @ rhs
    correction = 0.0
    if solved is not None:
        projected = adj(solved) @ projected
        correction = adj(projected) @ (projected / eig[:, None]) / alpha**2
    spectrum = eigh(herm(h0 - adj(rhs) @ rhs / alpha + correction), eigvals_only=True)[::-1]
    return float(np.maximum(spectrum[max(0, nu):], 0).sum()**0.5), nu, spectrum.tolist(), dict(
        keep=int(len(keep)), proof_dimension=int(len(d)))
