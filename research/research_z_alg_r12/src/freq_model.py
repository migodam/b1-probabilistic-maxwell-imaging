"""Fixed-geometry pulse TM map and derivatives for the signed-frequency test."""

from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.special import hankel1, h1vp, jv, jvp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "research_z_alg_r2" / "src"))
from physics import Model

BOX = (.5, 30.)
SCALE = .05


def ri(values):
    values = np.asarray(values)
    return np.stack([values.real, values.imag], axis=1).reshape(
        (2 * values.shape[0],) + values.shape[1:]
    )


class FrequencyTM:
    def __init__(self, geometry):
        self.geometry = geometry
        self.k0 = geometry["k"]
        self.factorizations = 0

    @lru_cache(maxsize=128)
    def kernels(self, delta):
        k = self.k0 * (1 + delta)
        model = Model(**dict(self.geometry, k=k))
        radius = model.radius
        d = np.linalg.norm(model.centers[:, None] - model.centers[None, :], axis=-1)
        np.fill_diagonal(d, 1.)
        dr = np.linalg.norm(model.receivers[:, None] - model.centers[None, :], axis=-1)
        ds = np.linalg.norm(model.centers[:, None] - model.sources[None, :], axis=-1)
        coeff = .5j * np.pi * k * radius * jv(1, k * radius)
        coeff_p = .5j * np.pi * radius * (jv(1, k * radius) + k * radius * jvp(1, k * radius))
        gp = coeff_p * hankel1(0, k * d) + coeff * d * h1vp(0, k * d)
        np.fill_diagonal(gp, .5j * np.pi * radius * (
            hankel1(1, k * radius) + k * radius * h1vp(1, k * radius)))
        sp = coeff_p * hankel1(0, k * dr) + coeff * dr * h1vp(0, k * dr)
        # For E=C iH0(kr)/4, free-background power per unit length is |C|²/(8 Z k).
        amplitude = np.sqrt(k / self.k0)
        b = amplitude * model.E0
        bp = amplitude * .25j * ds * h1vp(0, k * ds) + b / (2 * k)
        return model.G, b, model.S, self.k0 * gp, self.k0 * bp, self.k0 * sp

    def response(self, m, delta=0., hessian=False, derivative=False):
        self.factorizations += 1
        g, b, s, gp, bp, sp = self.kernels(float(delta))
        resolvent = np.linalg.inv(np.eye(m.size) - m[:, None] * g)
        current = resolvent @ (m[:, None] * b)
        field = b + g @ current
        transfer = s @ resolvent
        y = (s @ current).T.reshape(-1)
        j = np.einsum("ri,is->sri", transfer, field).reshape(y.size, m.size)
        if derivative:
            current_p = resolvent @ (m[:, None] * (gp @ current + bp))
            field_p = bp + gp @ current + g @ current_p
            transfer_p = (sp + transfer @ (m[:, None] * gp)) @ resolvent
            yp = (sp @ current + s @ current_p).T.reshape(-1)
            jp = (np.einsum("ri,is->sri", transfer_p, field) +
                  np.einsum("ri,is->sri", transfer, field_p)).reshape(j.shape)
            return y, j, yp, jp
        if hessian:
            shared = np.einsum("ri,ij,js->srij", transfer, g @ resolvent, field)
            h = (shared + shared.swapaxes(-1, -2)).reshape(y.size, m.size, m.size)
            return y, j, h
        return y, j


class Fit:
    def __init__(self, model, deltas, data, weights=None, gain_calibration=None):
        self.model = model
        self.deltas = deltas
        self.weights = tuple(np.ones(len(deltas))) if weights is None else weights
        self.data = np.concatenate([w * y for w, y in zip(self.weights, data, strict=True)])
        self.calibration = gain_calibration
        self.last = None

    def evaluate(self, m, second=False):
        if self.last is not None and np.array_equal(self.last[0], m) and (
            self.last[1] is not None or not second
        ):
            return self.last[2]
        outputs = [self.model.response(m, d, hessian=second) for d in self.deltas]
        y, j = [np.concatenate([w * row[i] for w, row in zip(self.weights, outputs, strict=True)])
                for i in range(2)]
        h = np.concatenate([w * row[2] for w, row in zip(self.weights, outputs, strict=True)]) if second else None
        if self.calibration is None:
            r, jr, gain = y - self.data, j, 1.
        else:
            lam, center = self.calibration
            denom = np.vdot(y, y).real + lam
            gain = (np.vdot(y, self.data) + lam * center) / denom
            r = np.r_[gain * y - self.data, np.sqrt(lam) * (gain - center)]
            numerator_p = j.conj().T @ self.data
            denominator_p = 2 * np.real(y.conj() @ j)
            gain_p = (numerator_p - gain * denominator_p) / denom
            jr = np.vstack([gain * j + y[:, None] * gain_p, np.sqrt(lam) * gain_p])
        grad = np.real(jr.conj().T @ r)
        hh = None
        if second:
            fixed_j = gain * j
            hh = np.real(fixed_j.conj().T @ fixed_j +
                         np.einsum("k,kij->ij", r[:y.size].conj(), gain * h))
            if self.calibration is not None:
                nuisance_j = np.column_stack((y, 1j * y))
                cross = np.real(fixed_j.conj().T @ nuisance_j +
                                np.column_stack((j.conj().T @ r[:y.size],
                                                 -1j * j.conj().T @ r[:y.size])))
                hh -= cross @ cross.T / denom
            hh = (hh + hh.T) / 2
        result = r, jr, grad, hh, gain
        self.last = m.copy(), h, result
        return result

    def solve(self, initial, limit=180):
        start = self.model.factorizations
        sol = least_squares(lambda m: ri(self.evaluate(m)[0]) / SCALE,
                            initial, jac=lambda m: ri(self.evaluate(m)[1]) / SCALE,
                            bounds=BOX, xtol=1e-13, ftol=1e-13, gtol=1e-14, max_nfev=limit)
        return sol.x, self.model.factorizations - start

    def describe(self, m):
        residual, _, gradient, hessian, gain = self.evaluate(m, True)
        active = (m <= BOX[0] + 1e-8) | (m >= BOX[1] - 1e-8)
        kkt = m - np.clip(m - gradient, *BOX)
        ev = np.linalg.eigvalsh(hessian[np.ix_(~active, ~active)])
        return dict(loss=float(np.vdot(residual, residual).real / 2),
                    gradient=float(np.linalg.norm(gradient, ord=np.inf)),
                    kkt=float(np.linalg.norm(kkt, ord=np.inf)),
                    active=np.flatnonzero(active).tolist(), eig=ev.tolist(),
                    gain=[float(np.real(gain)), float(np.imag(gain))])
