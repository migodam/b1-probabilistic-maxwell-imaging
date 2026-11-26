"""Sourcewise Dual-CIE continuation on the exact fixed dressed Green.

The alternating update and PR quadratic line search follow SRAP engine.py.
The Fourier basis and shared real readout call the unchanged R8 dependencies.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import svd
from scipy.ndimage import gaussian_filter, label
from scipy.sparse.linalg import LinearOperator, eigsh

from shared import FREEZE, NativeFourierAPCS, _readout, col_energy, elapsed, energy, now

BASE = FREEZE["baseline"]
UPPER = float(BASE["material_box"][1])


@dataclass
class CurrentState:
    x: np.ndarray
    coefficients: np.ndarray
    frequencies: np.ndarray
    stage: int
    common_beta: float | None = None


class CIE:
    def __init__(self, model):
        self.model = model
        started = now()
        u, sv, vh = svd(model.S, full_matrices=False, check_finite=False)
        rank = int(np.count_nonzero(sv > BASE["dpcs_relative_singular_cutoff"] * sv[0]))
        self.singular_values = sv
        self.Vd = vh[:rank].conj().T.copy()
        self.JD = self.Vd @ ((u[:, :rank].conj().T @ model.Y) / sv[:rank, None])
        self.D0 = model.S @ self.JD - model.Y
        self.Q0 = model.B + model.green.apply(self.JD)
        self.wd = 1 / np.maximum(col_energy(model.Y), 1e-300)
        self.ws = 1 / np.maximum(col_energy(model.B), 1e-300)
        self.setup = dict(receiver_rank=rank, singular_values=sv.tolist(), **elapsed(started))
        self.stage_basis = {}

    def stage(self, k):
        if k not in self.stage_basis:
            self.stage_basis[k] = NativeFourierAPCS(self.model.shape, BASE["fourier_corner_width"][k],
                                                   self.Vd, self.model.cost)
        return self.stage_basis[k]

    def auxiliary_current(self, state):
        a = self.stage(state.stage).apply(state.coefficients)
        j = self.JD + a
        q = self.Q0 + self.model.green.apply(a)
        return j, q

    def advance(self, initial, stage, cap, *, box_upper=None, stop=False,
                penalty=None, audit_stride=None, coefficient_rms_tol=None):
        model = self.model
        counts = model.cost
        before = counts.copy()
        started = now()
        upper = np.full(model.n, UPPER) if box_upper is None else box_upper
        lower = np.zeros_like(upper)
        x = np.clip(initial.x, lower, upper)
        U = self.stage(stage)
        if initial.frequencies.size:
            a = U.transport(initial.coefficients, initial.frequencies)
        else:
            a = np.zeros((U.k, model.ns), complex)
        if initial.common_beta is None:
            b_j, b_m = BASE["beta_current"][stage], BASE["beta_material"][stage]
        else:
            b_j = b_m = initial.common_beta
        cur = U.apply(a)
        j, q = self.JD + cur, self.Q0 + model.green.apply(cur)
        d = self.D0 + model.S @ cur
        previous_gradient = previous_direction = None
        trace, audits = [], []
        reason = "MAX_UPDATES" if stop else "MATCHED_UPDATES"
        tolerance = (BASE["relative_coefficient_rms_stop_strictly_below"]
                     if coefficient_rms_tol is None else coefficient_rms_tol)
        rms = math.inf
        cfg = dict(material_inner_maxiter=BASE["material_inner_maxiter"])
        def reg(value):
            return model.reg.evaluate(value)[:2]

        for step in range(1, cap + 1):
            p = (b_j / (1 + b_j * x))[:, None]
            r = p * (j - x[:, None] * q)
            v = p * r * self.ws
            back = model.SH @ (d * self.wd) + v - model.green.adjoint(x[:, None] * v)
            if penalty is not None:
                back += penalty.adjoint(j) * self.ws
            gradient = U.adjoint(back)
            direction = -gradient
            pr = np.zeros(model.ns)
            restart_count = 0
            if previous_gradient is not None:
                numerator = np.sum((gradient.conj() * (gradient - previous_gradient)).real, axis=0)
                pr = np.maximum(0., numerator / np.maximum(col_energy(previous_gradient), 1e-300))
                direction += previous_direction * pr
                restart = np.sum((gradient.conj() * direction).real, axis=0) >= 0.
                direction[:, restart] = -gradient[:, restart]
                restart_count = int(restart.sum())
            change_j = U.apply(direction)
            change_q = model.green.apply(change_j)
            change_d = model.S @ change_j
            change_r = p * (change_j - x[:, None] * change_q)
            denominator = self.wd * col_energy(change_d) + self.ws * col_energy(change_r)
            if penalty is not None:
                denominator += penalty.quadratic(change_j) * self.ws
            alpha = -np.sum((gradient.conj() * direction).real, axis=0) / np.maximum(denominator, 1e-300)
            da = direction * alpha
            rms = float(np.sqrt(np.mean(col_energy(da) / np.maximum(col_energy(a), 1e-300))))
            a += da
            j += change_j * alpha
            q += change_q * alpha
            d += change_d * alpha
            x, inner = _readout(j, q, x, b_m, self.ws, 1., (lower, upper), reg, True, cfg, counts)
            previous_gradient, previous_direction = gradient, direction
            after = (b_m / (1 + b_m * x))[:, None] * (j - x[:, None] * q)
            record = dict(step=step, coefficient_rms=rms, data_term=.5 * float(self.wd @ col_energy(d)),
                          state_term=.5 * float(self.ws @ col_energy(after)),
                          current_penalty=0. if penalty is None else
                          .5 * float(self.ws @ penalty.quadratic(j)),
                          pr=pr.tolist(), step_lengths=alpha.tolist(), source_restarts=restart_count,
                          **inner, **elapsed(started))
            trace.append(record)
            if audit_stride and step % audit_stride == 0:
                audit = model.audit(x, step=step)
                audits.append((x.copy(), dict(audit, **elapsed(started))))
            if stop and step >= BASE["minimum_updates_per_stage"] and rms < tolerance:
                reason = "COEFFICIENT_RMS"
                break
        out = CurrentState(x, a, U.frequencies.copy(), stage, initial.common_beta)
        if not audits or audits[-1][1]["step"] != step:
            audit = model.audit(x, step=step)
            audits.append((x.copy(), dict(audit, **elapsed(started))))
        run = dict(stage=stage, fourier_columns=U.k, beta_current=b_j, beta_material=b_m,
                   stop=reason, last_coefficient_rms=rms, updates=step,
                   coefficient_rms_tolerance=tolerance,
                   audits=[info for _, info in audits], trace=trace,
                   counts={k: counts[k] - before[k] for k in counts if counts[k] != before[k]},
                   **elapsed(started))
        return out, run, audits


def run_baseline(model):
    started = now()
    cie = CIE(model)
    initial = CurrentState(np.full(model.n, BASE["initial_contrast"]),
                           np.empty((0, model.ns), complex), np.empty((0, 2), int), 0)
    incumbent = model.audit(initial.x, stage=-1, **elapsed(started))
    states = []
    audits = [(initial.x, incumbent)]
    reports = []
    for stage in range(len(BASE["fourier_corner_width"])):
        initial, report, evaluated = cie.advance(initial, stage, BASE["maximum_updates_per_stage"], stop=True)
        states.append(initial)
        reports.append(report)
        audited_x, examined = evaluated[-1]
        if report["stop"] == "COEFFICIENT_RMS":
            audits.append((audited_x, dict(examined, stage=stage)))
    selected, selected_report = min(audits, key=lambda item: item[1]["physical_objective"])
    converged = all(r["stop"] == "COEFFICIENT_RMS" for r in reports)
    metadata = dict(stages=reports, converged=converged, selected=selected_report,
                    cie_setup=cie.setup, counts=dict(model.cost), **elapsed(started))
    return selected, states, metadata


def segment(x):
    side = x.shape[0]
    smooth = gaussian_filter(x, sigma=side/32.)
    threshold = max(.25, .4*float(np.max(smooth)))
    labels, total = label(smooth >= threshold)
    all_parts = [labels == k for k in range(1, total+1)]
    all_parts.sort(key=lambda p: -float(np.sum(smooth[p])))
    minimum = math.ceil(4 * (side/32.)**2)
    major = [p for p in all_parts if p.sum() >= minimum][:3]
    if not major:
        peak = np.zeros_like(x, bool)
        peak.flat[np.argmax(smooth)] = True
        major = all_parts[:1] or [peak]
    masks = np.stack(major)
    y, xx = np.indices(x.shape)
    low = (smooth < .5*threshold) & ~masks.any(axis=0)
    tiles = np.stack([low & ((y < side//2) == u) & ((xx < side//2) == v)
                      for u, v in ((True, True), (True, False), (False, True), (False, False))])
    info = dict(threshold=threshold, contrast_peak=float(smooth.max()), minimum_cells=minimum,
                components=[dict(id=f"roi_{i+1:02d}", cells=int(p.sum()),
                                 mass=float(np.sum(x[p])), peak=float(np.max(x[p])))
                            for i, p in enumerate(masks)],
                background_cells=[int(z.sum()) for z in tiles])
    return masks, tiles, smooth, info


class Penalty:
    def __init__(self, weight, vectors=None):
        self.weight = weight
        self.vectors = vectors
        self.vectors_h = None if vectors is None else vectors.conj().T.copy()

    def quadratic(self, j):
        return self.weight * col_energy(j if self.vectors is None else self.vectors_h @ j)

    def adjoint(self, j):
        projected = j if self.vectors is None else self.vectors @ (self.vectors_h @ j)
        return self.weight * projected


def compute_k_modes(model, max_rank=8):
    """Uniform-sector K=S* S+I-U Herm(G) on the exact receiver kernel."""
    started = now()
    _, sv, vh = svd(model.S, full_matrices=False, check_finite=False)
    rank = int(np.count_nonzero(sv > np.finfo(float).eps * max(model.S.shape) * sv[0]))
    rows = vh[:rank].conj().T.copy()
    rows_h = rows.conj().T.copy()

    def null(z):
        return z - rows @ (rows_h @ z)

    def k_apply(z):
        p = null(z)
        herm = .5*(model.green.apply(p) + model.green.adjoint(p))
        # Lift the eliminated row directions to a positive barrier.
        return p - UPPER*null(herm) + (z-p)

    n = model.n
    rng = np.random.default_rng(FREEZE["B2"]["random_seed"])
    initial = null((rng.standard_normal(n) + 1j*rng.standard_normal(n))[:, None])[:, 0]
    operator = LinearOperator((n, n), matvec=lambda v: k_apply(v[:, None])[:, 0], dtype=complex)
    _, modes = eigsh(operator, k=max_rank, which="SA", tol=2e-8, maxiter=600, v0=initial)
    order = np.argsort(np.real(np.sum(modes.conj()*k_apply(modes), axis=0)))
    modes = null(modes[:, order])
    # eigsh returns normalized eigenvectors to its tolerance; rotate QR only
    # in each degenerate cluster by preserving the returned sorted order.
    modes /= np.linalg.norm(modes, axis=0)
    values = np.sum(modes.conj()*k_apply(modes), axis=0).real
    negative = values < -1e-7
    dark = modes[:, negative]
    rnd = null(rng.standard_normal(dark.shape) + 1j*rng.standard_normal(dark.shape))
    comparison = np.linalg.qr(rnd, mode="reduced")[0]
    true_defects = np.linalg.norm(k_apply(modes) - modes * values, axis=0)
    info = dict(restricted_values=values.tolist(), negative_kept=int(negative.sum()),
                receiver_row_rank=rank, receiver_row_singular_values=sv.tolist(),
                receiver_leak_absolute=float(np.linalg.norm(model.S @ modes)),
                k_eigenpair_defects=true_defects.tolist(),
                counts=dict(model.cost), **elapsed(started))
    return dark, comparison, info
