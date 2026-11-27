"""Certified finite-area trajectories for the frozen loaded 400 MHz scalar TM model."""
from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True

import numpy as np
from scipy.linalg import lu_factor, lu_solve, qr, solve, svdvals
from scipy.optimize import linprog
from scipy.ndimage import label
from scipy import sparse
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
FROZEN = HERE.parents[1] / "z_alg_optical_inertia" / "blind_uq_mimo_theory"
sys.path.insert(0, str(FROZEN / "imaging"))
sys.path.insert(0, str(FROZEN / "diagnostics"))
from shared import Model, write_arrays, write_json
from verified_wave import VerifiedWave, fro, opnorm, rd, up

CASES = {
    "n32_c01": ("results/stable_profiles/n32_c01/roi_01.npz", 256),
    "n96_c01": ("results/gn_profiles/n96_c01/roi_02.npz", 256),
}
BOX = 4.5


def timer():
    return perf_counter(), process_time()


def elapsed(start):
    return {"wall_s": perf_counter() - start[0], "cpu_s": process_time() - start[1]}


def realify(matrix):
    return np.concatenate((matrix.real, matrix.imag), axis=0)


def spectral_upper(matrix):
    """Gershgorin on a Gram matrix plus an outward multiplication allowance."""
    if matrix.shape[0] < matrix.shape[1]:
        matrix = matrix.conj().T
    gram = matrix.conj().T @ matrix
    norm = np.sum(np.abs(gram), axis=1).max()
    return min(opnorm(matrix), up(np.sqrt(up(norm + rd(matrix.shape[0]) * fro(matrix)**2))))


def adjoint_green(wave, v):
    padded = np.zeros((v.shape[1], wave.p, wave.p), complex)
    padded[:, :wave.side, :wave.side] = v.T.reshape(-1, wave.side, wave.side)
    fourier = wave.transform(wave.transform(padded) * wave.sym.conj(), True)
    return (fourier[:, :wave.side, :wave.side].reshape(v.shape[1], -1).T
            + wave.R.conj().T @ (wave.L.conj().T @ v))


def reciprocity_defect(wave):
    """Separate the exact periodic transpose and the condensed load transpose."""
    reflect = (-np.arange(wave.p)) % wave.p
    flipped = wave.sym[np.ix_(reflect, reflect)]
    symbol_defect = up(np.max(abs(wave.sym - flipped)) * (1 + rd(4)))
    lt = wave.L.T
    # M need only be an approximate representation of R in the row space of L^T.
    center = np.linalg.lstsq(lt.T, wave.R.T, rcond=None)[0].T
    mismatch = fro(wave.R - center @ lt)
    mismatch = up(mismatch + rd(center.shape[1]) * fro(center) * fro(lt)
                  + rd(2) * (fro(center @ lt) + fro(wave.R)))
    skew = up(fro(center - center.T) + rd(2) * 2 * fro(center))
    return up(symbol_defect + wave.lnorm**2 * skew + 2 * wave.lnorm * mismatch)


def select_exterior(m, x, field, receiver, roi, count, reference=None):
    ext = np.flatnonzero(~roi)
    jac = realify((receiver[:, None, :] * field.T[None, :, :]).reshape(m.Y.size, m.n))
    h = jac[:, ext]
    margin = np.minimum(x[ext], BOX - x[ext])
    priorities = (.25 + np.sqrt(margin / BOX)) * np.linalg.norm(h, axis=0)
    nonzero = np.maximum(np.linalg.norm(h, axis=0), 1e-30)
    weighted = h / nonzero[None, :] * priorities
    _, _, pivots = qr(weighted, pivoting=True, mode="economic")
    # All pivots needed for the real data dimension. Spread the other choices
    # over a regular permutation of the available domain, with an interior-room preference.
    dimension = h.shape[0]
    picked = pivots[:dimension]
    spare = np.setdiff1d(np.arange(ext.size), picked)
    room = spare[margin[spare] > .15]
    stride = np.linspace(0, room.size - 1, min(count - dimension, room.size), dtype=int)
    if reference is not None and count > 256:
        movement = np.abs(reference[ext] - x[ext])
        candidates = spare[np.argsort(movement[spare] * (.35 + np.sqrt(margin[spare])))[::-1]]
        prior = candidates[:(count - dimension) // 2]
        room = np.setdiff1d(room, prior)
        stride = np.linspace(0, room.size - 1,
                             min(count - dimension - len(prior), room.size), dtype=int)
        picked = np.r_[picked, prior, room[stride]]
    else:
        picked = np.r_[picked, room[stride]]
    if picked.size < count:
        picked = np.r_[picked, np.setdiff1d(spare, picked)[:count - picked.size]]
    return np.sort(ext[picked]), svdvals(h)


@dataclass
class Reduction:
    m: Model
    x0: np.ndarray
    indices: np.ndarray
    initial: np.ndarray
    transfer: np.ndarray
    kernel: np.ndarray
    field: np.ndarray
    kernel_error: float
    field_error: float
    transfer_error: float
    initial_error: float
    receipt: dict

    @classmethod
    def build(cls, m, wave, x0, indices, j0, receiver):
        start = timer()
        inv, info = wave.inverse(x0)
        print("verified native anchor inverse", inv, elapsed(start), flush=True)
        jgreen = wave.apply(j0, True)
        b0 = x0[:, None] * m.B
        max_x = float(x0.max())
        jerr = up(fro(j0 - x0[:, None] * jgreen - b0)
                  + max_x * wave.full_green_err * fro(j0)
                  + rd(5) * (fro(j0) + max_x * (fro(jgreen) + fro(m.B))))
        u = receiver.conj().T
        u_defect = adjoint_green(wave, x0[:, None] * u)
        uerr = up(fro(u - u_defect - m.SH) + max_x * wave.full_green_err * fro(u)
                  + rd(5) * (fro(u) + max_x * wave.gnorm * fro(u) + fro(m.SH)))
        gain = up(spectral_upper(receiver) + inv * uerr)

        width = len(indices)
        w = np.empty((m.n, width), complex)
        gw = np.empty_like(w)
        rsq = 0.0
        gsq = 0.0
        step = 16
        for offset in range(0, width, step):
            ids = indices[offset:offset + step]
            col = slice(offset, offset + len(ids))
            rhs = np.zeros((m.n, len(ids)), complex)
            rhs[ids, np.arange(len(ids))] = 1
            if m.green.dense is None:
                wp = m.solve(x0, rhs, rtol=2e-12)
            else:
                wp = lu_solve(m.factor, rhs, check_finite=False)
            gq = wave.apply(wp, True)
            green_bound = up(wave.full_green_err * fro(wp))
            rs = up(fro(wp - x0[:, None] * gq - rhs) + max_x * green_bound
                    + rd(4) * (fro(wp) + max_x * fro(gq) + fro(rhs)))
            rsq = up(rsq + up(rs**2))
            gsq = up(gsq + up(green_bound**2))
            w[:, col] = wp
            gw[:, col] = gq
            if offset % 128 == 0 or offset + len(ids) == width:
                print("audited support solves", offset + len(ids), "/", width, elapsed(start), flush=True)

        solve_defect = up(np.sqrt(rsq))
        green_error = up(np.sqrt(gsq))
        skew = reciprocity_defect(wave)
        support_gain = up(spectral_upper(gw) + green_error
                          + wave.gnorm * inv * solve_defect + inv**2 * skew)
        support_gain = min(support_gain, up(wave.gnorm * inv))
        kernel_error = up(green_error + support_gain * solve_defect)
        field_error = up(wave.full_green_err * fro(j0) + support_gain * jerr
                         + rd(2) * (fro(m.B) + fro(jgreen)))
        transfer_error = up(gain * solve_defect + rd(m.n) * wave.snorm * fro(w))
        initial_error = up(gain * jerr + rd(m.n) * wave.snorm * fro(j0))
        initial = m.S @ j0
        report = dict(inverse=info, full_green_error_upper=wave.full_green_err,
                      anchor_original_state_defect_upper=jerr,
                      anchor_adjoint_defect_upper=uerr, anchor_gain_upper=gain,
                      native_support_state_defect_upper=solve_defect,
                      native_support_G_application_error_upper=green_error,
                      exact_transpose_defect_upper=skew, support_green_gain_upper=support_gain,
                      kernel_error_upper=kernel_error, incident_field_error_upper=field_error,
                      observation_transfer_error_upper=transfer_error,
                      initial_prediction_error_upper=initial_error,
                      anchor_residual_upper_over_r=up((fro(initial - m.Y) + initial_error
                        + rd(3) * (fro(initial) + fro(m.Y))) / m.radius),
                      support_size=width, **elapsed(start))
        return cls(m, x0, indices, initial, m.S @ w, gw[indices],
                   (m.B + jgreen)[indices], kernel_error, field_error,
                   transfer_error, initial_error, report)

    def evaluation(self, x, derivative=False, audit=False):
        d = x - self.x0[self.indices]
        operator = np.eye(len(d)) - d[:, None] * self.kernel
        factor = lu_factor(operator, check_finite=False)
        q = lu_solve(factor, d[:, None] * self.field, check_finite=False)
        pred = self.initial + self.transfer @ q
        if derivative or audit:
            adjoint = lu_solve(factor, self.transfer.T, trans=1, check_finite=False).T
        if derivative:
            field = self.field + self.kernel @ q
            jacobian = (adjoint[:, None, :] * field.T[None, :, :]).reshape(self.m.Y.size, len(d))
            return pred, realify(jacobian)
        if audit:
            inverse = lu_solve(factor, np.eye(len(d)), check_finite=False)
            return Knot(x.copy(), d, operator, q, pred, inverse, adjoint)
        return pred


@dataclass
class Knot:
    x: np.ndarray
    d: np.ndarray
    operator: np.ndarray
    current: np.ndarray
    prediction: np.ndarray
    inverse: np.ndarray
    adjoint: np.ndarray


def certify_segment(reduction, a, b):
    """Bound all response states of the polyline segment with signed Bernstein terms."""
    if min(a.x.min(), b.x.min()) < 0. or max(a.x.max(), b.x.max()) > BOX:
        return dict(strict=False, cause="material box at an endpoint")
    length = len(a.d)
    q0, q1 = a.current, b.current
    c0, c1 = a.adjoint, b.adjoint
    p0, p1 = a.inverse, b.inverse
    dr = b.d - a.d
    qa = q1 - q0
    pa = p1 - p0
    ca = c1 - c0
    shift = max(np.max(np.abs(a.d)), np.max(np.abs(b.d)))
    ds = float(np.max(abs(dr)))
    identity_norm = up(np.sqrt(length))
    roundshift = rd(4) * (2 * BOX)
    knorm = fro(reduction.kernel)
    enorm = fro(reduction.field)
    aerr = up((shift + roundshift) * reduction.kernel_error + roundshift * knorm
              + rd(3) * (identity_norm + shift * knorm))
    berr = up((shift + roundshift) * reduction.field_error
              + roundshift * enorm + rd(2) * shift * enorm)
    da = -dr[:, None] * reduction.kernel
    daerr = up((ds + 4 * roundshift) * reduction.kernel_error
               + 4 * roundshift * knorm + rd(2) * ds * knorm)
    vbound = max(spectral_upper(p0), spectral_upper(p1))
    inv_end = []
    adj_end = []
    state_end = []
    state_error = []
    for k in (a, b):
        afrob = fro(k.operator)
        inv_end.append(up(fro(np.eye(length) - k.operator @ k.inverse)
                    + aerr * spectral_upper(k.inverse)
                    + rd(length) * (identity_norm + afrob * fro(k.inverse))))
        adj_end.append(up(fro(reduction.transfer - k.adjoint @ k.operator)
                    + reduction.transfer_error + aerr * spectral_upper(k.adjoint)
                    + rd(length) * (fro(reduction.transfer) + fro(k.adjoint) * afrob)))
        state = k.d[:, None] * reduction.field - k.operator @ k.current
        state_end.append(state)
        state_error.append(up(berr + aerr * fro(k.current)
                    + rd(length) * (shift * enorm + afrob * fro(k.current))))
    inv_cross = da @ pa
    inv_cross_err = up(daerr * spectral_upper(pa)
                 + rd(length) * fro(da) * (fro(p0) + fro(p1)))
    inv_defect = up(max(inv_end) + .25 * up(fro(inv_cross) + inv_cross_err))
    result = dict(inverse_defect_upper=inv_defect,
                  inverse_endpoint_defects_upper=inv_end,
                  inverse_cross_defect_upper=up(fro(inv_cross) + inv_cross_err))
    if inv_defect >= 1:
        return dict(**result, strict=False, cause="reduced state inverse defect at least one")
    inverse = up(vbound / (1 - inv_defect))
    adj_cross = ca @ da
    adj_cross_err = up(daerr * spectral_upper(ca)
                    + rd(length) * (fro(c0) + fro(c1)) * fro(da))
    adjoint_residual = up(max(adj_end) + .25 * up(fro(adj_cross) + adj_cross_err))
    adjoint_gain = max(spectral_upper(c0), spectral_upper(c1))
    state_cross = da @ qa
    state_cross_error = up(daerr * fro(qa)
                 + rd(length) * fro(da) * (fro(q0) + fro(q1)))
    state_bound = up(max(up(fro(r) + er) for r, er in zip(state_end, state_error))
                     + .25 * up(fro(state_cross) + state_cross_error))
    state_nominal_error = up(max(state_error) + .25 * state_cross_error
                         + rd(6) * (fro(state_end[0]) + fro(state_end[1]) + fro(state_cross)))
    # D(t) times the quadratic original-state defect has a signed cubic control polygon.
    r0, r1 = state_end
    rm = (r0 + r1 + state_cross) / 2
    f0 = a.prediction - reduction.m.Y
    f1 = b.prediction - reduction.m.Y
    signed = np.array([
        f0 + c0 @ r0,
        (2 * f0 + f1 + c1 @ r0 + 2 * (c0 @ rm)) / 3,
        (f0 + 2 * f1 + c0 @ r1 + 2 * (c1 @ rm)) / 3,
        f1 + c1 @ r1,
    ])
    prediction_error = up(reduction.initial_error
                     + reduction.transfer_error * max(fro(q0), fro(q1))
                     + rd(length) * (fro(reduction.initial) + fro(reduction.m.Y)
                       + fro(reduction.transfer) * max(fro(q0), fro(q1))))
    polynomial_round = up(rd(length + 12) * (fro(f0) + fro(f1)
                     + (fro(c0) + fro(c1)) * (fro(r0) + fro(r1) + fro(rm))))
    outward = up(prediction_error + polynomial_round + adjoint_gain * state_nominal_error
               + adjoint_residual * inverse * state_bound)
    signed_max = max(fro(data) for data in signed)
    ratio = up((signed_max + outward) / reduction.m.radius)
    return dict(**result, strict=bool(ratio < 1), residual_upper_over_r=ratio,
                inverse_norm_upper=inverse, adjoint_defect_upper=adjoint_residual,
                original_state_defect_upper=state_bound,
                original_state_cross_upper=up(fro(state_cross) + state_cross_error),
                signed_bernstein_control_norms_over_r=[up(fro(z) / reduction.m.radius) for z in signed],
                unsigned_output_error_over_r=up(outward / reduction.m.radius),
                cause="accepted" if ratio < 1 else "forward data-ball upper bound at least one")


def min_margin_step(jac, target, x):
    """Minimum maximum direction-specific fraction of remaining material margin."""
    n = len(x)
    eye = sparse.eye(n, format="csr")
    upper = sparse.hstack((eye, -np.maximum(0., BOX - x)[:, None]), format="csr")
    lower = sparse.hstack((-eye, -np.maximum(0., x)[:, None]), format="csr")
    scale = np.maximum(np.linalg.norm(jac, axis=1), 1e-20)
    ans = linprog(np.r_[np.zeros(n), 1.],
                  A_ub=sparse.vstack((upper, lower), format="csr"), b_ub=np.zeros(2*n),
                  A_eq=np.column_stack((jac / scale[:, None], np.zeros(len(scale)))),
                  b_eq=target / scale, bounds=[(None, None)] * n + [(0., None)],
                  method="highs", options={"primal_feasibility_tolerance": 1e-9,
                                           "dual_feasibility_tolerance": 1e-9})
    return None if not ans.success else (ans.x[:-1], float(ans.x[-1]))


def min_l2_step(jac, target, x):
    n = len(x)
    active = np.ones(n, bool)
    for _ in range(8):
        d = np.linalg.lstsq(jac[:, active], target, rcond=1e-13)[0]
        blocked = ((x[active] <= 1e-10) & (d < 0)) | ((x[active] >= BOX - 1e-10) & (d > 0))
        if not np.any(blocked):
            break
        active[np.flatnonzero(active)[blocked]] = False
        if np.count_nonzero(active) < jac.shape[0]:
            return None
    direction = np.zeros(n)
    direction[active] = d
    margin = np.where(direction >= 0., BOX - x, x)
    ratio = np.max(np.abs(direction) / np.maximum(margin, 1e-12))
    return direction, float(ratio)


def project_sum(values, total):
    """Shift then clip onto the box with a prescribed sum."""
    if total == 0:
        return np.zeros_like(values)
    lo, hi = -values.max(), BOX - values.min()
    for _ in range(46):
        mid = (lo + hi) / 2
        if np.clip(values + mid, 0., BOX).sum() < total:
            lo = mid
        else:
            hi = mid
    answer = np.clip(values + (lo + hi) / 2, 0., BOX)
    room = np.minimum(answer, BOX - answer)
    k = int(np.argmax(room))
    answer[k] = np.clip(answer[k] + total - answer.sum(), 0, BOX)
    return answer


def propose_ball(reduction, last, target, is_patch, method, max_correct=22):
    """Use the original noise ball and regularize, retaining all real readout coordinates."""
    base = last.x.copy()
    candidate = base.copy()
    candidate[is_patch] = (project_sum(base[is_patch], target * is_patch.sum())
                           if method == "ball_shape_free" else target)
    movable = np.ones(len(base), bool) if method == "ball_shape_free" else ~is_patch
    history = []
    damping = .012
    spectrum = None
    best = None
    for iteration in range(max_correct + 1):
        data, jac = reduction.evaluation(candidate, derivative=True)
        residual = (reduction.m.Y - data).ravel()
        value = float(np.linalg.norm(residual) / reduction.m.radius)
        err = realify(residual)
        history.append(dict(iteration=iteration, residual_over_r=value,
                            same_field_defect_over_r=float(np.linalg.norm(data-reduction.initial)/reduction.m.radius),
                            damping=damping))
        if best is None or value < best[0]:
            best = (value, candidate.copy())
        if spectrum is None:
            spectrum = svdvals(jac[:, ~is_patch]).tolist()
        if value < .43 or iteration == max_correct:
            break

        # The option to leave a bound remains available if the current gradient
        # points into the box. Small directions are damped but never discarded.
        gradient = jac.T @ err
        can_move = movable & (((candidate > 1e-10) & (candidate < BOX-1e-10))
                    | ((candidate <= 1e-10) & (gradient > 0))
                    | ((candidate >= BOX-1e-10) & (gradient < 0)))
        if method == "ball_shape_free":
            can_move |= is_patch
        room = np.minimum(candidate, BOX - candidate)
        if method == "ball_min_l2":
            weights = can_move.astype(float)
        else:
            weights = can_move * (.22 + .78 * np.sqrt(room / (BOX / 2)))
        scaled = jac * weights[None, :]
        if method == "ball_shape_free":
            # The weighted update has zero sum on the patch before box projection.
            # This puts the material mean into the constraint, leaving the shape free.
            normal = weights * is_patch
            normal /= np.linalg.norm(normal)
            scaled -= (scaled @ normal)[:, None] * normal[None, :]
        gram = scaled @ scaled.T
        regularizer = damping * max(np.trace(gram) / len(err), 1e-18)
        direction = weights * (scaled.T @ solve(gram + regularizer * np.eye(len(err)),
                                                err, assume_a="pos"))
        direction *= min(1., .55 / max(np.max(abs(direction)), 1e-30))
        for fraction in (1., .5, .25, .125, .0625):
            trial = np.clip(candidate + fraction * direction, 0., BOX)
            if method == "ball_shape_free":
                trial[is_patch] = project_sum(trial[is_patch], target * is_patch.sum())
            if np.linalg.norm(reduction.evaluation(trial) - reduction.m.Y) < value * reduction.m.radius:
                candidate = trial
                damping = max(2e-5, damping * .58)
                break
        else:
            damping *= 3
            if damping > 2.:
                break
    diag = dict(first_linear_relative_defect=None, predictor_margin_fraction=None,
                exterior_real_jacobian_singular_values=spectrum, correctors=history)
    if best[0] >= .955:
        return None, dict(**diag, cause="proposed endpoint lacks residual slack")
    return reduction.evaluation(best[1], audit=True), diag


def propose(reduction, last, target, is_patch, method, max_correct=6):
    base = last.x.copy()
    outside = ~is_patch
    goal = reduction.initial
    _, jac = reduction.evaluation(base, derivative=True)
    new_patch = np.full(is_patch.sum(), target)
    if method == "shape_free":
        # Legacy comparison: the old run froze a nonuniform checker shape.
        trend = np.where(np.arange(is_patch.sum()) % 2 == 0, .4, 1.6)
        low, high = 0., BOX / trend.min()
        for _ in range(44):
            mid = (low + high) / 2
            if np.clip(mid * trend, 0., BOX).mean() < target:
                low = mid
            else:
                high = mid
        new_patch = np.clip((low+high)/2 * trend, 0., BOX)
    d_roi = new_patch - base[is_patch]
    optimize = min_l2_step if method == "min_l2" else min_margin_step
    solve = optimize(jac[:, outside], -jac[:, is_patch] @ d_roi, base[outside])
    diagnostic = {"first_linear_relative_defect": None,
                  "predictor_margin_fraction": None, "correctors": []}
    sv = svdvals(jac[:, outside])
    diagnostic["exterior_real_jacobian_singular_values"] = [float(v) for v in sv]
    if solve is None:
        return None, dict(**diagnostic, cause="linear compensation solver failed")
    d_ext, lam = solve
    diagnostic["predictor_margin_fraction"] = lam
    diagnostic["first_linear_relative_defect"] = float(np.linalg.norm(
        jac[:, outside] @ d_ext + jac[:, is_patch] @ d_roi) / reduction.m.radius)
    candidate = base.copy()
    candidate[is_patch] = new_patch
    candidate[outside] = np.clip(base[outside] + d_ext, 0., BOX)
    for iteration in range(max_correct + 1):
        data, jac = reduction.evaluation(candidate, derivative=True)
        value = float(np.linalg.norm(data - goal) / reduction.m.radius)
        ball = float(np.linalg.norm(data - reduction.m.Y) / reduction.m.radius)
        diagnostic["correctors"].append(dict(iteration=iteration, same_field_defect_over_r=value,
                                              residual_over_r=ball))
        if value < 3e-5 or iteration == max_correct:
            break
        solve = optimize(jac[:, outside], realify((goal - data).ravel()), candidate[outside])
        if solve is None:
            break
        d_ext, lam = solve
        for fraction in (1., .5, .25, .125):
            trial = candidate.copy()
            trial[outside] = np.clip(trial[outside] + fraction * d_ext, 0., BOX)
            if np.linalg.norm(reduction.evaluation(trial) - goal) < value * reduction.m.radius:
                candidate = trial
                break
        else:
            break
    diag = diagnostic["correctors"][-1]
    if diag["residual_over_r"] >= .9:
        return None, dict(**diagnostic, cause="proposed endpoint lacks residual slack")
    return reduction.evaluation(candidate, audit=True), diagnostic


def exact_chart(reduction, is_patch, *, free_shape=False):
    """One-sided Banach certificate for exactly the anchor's complex response."""
    start = timer()
    x = reduction.x0[reduction.indices]
    _, jac = reduction.evaluation(x, derivative=True)
    margin = np.minimum(x, BOX - x)
    free = ~is_patch & (margin > .05)
    h = jac[:, free]
    weights = margin[free]
    u, values, vh = np.linalg.svd(h * weights[None, :], full_matrices=False)
    if len(values) != jac.shape[0] or values[-1] <= 0:
        return dict(strict=False, stop="insufficient free exterior directions", **elapsed(start)), {}
    lift = np.zeros((len(x), jac.shape[0]))
    lift[free] = weights[:, None] * ((vh.T / values) @ u.T)
    tangent = np.zeros_like(x)
    tangent[is_patch] = 1.
    if free_shape:
        # Fix only the regional mean. Diffuse positive bounds keep the
        # optimized direction physically extended on the whole task area.
        n_free, n_region = free.sum(), is_patch.sum()
        eye = sparse.eye(n_free, format="csr")
        zeros = sparse.csr_matrix((n_free,n_region))
        upper = sparse.hstack((zeros,eye,-(BOX-x[free])[:,None]),format="csr")
        lower = sparse.hstack((zeros,-eye,-x[free,None]),format="csr")
        block = np.column_stack((jac[:,is_patch],h))
        row_scale = np.maximum(np.linalg.norm(block,axis=1),1e-20)
        data_rows = np.column_stack((block/row_scale[:,None],np.zeros(len(row_scale))))
        mean_row = np.r_[np.full(n_region,1./n_region),np.zeros(n_free+1)]
        opt = linprog(np.r_[np.zeros(n_region+n_free),1.],
            A_eq=np.vstack((data_rows,mean_row)),b_eq=np.r_[np.zeros(len(row_scale)),1.],
            A_ub=sparse.vstack((upper,lower),format="csr"),b_ub=np.zeros(2*n_free),
            bounds=[(.25,2.5)]*n_region+[(None,None)]*n_free+[(0.,None)],method="highs")
        if opt.success:
            tangent[is_patch]=opt.x[:n_region]
            tangent[free]=opt.x[n_region:n_region+n_free]
        else:
            tangent[free]=lift[free] @ -(jac @ tangent)
    else:
        rhs = -(jac @ tangent)
        proposal = min_margin_step(h, rhs, x[free])
        tangent[free] = lift[free] @ rhs if proposal is None else proposal[0]
    arrays = dict(exact_lift=lift, exact_tangent=tangent)
    dnorm = spectral_upper(reduction.transfer)
    enorm = fro(reduction.field)

    def derivative_uncertainty(v):
        vinf = up(np.max(abs(v)))
        dd = reduction.transfer * v[None, :]
        de = v[:, None] * reduction.field
        return up(reduction.transfer_error * (fro(de) + rd(2)*vinf*enorm)
                  + (spectral_upper(dd) + rd(2)*vinf*dnorm
                     + reduction.transfer_error * vinf) * reduction.field_error)

    lift_error = up(np.sqrt(up(sum(up(derivative_uncertainty(col)**2) for col in lift.T))))
    # The entries of the realified nominal DF are products of two complex
    # stored coefficients. Its formation error enters each linear proof.
    jacerr = up(rd(3) * fro(reduction.transfer) * enorm)
    etal = up(fro(np.eye(jac.shape[0]) - jac @ lift) + lift_error + jacerr * fro(lift)
               + rd(len(x)) * (np.sqrt(jac.shape[0]) + fro(jac) * fro(lift)))
    zeta = up(fro(jac @ tangent) + derivative_uncertainty(tangent)
              + (jacerr + rd(len(x))*fro(jac)) * fro(tangent))
    k = up(spectral_upper(reduction.kernel) + reduction.kernel_error)
    d0 = up(up(dnorm + reduction.transfer_error) * up(enorm + reduction.field_error))
    def quadratic_coefficient(left, right):
        left_norm = up(np.max(abs(left)))
        right_norm = up(np.max(abs(right)))
        forced = right[:, None] * reduction.field
        forced_error = up(right_norm * (reduction.field_error + rd(2) * enorm))
        feedback = reduction.kernel @ forced
        feedback_error = up(k * forced_error + reduction.kernel_error * fro(forced)
                           + rd(len(x)) * fro(reduction.kernel) * fro(forced))
        weighted = left[:, None] * feedback
        weighted_error = up(left_norm * (feedback_error + rd(2) * fro(feedback)))
        observed = reduction.transfer @ weighted
        error = up((dnorm + reduction.transfer_error) * weighted_error
                   + reduction.transfer_error * fro(weighted)
                   + rd(len(x)) * fro(reduction.transfer) * fro(weighted))
        return observed, error

    g2, g2_error = quadratic_coefficient(tangent, tangent)
    quadratic_center_bound = up(fro(g2) + g2_error)
    directional = []
    directional_errors = []
    for col in lift.T:
        one, err_one = quadratic_coefficient(tangent, col)
        two, err_two = quadratic_coefficient(col, tangent)
        directional.append(realify((one + two).ravel()))
        directional_errors.append(up(err_one + err_two + rd(2)*(fro(one) + fro(two))))
    directional_bound = up(spectral_upper(np.asarray(directional).T) + fro(directional_errors))
    rownorm = np.asarray([fro(row) for row in lift])
    lnorm = up(rownorm.max())
    vnorm = up(np.max(abs(tangent)))
    receipt = dict(eta_upper=etal, zeta_upper=zeta, kernel_norm_upper=k,
                   observation_incident_product_upper=d0, tangent_inf_upper=vnorm,
                   lift_two_to_inf_upper=lnorm, free_exterior_pixels=int(free.sum()),
                   center_quadratic_coefficient_upper=quadratic_center_bound,
                   directional_lift_jacobian_first_coefficient_upper=directional_bound,
                   weighted_exterior_jacobian_singular_values=values.tolist(),
                   target="exact F(x0) with stored nominal instrument",
                   tangent_task_shape="shape-free mean constrained tangent" if free_shape else "uniform",
                   direction_on_task_min=float(tangent[is_patch].min()),
                   direction_on_task_max=float(tangent[is_patch].max()),
                   signed_parameter_domain="positive amplitude from originally vacuum patch",
                   anchor_residual_upper_over_r=reduction.receipt["anchor_residual_upper_over_r"])
    if etal >= 1:
        return dict(**receipt, strict=False, stop="linear lift defect at least one", **elapsed(start)), arrays

    def evaluate(a):
        dlead = up(a * vnorm)
        lead_q = up(k * dlead)
        if lead_q >= 1:
            return None
        gamma = up(a * zeta + a**2 * quadratic_center_bound
                   + up(d0 * dlead * lead_q**2 / ((1-lead_q) * (1-rd(2)))))
        b = up(2 * gamma / ((1-etal) * (1-rd(2))))
        depth = up(dlead + b * lnorm)
        q = up(k * depth)
        if q >= 1:
            return None
        # Keep the first signed term of DF(tv)L. Bound its remaining
        # resolvent series and the perturbation off the line separately.
        tail = up(3 * lead_q**2 / ((1-lead_q) * (1-rd(2)))**2)
        transverse = up(2 * k * lnorm * b / ((1-q) * (1-rd(2)))**3)
        kappa = up(etal + a * directional_bound + lnorm * d0 * (tail + transverse))
        containment = up(gamma + kappa * b)
        collision = up(np.max((a*abs(tangent[free]) + b*rownorm[free])
                            / np.maximum(1e-30, margin[free] * (1-rd(3)))))
        if not (kappa < 1. and containment <= b and collision < 1.
                and up(a * np.max(abs(tangent[is_patch]))) <= BOX):
            return None
        change = a * abs(tangent[free]) + b * rownorm[free]
        area = (2 / reduction.m.side)**2
        return dict(strict=True, amplitude=a, tube_radius=b, k_dot_d_upper=q,
                    contraction_upper=kappa, center_error_upper=gamma,
                    image_radius_upper=containment, exterior_margin_fraction_upper=collision,
                    exterior_cost_upper=dict(L1_m2=up(area * fro(np.ones(len(change))) * fro(change)),
                       L2_m=up(np.sqrt(area) * fro(change)), Linf=up(np.max(change))))

    candidates = []
    for n in range(3, 43):
        a = 2.**-n
        result = evaluate(a)
        candidates.append(a)
        if result is not None:
            return dict(**receipt, **result, tested_dyadic_amplitudes=candidates,
                        patch_mean_interval=[0., a if not free_shape else
                           float(a * np.mean(tangent[is_patch]) * (1-rd(len(tangent))))],
                        **elapsed(start)), arrays
    return dict(**receipt, strict=False, stop="all dyadic widths outside contraction or material tube",
                tested_dyadic_amplitudes=candidates, **elapsed(start)), arrays


def compensation(change, cell_area):
    return dict(L1_m2=float(cell_area * np.sum(abs(change))),
                L2_m=float(np.sqrt(cell_area) * np.linalg.norm(change)),
                Linf=float(np.max(abs(change))), changed_cells=int(np.count_nonzero(change)))


def run_path(reduction, is_patch, roi, method, max_target, minimum_step, cell_area,
             *, seed=None, start_target=0., start_step=.125, attempt_budget=140):
    start = timer()
    origin = reduction.x0 if seed is None else seed
    a = reduction.evaluation(origin[reduction.indices], audit=True)
    knots = [a]
    attempts = []
    target = start_target
    step = start_step
    stop = "target reached"
    while target < max_target:
        proposed = min(target + step, max_target)
        proposer = propose_ball if method.startswith("ball_") else propose
        b, detail = proposer(reduction, a, proposed, is_patch, method)
        row = dict(t0=target, t1=proposed, attempted_step=proposed-target, **detail)
        if b is not None:
            proof = certify_segment(reduction, a, b)
            row["verification"] = proof
            row["strict"] = proof["strict"]
            if proof["strict"]:
                local = b.x[~is_patch] - a.x[~is_patch]
                row["active_switches"] = int(np.count_nonzero(((a.x[~is_patch] == 0) |
                    (a.x[~is_patch] == BOX)) != ((b.x[~is_patch] == 0) | (b.x[~is_patch] == BOX))))
                row["local_exterior_cost"] = compensation(local, cell_area)
                row["cumulative_exterior_cost"] = compensation(
                    b.x[~is_patch] - reduction.x0[reduction.indices[~is_patch]], cell_area)
                row["box_margin_min"] = float(np.minimum(b.x, BOX - b.x).min())
                row["box_margin_quantiles_exterior"] = [float(z) for z in
                    np.quantile(np.minimum(b.x[~is_patch], BOX - b.x[~is_patch]), [0, .1, .5, .9, 1])]
                knots.append(b)
                a = b
                target = proposed
                step = min(.25, step * 2)
            else:
                step /= 2
        else:
            row["strict"] = False
            step /= 2
        attempts.append(row)
        if row["strict"] or step < minimum_step:
            detail_out = row.get("verification", {})
            print(method, "accepted" if row["strict"] else "stopped", target,
                  detail_out.get("residual_upper_over_r"), "attempts", len(attempts),
                  elapsed(start), flush=True)
        if step < minimum_step:
            stop = "minimum step reached without a verified next segment"
            break
        if len(attempts) >= attempt_budget and target < max_target:
            stop = f"{attempt_budget} trial budget exhausted"
            break
    stack = np.repeat(reduction.x0[None, :], len(knots), axis=0)
    stack[:, reduction.indices] = np.stack([k.x for k in knots])
    ratios = [row["verification"]["residual_upper_over_r"] for row in attempts if row["strict"]]
    worst = max(ratios) if ratios else None
    return dict(method=method,
                certified_patch_mean_interval=[float(stack[i, reduction.indices[is_patch]].mean()) for i in (0,-1)],
                certified_registered_region_mean_interval=[float(stack[i, roi].mean()) for i in (0,-1)],
                strict_segments=len(knots)-1, peak=float(stack[-1, roi].max()),
                accepted_ball_upper_over_r=worst,
                exterior_change=compensation(stack[-1, ~roi] - reduction.x0[~roi], cell_area),
                attempts=attempts, stop=stop, **elapsed(start)), stack, knots


def run(args):
    start = timer()
    case = args.case
    old, count = CASES[case]
    with np.load(FROZEN / old, allow_pickle=False) as arr:
        original = arr["material"].ravel()
        roi = arr["roi"].astype(bool).ravel()
    if np.any(original[roi] != 0.) or np.any(original < 0.) or np.any(original > BOX):
        raise ValueError("Seed is not a material-box vacuum ROI map")
    extended_path = FROZEN / "results" / "extended" / f"{case}.npz"
    with np.load(extended_path, allow_pickle=False) as arr:
        patch = np.sort(arr["patch"].astype(int))
        ext_reference = arr["material"].ravel()
        if not np.array_equal(roi, arr["roi"].ravel()) or not np.all(roi[patch]):
            raise ValueError("Published patch and ROI have changed")
    m = Model.read(case)
    cells = np.zeros(m.n, bool)
    cells[patch] = True
    if label(cells.reshape(m.shape))[1] != 1:
        raise ValueError("Published patch is not four-connected")
    with np.load(FROZEN / "ports" / "public" / "cases" / f"{case}.npz", allow_pickle=False) as arr:
        obs = str(arr["observation_id"].item())
    with np.load(FROZEN / "ports" / "public" / "operators" / f"{m.operator_id}.npz", allow_pickle=False) as arr:
        centers = arr["centers"]
        pitch = float(np.max(centers[:, 0]) - np.min(centers[:, 0])) / (m.side - 1)
        cell_area = (2 / m.side)**2
        if abs(pitch - 2 / m.side) > 1e-12:
            raise ValueError("Frozen grid pitch no longer matches the 2 m domain")
        wave = VerifiedWave(arr)
    anchor = m.forward(original, fresh=True)
    receiver = m.solve(original, m.SH, adjoint=True, rtol=2e-12).conj().T
    count = args.exterior_count or count
    exterior, spectrum = select_exterior(m, original, anchor["field"], receiver,
                                        roi, min(count, int((~roi).sum())), ext_reference)
    task_cells = np.flatnonzero(roi) if args.whole_region else patch
    indices = np.r_[task_cells, exterior]
    is_patch = np.isin(indices, task_cells)
    reused_reduction = args.cache.exists() and not args.recompute
    if reused_reduction:
        with np.load(args.cache, allow_pickle=False) as cache:
            if not np.array_equal(indices, cache["indices"]) or not np.array_equal(original, cache["x0"]):
                raise ValueError("Reduction cache references a different frozen support")
            rec = json.loads(str(cache["receipt"].item()))
            reduced = Reduction(m, original, indices, cache["initial"], cache["transfer"],
                                cache["kernel"], cache["field"], rec["kernel_error_upper"],
                                rec["incident_field_error_upper"], rec["observation_transfer_error_upper"],
                                rec["initial_prediction_error_upper"], rec)
    else:
        reduced = Reduction.build(m, wave, original, indices, anchor["j"], receiver)
        write_arrays(args.cache, indices=indices, x0=original, initial=reduced.initial,
                    transfer=reduced.transfer, kernel=reduced.kernel, field=reduced.field,
                    receipt=json.dumps(reduced.receipt))
    print("reduction receipt", json.dumps(reduced.receipt), flush=True)
    ext_replay = m.forward(ext_reference, fresh=True)
    report = dict(case=case, observation_id=obs, operator_id=m.operator_id, grid=m.side,
                  frequency_MHz=400, reads_complex=64, pulses=8,
                  power_available_per_pulse_W=.125, total_available_power_over_all_source_columns_W=1.,
                  acquisition_pulse_duration_s=None, radius=m.radius, box=[0, BOX],
                  anchor_file=old, patch_source=str(extended_path.relative_to(FROZEN)),
                  roi_cells=int(roi.sum()), patch_cells=len(patch), patch_area_m2=cell_area*len(patch),
                  roi_area_m2=cell_area*int(roi.sum()), cell_area_m2=cell_area,
                  task_cells=len(task_cells), task_area_m2=cell_area*len(task_cells),
                  task_mode="registered whole regional mean" if args.whole_region else "registered connected patch",
                  ext_support_cells=len(exterior), reference_data_exterior_jacobian_singular_values=spectrum.tolist(),
                  reused_reduction_cache=reused_reduction,
                  strict_reduction=reduced.receipt,
                  old_extended_map_numerical_comparison=dict(residual_over_r=ext_replay["residual_ratio"],
                    region_mean=float(ext_reference[roi].mean()), patch_mean=float(ext_reference[patch].mean()),
                    exterior_change=compensation(ext_reference[~roi] - original[~roi], cell_area)),
                  fixed_operator_and_radius="exact real arithmetic on the stored Fourier symbol, readout, incident field, fixed loading, noisy A observation and deterministic radius",
                  paths={})
    arrays = dict(material_anchor=original, material_old_extended=ext_reference, roi_mask=roi,
                  patch=patch, task_cells=task_cells, exterior_support=exterior, Y_observed=m.Y,
                  Y_anchor_nominal=reduced.initial, radius=m.radius)
    if not args.whole_region:
        report["exact_constant_field_chart"], chart = exact_chart(reduced, is_patch)
        arrays.update(chart)
        print("strict exact chart", report["exact_constant_field_chart"].get("amplitude"),
              report["exact_constant_field_chart"].get("stop"), flush=True)
    if args.exact_only:
        write_json(HERE / f"{case}_exact_chart{args.run_label}.json", report)
        write_arrays(HERE / f"{case}_exact_chart{args.run_label}.npz", **arrays)
        return
    for method in args.methods:
        path, mats, knots = run_path(reduced, is_patch, roi, method, args.target, args.min_step, cell_area)
        report["paths"][method] = path
        arrays[method + "_materials"] = mats
        arrays[method + "_predictions"] = np.stack([k.prediction for k in knots])
        arrays[method + "_support_currents"] = np.stack([k.current for k in knots])
        replay = m.forward(mats[-1], fresh=True)
        path["terminal_independent_original_replay"] = dict(residual_over_r=replay["residual_ratio"],
                    reduced_vs_original_over_r=float(np.linalg.norm(replay["pred"] - knots[-1].prediction) / m.radius),
                    state_defect=replay["state_relative_defect"])
        arrays[method + "_terminal_original_current"] = replay["j"]
        arrays[method + "_terminal_original_prediction"] = replay["pred"]
        report["total"] = dict(**elapsed(start), original_native_counts=dict(m.cost))
        report["with_charged_reduction"] = {k: report["total"][k] +
                                  (reduced.receipt[k] if reused_reduction else 0.) for k in ("wall_s", "cpu_s")}
        write_json(HERE / f"{case}_paths{args.run_label}.json", report)
        write_arrays(HERE / f"{case}_paths{args.run_label}.npz", **arrays)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--target", type=float, default=2.)
    parser.add_argument("--min-step", type=float, default=2.**-10)
    parser.add_argument("--methods", nargs="+", choices=["min_margin", "min_l2", "shape_free", "ball_min_margin",
                        "ball_min_l2", "ball_shape_free"], default=["ball_min_margin", "ball_min_l2", "ball_shape_free"])
    parser.add_argument("--run-label", default="")
    parser.add_argument("--exterior-count", type=int)
    parser.add_argument("--whole-region", action="store_true")
    parser.add_argument("--exact-only", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    suffix = f"_ext{args.exterior_count}" if args.exterior_count else ""
    suffix += "_whole_roi" if args.whole_region else ""
    args.cache = HERE / f"{args.case}_reduction{suffix}.npz"
    with threadpool_limits(1):
        run(args)
