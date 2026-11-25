"""Entire-field signed receiver separation with a full-domain exterior."""
import argparse
from collections import defaultdict
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.linalg import cholesky, eigvalsh, solve
from threadpoolctl import threadpool_info, threadpool_limits

from common import (FREEZE, HERE, Geometry, clock, elapsed, geom_path, pattern,
                    save_json)


# Binary64 accumulation bound. Every dot length is <= 9216; eight
# real operations per complex term is sufficient. The multiplier leaves
# more than 800-fold room over gamma_(8*9216) and all scalar chains.
ARITH = 2.0**-27


def fnorm(x, axis=None):
    return np.linalg.norm(x, axis=axis)


def add_bound(x):
    return x + ARITH * np.maximum(1.0, x)


def verify_eigen_bound(a, upper, entry_error=0):
    """A proposed Cholesky factor proves a spectral bound by residual."""
    norm = fnorm(a)
    h = (a + a.conj().T) / 2
    dim = len(h)
    i = np.eye(dim)
    endpoint = eigvalsh(h, subset_by_index=[dim - 1, dim - 1] if upper else [0, 0],
                        check_finite=False)[0]
    slack = ARITH * max(1., norm)
    trial = endpoint + slack if upper else endpoint - slack
    sign = 1 if upper else -1
    positive = sign * (trial * i - h)
    try:
        L = cholesky(positive, lower=True, check_finite=False)
    except np.linalg.LinAlgError:
        return add_bound(norm + entry_error) if upper else -add_bound(norm + entry_error)
    approx = L @ L.conj().T
    arithmetic = ARITH * max(1., abs(trial) * np.sqrt(dim) + norm + fnorm(L)**2)
    error = fnorm(positive - approx) + arithmetic + entry_error
    return trial + sign * error


def op_upper(a):
    gram = a.conj().T @ a
    return np.sqrt(max(0., verify_eigen_bound(gram, True, ARITH * fnorm(a)**2)))


def local_models(g, shape):
    offsets, profile = pattern(g.side, shape)
    nodes = g.xy[:, None, :] + offsets[None, :, :]
    valid = np.all((nodes >= 0) & (nodes < g.side), axis=2)
    masks, order = np.unique(valid, axis=0, return_inverse=True)
    lut = nodes[..., 1] * g.side + nodes[..., 0]
    for mask_id, mask in enumerate(masks):
        ids = np.flatnonzero(order == mask_id)
        cells = lut[ids][:, mask]
        assert cells.shape[1] > 0
        yield ids, cells, profile[mask]


class ReceiverGeometry:
    def __init__(self, g):
        self.g = g
        self.bnorm = add_bound(fnorm(g.B))
        self.sgram = g.S @ g.S.conj().T
        self.snorm = np.sqrt(verify_eigen_bound(self.sgram, True, ARITH * fnorm(g.S)**2))
        first = np.einsum("ri,is->rsi", g.S, g.B).reshape(g.S.shape[0] * g.B.shape[1], -1)
        raw = np.r_[first.real, first.imag]
        self.bg_first_sum = raw.sum(axis=1)
        self.bg_first_gram = raw @ raw.T
        self.bg_first_scale = g.n * fnorm(raw)**2
        self.bg_total_bound = add_bound(np.sqrt(g.n) * fnorm(raw))

    def background_linear(self, w):
        flat = w.reshape(len(w), -1)
        wf = np.c_[flat.real, flat.imag]
        sq = np.einsum("bi,bi->b", wf @ self.bg_first_gram, wf)
        sq = np.maximum(0., self.g.n * sq) + ARITH * self.bg_first_scale
        return np.maximum(0., 0.5 * (wf @ self.bg_first_sum + np.sqrt(sq))) + (
            ARITH * self.bg_total_bound
        )

    def sensor_direction(self, w):
        t = np.einsum("brs,brs->b", w.conj(), self.sgram @ w).real
        return np.sqrt(np.maximum(t, 0.) + ARITH * max(1., fnorm(self.g.S)**2))


def active_coercivity(g_full, upper, entry_error):
    scales = np.sqrt(upper)
    scaled_green = scales[:, None] * g_full * scales
    a = np.eye(len(upper)) - scaled_green

    def hermitian(theta):
        rot = np.exp(1j * theta) * a
        return (rot + rot.conj().T) / 2

    angles = np.deg2rad([0, 30, 50, 65, 75, 82, 87, 89])
    estimates = [eigvalsh(hermitian(theta), subset_by_index=[0, 0], check_finite=False)[0]
                 for theta in angles]
    theta = angles[int(np.argmax(estimates))]
    error = len(upper) * np.max(upper) * entry_error + ARITH * fnorm(a)
    floor = verify_eigen_bound(hermitian(theta), False, error)
    return floor, float(theta)


def local_inverse(g, cells, p, lo, hi, full_floor, floor_theta, op):
    k = len(p)
    green = g.restrict(cells)
    gram = g.restrict(cells, g.gram_lag)
    err_g = k * op["kernel_entry_error"]
    err_gram = k * op["gram_kernel_entry_error"]
    gram12 = green.conj().T @ green
    gram21 = green @ green.conj().T
    gram_error = err_gram + err_g * (2 * fnorm(green) + err_g) + ARITH * (
        fnorm(gram) + fnorm(green)**2
    )
    coupling = np.sqrt(max(0.,
        verify_eigen_bound(gram - gram12, True, gram_error),
        verify_eigen_bound(gram - gram21, True, gram_error),
    ))
    low, high = p * lo, p * hi
    m = (low + high) / 2
    I = np.eye(k)
    state = np.diag(1.0 / m) - green
    T = solve(state, I, check_finite=False)
    tnorm = op_upper(T)
    delta = fnorm(I - T @ state) + err_g * tnorm + ARITH * max(
        1., np.sqrt(k) + fnorm(T) * fnorm(state)
    )
    inv_gap = delta * tnorm / (1.0 - delta) if delta < 1 else np.inf
    inv_norm = tnorm + inv_gap
    dmin, dmax = 1 / m - 1 / low, 1 / m - 1 / high
    max_d = np.maximum(-dmin, dmax)
    d_op = op_upper(T * max_d) + inv_gap * np.max(max_d) + ARITH
    # If |D| <= D_max diagonal, ||T D|| <= ||T D_max||. The Hermitian
    # matrix with diag(1) minus sqrt(upper) G sqrt(upper) restricts to
    # a principal submatrix when a raster footprint is clipped.
    coercive_norm = np.max(high) / full_floor if full_floor > 0 else np.inf
    return dict(T=T, inv_gap=inv_gap, inv_norm=inv_norm, d_op=d_op, coupling=coupling,
                coercive_norm=coercive_norm, dmin=dmin, dmax=dmax,
                max_d=max_d, theta=floor_theta, floor=full_floor)


def bound_currents(local, J0norm, d_norm, exterior0, b, Lambda, use_coercivity):
    L = local["coupling"]
    if use_coercivity:
        t = local["coercive_norm"]
        remaining = 1.0
    else:
        t = local["inv_norm"]
        remaining = 1.0 - local["d_op"]
    available = 1 - b * Lambda
    numerator = b * L
    if not np.isfinite(t) or available <= 0:
        return None
    feedback = t * L * numerator / available
    gap = remaining - feedback
    if gap <= 100 * ARITH:
        return None
    active = (t*d_norm + t*L*b*exterior0/available) / gap
    exterior = b * (exterior0 + L*active) / available
    return add_bound(active), add_bound(exterior)


def apply_bounds(ports, local, Y, r, B_A, S_A, J0, P_A, F0, b, Lambda):
    delta = Y - F0
    residual = fnorm(delta, axis=(1, 2))
    W = delta / np.maximum(residual[:, None, None], 1e-100)
    Wnorm = add_bound(fnorm(W, axis=(1, 2)))
    separation = (W.conj() * delta).sum(axis=(1, 2)).real

    norm_B = add_bound(fnorm(B_A, axis=(1, 2)))
    norm_S = np.broadcast_to(ports.snorm, norm_B.shape)
    norm_Z = add_bound(fnorm(J0, axis=(1, 2)))
    norm_T_error = local["inv_gap"]
    loss_J = norm_T_error * norm_B + ARITH * (
        norm_Z + fnorm(local["T"]) * norm_B
    )
    loss_P = norm_T_error * norm_S + ARITH * (
        norm_S * fnorm(local["T"])
    )
    loss_F = (norm_S * loss_J + ARITH * (
        (norm_Z + loss_J) * norm_S + fnorm(F0, axis=(1, 2)) + fnorm(Y)
    ))

    WA = W.conj().transpose(0, 2, 1) @ P_A
    norm_WA = add_bound(fnorm(WA, axis=(1, 2)))
    wa_error = Wnorm * loss_P + ARITH * Wnorm * fnorm(P_A, axis=(1, 2))
    norm_WA += wa_error
    gradient = np.einsum("bsi,bis->bi", WA, J0).real
    dabs = np.max(local["max_d"])
    remain = np.maximum(gradient * local["dmin"], gradient * local["dmax"]).sum(axis=1)
    remain += dabs * (wa_error * norm_Z + norm_WA * loss_J + ARITH * norm_Z * norm_WA)

    d_norm = add_bound(fnorm(local["max_d"][None, :, None] * J0, axis=(1, 2))) + dabs*loss_J
    L = local["coupling"]
    background_term = np.zeros_like(residual)
    exterior0 = np.zeros_like(residual)
    direction_E = np.zeros_like(residual)
    if b:
        # Keep the positive-part support of the exact free-background
        # derivative. Corrections due to active interaction are enclosed
        # separately by their projected norm.
        norm_WSE = ports.sensor_direction(W)
        direction_E = norm_WSE + L * norm_WA
        scattered_background = L * (norm_Z + loss_J)
        direct_B = np.sqrt(np.maximum(0.,
            ports.bnorm**2 - fnorm(B_A, axis=(1, 2))**2
        )) + ARITH * ports.bnorm
        exterior0 = direct_B + scattered_background
        background_term = b * (ports.background_linear(W) +
                               direction_E * scattered_background + L * norm_WA * direct_B)

    passed = []
    for coercive in (False, True):
        bounds = bound_currents(local, norm_Z, d_norm, exterior0, b, Lambda, coercive)
        if bounds is None:
            passed.append(np.full(len(delta), -np.inf))
            continue
        pA, pE = bounds
        remainder = dabs * norm_WA * pA + b * direction_E * (L*pA + Lambda*pE)
        tolerance = ARITH * (
            1. + np.abs(separation) + remain + background_term + remainder + fnorm(Y)
        )
        passed.append(separation - r*Wnorm - Wnorm*loss_F - remain -
                      background_term - remainder - tolerance)
    return residual, passed


def osm(g, Y):
    # Exact sampled Born atoms with their column-specific norms.
    n = g.side
    S = g.S / np.linalg.norm(g.S, axis=0)
    B = g.B / np.linalg.norm(g.B, axis=1)[:, None]
    return np.abs(np.einsum("ri,rs,is->i", S.conj(), Y, B.conj())).reshape(n, n)


def write_shape_summary(g, shape, observations, maps, costs, side, label):
    result = {}
    prep = np.asarray(costs)
    for oid, data in observations.items():
        nom = maps[oid]
        result[oid] = {}
        for b, registered in nom.items():
            arrays = {}
            for k, v in registered.items():
                arrays[k] = v
            Y, r = data
            arrays["osm"] = osm(g, Y)
            prefix = HERE / "results" / f"{oid}__b{b:g}__{label}"
            np.savez_compressed(prefix.with_suffix(".npz") if "." not in prefix.name
                                else prefix.parent / (prefix.name + ".npz"), **arrays)
            maxbound = np.maximum(arrays["margin_neumann"], arrays["margin_coercive"])
            result[oid][str(b)] = dict(
                centres=g.n,
                excluded=int(np.count_nonzero(maxbound > 0)),
                neumann_excluded=int(np.count_nonzero(arrays["margin_neumann"] > 0)),
                coercivity_excluded=int(np.count_nonzero(arrays["margin_coercive"] > 0)),
                coercivity_extra=int(np.count_nonzero((arrays["margin_coercive"] > 0) &
                                                      ~(arrays["margin_neumann"] > 0))),
                coercivity_bound_available=int(np.count_nonzero(np.isfinite(arrays["margin_coercive"]))),
                neumann_bound_available=int(np.count_nonzero(np.isfinite(arrays["margin_neumann"]))),
            )
    save_json(HERE / "results" / f"map_{side}_{shape}_{label}.json",
              dict(summary=result, costs=dict(unique_clips=len(costs),
                                             wall_s=float(prep[:, 0].sum()),
                                             cpu_s=float(prep[:, 1].sum())),
                   interval=label))


def main(side, shape, label, lo, hi):
    start = clock()
    g = Geometry(geom_path(side))
    ports = ReceiverGeometry(g)
    op = json.loads((HERE / "results" / "operator_bounds.json").read_text())[str(side)]
    acq = json.loads((HERE / "public" / "acquisition.json").read_text())
    observations = {}
    for key, meta in acq["observations"].items():
        if meta["side"] != side or meta["shape"] != shape:
            continue
        with np.load(HERE / "public" / f"{key}.npz", allow_pickle=False) as f:
            observations[key] = f["Y"], float(f["r"])
    maps = defaultdict(dict)
    b_pos = FREEZE["model"]["background_upper"]
    for scene in observations:
        values = (b_pos,) if "__weak__" in scene else (0.0, b_pos)
        for b in values:
            maps[scene][b] = {k: np.full(g.n, -np.inf)
                             for k in ("margin_neumann", "margin_coercive", "reference_error")}
    Lambda = op["full_operator_norm_upper"]
    offsets, p = pattern(side, shape)
    loc = offsets + side//2
    green = g.restrict(loc[:, 1] * side + loc[:, 0])
    floor, theta = active_coercivity(green, p * hi, op["kernel_entry_error"])
    costs = []
    for ids, cells, profile in local_models(g, shape):
        one = clock()
        local = local_inverse(g, cells[0], profile, lo, hi, floor, theta, op)
        if not np.isfinite(local["inv_gap"]):
            continue
        for chunk in range(0, len(ids), 128):
            nodes = ids[chunk:chunk+128]
            chunk_cells = cells[chunk:chunk+128]
            B_A = g.B[chunk_cells]
            S_A = g.S[:, chunk_cells].transpose(1, 0, 2)
            J0 = local["T"] @ B_A
            P_A = S_A @ local["T"]
            F0 = S_A @ J0
            for key, (Y, r) in observations.items():
                for b, arrays in maps[key].items():
                    err, (neumann, coercive) = apply_bounds(
                        ports, local, Y, r, B_A, S_A, J0, P_A, F0, b, Lambda
                    )
                    arrays["margin_neumann"][nodes] = neumann
                    arrays["margin_coercive"][nodes] = coercive
                    arrays["reference_error"][nodes] = err
        ends = elapsed(one)
        costs.append((ends["wall_s"], ends["cpu_s"]))
    write_shape_summary(g, shape, observations, maps, costs, side, label)
    info = json.loads((HERE / "results" / f"map_{side}_{shape}_{label}.json").read_text())
    print(side, shape, label, "coercive", floor, f"angle={np.rad2deg(theta):g}", elapsed(start),
          {k.split("__")[1] + k.split("__")[2]: [x["excluded"] for x in v.values()]
           for k, v in info["summary"].items()}, flush=True)


if __name__ == "__main__":
    arg = argparse.ArgumentParser()
    arg.add_argument("side", type=int)
    arg.add_argument("shape")
    arg.add_argument("label")
    arg.add_argument("lower", type=float)
    arg.add_argument("upper", type=float)
    parsed = arg.parse_args()
    with threadpool_limits(limits=1, user_api="blas"):
        main(parsed.side, parsed.shape, parsed.label, parsed.lower, parsed.upper)
