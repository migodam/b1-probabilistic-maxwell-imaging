"""Global old-data enclosures from rank-one material rows and integer Bernstein bounds.

Directions from numerical linear algebra only propose scalar inequalities.
All tests that remove material use integer controls and integer inequalities.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import json
import math
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from numba import njit

from data14 import ROOT, FrequencyTM, load_cases

ctx.prec = 128
K = 6
GRID = 1 << 20
MASS = 1 << 42
WEIGHT = 1 << 16
SIGNS = np.array([[(v >> i) & 1 for i in range(K)] for v in range(1 << K)], dtype=np.int64)
LO, HI = .5, 30.
RESULTS = ROOT / "results"


def complex_ball(x):
    return acb(float(np.real(x)), float(np.imag(x)))


@lru_cache(maxsize=160)
def wave_table(geometry_json, delta=0.):
    """An exactly rounded determinant and block numerator at all 64 corners."""
    model = FrequencyTM(json.loads(geometry_json))
    green, incident, receiver = model.kernels(float(delta))[:3]
    k, s = incident.shape
    assert k == K and receiver.shape == (1, K)
    G = [[complex_ball(z) for z in row] for row in green]
    B = [[complex_ball(z) for z in row] for row in incident]
    S = [complex_ball(z) for z in receiver[0]]
    values = []
    for bits in SIGNS:
        m = [arb(HI if x else LO) for x in bits]
        rows = [[int(i == j) - m[i] * G[i][j] for j in range(k)] for i in range(k)]
        q = acb_mat(rows).det()
        nums = [acb_mat([line + [-m[i] * B[i][l]] for i, line in enumerate(rows)]
                        + [S + [0]]).det() for l in range(s)]
        values.append((q, nums))
    return values


def quantize(rows):
    """Represent each ball as an integer times 2^e with absolute error <= 1 unit."""
    peak = max(float(abs(a).upper()) for row in rows for a in row)
    exponent = math.frexp(peak or 1.)[1] - 36
    scale = arb(2) ** (-exponent)
    rows_int = []
    for row in rows:
        out = []
        for value in row:
            v = value * scale
            n = int((v.mid() + arb(.5)).floor().unique_fmpz())
            if not abs(v - n) < 1:
                raise ArithmeticError("Arb enclosure too wide for quantization")
            out.append(n)
        rows_int.append(out)
    return np.array(rows_int, dtype=np.int64), exponent


@njit(cache=True)
def restrict(table, lo, hi):
    """Each tensor stage rounds down by <1 integer. Absolute error <= 1+k."""
    out = table.copy()
    for axis in range(K):
        stride = 1 << axis
        for v in range(1 << K):
            if v & stride:
                continue
            j = v | stride
            for h in range(table.shape[1]):
                a = out[v, h]
                change = out[j, h] - a
                out[v, h] = a + (change * lo[axis]) // GRID
                out[j, h] = a + (change * hi[axis]) // GRID
    return out


def ceilroot(t):
    i = math.isqrt(t)
    return i + (i * i < t)


@dataclass
class OldTable:
    controls: np.ndarray
    quant_eps: int
    gain_re: np.ndarray
    gain_im: np.ndarray
    p_exponent: int
    q_exponent: int


def case_table(case):
    tab = wave_table(json.dumps(case.geometry, sort_keys=True))
    D = [complex_ball(x) for x in case.old]
    residuals = [[n - z * q for n, z in zip(nums, D)] for q, nums in tab]
    p, ep = quantize([[v.real for v in row] + [v.imag for v in row] for row in residuals])
    den, eq = quantize([[q.real, q.imag] for q, _ in tab])
    # Radii at the public decimal endpoints are rounded outward as well.
    gamma = arb(float(case.gain_radius)) + arb(2) ** -54
    eps = arb(float(case.additive_radius)) + arb(2) ** -54
    scale = arb(MASS) * arb(2) ** (eq - ep) / (1 - gamma)
    eps_upper = int((eps * scale).upper().ceil().unique_fmpz())
    g_vec = [z * (gamma * scale) for z in D]
    def fixed(x):
        a = int((x.mid() + arb(.5)).floor().unique_fmpz())
        assert abs(x-a) < 1
        return a
    return OldTable(np.column_stack([p, den]), eps_upper,
                    np.array([fixed(x.real) for x in g_vec], dtype=np.int64),
                    np.array([fixed(x.imag) for x in g_vec], dtype=np.int64), ep, eq)


def int_directions(p):
    """Any directions work. A small SVD identifies the most informative directions locally."""
    mid = p.mean(axis=0)
    j = p.T.astype(float) @ (2. * SIGNS - 1.) / len(p)
    u, s, vh = np.linalg.svd(j, full_matrices=False)
    inv = (vh.T / np.maximum(s, max(s[0] * 1e-11, 1e-15))) @ u.T
    rows = np.vstack([inv, u.T, mid[None]])
    sizes = np.max(abs(rows), axis=1)
    return np.rint(WEIGHT * rows / np.where(sizes > 0, sizes, 1.)[:, None]).astype(np.int64)


def projected_radii(old, weights, qmax):
    """Upper integer bounds for the common-gain real projection, incl. table interpolation."""
    d = old.gain_re.size
    a, b = weights[:, :d], weights[:, d:]
    # Python integers for the products with 42-bit quantized gains, to avoid overflow.
    sr = a.astype(object) @ old.gain_re.astype(object) + b.astype(object) @ old.gain_im.astype(object)
    si = a.astype(object) @ old.gain_im.astype(object) - b.astype(object) @ old.gain_re.astype(object)
    n1 = np.abs(weights).sum(axis=1)
    n2 = (weights * weights).sum(axis=1)
    t = []
    for real, imag, l1, l2 in zip(sr, si, n1.tolist(), n2.tolist()):
        upper = ceilroot(old.quant_eps ** 2 * l2 + (abs(real) + l1)**2 + (abs(imag) + l1)**2)
        t.append((qmax * upper + MASS - 1) // MASS + (1 + K) * l1)
    return np.array(t, dtype=np.int64)


def contract(old, lo, hi, precondition=True):
    tab = restrict(old.controls, lo, hi)
    p, q = tab[:, :-2], tab[:, -2:]
    qa, qb = (np.max(np.abs(q), axis=0) + K + 1).tolist()
    qmax = ceilroot(qa * qa + qb * qb)
    weights = int_directions(p) if precondition else np.eye(p.shape[1], dtype=np.int64) * WEIGHT
    proj = p @ weights.T
    radii = projected_radii(old, weights, qmax)
    if np.any(proj.min(axis=0) > radii) or np.any(proj.max(axis=0) < -radii):
        return None
    if not precondition:
        return lo, hi
    return material_ranges(lo, hi, proj, radii)


def material_ranges(lo, hi, proj, radii):
    """Round every material coordinate outwards from its six scalar dual bounds."""
    width = hi - lo
    # Use a fully integer certificate for each proposed affine material estimate.
    slopes = np.rint(2.*np.sum((2*SIGNS - 1) * proj[:, :K].astype(float), axis=0)
                     / (len(proj) * np.maximum(1, width))).astype(np.int64)
    # Bound before multiplying. This avoids overflow even on nearly singular trial Jacobians.
    limit = (1 << 58) // np.maximum(width, 1)
    slopes = np.clip(slopes, -limit, limit)
    vals = SIGNS * (slopes * width) - proj[:, :K]
    lower = vals.min(axis=0) - radii[:K]
    upper = vals.max(axis=0) + radii[:K]
    negative = slopes < 0
    lower, upper = np.where(negative, -upper, lower), np.where(negative, -lower, upper)
    slope = abs(slopes)
    # Output boxes represent real continuous material coordinates, so rounding is OUTWARD.
    a = lo + lower // np.maximum(slope, 1)
    b = lo - ((-upper) // np.maximum(slope, 1))
    new_lo = np.where(slope, np.maximum(lo, a), lo)
    new_hi = np.where(slope, np.minimum(hi, b), hi)
    if np.any(new_lo > new_hi):
        return None
    return new_lo, new_hi


def cover(case, leaf_width=.40, max_nodes=18000, precondition=True, start_boxes=None):
    from collections import deque
    wall, cpu = perf_counter(), process_time()
    old = case_table(case)
    setup = perf_counter() - wall
    if start_boxes is None:
        queue = deque([(np.zeros(K, dtype=np.int64), np.full(K, GRID, dtype=np.int64))])
    else:
        # These data are already outward bounds. Integers, a power-of-two grid and
        # an extra grid unit on each side ensure the change of grid encloses them.
        cell = (np.asarray(start_boxes) - LO) * (GRID / (HI - LO))
        lower = np.maximum(0, np.floor(cell[:, 0]).astype(np.int64) - 1)
        upper = np.minimum(GRID, np.ceil(cell[:, 1]).astype(np.int64) + 1)
        queue = deque(zip(lower, upper))
    leaves = []
    visited = rejected = trimmed = 0
    target = max(1, math.ceil(leaf_width * GRID / (HI - LO)))
    while queue and visited < max_nodes:
        lo, hi = queue.popleft()
        out = contract(old, lo, hi, precondition)
        visited += 1
        if out is None:
            rejected += 1
            continue
        new_lo, new_hi = out
        before, after = hi - lo, new_hi - new_lo
        if np.max(after) <= target:
            leaves.append(out)
        elif after.sum() < .77 * before.sum():
            trimmed += 1
            queue.append(out)
        else:
            axis = np.argmax(after)
            mid = (new_lo[axis] + new_hi[axis]) // 2
            a, b = new_hi.copy(), new_lo.copy()
            a[axis] = b[axis] = mid
            queue.extend(((new_lo, a), (b, new_hi)))
    # Every pending box is part of the returned outer set.
    unresolved = len(queue)
    leaves.extend(queue)
    boxes = np.array(leaves, dtype=np.int64).reshape(-1, 2, K)
    f = (HI - LO) / GRID
    widths = (boxes[:, 1] - boxes[:, 0]) * f
    span = widths.max(axis=1) if len(boxes) else np.array([])
    envelope = (boxes[:, 1].max(axis=0) - boxes[:, 0].min(axis=0)) * f if len(boxes) else np.zeros(K)
    rms_hull = (sum((arb(float(w))**2 for w in envelope), arb(0)) / K).sqrt()
    summary = dict(case_id=case.case_id, cohort=case.cohort, visited=visited, rejected=rejected,
        start_boxes=len(start_boxes) if start_boxes is not None else 1,
        contracted_again=trimmed, returned=len(boxes), unfinished=unresolved,
        outer_volume_fraction=float(np.prod(widths/(HI-LO), axis=1).sum()),
        coordinate_hull_width=envelope.tolist(), rms_hull_upper=math.nextafter(float(rms_hull.upper()), math.inf),
        max_leaf_width=float(np.max(span,initial=0.)), leaf_target=leaf_width, max_nodes=max_nodes,
        setup_wall_s=setup, wall_s=perf_counter()-wall, cpu_s=process_time()-cpu)
    return boxes, summary


def run(names=None, leaf=.4, nodes=18000, tag="", seed=None):
    target = RESULTS / ("material_covers" + tag)
    target.mkdir(exist_ok=True)
    regions = np.load(RESULTS / "pivot_old_regions.npz") if seed else None
    for case in load_cases():
        if names is not None and case.case_id not in names:
            continue
        boxes, row = cover(case, leaf, nodes, start_boxes=regions[case.case_id+"/"+seed] if seed else None)
        np.savez_compressed(target/(case.case_id+".npz"), boxes=boxes)
        (target/(case.case_id+".json")).write_text(json.dumps(row, indent=2)+"\n")
        print(case.case_id, row["visited"], row["returned"], row["unfinished"], row["rms_hull_upper"],
              round(row["wall_s"],2), flush=True)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append")
    ap.add_argument("--leaf", type=float, default=.4)
    ap.add_argument("--nodes", type=int, default=18000)
    ap.add_argument("--tag", default="")
    ap.add_argument("--seed", choices=["r10_scalar", "r13", "r10_shared"])
    a = ap.parse_args()
    run(a.case, a.leaf, a.nodes, a.tag, a.seed)
