"""R10 cofactor and direct shared-gain material-box controls."""
from itertools import combinations
from math import frexp

import numpy as np
from flint import acb, acb_mat, arb

from pivot_bounds import ball, midpoint, norm, pick_cols, sqnorm


def determinant(columns):
    return acb_mat(columns).det()


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), acb(0))


def bernstein(values):
    """Round downward onto a fixed-point lattice for convex subdivision."""
    for axis in range(values.ndim):
        v = np.moveaxis(values, axis, 0)
        v[1] = 2 * v[1] - (v[0] + v[2]) / 2
    peak = max(float(abs(v).upper()) for v in values.flat)
    scale = arb(2) ** (48 - frexp(peak)[1]) if peak else arb(1)
    controls = [int((v * scale).lower().floor().unique_fmpz()) for v in values.flat]
    assert max(map(abs, controls)) < 2**49
    return np.asarray(controls, dtype=np.int64).reshape(values.shape)


def split(C, axis):
    a, b, c = np.moveaxis(C, axis, 0)
    ab, bc = (a + b) // 2, (b + c) // 2
    m = (ab + bc) // 2
    return np.stack((a, ab, m), axis=axis), np.stack((m, bc, c), axis=axis)


def quadratic_positive(C):
    """Strict positivity in one Bernstein axis, after downward integer scaling."""
    a, b, c = np.moveaxis(C, -1, 0)
    largest = np.maximum(np.maximum(abs(a), abs(b)), abs(c))
    shift = np.maximum(np.frexp(largest.astype(float))[1].astype(np.int64) - 29, 0)
    neg_middle_up = (-b + (1 << shift) - 1) >> shift
    end_product_down = (a >> shift) * (c >> shift)
    return (a > 0) & (c > 0) & ((b >= 0) | (end_product_down > neg_middle_up**2))


def tail_split(controls, depth):
    bins = [controls]
    for _ in range(depth):
        bins = [child for parent in bins for child in split(parent, -1)]
    return np.stack(bins)


class Cofactors:
    """Shared pivot chart with independent and common-error R10 tests."""
    def __init__(self, case, tail_depth=4):
        self.case = case
        self.pivots, self.tail = case.pivots, case.tail
        self.tail_axis = {j: i for i, j in enumerate(self.tail)}
        self.tail_bins = 2**tail_depth
        table = []
        for j, (a, b) in enumerate(zip(case.lower, case.upper)):
            u0 = self.cell(j, arb(a))
            u1 = self.cell(j, (arb(a) + arb(b)) / 2)
            u2 = self.cell(j, arb(b))
            table.append([u0, u1, u2])
        controls = []
        for t in self.tail:
            ids = (*self.pivots, t)
            values = np.empty((3,) * len(ids), dtype=object)
            for index in np.ndindex(values.shape):
                G, r = self.circuit([table[j][i] for j, i in zip(ids, index)])
                values[index] = sqnorm([r]) - sqnorm(G)
            controls.append(tail_split(bernstein(values), tail_depth))
        self.controls = np.stack(controls)

    def cell(self, j, q):
        case = self.case
        U = pick_cols(case.eye, [j]) - pick_cols(case.H, [j]) * q
        columns = (case.N * U).entries()
        x = case.left * U
        c = (case.host * x)[0, 0] + case.R[0, j] * q
        r = (case.D * x)[0, 0] - c
        return columns, x.entries(), c, r

    def circuit(self, cells):
        case = self.case
        cols = [c[0] for c in cells]
        weights = [(-1)**i * determinant(cols[:i] + cols[i+1:]) for i in range(len(cols))]
        x = [dot(weights, xi) for xi in zip(*(c[1] for c in cells))]
        c = dot(weights, [c[2] for c in cells])
        r = dot(weights, [c[3] for c in cells])
        return [case.gamma * c, *(case.eps * xi for xi in x)], r

    def bound(self, lo, hi, controls):
        passed = ~np.all(quadratic_positive(controls), axis=tuple(range(2, controls.ndim - 1)))
        lo, hi = lo.copy(), hi.copy()
        for j, keep in zip(self.tail, passed):
            ids = np.flatnonzero(keep)
            if not len(ids):
                return None
            origin = self.case.lower[j]
            step = (self.case.upper[j] - origin) / self.tail_bins
            lo[j] = max(lo[j], origin + step * ids[0])
            hi[j] = min(hi[j], origin + step * (ids[-1] + 1))
            if lo[j] > hi[j]:
                return None
        return lo, hi

    def evaluate(self, q):
        cells = [self.cell(j, qj) for j, qj in enumerate(q)]
        out = [self.circuit([cells[j] for j in (*self.pivots, t)]) for t in self.tail]
        return [v[0] for v in out], [v[1] for v in out]

    def joint_bound(self, lo, hi):
        """Use one fixed separating direction for all materials in this box."""
        m = midpoint(lo, hi)
        gg, rr = self.evaluate([arb(float(x)) for x in m])
        g = np.asarray([[complex(float(a.real.mid()), float(a.imag.mid())) for a in v] for v in gg]).T
        r = np.asarray([complex(float(a.real.mid()), float(a.imag.mid())) for a in rr])
        try:
            q, T = np.linalg.qr(g, mode='reduced')
            u = np.linalg.solve(T.conj().T, r.conj())
            if np.linalg.norm(u) <= 1:
                return False
            dual = np.linalg.solve(T, u)
        except np.linalg.LinAlgError:
            return False
        scale = max(abs(dual))
        if not np.isfinite(scale) or not scale > 0:
            return False
        direction = [ball(v / scale) for v in dual]
        coords = [arb(float(x)).union(arb(float(y))) for x, y in zip(lo, hi)]
        G, r = self.evaluate(coords)
        return abs(dot(direction, r)).lower() > norm([dot(direction, values) for values in zip(*G)])


def direct_bernstein(case):
    """Six-dimensional interval profile with the exact common complex gain."""
    k, s = case.k, case.s
    H, B, R = case.H, case.B, case.R
    hy, obs = case.host.entries(), case.D.entries()
    eps, gamma = case.eps, case.gamma
    ranges = [(arb(a), (arb(a) + arb(b)) / 2, arb(b))
              for a, b in zip(case.lower, case.upper)]
    pairs = tuple(combinations(range(s), 2))
    g = (eps / gamma)**2
    values = np.empty((3,) * k, dtype=object)
    for ids in np.ndindex(values.shape):
        m = [ranges[j][idx] for j, idx in enumerate(ids)]
        U = [[int(i == j) - H[i, j] * m[j] for j in range(k)] for i in range(k)]
        det = acb_mat(U).det()
        border = [-R[0, j] * m[j] for j in range(k)]
        Y = [acb_mat([U[i] + [B[i, a]] for i in range(k)] + [border + [hy[a]]]).det()
             for a in range(s)]
        residual = [det * obs[a] - Y[a] for a in range(s)]
        wedge = [obs[a] * Y[b] - obs[b] * Y[a] for a, b in pairs]
        values[ids] = sqnorm(wedge) + g * sqnorm(residual) - eps**2 * (sqnorm(Y) + g * sqnorm([det]))
    return bernstein(values)


def profile_witness(case, m):
    """Score a revealed complete material with exact-operator Arb forward and shared error."""
    k = case.k
    H, B, R = case.H, case.B, case.R
    T = acb_mat([[(1 / ball(m[i]) if i == j else 0) - H[i, j] for j in range(k)] for i in range(k)])
    f = (case.host + R * T.solve(B)).entries()
    residual = [d - y for d, y in zip(case.D.entries(), f)]
    g = dot([z.conjugate() for z in f], residual) / (sqnorm(f) + (case.eps / case.gamma)**2)
    return (sqnorm([g]) / case.gamma**2
            + sqnorm([r - g * y for r, y in zip(residual, f)]) / case.eps**2)
