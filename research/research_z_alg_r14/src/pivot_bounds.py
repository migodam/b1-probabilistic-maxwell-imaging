"""Rounded-outward old-data bounds for the original R13 pivot construction."""
from itertools import combinations
from math import inf, nextafter

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import qr

ctx.prec = 96


def ball(z):
    return acb(float(z.real), float(z.imag))


def mat(a):
    return acb_mat([[ball(z) for z in line] for line in a])


def row(a):
    return mat([a])


def sqnorm(values):
    return sum((a.real**2 + a.imag**2 for a in values), arb(0))


def norm(values):
    return sqnorm(values).upper().sqrt().upper()


def pick_cols(A, cols):
    return acb_mat([[A[i, j] for j in cols] for i in range(A.nrows())])


def scale_cols(A, scale):
    return acb_mat([[A[i, j] * scale[j] for j in range(A.ncols())]
                    for i in range(A.nrows())])


def opnorm(A):
    """Upper bound on the spectral norm from Frobenius and Gram row sums."""
    frobenius = sqnorm(A.entries()).upper()
    sums = []
    for i in range(A.ncols()):
        sums.append(sum((abs(sum((A[k, i].conjugate() * A[k, j]
                                   for k in range(A.nrows())), acb(0))).upper()
                         for j in range(A.ncols())), arb(0)).upper())
    gershgorin = max(sums)
    return min(frobenius, gershgorin).sqrt().upper()


def lower_float(x):
    return nextafter(float(x.lower()), -inf)


def upper_float(x):
    return nextafter(float(x.upper()), inf)


def midpoint(lo, hi):
    return lo + (hi - lo) * .5


class OldCase:
    def __init__(self, H, R, B, host, D, eps, gamma, lower, upper):
        self.k, self.s = B.shape
        self.d = self.k - self.s
        self.H, self.R, self.B, self.host, self.D = mat(H), row(np.ravel(R)), mat(B), row(host), row(D)
        self.eps, self.gamma = arb(float(eps)), arb(float(gamma))
        self.lower, self.upper = np.asarray(lower), np.asarray(upper)
        self.eye = mat(np.eye(self.k))
        idx = qr(B.T, pivoting=True, mode='economic')[2][:self.s].tolist()
        free = [j for j in range(self.k) if j not in idx]
        inverse = mat(B[idx]).inv()
        self.left = acb_mat([[inverse[i, idx.index(j)] if j in idx else acb(0)
                              for j in range(self.k)] for i in range(self.s)])
        Nf = mat(B[free]) * inverse
        self.N = acb_mat([[acb(int(i == free.index(j))) if j in free else -Nf[i, idx.index(j)]
                           for j in range(self.k)] for i in range(self.d)])
        self.NH = self.N * self.H
        self.wD = (self.D - self.host) * self.left
        self.pivots, self.chart_q = self.choose_chart(H)
        self.tail = tuple(j for j in range(self.k) if j not in self.pivots)
        self.NHI = pick_cols(self.NH, self.pivots)
        self.HI = pick_cols(self.H, self.pivots)
        self.NI = pick_cols(self.N, self.pivots)

    def choose_chart(self, H):
        """Use a midpoint conditioning proxy for a fixed chart, without measured fits."""
        N = np.array([[complex(float(self.N[i,j].real.mid()), float(self.N[i,j].imag.mid()))
                       for j in range(self.k)] for i in range(self.d)])
        x, radius = (self.lower + self.upper) / 2, (self.upper - self.lower) / 2
        U = np.eye(self.k) - H * x
        q = {}
        for p in combinations(range(self.k), self.d):
            ix = list(p)
            try:
                update = np.linalg.solve(N @ U[:, ix], (N @ H[:, ix]) * radius[ix])
                q[p] = float(np.linalg.norm(update, 2))
            except np.linalg.LinAlgError:
                q[p] = np.inf
        p = min(q, key=q.get)
        return p, q[p]

    def center(self, p):
        A = self.NI - scale_cols(self.NHI, p)
        Ainv = A.inv()
        U = pick_cols(self.eye, self.pivots) - scale_cols(self.HI, p)
        Q = Ainv * self.N
        w = self.wD + (scale_cols(pick_cols(self.R, self.pivots), p) - self.wD * U) * Q
        K = self.left - (self.left * U) * Q
        v = self.R + w * self.H
        return Ainv, w, K, v

    def weighted(self, X):
        gain = self.D * X
        return acb_mat([[self.eps * X[i,j] for j in range(X.ncols())] for i in range(X.nrows())]
                        + [[self.gamma * z for z in gain.entries()]])

    def cover(self, lo, hi):
        """Return contracted material bounds on an entire pivot box and its status."""
        before = (lo, hi)
        p = [arb(float(midpoint(lo[j], hi[j]))) for j in self.pivots]
        r = [max(p0 - arb(float(lo[j])), arb(float(hi[j])) - p0).upper()
             for p0, j in zip(p, self.pivots)]
        try:
            inv, w, K, v = self.center(p)
        except (ZeroDivisionError, ValueError):
            return before, 'center_singular', 0, None
        if any(not z.is_finite() for z in inv.entries()):
            return before, 'center_singular', 0, None
        q = opnorm(inv * scale_cols(self.NHI, r))
        if not q < 1:
            return before, 'neumann_open', 0, float(q)
        K_H = K * self.H
        C = self.weighted(K)
        C_H = C * self.H
        gain_factor = 1 - self.gamma
        delta = ((norm(scale_cols(pick_cols(v, self.pivots), r).entries())
                  + opnorm(scale_cols(pick_cols(C_H, self.pivots), r)) / gain_factor)
                 * opnorm(inv) / (1 - q)).upper()
        lo, hi = lo.copy(), hi.copy()
        dark = 0
        b_min = None
        for j in self.tail:
            rho_H = (norm(pick_cols(C_H, [j]).entries()) / gain_factor
                     + delta * norm(pick_cols(self.NH, [j]).entries())).upper()
            b = (abs(v[0,j]).lower() - rho_H).lower()
            b_min = float(b) if b_min is None else min(b_min, float(b))
            if not b > 0:
                dark += 1
                continue
            mz = w[0,j] / v[0,j]
            C_l = pick_cols(C, [j]) - pick_cols(C_H, [j]) * mz
            N_l = pick_cols(self.N, [j]) - pick_cols(self.NH, [j]) * mz
            rho_l = (norm(C_l.entries()) / gain_factor + delta * norm(N_l.entries())).upper()
            t = (rho_l / b).upper()
            imag = abs(mz.imag).lower()
            if imag > t:
                return None, 'imaginary_exclusion', dark, float(q)
            width = max(arb(0), (t**2 - imag**2).upper()).sqrt().upper()
            lo[j] = max(lo[j], lower_float(mz.real - width))
            hi[j] = min(hi[j], upper_float(mz.real + width))
            if lo[j] > hi[j]:
                return None, 'material_exclusion', dark, float(q)
        return (lo, hi), 'dark_open' if dark else 'bounded', dark, float(q)
