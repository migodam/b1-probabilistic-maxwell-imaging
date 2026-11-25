"""Fixed trace congruence, full bordered reference, and verified VIE witness."""
from dataclasses import dataclass
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, eigvalsh, solve
from scipy.optimize import least_squares


def clock():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(('wall_seconds', 'cpu_seconds'),
                    (a-b for a, b in zip(clock(), start)), strict=True))


def star(x):
    return x.conj().T


def herm(x):
    return (x + star(x)) * .5


def ball(x):
    return acb_mat(*x.shape, [acb(float(v.real), float(v.imag)) for v in x.flat])


def dagger(a):
    return a.conjugate().transpose()


def identity(n):
    return acb_mat(n, n, [int(i == j) for i, j in np.ndindex(n, n)])


def sub(a, rows, cols):
    return acb_mat(len(rows), len(cols), [a[int(i), int(j)] for i in rows for j in cols])


def block(a, b, c, d):
    entries = [v for left, right in ((a, b), (c, d))
               for i in range(left.nrows()) for matrix in (left, right)
               for v in [matrix[i, j] for j in range(matrix.ncols())]]
    return acb_mat(a.nrows()+c.nrows(), a.ncols()+b.ncols(), entries)


def diagonal_times(diag, a):
    return acb_mat(a.nrows(), a.ncols(),
                   [diag[i] * a[i, j] for i, j in np.ndindex(a.nrows(), a.ncols())])


def midpoint(a):
    return np.array([complex(v) for v in a.entries()]).reshape(a.nrows(), a.ncols())


def norm(a):
    return sum((v.real**2 + v.imag**2 for v in a.entries()), arb(0)).sqrt()


def lower(a):
    return float(np.nextafter(float(a.lower()), -np.inf))


def upper(a):
    return float(np.nextafter(float(a.upper()), np.inf))


def gershgorin_floor(a):
    return min((a[i, i].real - sum((abs(a[i, j]) for j in range(a.ncols()) if i != j),
                                  arb(0))).lower() for i in range(a.nrows()))


def sign_floor(a, vectors, sign):
    """Check a proposed frame in the exact interval matrix."""
    if not vectors.shape[1]:
        return None, None
    q = ball(vectors)
    metric = gershgorin_floor(dagger(q) * q)
    margin = gershgorin_floor(sign * (dagger(q) * a * q))
    if not metric > 0 or not margin > 0:
        return None, None
    return lower(margin), lower(metric)


def interval_inertia(a, predicted=None, full_count=False, frames=None):
    values, directions = eigh(herm(midpoint(a)), check_finite=False, driver='evd')
    count = int(np.count_nonzero(values > 0))
    k = count if full_count else predicted
    if k > count:
        return None
    p, metric = sign_floor(a, directions[:, -k:] if k else directions[:, :0], 1)
    if k and p is None:
        return None
    nfloor, nmetric = None, None
    negative = len(values)-count
    if full_count:
        if not np.all(values != 0):
            return None
        nfloor, nmetric = sign_floor(a, directions[:, :negative], -1)
        if negative and nfloor is None:
            return None
    if frames is not None:
        frames['positive'] = directions[:, -k:] if k else directions[:, :0]
        if full_count:
            frames['negative'] = directions[:, :negative]
    return dict(verified_positive=int(k), positive_gap_lower=p, positive_frame_gram_lower=metric,
                verified_negative=negative if full_count else None, negative_gap_lower=nfloor,
                negative_frame_gram_lower=nmetric, proposal_positive=count, dimension=len(values))


def full_float(G, S, B, Y, ell, cap, w):
    m = .5 * (ell + cap)
    U = np.eye(len(ell)) - m[:, None] * G
    W = (ell * cap - m*m)
    F = m[:, None] * B
    K = herm(star(S) @ S + w * (star(U) @ U + star(G) @ (W[:, None] * G)))
    L = star(S) @ Y + w * (star(U) @ F - star(G) @ (W[:, None] * B))
    C = herm(star(Y) @ Y + w * star(B) @ ((ell*cap)[:, None] * B))
    return K, L, C


def full_ball(G, S, B, Y, ell, cap, w):
    l, u = ([arb(float(v)) for v in x] for x in (ell, cap))
    m = [(a+b)/2 for a, b in zip(l, u, strict=True)]
    delta2 = [-((b-a)/2)**2 for a, b in zip(l, u, strict=True)]
    product = [a*b for a, b in zip(l, u, strict=True)]
    U = identity(len(l)) - diagonal_times(m, G)
    F = diagonal_times(m, B)
    K = dagger(S)*S + w*(dagger(U)*U + dagger(G)*diagonal_times(delta2, G))
    L = dagger(S)*Y + w*(dagger(U)*F - dagger(G)*diagonal_times(delta2, B))
    C = dagger(Y)*Y + w*dagger(B)*diagonal_times(product, B)
    return K, L, C


def residual_proposal(state, coupling, corner, adjusted_negative=0):
    v, q = eigh(state, check_finite=False, driver='evd')
    if not len(v):
        return 0, eigvalsh(herm(corner), check_finite=False)
    if np.min(abs(v)) <= 1e-12 * np.max(abs(v)):
        return None
    nu = int(np.count_nonzero(v < 0)) - adjusted_negative
    if nu < 0:
        return None
    schur = corner - star(coupling) @ (q @ ((star(q) @ coupling) / v[:, None]))
    return nu, eigvalsh(herm(schur), check_finite=False)


def threshold(proposal, radius):
    if proposal is None:
        return None, None, None, None
    nu, lam = proposal
    avail = lam[:max(0, len(lam)-nu)][::-1]
    if not len(avail):
        return None, None, 0., nu
    jmax = np.arange(1, len(avail)+1)
    best = int(np.argmax(jmax*avail))
    numeric = float(np.sqrt(np.maximum(avail, 0).sum()) / radius)
    j = best + 1
    t = float(.5 * (radius**2/j + avail[best]))
    if not j*t > radius**2:
        return None, None, numeric, nu
    return j, t, numeric, nu


@dataclass
class Geometry:
    side: int
    sensor: np.ndarray
    incident: np.ndarray
    green: np.ndarray | None = None
    spectrum: np.ndarray | None = None

    def __post_init__(self):
        if self.spectrum is not None:
            self.kernel = np.fft.ifft2(self.spectrum)
        self.ball_kernel = None

    @property
    def n(self):
        return self.side**2

    def _index(self, cells):
        y, x = np.divmod(cells, self.side)
        pad = self.spectrum.shape[0]
        return (y[:, None]-y[None, :]) % pad, (x[:, None]-x[None, :]) % pad

    def restrict(self, cells):
        g = (self.kernel[self._index(cells)] if self.green is None
             else self.green[np.ix_(cells, cells)])
        return g, self.sensor[:, cells], self.incident[cells]

    def prepare_arb(self):
        if self.spectrum is not None and self.ball_kernel is None:
            horizontal = [acb.dft([acb(float(z.real), float(z.imag)) for z in row], inverse=True)
                          for row in self.spectrum]
            self.ball_kernel = [acb.dft(list(col), inverse=True) for col in zip(*horizontal)]

    def restrict_arb(self, cells):
        if self.green is not None:
            return tuple(ball(x) for x in self.restrict(cells))
        self.prepare_arb()
        iy, ix = self._index(cells)
        g = acb_mat(len(cells), len(cells), [self.ball_kernel[int(x)][int(y)]
                                           for y, x in zip(iy.flat, ix.flat)])
        return g, ball(self.sensor[:, cells]), ball(self.incident[cells])


class FixedKernel:
    """Auxiliary congruence of a full n-state bordered sector certificate.

    Kb is the common [0,b] state block. If all other endpoints change only on
    R, the 2|R| trace auxiliary uses a nonsingular 2×2 D difference per row.
    Inertia(full threshold border) = inertia(Reduced) + (n,0,0) - (|R|,|R|,0).
    """

    def __init__(self, geometry, Y, b, weight):
        self.n = geometry.n
        self.Y = Y
        self.b = b
        self.weight = weight
        self.g = geometry.green
        self.s = geometry.sensor
        self.e = geometry.incident
        n = self.n
        G, S, B = self.g, self.s, self.e
        K = herm(star(S) @ S + weight*(np.eye(n)-.5*b*(G+star(G))))
        L = star(S) @ Y + .5*b*weight*B
        C = star(Y) @ Y
        T = np.vstack((np.eye(n), G))
        E = np.vstack((np.zeros_like(B), B))
        X = solve(K, np.column_stack((star(T), L)), assume_a='pos', check_finite=False)
        self.gram = herm(T @ X[:, :2*n])
        self.coupling = E + T @ X[:, 2*n:]
        self.corner = herm(C - star(L) @ X[:, 2*n:])
        self.strict = None

    def prepare_arb(self):
        if self.strict is not None:
            return
        n = self.n
        G, S, B, Y = (ball(x) for x in (self.g, self.s, self.e, self.Y))
        b, w = arb(self.b), arb(self.weight)
        K = dagger(S)*S + w*(identity(n)-(b/2)*(G+dagger(G)))
        L = dagger(S)*Y + (b*w/2)*B
        C = dagger(Y)*Y
        bound = sum((abs(z)**2 for z in G.entries()), arb(0)).sqrt()
        floor = w*(1-b*bound)
        if not floor > 0:
            raise ArithmeticError("The common state block lacks the declared norm-based floor.")
        T = acb_mat(2*n, n, identity(n).entries()+G.entries())
        E = acb_mat(2*n, B.ncols(), acb_mat(n, B.ncols()).entries()+B.entries())
        rhs = acb_mat(n, 2*n+B.ncols(), [x for i in range(n) for M in (dagger(T), L)
                                        for x in [M[i, j] for j in range(M.ncols())]])
        X = K.solve(rhs)
        gram = T * sub(X, range(n), range(2*n))
        tail_x = sub(X, range(n), range(2*n, 2*n+B.ncols()))
        N = E + T*tail_x
        H = C - dagger(L)*tail_x
        self.strict = gram, N, H
        self.positive_floor = lower(floor)

    def _changed(self, ell, cap):
        where = np.flatnonzero((ell != 0.) | (cap != self.b))
        return where, np.r_[where, self.n+where]

    def float_reduced(self, ell, cap):
        ids, traces = self._changed(ell, cap)
        k = len(ids)
        R = -self.gram[np.ix_(traces, traces)]
        w = self.weight
        a = -w * .5 * (ell[ids]+cap[ids]-self.b)
        if np.any(a == 0.):
            raise ArithmeticError("Singular local D difference.")
        ix = np.arange(k)
        R[ix, ix] += w * ell[ids] * cap[ids] / (a*a)
        R[ix, k+ix] -= 1/a
        R[k+ix, ix] -= 1/a
        return R, self.coupling[traces], self.corner, k

    def interval_reduced(self, ell, cap, t):
        self.prepare_arb()
        ids, traces = self._changed(ell, cap)
        k = len(ids)
        m, N, H = self.strict
        R = -sub(m, traces, traces)
        w = arb(self.weight)
        for i, loc in enumerate(ids):
            l, u = arb(float(ell[loc])), arb(float(cap[loc]))
            a = -w * (l+u-arb(self.b))/2
            if not abs(a) > 0:
                raise ArithmeticError("Singular local D difference.")
            R[i, i] += w*l*u/(a*a)
            R[i, k+i] -= 1/a
            R[k+i, i] -= 1/a
        N = sub(N, traces, range(H.ncols()))
        A = block(R, N, dagger(N), H-arb(t)*identity(H.nrows()))
        return A, k

    def certify(self, ell, cap, radius):
        R, N, H, k = self.float_reduced(ell, cap)
        j, t, numeric, nu = threshold(residual_proposal(R, N, H, k), radius)
        row = dict(numeric_tail_over_radius=numeric, numeric_state_negative=nu,
                   changed_rows=k, state_dimension=self.n, certificate=None)
        if j is None or not arb(j)*arb(t) >= arb(radius)**2:
            return row
        A, count = self.interval_reduced(ell, cap, t)
        proof = interval_inertia(A, count+j)
        if proof is not None:
            row["certificate"] = dict(**proof, j=j, t=t, weight=self.weight,
                   precision_bits=ctx.prec, common_state_floor_lower=self.positive_floor,
                   residual_over_radius_lower=lower((arb(j)*arb(t)).sqrt()/arb(radius)),
                   verified_full_positive_at_least=self.n+j,
                   algorithm="exact_fixed_trace_augmented_border")
        return row


def full_border(K, L, C, t):
    return np.block([[K, -L], [-star(L), C-t*np.eye(C.shape[0])]])


def full_interval_border(K, L, C, t):
    return block(K, -L, -dagger(L), C-arb(t)*identity(C.ncols()))


def certify_direct(geom, Y, ell, cap, radius, weights, max_cells):
    ids = np.flatnonzero(cap)
    row = dict(state_dimension=len(ids), weights_tried=0, certificate=None)
    if len(ids) > max_cells:
        row["reason"] = "retained vacuum-support dimension exceeds node budget"
        return row
    G, S, B = geom.restrict(ids)
    best = None
    for w in weights:
        row["weights_tried"] += 1
        K, L, C = full_float(G, S, B, Y, ell[ids], cap[ids], w)
        j, t, numeric, nu = threshold(residual_proposal(K, L, C), radius)
        if numeric is not None and (best is None or numeric > best):
            best = numeric
        if j is None or not arb(j)*arb(t) >= arb(radius)**2:
            continue
        g, s, b = geom.restrict_arb(ids)
        A = full_interval_border(*full_ball(g, s, b, ball(Y), ell[ids], cap[ids], arb(w)), t)
        proof = interval_inertia(A, len(ids)+j)
        if proof is not None:
            row["certificate"] = dict(**proof, j=j, t=t, weight=w, precision_bits=ctx.prec,
                    residual_over_radius_lower=lower((arb(j)*arb(t)).sqrt()/arb(radius)),
                    verified_full_positive_at_least=len(ids)+j,
                    algorithm="full_border_after_exact_exterior_vacuum_elimination")
            break
    row["numeric_tail_over_radius"] = best
    return row


def strict_witness(geom, Y, ell, cap, radius, material):
    if not (np.all(ell <= material) and np.all(material <= cap)):
        return None
    ids = np.flatnonzero(material)
    if len(ids):
        g, s, b = geom.restrict_arb(ids)
        chi = list(map(arb, map(float, material[ids])))
        try:
            j = (identity(len(ids))-diagonal_times(chi, g)).solve(diagonal_times(chi, b))
        except (ZeroDivisionError, ValueError):
            return None
        residual = norm(s*j-ball(Y))
    else:
        residual = norm(ball(Y))
    if not residual <= arb(radius):
        return None
    return dict(nonzero_rows=len(ids), verified_invertible=True, precision_bits=ctx.prec,
                residual_upper=upper(residual), residual_over_radius_upper=upper(residual/arb(radius)))


def float_forward(geom, material):
    cells = np.flatnonzero(material)
    G, S, B = geom.restrict(cells)
    j = solve(np.eye(len(cells))-material[cells, None]*G, material[cells, None]*B,
              check_finite=False)
    return S @ j


def witness_search(geom, Y, ell, cap, radius, max_evaluations=9):
    ids = np.flatnonzero(cap)
    G, S, B = geom.restrict(ids)
    low, high = ell[ids], cap[ids]
    varying = np.flatnonzero(low < high)
    eye = np.eye(len(ids))
    scale = np.linalg.norm(Y)
    cache = [None, None, None]
    count = dict(forward_evaluations=0, jacobian_evaluations=0)
    # Keep all rows with an upper bound free. Start without an image or truth.
    start = low[varying] + .35*(high-low)[varying]

    def state(x):
        if cache[0] is None or not np.array_equal(cache[0], x):
            m = low.copy()
            m[varying] = x
            a = eye-m[:, None]*G
            J = solve(a, m[:, None]*B, check_finite=False)
            cache[:] = [x.copy(), a, J]
            count["forward_evaluations"] += 1
        return cache[1], cache[2]

    def fun(x):
        _, J = state(x)
        v = (S @ J-Y).ravel() / scale
        return np.r_[v.real, v.imag]

    def jac(x):
        a, J = state(x)
        se = solve(a.T, S.T, check_finite=False).T[:, varying]
        field = (B+G@J)[varying]
        d = np.einsum('rk,ks->rsk', se, field).reshape(-1, len(varying)) / scale
        count["jacobian_evaluations"] += 1
        return np.vstack((d.real, d.imag))

    if len(varying) == 0:
        fit = start
    else:
        result = least_squares(fun, start, jac=jac, bounds=(low[varying], high[varying]),
                               max_nfev=max_evaluations, ftol=1e-10, gtol=1e-10, xtol=1e-10)
        fit = np.clip(result.x, low[varying], high[varying])
    candidate = np.zeros(geom.n)
    candidate[ids] = low
    candidate[ids[varying]] = fit
    numeric = float(np.linalg.norm(fun(fit))*scale/radius)
    proof = strict_witness(geom, Y, ell, cap, radius, candidate) if numeric <= 1.001 else None
    return candidate, proof, dict(**count, residual_over_radius_numeric=numeric)
