"""Outward-rounded finite-operator bounds for candidate boxes.

The supplied matrices and endpoints are exact stored binary64 numbers. Float
factorizations propose certificates, while Arb re-evaluates all decisive bounds.
"""

import numpy as np
from flint import acb, acb_mat, arb
from scipy.linalg import eigh, solve, svd

from pool_common import HERE, clock, elapsed
import sys

sys.path.insert(0, str(HERE.parents[1]))
from core.certificate import ball, frobenius_squared, gershgorin, safe_float, star


def herm(a):
    return (a + a.conj().T) * .5


def norm(a):
    return frobenius_squared(a).sqrt()


def midpoint(a):
    return np.array([[complex(v) for v in row] for row in a.tolist()])


def lower(a):
    return safe_float(a, upper=False)


def upper(a):
    return safe_float(a, upper=True)


def identity(n):
    return acb_mat(n, n, [1 if i == j else 0 for i, j in np.ndindex(n, n)])


def diagonal(values):
    n = len(values)
    return acb_mat(n, n, [values[i] if i == j else 0 for i, j in np.ndindex(n, n)])


def exact_endpoints(lo, hi):
    pair = [(arb(float(a)), arb(float(b))) for a, b in zip(lo, hi, strict=True)]
    return diagonal([(a + b) / 2 for a, b in pair]), diagonal([(b - a) / 2 for a, b in pair])


def spectral_upper(a, reference):
    """Suggest a diagonalization, then bound the generalized Rayleigh quotient."""
    gram = star(a) * a
    reference = herm(reference.conj().T @ reference)
    basis = ball(eigh(reference, check_finite=False)[1])
    diagonalized = star(basis) * gram * basis
    metric = star(basis) * basis
    metric_lower = gershgorin(metric, upper=False)
    return (gershgorin(diagonalized, upper=True) / metric_lower).sqrt().upper()


def quadratic(green, sensor, incident, observed, lo, hi):
    center = .5 * (lo + hi)
    delta = .5 * (hi - lo)
    p = np.eye(green.shape[0]) - center[:, None] * green
    gr = delta[:, None] * green
    br = delta[:, None] * incident
    bm = center[:, None] * incident
    return (herm(sensor.conj().T @ sensor),
            herm(p.conj().T @ p - gr.conj().T @ gr), herm(green / 1j),
            sensor.conj().T @ observed,
            p.conj().T @ bm + gr.conj().T @ br, incident / (2j),
            observed.conj().T @ observed, herm(bm.conj().T @ bm - br.conj().T @ br))


def search(terms):
    g0, g1, g2, l0, l1, l2, c0, c1 = terms
    scale = max(1e-7, float(eigh(g0, eigvals_only=True)[-1]))
    proposals = []
    counts = dict(proposals=0, numeric_only_usable=0, max_negative_inertia=0)
    for optical in (0., -2., 2., -8., 8., -32., 32.):
        for tau in scale * (2. ** np.arange(-13, 4, 2)):
            counts["proposals"] += 1
            a = g0 + tau * (g1 + optical * g2)
            lin = l0 + tau * (l1 - optical * l2)
            const = c0 + tau * c1
            vals, vec = eigh(a, check_finite=False)
            nu = np.count_nonzero(vals < 0)
            counts["max_negative_inertia"] = max(counts["max_negative_inertia"], int(nu))
            if nu >= l0.shape[1] or np.min(np.abs(vals)) < 1e-12:
                continue
            x = vec @ ((vec.conj().T @ lin) / vals[:, None])
            source = eigh(herm(const - lin.conj().T @ x), eigvals_only=True, check_finite=False)[::-1]
            estimate = float(np.sqrt(np.maximum(source[nu:], 0).sum()))
            counts["numeric_only_usable"] += 1
            proposals.append((estimate, float(tau), optical, int(nu), vals, vec))
    return sorted(proposals, key=lambda p: p[0], reverse=True), counts


class Bounds:
    def __init__(self, g, s, b, y):
        self.g, self.s, self.b, self.y = g, s, b, y
        self._exact = None

    def exact(self):
        if self._exact is None:
            g, s, b, y = map(ball, (self.g, self.s, self.b, self.y))
            self._exact = (g, s, b, y, identity(g.nrows()), star(s)*s, star(s)*y, star(y)*y)
        return self._exact

    def shifted(self, lo, hi, radius):
        start = clock()
        proposals, counts = search(quadratic(self.g, self.s, self.b, self.y, lo, hi))
        best_num = proposals[0][0] if proposals else 0.
        output = dict(strict_lower=0., numeric_lower=best_num, proof=None, errors=[], **counts)
        for entry in proposals[:3]:
            if entry[0] <= radius:
                break
            try:
                result = self.strict_shift(lo, hi, *entry[1:])
                output["strict_lower"] = max(result["strict_lower"], output["strict_lower"])
                if result["strict_lower"] > radius:
                    output["proof"] = result
                    break
            except (ArithmeticError, ZeroDivisionError, ValueError) as exc:
                output["errors"].append(str(exc))
        return {**output, **elapsed(start)}

    def strict_shift(self, lo, hi, tau, gamma, nu, vals, vec):
        g, _, b, _, identity_n, s2, sy, y2 = self.exact()
        m, d = exact_endpoints(lo, hi)
        p = identity_n - m*g
        gr, br, bm = d*g, d*b, m*b
        t, optical = arb(tau), arb(gamma)
        k = s2 + t * (star(p)*p - star(gr)*gr + optical * acb(0, -.5) * (g-star(g)))
        lin = sy + t * (star(p)*bm + star(gr)*br + optical * acb(0, .5)*b)
        constant = y2 + t * (star(bm)*bm - star(br)*br)
        vp = ball(vec[:, nu:])
        vpk = star(vp) * k * vp
        pos = gershgorin(vpk, upper=False)
        if not pos > 0:
            raise ArithmeticError("Positive subspace inertia could not be verified.")
        neg = None
        if nu:
            vn = ball(vec[:, :nu])
            neg = gershgorin(star(vn)*k*vn, upper=True)
            if not neg < 0:
                raise ArithmeticError("Negative subspace inertia could not be verified.")
        h = constant - star(lin) * k.solve(lin)
        _, source_vec = eigh(herm(midpoint(h)), check_finite=False)
        floors = []
        for j in range(nu + 1, b.ncols() + 1):
            v = ball(source_vec[:, -j:])
            floor = gershgorin(star(v) * h * v, upper=False)
            if not floor > 0:
                break
            metric = star(v)*v
            floors.append((floor / gershgorin(metric, upper=True)).lower())
        total = sum(floors, arb(0))
        return dict(strict_lower=max(0., lower(total.sqrt())), negative_inertia=nu,
                    tau=tau, optical_multiplier=gamma, positive_subspace_gap=lower(pos),
                    negative_subspace_upper=upper(neg) if neg is not None else None,
                    source_tail_eigenvalue_lowers=list(map(lower, floors)))

    def neumann(self, lo, hi, radius):
        start = clock()
        center = .5*(lo + hi)
        d = .5*(hi - lo)
        p = np.eye(len(center)) - center[:, None] * self.g
        combined_rhs = np.column_stack((center[:, None] * self.b, np.diag(d)))
        solution = solve(p, combined_rhs, check_finite=False)
        j, transfer = solution[:, :self.b.shape[1]], solution[:, self.b.shape[1]:]
        ec = self.b + self.g @ j
        q = self.g @ transfer
        s = self.s @ transfer
        yc = self.y - self.s @ j
        _, qsv, qvh = svd(q, check_finite=False)
        k = float(qsv[0])
        result = dict(numeric_contraction=k, center_residual=float(np.linalg.norm(yc)),
                      numeric_standard_lower=0., standard_lower=0., directional_lower=0.,
                      center_factorizations=1, center_incident_rhs=self.b.shape[1],
                      center_transfer_rhs=len(center), verified_inverse_factors=0)
        if k >= 1.:
            g, _, _, _, i, _, _, _ = self.exact()
            m, delta = exact_endpoints(lo, hi)
            q_exact = g * (i - m*g).solve(delta)
            witness = ball(qvh[0].conj()[:, None])
            result["verified_inverse_factors"] = 1
            result["strict_contraction_lower"] = lower((norm(q_exact*witness) / norm(witness)).lower())
            return {**result, **elapsed(start)}

        incident_norm = np.linalg.norm(ec) / (1 - k)
        left, _, _ = svd(s, full_matrices=True, check_finite=False)
        all_port = float(np.linalg.norm(yc) - np.linalg.norm(s, 2) * incident_norm)
        result["numeric_standard_lower"] = max(0., all_port)
        proposals = []
        for keep in range(1, left.shape[1] + 1):
            projection = left[:, -keep:].conj().T
            val = np.linalg.norm(projection @ yc) - np.linalg.norm(projection @ s, 2) * incident_norm
            proposals.append(float(val))
        result["numeric_directional_lower"] = max(0., max(proposals))
        if max(all_port, max(proposals)) <= radius:
            return {**result, **elapsed(start)}

        g, sensor, b, y, i, _, _, _ = self.exact()
        m, delta = exact_endpoints(lo, hi)
        pe = i - m*g
        j_exact = pe.solve(m*b)
        t_exact = pe.solve(delta)
        result["verified_inverse_factors"] = 2
        q_exact = g*t_exact
        k_hi = spectral_upper(q_exact, q)
        result["strict_contraction_upper"] = upper(k_hi)
        if not k_hi < 1:
            return {**result, **elapsed(start)}
        gain = (norm(b + g*j_exact).upper() / (1 - k_hi)).upper()
        test, dshift = sensor * t_exact, y - sensor*j_exact
        standard = norm(dshift).lower() - spectral_upper(test, s) * gain
        result["standard_lower"] = max(0., lower(standard))
        j_best = int(np.argmax(proposals)) + 1
        projection = left[:, -j_best:].conj().T
        w = ball(projection)
        w_norm = gershgorin(w*star(w), upper=True).sqrt()
        projected = (norm(w*dshift).lower() - spectral_upper(w*test, projection @ s)*gain) / w_norm
        result["directional_lower"] = max(0., lower(projected), result["standard_lower"])
        return {**result, **elapsed(start)}

    def point_upper(self, material):
        """Strictly verify a numerically selected feasible full-wave material."""
        g, sensor, b, y, i, _, _, _ = self.exact()
        x = diagonal([arb(float(p)) for p in material])
        data = sensor * (i - x*g).solve(x*b)
        return upper(norm(data-y))
