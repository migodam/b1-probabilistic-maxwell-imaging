"""Physical quadratic assembly, inertia proposals, and Arb enclosures."""
from dataclasses import dataclass

from common import REGION, herm, star
import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh, solve

ctx.prec = 160
ctx.threads = 1
BASE_WEIGHTS = 2.0 ** np.arange(-13, 2, 2)
OPTICAL = (0., -.5, .5)
LOCAL = (1., .25, 4.)


def ball(array):
    array = np.asarray(array)
    return acb_mat(*array.shape, [acb(float(v.real), float(v.imag)) for v in array.flat])


def bt(matrix):
    return matrix.conjugate().transpose()


def frob(matrix):
    return sum((abs(v)**2 for v in matrix.entries()), arb(0)).sqrt()


def low(value):
    return float(np.nextafter(float(value.lower()), -np.inf))


def high(value):
    return float(np.nextafter(float(value.upper()), np.inf))


def midpoint(matrix):
    return np.array([complex(v) for v in matrix.entries()]).reshape(matrix.nrows(), matrix.ncols())


def bdiag(diagonal, matrix):
    return acb_mat(matrix.nrows(), matrix.ncols(),
                   [a * matrix[i, j] for i, a in enumerate(diagonal) for j in range(matrix.ncols())])


def gersh(matrix):
    return min((matrix[i, i].real - sum((abs(matrix[i, j]) for j in range(matrix.ncols())
                                        if i != j), arb(0))).lower()
               for i in range(matrix.nrows()))


def propose_frame(matrix, count):
    if count == 0:
        return np.zeros((len(matrix), 0), dtype=complex), None
    eig, vec = eigh(herm(matrix), check_finite=False, driver="evd")
    edge = float(eig[-count])
    if edge <= 0.:
        return None, edge
    return vec[:, -count:] / np.sqrt(eig[-count:]), edge


def validate_frame(matrix, frame):
    p = ball(frame)
    gram = bt(p) * matrix * p
    gap = gersh(gram)
    return gap, dict(positive_count=frame.shape[1], stretched_positive_gap_lower=low(gap),
                     frame_rank_proved_by_positive_pullback=bool(gap > 0))


def strict_pd(matrix):
    frame, edge = propose_frame(midpoint(matrix), matrix.nrows())
    if frame is None:
        return False, dict(eigenvalue_proposal=edge)
    gap, record = validate_frame(matrix, frame)
    return bool(gap > 0), {**record, "eigenvalue_proposal": edge}


def shift(matrix, n, v):
    if isinstance(matrix, np.ndarray):
        matrix = matrix.copy()
        i = np.arange(n, len(matrix))
        matrix[i, i] -= v
        return matrix
    matrix = acb_mat(matrix)
    for i in range(n, matrix.nrows()):
        matrix[i, i] -= arb(v)
    return matrix


def block(aa, ab, bb):
    if isinstance(aa, np.ndarray):
        return np.block([[aa, ab], [star(ab), bb]])
    n, s = aa.nrows(), bb.nrows()
    return acb_mat(n + s, n + s, [aa[i, j] if i < n and j < n else
                                 ab[i, j - n] if i < n else
                                 ab[j, i - n].conjugate() if j < n else bb[i - n, j - n]
                                 for i in range(n+s) for j in range(n+s)])


def select(matrix, rows, columns):
    if isinstance(matrix, np.ndarray):
        return matrix[np.ix_(rows, columns)]
    return acb_mat(len(rows), len(columns), [matrix[int(i), int(j)] for i in rows for j in columns])


def rows(matrix, ids):
    n = matrix.ncols() if isinstance(matrix, acb_mat) else matrix.shape[1]
    return select(matrix, ids, range(n))


@dataclass
class Acquisition:
    G: np.ndarray
    S: np.ndarray
    B: np.ndarray
    Y: np.ndarray
    radius: float

    def __post_init__(self):
        self.n, self.s = self.B.shape
        self.J = np.hstack((np.eye(self.n), np.zeros((self.n, self.s))))
        self.E = np.hstack((self.G, self.B))
        self.O = np.hstack((self.S, -self.Y))
        self.receiver = herm(star(self.O) @ self.O)
        self.Gb, self.Sb, self.Bb, self.Yb = [ball(x) for x in (self.G, self.S, self.B, self.Y)]
        self.Jb, self.Eb, self.Ob = [ball(x) for x in (self.J, self.E, self.O)]
        self.receiver_ball = bt(self.Ob) * self.Ob

    def trace_update(self, ell, upp, mu, local, optical, selected=()):
        """Return numeric/rigorous local block, in [j_R, (GJ)_R, Z] coordinates."""
        ell = ell[REGION]
        upp = upp[REGION]
        w = np.full(len(REGION), mu * local)
        lam = np.full(len(REGION), mu * optical)
        for i in selected:
            w[i], lam[i] = 0., 0.
        dtype = complex
        q = len(REGION)
        ident = np.eye(q)
        j = np.hstack((ident, np.zeros((q, q + self.s))))
        e = np.hstack((np.zeros((q, q)), ident, self.B[REGION]))
        return sector(j, e, ell, upp, w, lam), sector_ball(ball(j), ball(e), ell, upp, w, lam)

    def strict_witness(self, material, ell, upp):
        if np.any(material < ell) or np.any(material > upp):
            return {"proved": False, "reason": "outside_material_box"}
        try:
            c = [arb(float(v)) for v in material]
            M = ball(np.eye(self.n)) - bdiag(c, self.Gb)
            exact = M.solve(bdiag(c, self.Bb))
            port = self.Sb * exact
            distance = frob(port - self.Yb)
            return dict(proved=bool(distance <= arb(self.radius)), material_in_box=True,
                        verified_original_square_solve=True, residual_lower=low(distance),
                        residual_upper=high(distance), radius=self.radius)
        except (ZeroDivisionError, ValueError) as exc:
            return dict(proved=False, reason=str(exc))


def sector(j, e, ell, upp, weight, optical):
    middle = (ell + upp) / 2
    halfwidth = (upp - ell) / 2
    u = j - middle[:, None] * e
    a = star(j) @ ((-0.5j * optical)[:, None] * e)
    return herm(star(u) @ (weight[:, None] * u) -
                star(e) @ ((weight * halfwidth**2)[:, None] * e) + a + star(a))


def sector_ball(j, e, ell, upp, weight, optical):
    l, u, w = [[arb(float(v)) for v in values] for values in (ell, upp, weight)]
    mid = [(x+y)/2 for x, y in zip(l, u)]
    hw2 = [((y-x)/2)**2 * z for x, y, z in zip(l, u, w)]
    uj = j - bdiag(mid, e)
    a = bt(j) * bdiag([acb(0, -arb(float(v))/2) for v in optical], e)
    return bt(uj) * bdiag(w, uj) - bt(e) * bdiag(hw2, e) + a + bt(a)


def weight_vector(n, mu, local, optical, skipped=()):
    w, v = np.full(n, mu), np.zeros(n)
    w[REGION], v[REGION] = mu * local, mu * optical
    if skipped:
        w[REGION[list(skipped)]] = 0
        v[REGION[list(skipped)]] = 0
    return w, v


def assemble(acq, ell, upp, mu, local, optical, strict=False, skipped=()):
    weight, opt = weight_vector(acq.n, mu, local, optical, skipped)
    if strict:
        return acq.receiver_ball + sector_ball(acq.Jb, acq.Eb, ell, upp, weight, opt)
    return acq.receiver + sector(acq.J, acq.E, ell, upp, weight, opt)


def tail_proposal(matrix, states, radius):
    try:
        d, p = eigh(herm(matrix[:states, :states]), check_finite=False, driver="evd")
        if min(abs(d)) < 1e-12 * max(1., max(abs(d))):
            return dict(strength=0., nu=None, terms=[], K_gap=float(min(abs(d))))
        nu = int(np.count_nonzero(d < 0))
        if nu >= len(matrix) - states:
            return dict(strength=0., nu=nu, terms=[], K_gap=float(min(abs(d))))
        l = star(p) @ matrix[:states, states:]
        schur = herm(matrix[states:, states:] - star(l) @ (l / d[:, None]))
        eigenvalues = eigh(schur, eigvals_only=True, check_finite=False)[::-1][nu:]
        terms = [(float(v) * (1-1e-5), j) for j, v in enumerate(eigenvalues, 1) if v > 0.]
        strength = np.sqrt(sum(t for t, j in terms)) / radius
        return dict(strength=float(strength), nu=nu, terms=terms, K_gap=float(min(abs(d))))
    except np.linalg.LinAlgError:
        return dict(strength=0., nu=None, terms=[], K_gap=None)


def threshold_options(candidate, radius):
    """Single jt proofs first. The full tail proof can certify stronger bounds."""
    terms = candidate["terms"]
    best = max(terms, key=lambda row: row[0]*row[1], default=None)
    if best is not None:
        t, j = best
        if t*j > radius**2:
            yield [(t, j, j)]
    if len(terms) > 1 and sum(t for t, _ in terms) > radius**2:
        yield [(t, j, 1) for t, j in terms]


def verify_inertia(matrix, trial, states, terms, radius, save_frame):
    records = []
    total = arb(0)
    for t, j, copies in terms:
        proposed, edge = propose_frame(shift(trial, states, t), states + j)
        if proposed is None:
            return None
        gap, record = validate_frame(shift(matrix, states, t), proposed)
        if not gap > 0:
            return None
        records.append(dict(threshold=t, j=j, term_multiplicity=copies,
                            frame=save_frame(proposed), floated_positive_edge=edge, **record))
        total += arb(t) * copies
    if not total > arb(radius)**2:
        return None
    return dict(terms=records, states=states, proof_lower=low(total.sqrt()),
                strict_distance_over_radius_lower=low(total.sqrt()/arb(radius)))
