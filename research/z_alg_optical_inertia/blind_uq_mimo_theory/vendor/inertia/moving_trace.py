"""Fixed-operator congruence with exact zero update directions and retained base pivots."""
from dataclasses import dataclass
from fractions import Fraction
from math import frexp, isfinite, ldexp

import numpy as np
from scipy.linalg import eigh, eigvalsh, get_lapack_funcs, ldl


def adj(a):
    return a.conj().T


def herm(a):
    return (a + adj(a)) * 0.5


def eig_counts(a):
    """Numerical sign proposals. Near-zero values have no inferred exact sign."""
    if not len(a):
        return (0, 0, 0)
    d = eigvalsh(herm(a), driver='evd', check_finite=False)
    tol = (32 * np.finfo(float).eps * len(a)) * max(np.max(abs(d)), 1.e-300)
    np_, nn = int(np.count_nonzero(d > tol)), int(np.count_nonzero(d < -tol))
    return np_, nn, len(d)-np_-nn


def bk_counts(a):
    """Bunch–Kaufman numerical proposals from 1- and 2-pivot blocks."""
    packed, pivots, info = get_lapack_funcs('hetrf', (a,))(np.array(a, order='F'), lower=1, overwrite_a=1)
    if info < 0:
        raise ValueError(f'hetrf rejected argument {-info}')
    threshold = 32 * np.finfo(float).eps * len(a) * max(np.max(abs(a)), 1.e-300)
    diag = packed.diagonal().real.copy()
    j = 0
    while j < len(pivots):
        if pivots[j] < 0:
            if pivots[j+1] != pivots[j]:
                raise ArithmeticError('invalid Hermitian pivot pair')
            off = packed[j+1, j]
            diag[j:j+2] = eigvalsh([[diag[j], off.conjugate()], [off, diag[j+1]]])
            j += 2
        else:
            j += 1
    positive = int(np.count_nonzero(diag > threshold))
    negative = int(np.count_nonzero(diag < -threshold))
    return positive, negative, len(diag)-positive-negative


def corrected_inertia(pivot_counts, reduced_counts, negative_aux_counts):
    return tuple(a + b - c for a, b, c in zip(pivot_counts, reduced_counts, negative_aux_counts, strict=True))


def propose_pivot(state, tolerance):
    """Propose an ordered nonsingular principal block without dropping state coordinates."""
    _, blocks, ordering = ldl(herm(state), check_finite=False)
    chosen = []
    first = 0
    while first < len(state):
        stop = first + (2 if first+1 < len(state) and blocks[first+1,first] != 0 else 1)
        if np.min(abs(eigvalsh(blocks[first:stop, first:stop], check_finite=False))) > tolerance:
            chosen.extend(ordering[first:stop])
        first = stop
    # A floating pivot pattern is only a proposal. Recheck the actual principal block.
    while chosen:
        principal = state[np.ix_(chosen,chosen)]
        if np.min(abs(eigvalsh(herm(principal), check_finite=False))) > tolerance:
            return chosen
        chosen.pop()
    return chosen


@dataclass(frozen=True)
class LocalUpdate:
    trace_rows: tuple
    scale: np.ndarray
    negative_scaled_inverse: np.ndarray
    negative_aux_counts: tuple
    full_block_dimension: int


def physical_update(ids, lower, upper, base_lower, base_upper, weight):
    """Rank-exact factorization for the two-row (current, field) TM sector update.

    The endpoint rank tests use the exact stored binary inputs. Optical weights
    are held fixed across updates.
    """
    rows, scales, pieces = [], [], []
    pos = neg = 0
    for i, l, u, l0, u0, w in zip(ids, lower, upper, base_lower, base_upper, weight, strict=True):
        if not (0 <= l <= u and 0 <= l0 <= u0 and w > 0 and all(map(isfinite,(u,u0,w)))):
            raise ValueError('real nonnegative endpoints and a positive sector weight are required')
        al,au,bl,bu = (Fraction(float(v)) for v in (l,u,l0,u0))
        dm = ((al+au)-(bl+bu))/2
        dp = al*au-bl*bu
        if not (dm or dp):
            continue
        a,b = -w*float(dm),w*float(dp)
        if (dm and not a) or (dp and not b) or not all(map(isfinite,(a,b))):
            raise ArithmeticError('stored exact update falls outside binary64 scaling')
        # The power of two is an exact trace rescaling; it balances the auxiliary pivots.
        s = ldexp(1., frexp(max(abs(a), abs(b)))[1] // 2)
        if a:
            rows.extend((2*i, 2*i+1))
            scales.extend((s, s))
            pieces.append((a, b, s, True))
            pos += 1
            neg += 1
        else:
            rows.append(2*i+1)
            scales.append(s)
            pieces.append((a, b, s, False))
            pos += b < 0
            neg += b > 0
    inverse = np.zeros((len(rows), len(rows)))
    j = 0
    for a, b, s, double in pieces:
        if double:
            inverse[j:j+2, j:j+2] = [[s*s*b/(a*a), -s*s/a], [-s*s/a, 0]]
            j += 2
        else:
            inverse[j, j] = -s*s/b
            j += 1
    return LocalUpdate(tuple(rows), np.array(scales), inverse, (int(pos), int(neg), 0), 2*len(ids))


@dataclass(frozen=True)
class TraceCache:
    gram: np.ndarray
    source_trace: np.ndarray
    source: np.ndarray
    retained_trace: np.ndarray
    retained_state: np.ndarray
    retained_source: np.ndarray
    pivot_counts: tuple
    pivot_indices: tuple

    @classmethod
    def from_dense(cls, state, coupling, corner, trace, affine, *, pivot=None):
        """Use a nonsingular principal subblock. Every other coordinate is retained.

        `coupling` is H in the original block [[K, -H], [-H*, C]]. With no
        explicit pivot, try the full K before proposing a smaller principal
        block from a pivoted LDL factorization. Uncertain coordinates stay in
        the retained block. Certification must still check the chosen signs.
        """
        n = len(state)
        automatic = pivot is None
        pivot = np.arange(n) if pivot is None else np.asarray(pivot, dtype=int)
        if len(np.unique(pivot)) != len(pivot) or np.any((pivot < 0) | (pivot >= n)):
            raise ValueError('invalid state pivot')
        retained = np.setdiff1d(np.arange(n), pivot)
        q = herm(state[np.ix_(pivot, pivot)])
        pz = state[np.ix_(pivot, retained)]
        lq = coupling[pivot]
        tq = trace[:, pivot]
        d, v = eigh(q, driver='evd', check_finite=False)
        tol = 32 * np.finfo(float).eps * len(q) * max(abs(d), default=0)
        if len(d) and np.any(abs(d) <= tol):
            if automatic:
                return cls.from_dense(state, coupling, corner, trace, affine, pivot=propose_pivot(state,tol))
            raise ArithmeticError('unresolved common pivot, retain its uncertain directions')
        k = lq.shape[1]
        solved = v @ ((adj(v) @ np.column_stack((lq, pz, adj(tq)))) / d[:, None])
        ql, qz, qt = np.split(solved, [k, k+len(retained)], axis=1)
        return cls(herm(tq @ qt), affine + tq @ ql,
                   herm(corner - adj(lq) @ ql), trace[:, retained] - tq @ qz,
                   herm(state[np.ix_(retained, retained)] - adj(pz) @ qz),
                   -coupling[retained] + adj(pz) @ ql,
                   (int(np.count_nonzero(d > 0)), int(np.count_nonzero(d < 0)), 0),
                   tuple(int(i) for i in pivot))

    def reduced(self, update, threshold):
        ids = np.array(update.trace_rows, dtype=int)
        s = update.scale
        a = update.negative_scaled_inverse - s[:, None]*self.gram[np.ix_(ids, ids)]*s[None, :]
        t = s[:, None]*self.retained_trace[ids]
        top = np.block([[a, t], [adj(t), self.retained_state]])
        if threshold is None:
            return herm(top)
        e = np.vstack((s[:, None]*self.source_trace[ids], self.retained_source))
        c = self.source - threshold*np.eye(len(self.source))
        return herm(np.block([[top, e], [adj(e), c]]))

    def counts(self, update, threshold):
        return corrected_inertia(self.pivot_counts, eig_counts(self.reduced(update, threshold)),
                                 update.negative_aux_counts)
