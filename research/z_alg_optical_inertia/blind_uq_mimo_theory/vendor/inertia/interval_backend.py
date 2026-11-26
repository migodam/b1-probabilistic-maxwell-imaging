"""128-bit one-sided verification of moving traces on the stored TM operator."""
from dataclasses import dataclass

import numpy as np
from flint import acb, acb_mat, arb
from scipy.linalg import eigh

from moving_trace import herm


def ball(a):
    a = np.asarray(a)
    return acb_mat(*a.shape, [acb(float(v.real), float(v.imag)) for v in a.flat])


def adj(a):
    return a.conjugate().transpose()


def sub(a, rows, cols):
    return acb_mat(len(rows), len(cols), [a[int(i), int(j)] for i in rows for j in cols])


def blocks(rows):
    height = sum(r[0].nrows() for r in rows)
    width = sum(x.ncols() for x in rows[0])
    entries = [a[i, j] for r in rows for i in range(r[0].nrows()) for a in r for j in range(a.ncols())]
    return acb_mat(height, width, entries)


def diagonal(values):
    out = acb_mat(len(values), len(values))
    for i, v in enumerate(values):
        out[i, i] = v
    return out


def row_scaled(a, scale):
    return acb_mat(a.nrows(), a.ncols(), [scale[i]*a[i,j] for i in range(a.nrows()) for j in range(a.ncols())])


def floor(a):
    return min((a[i,i].real - sum((abs(a[i,j]) for j in range(a.ncols()) if j != i), arb(0))).lower()
               for i in range(a.nrows()))


def lower(a):
    return float(np.nextafter(float(a.lower()), -np.inf))


def norm(a):
    return sum((z.real**2 + z.imag**2 for z in a.entries()), arb(0)).sqrt()


def verify_signs(a, approximate=None):
    """Definite projected forms prove lower sign counts via interval Gershgorin."""
    n = a.nrows()
    if n == 0:
        return {'positive_at_least': 0, 'negative_at_least': 0, 'unresolved': 0, 'floors': [None, None]}
    if approximate is None:
        approximate = np.array([complex(z) for z in a.entries()]).reshape(n, n)
    eig, vectors = eigh(herm(approximate), driver='evd', check_finite=False)
    counts, margins = [], []
    for sign in [1, -1]:
        ids = np.argsort(-sign*eig)
        ids = ids[sign*eig[ids] > 0]
        if not len(ids):
            counts.append(0)
            margins.append(None)
            continue
        v = ball(vectors[:, ids])
        restricted = sign*(adj(v)*(a*v))
        edge = floor(restricted)
        count = len(ids)
        # Check prefixes only if the smallest proposed sign is numerically difficult.
        while count and not edge > 0:
            count -= 1
            if count:
                edge = floor(sub(restricted, range(count), range(count)))
        counts.append(count)
        margins.append(lower(edge) if count else None)
    return {'positive_at_least': counts[0], 'negative_at_least': counts[1],
            'unresolved': n-sum(counts), 'floors': margins}


def pivot_ball(g, s, b, y, lower, upper, weight):
    """Original operator, independently assembled from all row endpoints."""
    l = [arb(float(v)) for v in lower]
    u = [arb(float(v)) for v in upper]
    w = [arb(float(v)) for v in weight]
    mid = [z*(a+b)/2 for z,a,b in zip(w,l,u,strict=True)]
    prod = [z*a*b for z,a,b in zip(w,l,u,strict=True)]
    mg = row_scaled(g, mid)
    k = adj(s)*s + diagonal(w) - mg - adj(mg)
    h = adj(s)*y + row_scaled(b, mid)
    c = adj(y)*y
    support = [i for i,v in enumerate(prod) if not v.is_zero()]
    if support:
        gr = sub(g, support, range(g.ncols()))
        br = sub(b, support, range(b.ncols()))
        pr = [prod[i] for i in support]
        gb = adj(gr)*row_scaled(br, pr)
        k += adj(gr)*row_scaled(gr, pr)
        h -= gb
        c += adj(br)*row_scaled(br, pr)
    return k, h, c


def trace_ball(g, b, sites):
    tt = [z for i in sites for row in ([int(j==i) for j in range(g.ncols())],
                                      [g[int(i),j] for j in range(g.ncols())]) for z in row]
    ee = [z for i in sites for row in ([0] * b.ncols(), [b[int(i),j] for j in range(b.ncols())]) for z in row]
    return acb_mat(2*len(sites), g.ncols(), tt), acb_mat(2*len(sites), b.ncols(), ee)


@dataclass(frozen=True)
class BallCache:
    gram: acb_mat
    source_trace: acb_mat
    source: acb_mat
    retained_trace: acb_mat
    retained_state: acb_mat
    retained_source: acb_mat
    pivot_proof: dict

    @classmethod
    def from_base(cls, state, h, c, trace, affine, *, numeric_pivot=None, pivot=None):
        n = state.nrows()
        pivot = np.arange(n) if pivot is None else np.asarray(pivot, dtype=int)
        retained = np.setdiff1d(np.arange(n), pivot)
        q = state if len(pivot) == n and np.all(pivot == np.arange(n)) else sub(state, pivot, pivot)
        pz = sub(state, pivot, retained)
        lq = sub(h, pivot, range(h.ncols()))
        tq = sub(trace, range(trace.nrows()), pivot)
        proof = verify_signs(q, numeric_pivot)
        if proof['unresolved']:
            raise ArithmeticError('uncertified common pivot signs')
        solved = q.solve(blocks([[lq, pz, adj(tq)]])) if len(pivot) else acb_mat(0, h.ncols()+len(retained)+trace.nrows())
        width, end = h.ncols(), h.ncols()+len(retained)
        ql, qz, qt = (sub(solved, range(len(pivot)), r) for r in
                      (range(width), range(width,end), range(end,solved.ncols())))
        return cls(tq*qt, affine+tq*ql, c-adj(lq)*ql,
                   sub(trace, range(trace.nrows()), retained)-tq*qz,
                   sub(state, retained, retained)-adj(pz)*qz,
                   -sub(h, retained, range(width))+adj(pz)*ql, proof)

    def reduced(self, update, exact_negative_inverse, threshold):
        s = [arb(float(v)) for v in update.scale]
        ids = update.trace_rows
        gram = sub(self.gram, ids, ids)
        a = exact_negative_inverse - row_scaled(adj(row_scaled(gram, s)), s)
        e = row_scaled(sub(self.retained_trace, ids, range(self.retained_trace.ncols())), s)
        top = blocks([[a,e],[adj(e),self.retained_state]])
        if threshold is None:
            return top
        c = self.source - diagonal([threshold] * self.source.nrows())
        h = blocks([[row_scaled(sub(self.source_trace, ids, range(c.ncols())), s)], [self.retained_source]])
        return blocks([[top,h],[adj(h),c]])

    def proof(self, update, exact_negative_inverse, threshold, proposal=None):
        reduced = self.reduced(update, exact_negative_inverse, threshold)
        proof = verify_signs(reduced, proposal)
        total_signs = ('positive_at_least','negative_at_least')
        size = reduced.nrows() - sum(update.negative_aux_counts)
        size += sum(self.pivot_proof[name] for name in total_signs)
        for idx, name in enumerate(total_signs):
            proof[name] = max(0,self.pivot_proof[name] + proof[name] - update.negative_aux_counts[idx])
        proof['unresolved'] = size - sum(proof[name] for name in total_signs)
        return proof


def exact_inverse(update, lower, upper, base_lower, base_upper, weight):
    """With exact binary64 endpoints, certify which physical rank case applies."""
    out = acb_mat(len(update.trace_rows), len(update.trace_rows))
    j = 0
    pos = neg = 0
    for vals in zip(lower, upper, base_lower, base_upper, weight, strict=True):
        l,u,l0,u0,w = [arb(float(v)) for v in vals]
        a = -w*((l+u)-(l0+u0))/2
        b = w*(l*u-l0*u0)
        if a.is_zero() and b.is_zero():
            continue
        scale = arb(float(update.scale[j]))**2
        if a.is_zero():
            if not (b > 0 or b < 0):
                raise ArithmeticError('unresolved rank-one update')
            out[j,j] = -scale/b
            pos += b < 0
            neg += b > 0
            j += 1
        elif a > 0 or a < 0:
            out[j,j] = scale*b/(a*a)
            out[j+1,j] = out[j,j+1] = -scale/a
            pos += 1
            neg += 1
            j += 2
        else:
            raise ArithmeticError('unresolved rank-two determinant')
    if j != out.nrows() or (pos,neg,0) != update.negative_aux_counts:
        raise ArithmeticError('exact update signs and numeric proposal differ')
    return out
