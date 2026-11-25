"""Reproduce exact singular-state and singular-update edge cases of the physical kernel."""
import json
from pathlib import Path

from flint import arb, ctx
import numpy as np
import sympy as sp

from interval_backend import BallCache, ball, exact_inverse
from moving_trace import TraceCache, physical_update


def exact_counts(matrix):
    """Exact rational symmetric congruence, including zero diagonal 2×2 pivots."""
    a = matrix
    p = n = 0
    while a.rows:
        indices = list(range(a.rows))
        diag = next((i for i in indices if a[i,i] != 0), None)
        if diag is not None:
            piv = [diag]
            p += int(bool(a[diag,diag] > 0))
            n += int(bool(a[diag,diag] < 0))
        else:
            off = next(((i,j) for i in indices for j in range(i) if a[i,j] != 0), None)
            if off is None:
                return p,n,a.rows
            piv = list(off)
            p += 1
            n += 1
        remain = [i for i in indices if i not in piv]
        cross = a.extract(piv,remain)
        a = a.extract(remain,remain) - cross.T*a.extract(piv,piv).inv()*cross
    return p,n,0


def numeric(a):
    return np.asarray(a, dtype=float).reshape(a.shape)


def run():
    result = []
    k = sp.Matrix([[2,2,0],[2,2,0],[0,0,-1]])
    h = sp.Matrix([[0],[1],[1]])
    c = sp.Matrix([[sp.Rational(1,4)]])
    e = sp.Matrix([[0],[2]])
    t = sp.Rational(1,8)
    ex = ([.5],[1.5],[0.],[2.],[1.])
    u = physical_update([0],*ex)
    d = sp.diag(0,sp.Rational(3,4))
    for label, field_row in (('surviving_state_zero', [1,1,1]), ('lifted_state_zero', [1,-1,1])):
        trace = sp.Matrix([[0,1,0],field_row])
        k1 = k + trace.T*d*trace
        full_trace = trace.row_join(e)
        a = k.row_join(-h).col_join((-h.T).row_join(c-sp.eye(1)*t)) + full_trace.T*d*full_trace
        original = [exact_counts(k1),exact_counts(a)]
        arrays = list(map(numeric,(k,h,c,trace,e)))
        move = TraceCache.from_dense(*arrays)
        pivot = move.pivot_indices
        moved = [move.counts(u,None),move.counts(u,float(t))]
        verified = BallCache.from_base(*(ball(x) for x in arrays),pivot=pivot)
        inv = exact_inverse(u,*ex)
        proofs = [verified.proof(u,inv,value) for value in (None,arb(float(t)))]
        computed = [(p['positive_at_least'],p['negative_at_least'],p['unresolved']) for p in proofs]
        assert original==moved==computed
        if label=='surviving_state_zero':
            assert original[0][2]==1 and original[1][2]==0
        assert exact_counts(k)==(1,1,1)
        result.append(dict(case=label,exact_base_inertia=exact_counts(k),avoided_singular_state_pivot=pivot,
                           original_exact_state=original[0],original_exact_bordered=original[1],
                           moved=moved,one_sided_proofs=proofs,
                           zeros_resolved_by_exact_rational_congruence=True,
                           local_update_dimension=2,local_update_rank=1))
    path = Path(__file__).resolve().parent/'pivot_regression.json'
    path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    ctx.prec = 128
    run()
