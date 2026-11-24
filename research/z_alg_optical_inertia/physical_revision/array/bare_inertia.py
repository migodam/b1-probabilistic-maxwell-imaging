"""Verify the common unobserved real-sector inertia at all contrast caps."""

from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

from screen import CAPS, HERE, STUDY, times, save_json
from verify_native import FourierEnvelope, Geometry, NativeState, Weights, ball, gershgorin, high, star, submatrix


def main():
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    with np.load(STUDY/"synthetic/public/fine_geometry.npz",allow_pickle=False) as f:
        spec=f["spectrum"]
        incident=f["PW360_10x20__B"][:,:1]
    g=Geometry(spec,np.zeros((1,len(incident)),complex),incident)
    proof=FourierEnvelope(g)
    max_cap=max(CAPS)
    prop=NativeState(g,Weights.global_(max_cap,1.,0.))
    modes=prop.modes(50)
    v=ball(prop.eig_basis[:,:prop.nu])
    gram=star(v)*v
    high_cap=star(v)*prop.exact_action(v)
    report=dict(scope="I-u Re_H G on all 9216 native cells, excluding receiver and optical terms",
                precision_bits=ctx.prec,numeric_high_cap=modes,rows=[])
    for u in CAPS:
        tw,tc=perf_counter(),process_time()
        vals=1+(u/max_cap)*(prop.eigs-1.)
        neg=int(np.count_nonzero(vals<0.))
        span=list(range(neg))
        factor=arb(u)/arb(max_cap)
        projected=submatrix(high_cap,span)*factor+submatrix(gram,span)*(1-factor)
        high0=gershgorin(projected,upper=True)
        if not high0<0:
            raise ArithmeticError("Outward negative subspace is unresolved.")
        row=dict(cap=u,native_nu_at_least=neg,negative_congruence_upper=high(high0),
                 numerical_smallest_scaled_K=vals.tolist())
        try:
            row["positive_gap_proof"]=proof.verify(Weights.global_(u,1.,0.),neg)
            row["exact_nu"]=neg
        except ArithmeticError as err:
            row["positive_gap_error"]=str(err)
        row.update(times(tw,tc))
        report["rows"].append(row)
        report.update(operations=dict(g.cost),**times(start,cpu))
        save_json(HERE/"bare_inertia.json",report)
        print(u,neg,row.get("exact_nu"),row.get("positive_gap_error"),flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
