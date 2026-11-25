"""Retrospective start/budget diagnosis of the one direct baseline failure."""

import json
from pathlib import Path

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from direct import geometry
from paired import CONFIRMATION, HERE, N_STARTS, RX, TX, proposals
from box_certificates import Bounds
from pool_common import FullWave, clock, elapsed, write_arrays, write_json

SEED = 20261125106
NAME = f"confirm49_{SEED}__cw4_rx16_n0125"
CONFIGS = ((42,"center"), (84,"best_sample"), (84,"center"))


def inspect(out):
    ctx.prec,ctx.threads=112,1
    g,s,b=geometry()
    with np.load(CONFIRMATION/"public"/f"{NAME}.npz") as p:
        y,radius=p["Y"],float(p["radius"])
    s,b=s[RX],b[:,TX]
    arrays, records = {}, {}
    for call_limit, initializer in CONFIGS:
        begin=clock()
        model=FullWave(g,s,b,y)
        samples,lower,upper,_=proposals(g,s,b,y,model,SEED+1007)
        screen=np.array([np.linalg.norm(model.evaluate(x)[0]-y) for x in samples]).reshape(-1,N_STARTS)
        best=screen.argmin(axis=1)
        visit=np.argsort(screen.min(axis=1),kind="stable")
        verifier=Bounds(g,s,b,y)
        audited=[]
        owner, best_residual = None, np.inf
        for i in visit:
            start=N_STARTS*i + (best[i] if initializer=="best_sample" else 0)
            candidate,predicted,fit=model.refine(samples[start],lower[i],upper[i],max_calls=call_limit)
            explicit=verifier.point_upper(candidate)
            audited.append(dict(branch=int(i),point_index=int(start),strict_upper=explicit,fit=fit))
            if explicit<best_residual:
                winner=candidate.copy(),predicted.copy()
                owner,best_residual=int(i),explicit
            if explicit<=radius:
                break
        key=f"{initializer}_{call_limit}"
        records[key]=dict(feasible=bool(best_residual<=radius),best_branch=owner,strict_ratio=best_residual/radius,
                          audit=audited,lu=dict(model.counts),**elapsed(begin))
        arrays[key+"__material"],arrays[key+"__train_prediction"]=winner
        print(key,len(audited),best_residual/radius,records[key]["wall_seconds"],flush=True)
    if not out.exists():
        out.mkdir(parents=True)
    write_arrays(out/"diagnostic_arrays.npz",**arrays)
    write_json(out/"diagnostic.json",dict(case=NAME,radius=radius,configs=records,
              protocol="Retrospective, only the one direct failure. Same all-64 rank; change initial sample or cap. "
                       "Strict endpoint audited for every tried branch (extra conservative cost)."))


if __name__=="__main__":
    with threadpool_limits(1):
        inspect(HERE/"direct_results")
