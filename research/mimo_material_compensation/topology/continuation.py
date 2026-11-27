"""Equal extra-budget firmware-free restart on each predeclared leading arm."""
from pathlib import Path
import argparse
import json

from ports import elapsed, load, now, save, dump
from loop import OUT,CFG,old_inputs,fit_family,row_fit
from fit import Attack,Observation
from threadpoolctl import threadpool_limits
import numpy as np

# New computation only, using each already-paid action's own new observed data.
ARMS=("B_interleaved","nuisance_Fisher","T_discrimination","compensation")
FAMILIES=("mean_low","mean_high","unrestricted")
MIX=(.15,.35,.60,.80,1.)
CAP=12


def run(case):
    start=now()
    engine,old,roi,patch=old_inputs(case)
    path=OUT/case
    des=json.loads((path/"design_frozen.json").read_text())
    d=load(path/"new_observed.npz")
    candidates=load(path/"design_candidates.npz")
    results={}
    out={}
    seen={}
    for arm in ARMS:
        s=now()
        cfg=des["actions"][arm]
        tx,rx=cfg["tx"],cfg["rx"]
        key=(tuple(tx),tuple(rx))
        if key in seen:
            results[arm]=dict(alias_of=seen[key])
            continue
        seen[key]=arm
        actual=Observation(tx,rx,d[arm+"_Y"],float(d[arm+"_r"]))
        table={}
        for family in FAMILIES:
            x=candidates[f"{family}_best"]
            if family=="unrestricted":
                # Best A-feasible class candidate is also a legal free seed.
                refs=["zero","patch","mean_low","mean_high","unrestricted"]
                x=min((candidates[f"{name}_best"] for name in refs),
                      key=lambda m:np.linalg.norm(old.residual(engine.state(m)["delta"])))
            initial=x.copy()
            old_proj=actual.project(engine.state(x)["delta"])
            r=[]
            for t in MIX:
                yt=old_proj+t*(actual.y-old_proj)
                trial=Observation(tx,rx,yt,actual.radius)
                alg=Attack(engine,[old,trial],CFG["fit_discrepancy_target"])
                st,log=alg.fit(x,fit_family(family,roi,patch),CAP)
                x=st["x"]
                r.append(dict(alpha=t,**log))
            out[f"{arm}_{family}"]=x
            out[f"{arm}_{family}_j"]=st["j"]
            out[f"{arm}_{family}_delta"]=st["delta"]
            table[family]=row_fit(x,st,dict(final=st["ratios"],steps=r),roi,patch,initial)
        results[arm]=dict(families=table,**elapsed(s))
        dump(path/"continuation_receipts.json",dict(case=case,arms=results,stages=MIX,max_GN_per_stage=CAP,**elapsed(start)))
        save(path/"continuation_materials.npz",**out)
        print(case,arm,{k:tuple(round(v,3) for v in row["final"]) for k,row in table.items()},flush=True)
    dump(path/"continuation_receipts.json",dict(case=case,arms=results,stages=MIX,max_GN_per_stage=CAP,
        total_counts=dict(engine.counters),stop="additional equal refit budget exhausted",**elapsed(start)))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("case")
    with threadpool_limits(limits=1,user_api="blas"):
        run(p.parse_args().case)
