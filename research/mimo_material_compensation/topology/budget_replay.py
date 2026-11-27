"""Scale-tested 0.5 mean attack under a preregistered exterior L2 cap."""
from argparse import ArgumentParser
from pathlib import Path
import json

from fit import Observation,box_project
from long_gn import NormConsistentAttack
from loop import CFG,OUT,old_inputs,fit_family,row_fit
from ports import A,PortEngine,dump,elapsed,load,now,save
from verify import VerifiedWave,audit_response,audit_projection
import numpy as np
from threadpoolctl import threadpool_limits

CASE="fresh_bg_b"
ARMS=("B_interleaved","compensation")
BETA=CFG["external_prior_weighted_l2_m"]
CAP=105
STAGES=(.1,.25,.4,.500001)


def run(arm):
    start=now()
    path=OUT/CASE
    engine,old,roi,patch=old_inputs(CASE)
    spec=json.loads((path/"design_frozen.json").read_text())["actions"][arm]
    acquired=load(path/"new_observed.npz")
    x_cie=load(path/f"dual_cie_{arm}.npz",("zero",))["zero"]
    assert np.all(x_cie[roi]==0)
    new=Observation(spec["tx"],spec["rx"],acquired[arm+"_Y"],float(acquired[arm+"_r"]))
    cell=2/engine.side
    budget=.9999999*BETA/cell
    v=VerifiedWave(dict(engine.op,E0=engine.b,S=engine.c))

    def cap_model(x,bounds):
        z=box_project(x,bounds)
        e=z[~roi]-x_cie[~roi]
        z[~roi]=x_cie[~roi]+e*min(1.,budget/max(np.linalg.norm(e),1e-30))
        z[~roi]=np.clip(z[~roi],0.,4.5)
        return z

    model=NormConsistentAttack(engine,[old,new],CFG["fit_discrepancy_target"],projection=cap_model)
    center=model.state(x_cie)
    if max(center["ratios"])>=1:
        raise RuntimeError("The preregistered budget is centered on a compatible material")
    def check(st):
        x=st["x"]
        q,e,proof=audit_response(v,x,st["j"])
        ext=x[~roi]-x_cie[~roi]
        return dict(old=audit_projection(q,e,A,np.arange(8,16),old.y,old.radius),
                    new=audit_projection(q,e,new.tx,new.rx,new.y,new.radius),
                    actual_ROI_mean=float(x[roi].mean()),
                    actual_exterior_L1_area=float(np.sum(abs(ext))*cell**2),
                    actual_exterior_L2_area=float(np.linalg.norm(ext)*cell),
                    actual_exterior_Linf=float(np.max(abs(ext))),
                    model=proof)
    rows={}
    data={"center":x_cie,"center_j":center["j"]}
    ref=check(center)
    for family in ("fixed_constant_ROI","shape_free_ROI_mean"):
        x=x_cie.copy()
        trace=[]
        for t in STAGES:
            p=fit_family(float(t),roi,patch)
            if family=="fixed_constant_ROI":
                lo,hi,_,_,_=p
                lo[roi]=hi[roi]=t
                p=lo,hi,roi,0.,4.5
            st,fit=model.fit(x,p,CAP)
            x=st["x"]
            receipt=check(st)
            accepted=receipt["old"]["strict_feasible"] and receipt["new"]["strict_feasible"]
            trace.append(dict(fit=fit,verified=receipt,threshold=t,
                              satisfies_data_balls=accepted,**elapsed(start)))
            key=f"{family}_{t}"
            data[key]=x
            data[key+"_j"]=st["j"]
            print(CASE,arm,family,t,st["ratios"],receipt["actual_exterior_L2_area"],accepted,flush=True)
        rows[family]=trace
        save(path/f"budget_replay_{arm}.npz",**data)
        dump(path/f"budget_replay_{arm}.json",dict(case=CASE,arm=arm,scope="Pointwise independently verified fixed nominal acquisition. Terminated layout geometry, center and exterior cap held fixed per action. Finite points do not certify a continuous path.",
            exterior_physical_area_L2_cap=BETA,maximum_major_steps_per_target=CAP,
            center=ref,registered_final_task=.5,results=rows,total_counts=dict(engine.counters),**elapsed(start)))


if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("arm",choices=ARMS)
    with threadpool_limits(limits=1,user_api="blas"):
        run(p.parse_args().arm)
