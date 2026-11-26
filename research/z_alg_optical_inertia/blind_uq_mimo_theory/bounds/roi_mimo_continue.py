"""Continue descending caps from a separately verified ceiling witness."""
import argparse
import json
import time

import numpy as np
from threadpoolctl import threadpool_limits

from roi_mimo import (OUT, ActiveModel, DiskGreen, ETAS, U, VerifiedWave, calibration_replay,
                      observation, op2, public_operator, try_box, write_arrays, write_json)


def cross_probe(case, q, reference_q):
    """Seed-only transfer from nearby verified q for a single measured case."""
    start=time.perf_counter(),time.process_time()
    base=OUT/"extended_budget"/f"n32_q{reference_q:02}"/f"{case}_continuation"
    tag,op=public_operator(32,q)
    disk=DiskGreen(32,400)
    bare=disk.dense()
    if "G" not in op:
        op["G"]=op2(disk,q,0,bare)["G"]
    with np.load(base.with_suffix(".npz")) as incoming:
        mask,lo_seed,hi_seed=incoming["roi"],incoming["low"],incoming["high"]
    check=VerifiedWave(op)
    y,r=observation(case,q,0)
    m=ActiveModel(op,y,r,tag,32)
    zero,ceil=np.zeros(m.n),np.full(m.n,U)
    low=zero
    low_cap=ceil.copy()
    low_cap[mask]=0
    pin=zero.copy()
    site=int(np.flatnonzero(mask)[np.argmax(hi_seed[mask])])
    pin[site]=U
    proposals={}
    trials={}
    for k,seed,ell,up in (("low",lo_seed,low,low_cap),("high",hi_seed,pin,ceil)):
        cand,rec=try_box(m,check,{"certified_donor":seed},ell,up,mask,160,.83)
        trials[k]=rec
        if rec["proof"]["verified_inside_noise"]:
            proposals[k]=cand
    builders={}
    for e in ETAS[1:]:
        other=op2(disk,q,e,bare)
        delta=np.linalg.norm(other["M_nominal"]-op["M_nominal"],2)
        builders[e]=other,VerifiedWave(other),float(delta)
    replay=calibration_replay(case,q,tag,builders,{k:(x,float(x[mask].max())) for k,x in proposals.items()})
    peaks=sorted(float(x[mask].max()) for x in proposals.values())
    dest=OUT/"extended_budget"/f"n32_q{q:02}"/f"{case}_warmfrom_q{reference_q:02}"
    write_arrays(dest.with_suffix(".npz"),roi=mask,**proposals)
    write_json(dest.with_suffix(".json"),dict(case=case,q=q,seed_only_donor_q=reference_q,
                  certified_outer=[0.,U],nominal_inner_peak_range=[peaks[0],peaks[-1]] if peaks else None,
                  nominal_inner_width=peaks[-1]-peaks[0] if peaks else None,
                  calibration_replay=replay,trials=trials,physical_cost=dict(m.cost),
                  wall_s=time.perf_counter()-start[0],cpu_s=time.process_time()-start[1]))
    print(case,q,peaks,{k:v["fit"]["residual_ratio"] for k,v in trials.items()},flush=True)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--case",required=True)
    a.add_argument("--q",required=True,type=int)
    a.add_argument("--from-q",type=int)
    args=a.parse_args()
    case,q=args.case,args.q
    if args.from_q:
        return cross_probe(case,q,args.from_q)
    base=OUT/"extended_budget"/f"n32_q{q:02}"/f"{case}.json"
    report=json.loads(base.read_text())
    if report["nominal_inner_peak_range"]!=[4.5,4.5]:
        raise ArithmeticError("This follow-up requires a verified ceiling")
    clock,cpu=time.perf_counter(),time.process_time()
    with np.load(base.with_suffix(".npz")) as arr:
        start,mask=arr["high"],arr["roi"]
    ident,op=public_operator(32,q)
    disk=DiskGreen(32,400)
    bare=disk.dense()
    if "G" not in op:
        op["G"]=op2(disk,q,0,bare)["G"]
    checker=VerifiedWave(op)
    y,r=observation(case,q,0)
    model=ActiveModel(op,y,r,ident,32)
    lo=np.zeros(model.n)
    limit=np.full(model.n,U)
    seed=start
    stages=[]
    best=report["nominal_endpoint_witnesses"]["high"]
    for cap in (4.25,4.,3.5,3.,2.5,2.,1.5,1.,.5,.25,0.):
        hi=limit.copy()
        hi[mask]=cap
        candidate,rec=try_box(model,checker,{"previous_certified_state":seed},lo,hi,mask,32,.83)
        stages.append(dict(cap=cap,**rec))
        if not rec["proof"]["verified_inside_noise"]:
            break
        seed,best=candidate,rec["proof"]
    nominal_cost=dict(model.cost)
    nominal_wall=time.perf_counter()-clock
    nominal_cpu=time.process_time()-cpu
    builders={}
    for e in ETAS[1:]:
        other=op2(disk,q,e,bare)
        delta=np.linalg.norm(other["M_nominal"]-op["M_nominal"],2)
        builders[e]=other,VerifiedWave(other),float(delta)
    replay=calibration_replay(case,q,ident,builders,{"low":(seed,best["peak"]),"high":(start,4.5)})
    write_arrays(base.with_name(f"{case}_continuation.npz"),roi=mask,low=seed,high=start)
    write_json(base.with_name(f"{case}_continuation.json"),
               dict(case=case,q=q,held_data_used=False,roi="roi_02",certified_outer=[0.,U],
                    nominal_inner_peak_range=[best["peak"],U],
                    nominal_inner_width=U-best["peak"],nominal_endpoint_witnesses={
                       "low":best,"high":report["nominal_endpoint_witnesses"]["high"]},
                    calibration_replay=replay,stages=stages,nominal_cost=nominal_cost,
                    nominal_wall_s=nominal_wall,nominal_cpu_s=nominal_cpu,
                    wall_s=time.perf_counter()-clock,cpu_s=time.process_time()-cpu))
    print(case,q,best["peak"],U,[r["inner_width"] for r in replay],flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
