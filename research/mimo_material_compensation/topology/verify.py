"""Strict fixed-operator witness replay and topology upper receipts."""
from pathlib import Path
import argparse
import json
import sys

from ports import (A, B, INCIDENCE, PRIOR, RX, TX, WATT, PortEngine, cut_id, dump, edges, elapsed,
                   load, now, project, save)
import numpy as np
from threadpoolctl import threadpool_limits

sys.path.insert(0,str(PRIOR/"diagnostics"))
from verified_wave import VerifiedWave, fro, rd, up

HERE=Path(__file__).resolve().parent
A_CUT=cut_id(A)


def audit_response(verifier,x,j):
    before=now()
    inv,receipt=verifier.inverse(x)
    response,response_error,curr_error=audit_matrix(verifier,x,j,verifier.B,verifier.S,inv)
    return response,response_error,dict(response_Frobenius_error_upper=response_error,
                              current_Frobenius_error_upper=curr_error,inverse=receipt,**elapsed(before))


def audit_matrix(verifier,x,j,excitation,readout,inverse_bound):
    rhs=x[:,None]*excitation
    ie=up(rd(3)*float(np.max(x))*fro(excitation))
    curr_error=verifier.solve_error(x,j,rhs,ie,inverse_bound)
    response=readout@j
    response_error=up(rd(len(x))*fro(readout)*fro(j)+fro(readout)*curr_error)
    return response,response_error,curr_error


def audit_published_a(v,engine,x,j,inv,observed,radius):
    p,e,_=audit_matrix(v,x,np.sqrt(WATT)*j[:,:8],engine.op["E0"],engine.op["S"],inv)
    return observed_receipt(p,e,observed,radius)


def audit_projection(model,e,tx,rx,observed,radius):
    p=project(model,tx,rx)
    return observed_receipt(p,up(e*np.sqrt(WATT)+rd(2)*fro(p)),observed,radius)


def two_norm_bounds(a,b,error):
    delta=a-b
    rounding=up(rd(2)*(fro(a)+fro(b)))
    e=up(error+rounding)
    lower=max(0.,np.linalg.norm(delta)/(1+rd(2*delta.size))-e)
    return float(np.nextafter(lower,0.)),up(fro(delta)+e)


def observed_receipt(a,e,observed,radius):
    low,upper=two_norm_bounds(a,observed,e)
    return dict(lower_over_r=float(np.nextafter(low/radius,0.)),upper_over_r=up(upper/radius),
                strict_feasible=bool(up(upper/radius)<1))


def fixed():
    stage=now()
    engine=PortEngine()
    verifier=VerifiedWave(dict(engine.op,S=engine.c,E0=engine.b))
    source=json.loads((HERE/"fixed_pair"/"summary.json").read_text())
    summary=[]
    upper_all=[]
    for case_row in source["cases"]:
        c=case_row["case"]
        q=load(HERE/"fixed_pair"/f"{c}.npz")
        radius=case_row["radius"]
        ps=[]
        signatures={}
        for name,x in zip(("zero","patch"),(q["original_zero"],q["original_patch"])):
            s=engine.state(x)
            model,e,r=audit_response(verifier,x,s["j"])
            ab=audit_published_a(verifier,engine,x,s["j"],r["inverse"]["inverse_norm_upper"],
                                 q["original_noisy_A"],radius)
            ps.append((model,e))
            signatures[name]=dict(physics=r,original_A=ab,real_box=bool(x.min()>=0 and x.max()<=4.5),
                                 ROI_mean=float(x[q["roi"]].mean()))
        z0,z1=(m for m,e in ps)
        diff=z1-z0
        fixed_raw=np.sum(abs(diff[RX[:,:,None],TX[:,None,:]])**2,axis=(1,2))*WATT
        e=up(np.sqrt(WATT)*(sum(n for m,n in ps)+rd(2)*(fro(z0)+fro(z1))))
        upper=np.nextafter((np.sqrt(fixed_raw)*(1+rd(256))+e)/(2*radius),np.inf)
        lower=np.zeros(len(upper))
        assert all(sig["original_A"]["strict_feasible"] and sig["real_box"] for sig in signatures.values())
        assert signatures["patch"]["ROI_mean"]-signatures["zero"]["ROI_mean"]>=.5
        upper_all.append(upper)
        summary.append(dict(case=c,task_threshold_roi_mean=.5,strict=signatures,boundary=dict(lower=0.,
                   incumbent_cut=case_row["best_cut"], certified_upper_max=float(upper.max()),
                   finite_action_upper_minus_incumbent_lower=float(upper.max()),cuts_pruned=0)))
    save(HERE/"fixed_pair"/"class_bounds.npz",lower=np.vstack([lower for _ in upper_all]),
         upper=np.vstack(upper_all),case_ids=np.asarray([r["case"] for r in summary]),cut_ids=np.arange(len(TX)))
    dump(HERE/"fixed_pair"/"verified_upper_receipts.json",dict(rows=summary,
         scope="Fixed calibrated stored discrete model. Each verified pair is compatible with original A and has at least 0.5 ROI mean gap. For all 6435 new layout actions it supplies an upper bound on the infimum over old-data-compatible legal pairs. Lower bound 0 retains full class; no novelty or optimality claim.",
         calls=dict(engine.counters),**elapsed(stage)))


def loop_audit(case):
    stage=now()
    path=HERE/"loop"/case
    from loop import CFG
    engine=PortEngine()
    v=VerifiedWave(dict(engine.op,S=engine.c,E0=engine.b))
    isolated=load(path/"new_observed.npz")
    old=load(path/"initial_A.npz")
    attack=json.loads((path/"attack_results.json").read_text())
    if "stop" not in attack:
        raise RuntimeError("Await the full attack receipts")
    design=json.loads((path/"design_frozen.json").read_text())
    mat=load(path/"attacked_materials.npz")
    roi=load(path/"design_candidates.npz",["roi"])["roi"]
    dedup=[]
    # Use a fixed cap of physical audits, prioritizing endpoints that maximize
    # observed interior task width. This ordering uses only candidate data.
    for name,row in attack["results"].items():
        if "alias_of" in row:
            continue
        feas={family: q["best"] for family,q in row["families"].items() if q["best"]["numerical_in_balls"]}
        if not feas:
            continue
        items=sorted(feas,key=lambda f:feas[f]["ROI_mean"])
        selected=[items[0],items[-1]] if len(items)>1 else [items[0]]
        gap=float(feas[items[-1]]["ROI_mean"]-feas[items[0]]["ROI_mean"])
        dedup.append((gap,name,selected))
    dedup.sort(key=lambda a:(-a[0],a[1]))
    outcomes={}
    cap=CFG["final_audit_budget_feasible_per_case"]
    candidates_audited=0
    for width,name,selected in dedup:
        if cap<len(selected):
            continue
        tx,rx=design["actions"][name]["tx"],design["actions"][name]["rx"]
        checks={}
        for f in selected:
            cap-=1
            candidates_audited+=1
            key=f"{name}_{f}"
            x,j=mat[key],mat[key+"_j"]
            try:
                q,e,receipt=audit_response(v,x,j)
            except ArithmeticError as err:
                checks[f]=dict(verified=False,reason=str(err))
                continue
            if case.startswith("n32_"):
                old_final=audit_published_a(v,engine,x,j,receipt["inverse"]["inverse_norm_upper"],
                                            old["Y"],float(old["radius"]))
            else:
                old_final=audit_projection(q,e,A,np.arange(8,16),old["Y"],float(old["radius"]))
            new_final=audit_projection(q,e,tx,rx,isolated[name+"_Y"],float(isolated[name+"_r"]))
            target=np.sum(x[roi])/roi.sum()
            error=up(rd(2*roi.sum())*np.sum(abs(x[roi]))/roi.sum())
            checks[f]=dict(verified=old_final["strict_feasible"] and new_final["strict_feasible"],
                     ROI_mean=float(target),ROI_mean_lower=float(np.nextafter(target-error,-np.inf)),
                     ROI_mean_upper=float(np.nextafter(target+error,np.inf)),
                     original_A=old_final,new=new_final,physical=receipt)
        verified=[q for q in checks.values() if q["verified"]]
        inner=sorted(q["ROI_mean"] for q in verified)
        width_lower=(max(q["ROI_mean_lower"] for q in verified)-min(q["ROI_mean_upper"] for q in verified)
                     if len(verified)>1 else 0.)
        outcomes[name]=dict(strict_candidates=checks,inner_ROI_means=inner,
                            ROI_width_lower=max(0.,width_lower),ROI_width_outer=4.5)
    dump(path/"strict_attacks.json",dict(case=case,arms=outcomes,physical_candidates_audited=candidates_audited,
        actual_candidate_count_budget=CFG["final_audit_budget_feasible_per_case"],
        scope="Strict nominal-calibration point residuals; unverified numerical material remains labeled. No full-class width bound beyond the material box.",
        **elapsed(stage)))
    print(case,{n:r["ROI_width_lower"] for n,r in outcomes.items()},flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=("fixed","loop"))
    p.add_argument("cases",nargs="*")
    a=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        if a.phase=="fixed":
            fixed()
        else:
            for case in a.cases:
                loop_audit(case)
