"""Strict receipts for the two-layout same-data solver follow-up."""
from pathlib import Path
import argparse
import json

from ports import A, PortEngine, dump, elapsed, load, now, project
from verify import VerifiedWave,audit_response,audit_projection
import numpy as np
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
CASE="fresh_bg_b"
ARMS=("B_interleaved","compensation")
FAMILIES=("unrestricted","zero","mean_low","mean_high","patch")


def run(method):
    p=HERE/"loop"/CASE
    start=now()
    spec=json.loads((p/"design_frozen.json").read_text())
    old=load(p/"initial_A.npz")
    new=load(p/"new_observed.npz")
    mask=load(p/"design_candidates.npz",("roi","patch"))
    truth=load(p/"sealed_scoring.npz",("full_calibrated_delta",))["full_calibrated_delta"]
    true_hold=project(truth,A,np.arange(16,24))
    engine=PortEngine()
    verifier=VerifiedWave(dict(engine.op,S=engine.c,E0=engine.b))
    rows={}
    for name in ARMS:
        arm_start=now()
        results=json.loads((p/f"{method}_{name}.json").read_text())
        arr=load(p/f"{method}_{name}.npz")
        families={}
        for f in FAMILIES:
            x,j=arr[f],arr[f+"_j"]
            q,e,verification=audit_response(verifier,x,j)
            acq=spec["actions"][name]
            old_audit=audit_projection(q,e,A,np.arange(8,16),old["Y"],float(old["radius"]))
            new_audit=audit_projection(q,e,acq["tx"],acq["rx"],new[name+"_Y"],float(new[name+"_r"]))
            response_ind,_,bare_cost=engine.independent_response(x)
            rel=float(np.linalg.norm(response_ind-q)/np.linalg.norm(response_ind))
            yhold=project(q,A,np.arange(16,24))
            families[f]=dict(original_A=old_audit,new_layout=new_audit,
               strict_joint_feasible=old_audit["strict_feasible"] and new_audit["strict_feasible"],
               strict_single_point_rejected=max(old_audit["lower_over_r"],new_audit["lower_over_r"])>1,
               ROI_mean=float(x[mask["roi"]].mean()),patch_min=float(x[mask["patch"]].min()),
               independent_bare_and_loaded_relative_difference=rel,
               offline_held_port_error=float(np.linalg.norm(yhold-true_hold)/np.linalg.norm(true_hold)),
               physics=verification,independent_bare_cost=bare_cost)
            assert rel<1e-8
        means=sorted((row["ROI_mean"],family) for family,row in families.items() if row["strict_joint_feasible"])
        pair=None
        if len(means)>1:
            lower,low=means[0]
            upper,high=means[-1]
            d=arr[high][~mask["roi"]]-arr[low][~mask["roi"]]
            cell=2/engine.side
            pair=dict(low_family=low,high_family=high,ROI_mean_gap=float(upper-lower),
                      exterior_L1_area=float(np.sum(abs(d))*cell**2),
                      exterior_L2_area=float(np.linalg.norm(d)*cell),
                      exterior_Linf=float(np.max(abs(d))),
                      changed_exterior_cells=int(np.count_nonzero(abs(d)>1e-12)))
        rows[name]=dict(families=families,verified_ROI_points=[q[0] for q in means],extreme_verified_pair=pair,
                        verified_inner_width=pair["ROI_mean_gap"] if pair is not None else 0.,
                        study_fit_cpu_s=results["cpu_s"],study_fit_wall_s=results["wall_s"],**elapsed(arm_start))
    dump(p/f"{method}_verified.json",dict(case=CASE,method=method,arms=rows,outer_width=4.5,
               independently_replayed_ground_truth="postseal_truth_verified.json",
               scope="Each value is for one fitted nominally calibrated sample. Rejected points leave the material class unresolved.",**elapsed(start)))


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("method",choices=("dual_cie","long_gn"))
    with threadpool_limits(limits=1,user_api="blas"):
        run(a.parse_args().method)
