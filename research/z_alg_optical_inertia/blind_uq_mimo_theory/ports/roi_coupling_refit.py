"""Exploratory rerun of one zero/high ROI pair on the exact 144 Ω data ball."""
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter, process_time

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from model import DiskGreen, network_terms, operator, port_point_matrices
from native_roi_path_verify import times
from native_verify import ExplicitDots
from roi_coupling import bare_certificate, candidates, circuit_models, observations, tighter_bare_transfer, trial
from sweep import tag, union_positions

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
CASE="n32_c02"
Q,SCALE=8,144
ITERS=64
TARGET_RATIO=.65
sys.path.insert(0,str(ROOT/"bounds"))
from roi_mimo import ActiveModel, box_start, gn, seed_files


def counts_since(model,start):
    return {k:v-start.get(k,0) for k,v in model.cost.items()
            if not k.startswith("max_") and v!=start.get(k,0)}


def optimize_endpoint(model,seed,lower,upper):
    before=dict(model.cost)
    started=perf_counter(),process_time()
    result,fit=gn(model,seed,lower,upper,ITERS,TARGET_RATIO)
    return result,fit|dict(paid=counts_since(model,before),**times(started))


def independent_audit(m,observation,G,E,S,mapping,net,dots):
    transfer,cost,state=bare_certificate(m,G,E,S,dots)
    primary=trial(transfer,mapping,net,observation)
    ambiguous=not (primary["train_strict_inside"] or primary["train_strict_outside"])
    result=primary|dict(bare_check=cost)
    if ambiguous:
        transfer,tight=tighter_bare_transfer(m,G,E,S,state)
        result=trial(transfer,mapping,net,observation)|dict(
            initial_bare_bound=primary,bare_check=cost,exact_Arb_defect=tight)
    for key in ("train_over_r","hold_over_r"):
        result[key][0]=max(0.,result[key][0])
    return result


def main():
    started=perf_counter(),process_time()
    ctx.prec,ctx.threads=128,1
    disk=DiskGreen(32,400)
    G=disk.dense()
    network=network_terms(disk,Q,SCALE)
    op=operator(disk,network,G)
    ident=tag(400,Q,SCALE,0)
    path=HERE/"public"/"cases"/f"{CASE}_mimo.npz"
    with np.load(path) as measured:
        Y=measured[ident+"_Y"]
        radius=float(measured[ident+"_radius"])
    model=ActiveModel(op,Y,radius,ident,32)
    fixed_mask,material_pool=candidates(CASE)
    online_mask,images=seed_files(CASE)
    if not np.array_equal(fixed_mask,online_mask) or fixed_mask.sum()!=6:
        raise ArithmeticError("ROI changed from frozen stationary second region")
    sources={f"archived_endpoint_{i+1:02}":point["sources"] for i,point in enumerate(material_pool)}
    seeds={f"archived_endpoint_{i+1:02}":point["x"] for i,point in enumerate(material_pool)}
    seeds["stationary_image"]=images["stationary_image"]
    sources["stationary_image"]=[dict(path=f"imaging/results/{CASE}/stationary/baseline.npz",key="chi")]
    if "unrestricted_public" in images:
        seeds["unrestricted_public"]=images["unrestricted_public"]
        sources["unrestricted_public"]=[dict(path=f"imaging/results/{CASE}/stationary/b1_controls.npz",key="physical_chi")]
    # One inexpensive new-data initializer. The historical point and image
    # candidates all remain available even if the Born map is a poor seed.
    receiver=op["S"]
    linear=(receiver[:,None,:]*op["E0"].T[None,:,:]).reshape(Y.size,len(fixed_mask))
    real=np.vstack((linear.real,linear.imag))
    y=Y.ravel()
    rhs=np.r_[y.real,y.imag]
    gram=real @ real.T
    ridge=max(1e-14,.025*np.max(np.diag(gram)))
    seeds["ridge_Born_screen_only"]=np.clip(real.T@np.linalg.solve(gram+ridge*np.eye(len(gram)),rhs),0,4.5)
    sources["ridge_Born_screen_only"]=[dict(analyzed_Y=ident,ridge=ridge)]
    start=perf_counter(),process_time()
    before=dict(model.cost)
    n=model.n
    lo,hi=np.zeros(n),np.full(n,4.5)
    cap=hi.copy()
    cap[fixed_mask]=0.
    x0,key0,initial0=box_start(model,seeds,lo,cap)
    screen0=dict(selected=key0,initial_ratios=initial0,paid=counts_since(model,before),**times(start))
    x_low,fit0=optimize_endpoint(model,x0,lo,cap)
    print("zero ROI",fit0["residual_ratio"],fit0["stop"],flush=True)

    # Test all six legal pin choices with only initializer forward solves.
    # After the scan there is one predeclared-budget optimization run.
    seeds["new_data_zero_ROI_fit"]=x_low
    sources["new_data_zero_ROI_fit"]=[dict(source="same data zero-ROI branch before high pin")]
    before=dict(model.cost)
    start=perf_counter(),process_time()
    pin_trials=[]
    for site in np.flatnonzero(fixed_mask):
        lower=lo.copy()
        lower[site]=4.5
        x0,key,rs=box_start(model,seeds,lower,hi)
        pin_trials.append((rs[key],int(site),key,rs,x0,lower))
    _,site,key1,_,x1,lower=min(pin_trials,key=lambda row:(row[0],row[1]))
    screen1=dict(selected=key1,pinned_cell_index=site,
                 pin_trials=[dict(site=pt[1],best_initial_ratio=pt[0],selected_initializer=pt[2],
                                  initial_ratios=pt[3]) for pt in pin_trials],
                 paid=counts_since(model,before),**times(start))
    x_high,fit1=optimize_endpoint(model,x1,lower,hi)
    print("peak 4.5 at cell",site,fit1["residual_ratio"],fit1["stop"],flush=True)
    assert np.all(x_low[fixed_mask]==0.) and x_high[site]==4.5
    assert all(np.all((x>=0)&(x<=4.5)) for x in (x_low,x_high))

    # Enter the holdout files only after both endpoints have been frozen.
    with np.load(path) as measured:
        obs=observations(measured,Q,SCALE)
    positions,ix=union_positions()
    E,S,_=port_point_matrices(disk,positions)
    net=circuit_models(disk)[Q,SCALE]
    audits={}
    with TemporaryDirectory() as tmp:
        dy=Path(tmp)/"explicit_dots.dylib"
        subprocess.run(["clang","-O3","-fno-fast-math","-ffp-contract=off","-dynamiclib",
                        "-o",str(dy),str(HERE/"certify_dot.c")],check=True)
        dots=ExplicitDots(dy,4)
        for name,m in (("low",x_low),("high",x_high)):
            audits[name]=independent_audit(m,obs,G,E,S,ix[Q],net,dots)
            print(name,"independent original circuit",audits[name]["train_over_r"],flush=True)
    before=dict(model.cost)
    st=[model.forward(x) for x in (x_low,x_high)]
    final_paid=counts_since(model,before)
    peaks=[float(np.max(x[fixed_mask])) for x in (x_low,x_high)]
    result=dict(case=CASE,roi="roi_02",region_cells=6,q=8,radiation_scale_ohms=144,frequency_MHz=400,
                total_available_W=1.0,calibration_eta=0.0,noise_radius_sqrtW=radius,
                data_npz=str(path.relative_to(ROOT)),data_Y_key=ident+"_Y",data_radius_key=ident+"_radius",
                exploratory=True,uses_only_legal_training_data=True,optimization_held_data_used=False,
                max_GN_steps_per_endpoint=ITERS,residual_target_over_r=TARGET_RATIO,
                initializer_sources=sources,low_screen=screen0,low_fit=fit0,
                high_screen=screen1,high_fit=fit1,independent_physics=audits,
                both_extreme_peaks_fit=all(v["train_strict_inside"] for v in audits.values()),
                endpoint_peaks=peaks,original_420_replay_preserved=True,
                optimization_and_final_forward_work=dict(model.cost),saving_original_forward_paid=final_paid,
                **times(started))
    out=HERE/"results"/"roi_coupling_refit"
    np.savez_compressed(out.with_suffix(".npz"),roi=fixed_mask,low=x_low,high=x_high,
                        current_low=st[0]["j"],current_high=st[1]["j"],
                        train_low=st[0]["pred"],train_high=st[1]["pred"],
                        hold_low=op["S_hold"]@st[0]["j"],hold_high=op["S_hold"]@st[1]["j"],
                        high_pinned_site=site,noise_radius=radius)
    out.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    print({k:result[k] for k in ("both_extreme_peaks_fit","endpoint_peaks","wall_seconds","CPU_seconds",
                                  "optimization_and_final_forward_work")},flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
