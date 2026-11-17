"""Validate clean anchor ambiguities and probe disconnected polynomial intervals."""
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from pencil_experiment import (
    CENTER, HIGH, LOW, ROUNDING, SCALE, ROOT, REGIONS, PatchResponse, encode,
    host_model, noisy_components, pack, subspaces,
)


def decoded(value):
    return np.array(value["real"])+1j*np.array(value["imag"])


def independent_response(model, ker, z, A):
    full = model.forward(ker.m0+ker.U@z, independent=True)
    d = (full["y"]-ker.y0).ravel() @ A
    return d, full["resolvent_sigma_min"]


def validate_clean_anchor(row, model, ker, data, feeds):
    forward = PatchResponse(ker,feeds,REGIONS["square"])
    results = []
    for trial in row["trf"]["starts"]:
        if trial["residual_l2"] > ROUNDING:
            continue
        z = forward.material(trial["x"])
        d,singular = independent_response(model,ker,z,feeds)
        residual=np.linalg.norm(d-data)
        if residual<=ROUNDING:
            results.append(dict(q=trial["x"][0],z=z,full_residual=residual,
                                reduced_residual=trial["residual_l2"],full_state_min_singular=singular))
    results.sort(key=lambda t:t["q"])
    spread=results[-1]["q"]-results[0]["q"] if len(results)>1 else None
    return dict(layout=row["layout"],case=row["case"],candidates=results,spread=spread,
                observed_truth=row["truth"],use="floating numerical feasibility at the declared radius")


def probe_component(row, model, ker, data, feeds, initial_points):
    target = REGIONS[row["region"]]
    forward = PatchResponse(ker,feeds,target)
    observation_scale=max(np.linalg.norm(data),1e-3)
    components=row["interval"]["components"]
    attempts = []
    for interval in components:
        # Interior representatives from the saved polynomial component.
        qs=np.linspace(*interval,5)[1:-1]
        for q in qs:
            candidates=[]
            for initialization in initial_points:
                def residual(nuis):
                    return pack(forward.predict(np.r_[q,nuis])-data)/observation_scale

                def jacobian(nuis):
                    return forward.predict(np.r_[q,nuis],derivative=True)[:,1:]/observation_scale

                sol=least_squares(residual,initialization[1:],jac=jacobian,bounds=(LOW,HIGH),max_nfev=250,
                                  ftol=1e-12,xtol=1e-12,gtol=1e-12)
                x=np.r_[q,sol.x]
                candidates.append(dict(x=x,residual=np.linalg.norm(sol.fun)*observation_scale,status=sol.status))
            best=min(candidates,key=lambda c:c["residual"])
            z=forward.material(best["x"])
            d,singular=independent_response(model,ker,z,feeds)
            residual=np.linalg.norm(d-data)
            attempts.append(dict(component=interval,q=q,candidates=candidates,z_best=z,
                                 full_residual=residual,full_state_min_singular=singular,
                                 independently_feasible=bool(residual<=row["noise_radius_with_rounding"])))
    return dict(layout=row["layout"],region=row["region"],case=row["case"],observed_truth=row["truth"],
                noise_radius_with_rounding=row["noise_radius_with_rounding"],attempts=attempts)


def run():
    start=time.process_time()
    directory=ROOT/"results"
    raw=json.loads((directory/"pencil_raw.json").read_text())
    saved=np.load(directory/"pencil_observations.npz")
    clean_validation=[]
    disjoint_probes=[]
    max_endpoint_diff=0.
    interval_checks=0
    for layout in raw["design"]:
        model,ker,Q,W=host_model(layout)
        designs,_=subspaces(layout,model,ker,Q,W)
        for row in raw["observations"]:
            if row["layout"]!=layout:
                continue
            region=row["region"]
            name=f"{layout}_{region}_{row['case']:02d}_{row['arm']}"
            level=raw["noise_norms"].index(row["epsilon"])
            if region=="square" and row["arm"]=="anchor" and level==0:
                clean_validation.append(validate_clean_anchor(row,model,ker,saved[name+"_noise0"],saved[name+"_feed"]))
            if row["arm"]=="pencil":
                inv=row["interval"]
                c=designs[region][0]["coord"]/np.sqrt(raw["common_total_energy"]/row["shots"])
                recomputed=[[CENTER+SCALE*l,CENTER+SCALE*u] for l,u in noisy_components(
                    decoded(inv["g"]),c,row["noise_radius_with_rounding"])]
                old=inv["components"]
                if len(recomputed)!=len(old):
                    raise RuntimeError(f"Number of connected components changed for {name}")
                max_endpoint_diff=max(max_endpoint_diff,np.max(abs(np.array(recomputed)-np.array(old)),initial=0.))
                interval_checks+=1
                if len(old)>1:
                    starts=saved[f"{layout}_{region}_initial_parameters"]
                    disjoint_probes.append(probe_component(row,model,ker,saved[f"{name}_noise{level}"],
                                                           saved[name+"_feed"],starts))
    result=dict(
        clean_anchor_validations=clean_validation,disjoint_probes=disjoint_probes,
        interval_solver_recheck=dict(count=interval_checks,max_endpoint_difference=max_endpoint_diff),
        cpu_seconds=time.process_time()-start,
    )
    (directory/"pencil_feasible_followup.json").write_text(json.dumps(result,default=encode,indent=2)+"\n")
    accepted=[x for x in clean_validation if x["spread"] is not None and x["spread"]>1e-5]
    print("Validated distinct clean-square anchor pairs",len(accepted),"of",len(clean_validation))
    for layout in raw["design"]:
        best=max((x for x in accepted if x["layout"]==layout),key=lambda x:x["spread"])
        print(layout, "widest observed clean-square witness",best["case"],"truth",best["observed_truth"],
              [(c["q"],c["full_residual"]) for c in (best["candidates"][0],best["candidates"][-1])])
    for probe in disjoint_probes:
        print("Disconnected-set probe",probe["layout"],probe["region"],probe["case"],
              [(v["q"],v["full_residual"],v["independently_feasible"]) for v in probe["attempts"]])


if __name__=="__main__":
    run()
