"""Public-input-only design and nonlinear post-acquisition re-attacks."""
from pathlib import Path
import argparse
import csv
import json
import subprocess
import sys

from ports import (A, B, EDGE_I, EDGE_J, INCIDENCE, PRIOR, RX, TX, WATT, PortEngine,
                   cut_id, dump, edges, elapsed, load, now, realify, save)
from rank import bistatic_geometry, cut_power, fisher_scores
from fit import Attack, Observation, born_initial, mask_family
import numpy as np
from scipy.ndimage import label
from scipy.linalg import cho_factor, cho_solve
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"campaign.json").read_text())
OUT=HERE/"loop"
FAMILIES=("zero","patch","mean_low","mean_high","unrestricted")
IA,IB=cut_id(A),cut_id(B)


def object_info(case):
    return next(x for x in CFG["objects"] if x["case"]==case)


def old_inputs(case):
    data=load(OUT/case/"initial_A.npz")
    engine=PortEngine(object_info(case)["side"])
    old=Observation(data["tx"],data["rx"],data["Y"],float(data["radius"]))
    coord=engine.op["centers"]
    distance=np.linalg.norm(coord-np.array(object_info(case)["roi_center_m"]),axis=1)
    roi,patch=distance<=CFG["roi_radius_m"],distance<=CFG["patch_radius_m"]
    if not roi.any() or not patch.any() or label(patch.reshape(engine.side,engine.side))[1]!=1:
        raise ArithmeticError("Registered finite-area patch must be nonempty and connected")
    return engine,old,roi,patch


def fit_family(family,roi,patch):
    return mask_family(family,roi,patch,CFG["mean_low"],CFG["mean_high"],CFG["patch_floor"])


def row_fit(x,st,log,roi,patch,old=None):
    q=dict(ROI_mean=float(x[roi].mean()),ROI_min=float(x[roi].min()),
           patch_min=float(x[patch].min()),patch_mean=float(x[patch].mean()),
           numerical_in_balls=bool(max(st["ratios"])<1))
    if old is not None:
        d=x[~roi]-old[~roi]
        side=np.sqrt(len(x))
        q.update(exterior_L1_area=float(np.sum(abs(d))*(2/side)**2),
                 exterior_L2_area=float(np.linalg.norm(d)*2/side),
                 exterior_Linf=float(np.max(abs(d))),changed_exterior_cells=int(np.sum(abs(d)>1e-10)))
    q.update(log)
    return q


def bounded_nuisance_score(j, d, cut, roi, xlow, xhigh, radius, side, cap):
    # Selected material pair with one endpoint's exterior refittable. Linear model
    # uses both old and new observations; each block is divided by its actual radius.
    e_mask=~roi
    cutoff=INCIDENCE[cut].astype(bool)
    old=INCIDENCE[IA].astype(bool)
    jp=edges(j)*np.sqrt(WATT)/radius
    delta=edges(d)*np.sqrt(WATT)/radius
    old_j=realify(jp[old][:,e_mask])
    new_j=realify(jp[cutoff][:,e_mask])
    # Old data are included as their actual fixed-candidate difference, not as a null
    # projection derived with knowledge of the target.
    matrix=np.vstack((old_j,new_j))
    target=realify(delta[old]),realify(delta[cutoff])
    y=np.concatenate(target)
    gram=matrix@matrix.T
    lips=float(np.max(np.sum(abs(gram),axis=1)))
    m=np.zeros(np.sum(e_mask))
    z=np.zeros_like(m)
    lo=-xhigh[e_mask]
    hi=4.5-xhigh[e_mask]
    t=1.
    beta=cap/(2/side)
    for _ in range(CFG["selection_projected_steps"]):
        m_next=np.clip(z-matrix.T@(matrix@z+y)/max(lips,1e-15),lo,hi)
        m_next*=min(1.,beta/max(np.linalg.norm(m_next),1e-30))
        t_next=(1+np.sqrt(1+4*t*t))/2
        z=m_next+(t-1)/t_next*(m_next-m)
        m,t=m_next,t_next
    residual=matrix@m+y
    direction=residual/max(np.linalg.norm(residual),1e-15)
    dual=matrix.T@direction
    box_support=np.sum(np.maximum(-dual*lo,-dual*hi))
    ball_support=beta*np.linalg.norm(dual)
    lower=float(max(0,direction@y-min(box_support,ball_support)))
    return dict(score=float(np.linalg.norm(residual)),linear_dual_lower_numerical=lower,
            old_predicted_difference=float(np.linalg.norm(residual[:128])),
            new_predicted_difference=float(np.linalg.norm(residual[128:])),
            exterior_L2_area=float(np.linalg.norm(m)*2/side),exterior_Linf=float(np.max(abs(m))))


def design(case):
    start=now()
    path=OUT/case
    if (path/"new_observed.npz").exists():
        raise FileExistsError("This acquisition has already run. Preserve the frozen design with its readings.")
    engine,old,roi,patch=old_inputs(case)
    attack=Attack(engine,[old],CFG["fit_discrepancy_target"])
    initial=born_initial(engine,old)
    uniform=np.full(engine.n,.85)
    states={}
    fits={}
    # c02 was acquired in the first execution, before the exact principal
    # restriction and uniform unrestricted restart were added. Replay its
    # already-frozen design with the actual run-time choices.
    engine.restrict_zero_cells=case!="n32_c02"
    unrestricted_seeds=(initial,uniform) if case!="n32_c02" else (initial,)
    state_options=[attack.fit(x,fit_family("unrestricted",roi,patch),CFG["old_fit_steps"])
                   for x in unrestricted_seeds]
    st=min((s for s,f in state_options), key=lambda s:max(s["ratios"]))
    states["unrestricted"]=[st]
    fits["unrestricted"]=[dict(seed=name,**row_fit(s["x"],s,fit,roi,patch))
                          for name,(s,fit) in zip(("Born","uniform"),state_options)]
    for family in FAMILIES[:-1]:
        states[family]=[]
        fits[family]=[]
        for seed_name,seed in [("A_unrestricted",st["x"]),("uniform",uniform)]:
            s,receipt=attack.fit(seed,fit_family(family,roi,patch),CFG["old_fit_steps"])
            states[family].append(s)
            fits[family].append(dict(seed=seed_name,**row_fit(s["x"],s,receipt,roi,patch)))
    selected={name:min(seq,key=lambda s:max(s["ratios"])) for name,seq in states.items()}
    # A bounded candidate-specific distinction. Any unfit input is identified in the receipt.
    pairs=[]
    for low,high in (("zero","patch"),("mean_low","mean_high"),("zero","mean_high")):
        for a,b in zip(states[low],states[high]):
            if b["x"][roi].mean()-a["x"][roi].mean()>=CFG["rho"] and max(a["ratios"]+b["ratios"])<1:
                pairs.append((low,high,a,b))
    fallback=not pairs
    if fallback:
        pairs=[(l,h,selected[l],selected[h]) for l,h in
                (("zero","patch"),("mean_low","mean_high"),("zero","mean_high"))]
    jac=engine.jacobian(selected["unrestricted"])
    local,nu=fisher_scores(jac,engine.op["centers"],roi,old.radius)
    geom,q=bistatic_geometry(engine.op["positions"],float(engine.op["k"]))
    signal=cut_power(selected["unrestricted"]["delta"])
    contrast=np.array([cut_power(b["delta"]-a["delta"])/(2*old.radius)**2 for _,_,a,b in pairs])
    t_score=np.min(contrast,axis=0)
    random=int(np.random.default_rng(3747292026+[r["case"] for r in CFG["objects"]].index(case)).integers(len(TX)))
    candidate_cuts=dict(A_repeat=IA, A_reverse=IA, B_interleaved=IB,
         random=random, Born_geometry=int(np.argmax(geom)), Fisher=int(np.argmax(local)),
         nuisance_Fisher=int(np.argmax(nu)), signal=int(np.argmax(signal)),
         T_discrimination=int(np.argmax(t_score)))
    # Shortlist includes the leading controls and independently scored finite-amplitude
    # alternatives. The common linearization is used for ranking only.
    pair_mats=[]
    for l,h,a,b in pairs:
        pair_mats.append((l,h,a,b,engine.jacobian(b)))
    finalists=list(dict.fromkeys([*candidate_cuts.values(),*np.argsort(t_score)[-CFG["selection_finalist_count"]:][::-1].tolist()]))
    scores={}
    for candidate in finalists:
        records=[]
        for l,h,a,b,jh in pair_mats:
            infos=bounded_nuisance_score(jh,b["delta"]-a["delta"],candidate,roi,a["x"],b["x"],
                    old.radius,engine.side,CFG["external_prior_weighted_l2_m"])
            records.append(dict(low=l, high=h,**infos))
        scores[candidate]=dict(per_pair=records,score=float(min(r["score"] for r in records)),
              minimum_candidate_linear_dual=float(min(r["linear_dual_lower_numerical"] for r in records)))
    compensation=max(finalists,key=lambda cut:(scores[cut]["score"],-cut))
    candidate_cuts["compensation"]=compensation
    for key, ix in candidate_cuts.items():
        assert len(TX[ix])==len(RX[ix])==8,key
    actions={}
    for name,c in candidate_cuts.items():
        tx,rx=TX[c],RX[c]
        if name=="A_reverse":
            tx,rx=rx,tx
        actions[name]=dict(cut=c,tx=tx.tolist(),rx=rx.tolist(),
                    predicted_pairwise_D_over_2r=[float(x) for x in np.sqrt(contrast[:,c])],
                    predicted_worst_pair_D_over_2r=float(np.sqrt(t_score[c])),
                    compensation_merit=scores[c],full_class_lower=0)
    array={}
    for family,values in states.items():
        for i,s in enumerate(values):
            array[f"{family}_{i}"]=s["x"]
            array[f"{family}_{i}_delta"]=s["delta"]
        array[f"{family}_best"]=selected[family]["x"]
    save(path/"design_candidates.npz",roi=roi,patch=patch,initial=initial,t_scores=t_score,
         candidate_contrast=contrast,geometry=geom,wavenumber_modes=q,fisher=local,
         nuisance_fisher=nu,signal=signal,**array)
    dump(path/"initial_fit_receipts.json",fits)
    dump(path/"design_frozen.json",dict(case=case,side=engine.side,
        provenance=object_info(case)["provenance"],case_source="initial_A.npz and campaign registration only",
        ROI_center_m=object_info(case)["roi_center_m"],ROI_cells=int(roi.sum()), ROI_area_m2=float(roi.sum()*(2/engine.side)**2),
        connected_patch_cells=int(patch.sum()),patch_area_m2=float(patch.sum()*(2/engine.side)**2),
        candidate_pairs=[dict(low=l,high=h,A_max_resid=max(a["ratios"]+b["ratios"]),
                              mean_difference=float((b["x"][roi]-a["x"][roi]).mean())) for l,h,a,b in pairs],
        candidate_pair_fallback=fallback,cut_count=len(TX),actions=actions,shortlist=scores,
        diversity_notes="B is the uniform-angle interleaved comparison on the fixed circle. A reversal is a separate acquisition and costs eight pulses. Duplicate design cuts are deduplicated during refit.",
        scope="all 6435 for edge baselines, shortlist only for compensation. Heuristic candidate-specific ranking with no global pruning",
        calls=dict(engine.counters),**elapsed(start)))
    print(case,"design", {n:int(c) for n,c in candidate_cuts.items()}, "initial", {k:round(max(v["ratios"]),3) for k,v in selected.items()}, "fallback", fallback, flush=True)


def attack(case):
    start=now()
    path=OUT/case
    engine,old,roi,patch=old_inputs(case)
    spec=json.loads((path/"design_frozen.json").read_text())
    data=load(path/"new_observed.npz")
    candidates=load(path/"design_candidates.npz")
    registered=tuple(dict.fromkeys([*FAMILIES]))
    materials={}
    result={}
    keyed_shots={}
    for name,action in spec["actions"].items():
        tx,rx=action["tx"],action["rx"]
        shot=(tuple(tx),tuple(rx))
        if shot in keyed_shots:
            original=keyed_shots[shot]
            np.testing.assert_array_equal(data[name+"_Y"],data[original+"_Y"])
            result[name]=dict(alias_of=original)
            continue
        keyed_shots[shot]=name
        s=now()
        before=engine.counters.copy()
        new=Observation(tx,rx,data[name+"_Y"],float(data[name+"_r"]))
        model=Attack(engine,[old,new],CFG["fit_discrepancy_target"])
        records={}
        replay={}
        new_seeds={}
        for f in registered:
            xs=tuple(candidates[f"{f}_{j}"] for j in range(1 if f=="unrestricted" else 2))
            init=candidates[f"{f}_best"]
            seeds=[("same_family_A_best",init)]
            if len(xs)>1:
                seeds.append(("same_family_A_other",next((x for x in xs if not np.array_equal(x,init)),xs[1])))
            else:
                seeds.append(("original_Born",candidates["initial"]))
            runs=[]
            out=[]
            for seed, x in seeds:
                state, log=model.fit(x,fit_family(f,roi,patch),CFG["attack_steps"])
                receipt=row_fit(state["x"],state,log,roi,patch,x)
                runs.append(dict(seed=seed,**receipt))
                out.append(state)
                if max(state["ratios"])<=CFG["fit_discrepancy_target"]:
                    break
            ix=min(range(len(out)),key=lambda q:max(out[q]["ratios"]))
            st=out[ix]
            fi=next((j for j,old_x in enumerate(xs) if np.array_equal(old_x,seeds[ix][1])),0)
            replay[f]=float(np.linalg.norm(new.project(candidates[f"{f}_{fi}_delta"])-new.y)/new.radius)
            records[f]=dict(chosen_start=ix,runs=runs,best=row_fit(st["x"],st,
                 {"residual_over_each_radius":st["ratios"]},roi,patch,seeds[ix][1]))
            materials[name+"_"+f]=st["x"]
            materials[name+"_"+f+"_delta"]=st["delta"]
            materials[name+"_"+f+"_j"]=st["j"]
            new_seeds[f]=st
        compatible=[r["best"]["ROI_mean"] for r in records.values() if r["best"]["numerical_in_balls"]]
        # Re-attack from both sides at intermediate mean values. Numerical points only.
        path_records=[]
        if len(compatible)>=2 and max(compatible)-min(compatible)>=CFG["rho"]:
            f0=min((f for f in records if records[f]["best"]["numerical_in_balls"]),key=lambda f:records[f]["best"]["ROI_mean"])
            f1=max((f for f in records if records[f]["best"]["numerical_in_balls"]),key=lambda f:records[f]["best"]["ROI_mean"])
            st0,st1=new_seeds[f0],new_seeds[f1]
            left=st0["x"].copy()
            for af in (.25,.50,.75):
                mean_target=float((1-af)*min(compatible)+af*max(compatible))
                guess=left*.55+.45*((1-af)*st0["x"]+af*st1["x"])
                state, log=model.fit(guess,fit_family(mean_target,roi,patch),min(CFG["attack_steps"],18))
                if max(state["ratios"])<1:
                    left=state["x"]
                materials[f"{name}_path_{af}"]=state["x"]
                path_records.append(dict(fraction=af,target_mean=mean_target,
                      **row_fit(state["x"],state,log,roi,patch,st0["x"])))
        result[name]=dict(families=records,intermediate_mean_refits=path_records,
                 fixed_old_material_replay_new_over_radius=replay,
                 numeric_inner_ROI_range=[min(compatible),max(compatible)] if compatible else [],
                 safe_outer_ROI_range=[0.,4.5],certificate_gap="whole-box unresolved",
                 counts=dict(engine.counters-before),**elapsed(s))
        print(case,name,{f:tuple(round(v,3) for v in r["best"]["residual_over_each_radius"]) for f,r in records.items()},flush=True)
        # Incremental receipts survive interruption. No scoring data are read here.
        dump(path/"attack_results.json",dict(case=case,results=result,counts=dict(engine.counters),**elapsed(start)))
        save(path/"attacked_materials.npz",**materials)
    dump(path/"attack_results.json",dict(case=case,results=result,counts=dict(engine.counters),
          stop="equal bounded numerical acquisition/refit campaign finished; no outer class width certificate",**elapsed(start)))


def run(case):
    path=OUT/case
    if not (path/"initial_A.npz").exists():
        subprocess.run([sys.executable,"-B",str(HERE/"simulator.py"),"initial",case],check=True)
    if not (path/"design_frozen.json").exists():
        design(case)
    if not (path/"new_observed.npz").exists():
        subprocess.run([sys.executable,"-B",str(HERE/"simulator.py"),"acquire",case],check=True)
    if not (path/"attack_results.json").exists() or not json.loads((path/"attack_results.json").read_text()).get("stop"):
        attack(case)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("phase",choices=("design","attack","run"))
    p.add_argument("cases",nargs="*",choices=tuple(x["case"] for x in CFG["objects"]))
    args=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        for case in args.cases or [obj["case"] for obj in CFG["objects"]]:
            dict(design=design,attack=attack,run=run)[args.phase](case)
