"""Matched wall-time two-dimensional SRAP study. Offline scoring lives in score_study."""
import os
for var in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[var] = "1"
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
import platform
from pathlib import Path
import random
import socket
import time
import numpy as np
import legacy_la as la

ROOT = Path(__file__).resolve().parent
HIST = Path(os.environ.get("SRAP_HIST", ROOT.parent/"limited_aperture_routed_dr_20261120"))
DATA = Path(os.environ.get("SRAP_DATA", ROOT/"data"))
OUTPUT = Path(os.environ.get("SRAP_OUTPUT", ROOT/"results"/"predictions"))
os.environ.setdefault("SRAP_LEGACY_FILE", str(ROOT/"legacy_la.py"))
DEV_IDS = ("la180_austria_rot90_eps3p5","la120_austria_eps3p5","la180_coupled_triplet_eps4p5",
           "fresnel_twin_8GHz_half180_tx")
METHODS = ("dual","srap")


def legacy_cases():
    out=[]
    for c in la.all_specs():
        table=la.stage_table(c)
        out.append(dict(case=c["case"],object_id=c.get("geom",c.get("object")),suite=c["suite"],
                        array=c.get("aperture",c.get("mode")),side=c["side"],width=c["width"],
                        public_schedule=dict(modes=table[0],betas=table[1],rank=table[2]),cohort="legacy"))
    return out


def cases():
    path=DATA/"cases.json"
    new=json.loads(path.read_text(encoding="utf8")) if path.is_file() else []
    return legacy_cases()+new


def public_path(c):
    folder=HIST/".runs" if c["cohort"]=="legacy" else DATA
    return folder/"public"/(c["case"]+".npz")


def load_public(c):
    with np.load(public_path(c),allow_pickle=False) as z:
        return {k:z[k] for k in z.files if k!="case_json"}


def output_path(case_id,method):
    return OUTPUT/f"{case_id}__{method}.npz"


def json_default(value):
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,np.generic): return value.item()
    raise TypeError(type(value).__name__)


def job(spec):
    c,method=spec
    out=output_path(c["case"],method)
    if out.is_file(): return "cached "+out.stem
    from scipy.optimize import minimize
    if method in ("dual","u11") and c["cohort"]=="legacy":
        name={"dual":"cie_dual","u11":"u11_routed_ap"}[method]
        # Legacy method uses eps only to dispatch between its two published schedules.
        setting=dict(c,eps=3.5 if c["public_schedule"]["modes"][1]==7 else 4.)
        iterate=lambda model,Y,m,delta: la.run_cie_arm(model,Y,setting,name,m,delta)
    else:
        if method=="srap_lazy":
            from engine_dev.lazy_srap import run_cie_arm
        else:
            from engine import run_cie_arm
        iterate=lambda model,Y,m,delta: run_cie_arm(model,Y,c,method,m,delta)
    public=load_public(c)
    Y=public["Y"]
    clock=la.Clock(30.)
    clock.watch(12.,15.,24.,30.)
    m0=np.zeros(c["side"]**2)
    hist=[]
    status="SCHEDULE_COMPLETE"
    step,stage=0,-1
    model=la.PaperModel(public,c,clock)
    # Air is an exact solution of the forward problem at zero contrast.
    def audit(m,step,stage):
        val,res,pred=model.objective(m,Y)
        pt=la.snapshot(clock,m,pred,val,res,step,stage)
        hist.append(pt)
    try:
        audit(m0,0,-1)
        due=clock.wall()
        last=None
        for mat,step,stage,at_event in iterate(model,Y,m0,float(public["delta"])):
            clock.check()
            last=mat
            if clock.wall()>=due:
                audit(mat,step,stage)
                due=clock.wall()+1.
        if last is not None and (not hist or hist[-1]["step"]!=step):
            audit(last,step,stage)
    except la.TimeLimit:
        status="WALL_LIMIT"
    except BaseException as e:
        if type(e).__name__=="TimeLimit": status="WALL_LIMIT"
        elif isinstance(e,(KeyboardInterrupt,SystemExit)): raise
        else:
            import traceback
            status="FAILED "+traceback.format_exc()
    clock.capture()
    end=(clock.wall(),clock.cpu(),dict(clock.counts))
    reached={"last_returned_step":int(step),"last_returned_stage":int(stage),"raw_spent":end}
    rows=[]; arrays={}; light_costs=[]
    for T in (15,30):
        cut=0.8*T
        raw=la.best_before(hist,T)
        seed=la.best_before(hist,cut)
        for view in ("raw","light"):
            p=raw if view=="raw" else seed
            used=clock.cuts.get(T,end) if view=="raw" else clock.cuts.get(cut,end)
            wall,cpu,counts=used
            selection_status=status if view=="raw" else "NO_ELIGIBLE_SEED"
            nit=0
            if view=="light" and seed is not None:
                extra=max(0.,min(T-cut,T-wall))
                lc=la.Clock(extra)
                # Avoid inheriting a dense forward fallback built after the prefix.
                # Light starts from the selected material; all new setup is charged.
                points=[seed]
                selection_status="NO_TIME"
                try:
                    lc.check()
                    lm=la.PaperModel(public,c,lc)
                    def callback(m):
                        nonlocal nit
                        nit+=1
                        u,r,y=lm.objective(m,Y)
                        q=la.snapshot(lc,m,y,u,r,nit,stage=seed["stage"],wall_off=wall,cpu_off=cpu,base_counts=counts)
                        if lc.wall()<=extra: points.append(q)
                        lc.check()
                    fit=minimize(lambda m:lm.data_value_gradient(m,Y),seed["m"],jac=True,method="L-BFGS-B",
                                bounds=[(0.,20.)]*lm.n,callback=callback,
                                options=dict(maxiter=20,maxfun=80,maxls=20,ftol=1e-12,gtol=1e-8))
                    selection_status=str(fit.message)
                except la.TimeLimit: selection_status="WALL_LIMIT"
                except Exception as e: selection_status="FAILED "+repr(e)
                p=la.best_before(points,T)
                extra_wall,extra_cpu=lc.wall(),lc.cpu()
                light_costs.append({"budget":T,"additional_wall":extra_wall,"additional_cpu":extra_cpu,"counts":dict(lc.counts)})
                wall+=extra_wall;cpu+=extra_cpu;counts=dict(Counter(counts)+lc.counts)
            row=dict(case=c["case"],object_id=c["object_id"],array=c["array"],cohort=c["cohort"],suite=c["suite"],
                     method=method,budget=T,view=view,status=selection_status,proposal_status=status,
                     host=socket.gethostname(),wall=wall,cpu=cpu,overrun=max(0.,wall-T),counts=counts,
                     iterations_corrector=nit,last_returned_step=step,last_returned_stage=stage)
            if p is not None:
                row.update(phi=p["phi"],residual=p["residual"],selected_stage=p["stage"],selected_step=p["step"],
                           available_wall=p["available_wall"])
                arrays[f"m_{T}_{view}"]=p["m"]
            else: row.update(phi=None,residual=None,selected_stage=None,selected_step=None,available_wall=None)
            rows.append(row)
    arrays.update(traj_wall=np.array([p["available_wall"] for p in hist]),traj_phi=np.array([p["phi"] for p in hist]),
                  traj_material=np.array([p["m"] for p in hist]),traj_res=np.array([p["residual"] for p in hist]),
                  traj_stage=np.array([p["stage"] for p in hist]),traj_step=np.array([p["step"] for p in hist]))
    arrays["rows_json"]=json.dumps(rows,default=json_default)
    arrays["cost_json"]=json.dumps(dict(raw=reached,light=light_costs,host=platform.platform()),default=json_default)
    arrays["route_json"]=json.dumps(getattr(model,"trace_route",[]),default=json_default)
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_suffix(".partial")
    with tmp.open("wb") as f: np.savez_compressed(f,**arrays)
    tmp.replace(out)
    return f"{out.stem} {status.splitlines()[0]} step={step} stage={stage} wall={end[0]:.1f}"


def run_pair(tasks):
    # Opposite method orders are balanced by a fixed pre-scoring seed.
    return [job(t) for t in tasks]

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cohort",default="legacy")
    p.add_argument("--arms",default="dual,srap")
    p.add_argument("--ids",default="")
    p.add_argument("--workers",type=int,default=4)
    p.add_argument("--partition",default="0/1",help="object index modulo divisor, keeps method pairs on one host")
    a=p.parse_args()
    cs=[c for c in cases() if c["cohort"] in a.cohort.split(",")]
    ids=a.ids.split(",") if a.ids else []
    if ids: cs=[c for c in cs if c["case"] in ids]
    objects=sorted({c["object_id"] for c in cs})
    parts,ns=a.partition.split("/"); which={int(x) for x in parts.split(",")};n=int(ns)
    cs=[c for c in cs if objects.index(c["object_id"])%n in which]
    rng=random.Random(20261121)
    batches=[]
    for c in cs:
        if a.arms=="remaining_legacy":
            arms=["srap"]
            if c["case"] in DEV_IDS: arms+=["reset","full_ap","graph"]
            if c["case"].startswith("la120_"): arms+=["dual_cap","srap_cap"]
        elif a.arms=="followup_all":
            arms=["dual","srap"]+(["dual_cap","srap_cap"] if c["array"].startswith("PW120") else [])
        else: arms=a.arms.split(",")
        group=[(c,m) for m in arms if not output_path(c["case"],m).exists()]
        rng.shuffle(group)
        if group: batches.append(group)
    rng.shuffle(batches)
    tasks=[t for group in batches for t in group]
    print(f"{len(tasks)} new runs; {len(cs)} settings; {socket.gethostname()}",flush=True)
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for messages in pool.map(run_pair,batches,chunksize=1): print("\n".join(messages),flush=True)

if __name__=="__main__": main()
