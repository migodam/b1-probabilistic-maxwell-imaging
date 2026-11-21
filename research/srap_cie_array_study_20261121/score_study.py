"""Open synthetic truth and held-out fields only after every requested prediction is sealed."""
import argparse
import csv
import json
from pathlib import Path
import time
import run_study as run
import numpy as np
from scipy import ndimage
from concurrent.futures import ProcessPoolExecutor

ROOT=run.ROOT
DEV_CROSS={"dual","srap","dual_cap","srap_cap"}
DIAGNOSTIC={"reset","full_ap","graph"}
LEGACY={x["case"]:x for x in run.la.all_specs()}


def expected(c):
    if c["cohort"]=="confirmation" and c["array"].startswith("LS_"):
        return sorted(DEV_CROSS)
    if c["cohort"]=="followup" and c["array"].startswith("PW120"):
        return sorted(DEV_CROSS)
    if c["cohort"]=="legacy":
        arms={"dual","srap","u11"}
        if c["case"] in run.DEV_IDS: arms.update(DIAGNOSTIC)
        if c["case"].startswith("la120_"): arms.update(("dual_cap","srap_cap"))
        return sorted(arms)
    return sorted(DEV_CROSS if c["cohort"]=="development" and c["suite"]=="LA_SYNTH" else run.METHODS)


def shape_metrics(m,truth):
    truth=truth.ravel()
    diff=m-truth
    norm=np.linalg.norm(truth)
    peak=truth.max()
    a=(m>=.5*peak).reshape(int(np.sqrt(m.size)),-1)
    b=(truth>=.5*peak).reshape(a.shape)
    def topo(mask):
        n=ndimage.label(mask)[1]
        holes=ndimage.binary_fill_holes(mask)&~mask
        return n,ndimage.label(holes)[1]
    cn,hn=topo(a);ct,ht=topo(b)
    return dict(rel_L2=np.linalg.norm(diff)/norm,amplitude_NRMSE=np.sqrt(np.mean(diff**2))/peak,
                half_IoU=np.count_nonzero(a&b)/max(np.count_nonzero(a|b),1),
                inside_L2=np.linalg.norm(diff[truth>0])/norm,background_L2=np.linalg.norm(diff[truth==0])/norm,
                foreground_mean_ratio=float(np.mean(m[b.ravel()])/np.mean(truth[b.ravel()])),
                holes=hn,true_holes=ht,components=cn,true_components=ct,peak_contrast=m.max())


def prediction_path(c,arm):
    if c["cohort"]=="confirmation" and arm in ("dual_cap","srap_cap"):
        folder="arc_cap" if c["array"].startswith("LS_ARC_") else "line_cap"
        return ROOT/"results"/folder/f"{c['case']}__{arm}.npz"
    return run.output_path(c["case"],arm)

def get_truth(c,pub):
    if c["cohort"]=="legacy": return run.la.coarse_truth(LEGACY[c["case"]],pub["centers"])
    with np.load(run.DATA/"truth"/(c["case"]+".npz")) as z: return z["m"].ravel()


def field_scores(model,m,y):
    model.clock=run.la.Clock(600)
    f=model.ports.S@model.state(m)["j"]
    ok=np.isfinite(y)
    return np.linalg.norm((f-y)[ok])/np.linalg.norm(y[ok])


def score_case(c, arms=None, directory=None):
    start=time.perf_counter(); pub=run.load_public(c)
    truth=None;held=[]
    if c["suite"]=="LA_SYNTH": truth=get_truth(c,pub)
    else:
        d=(run.HIST/".runs" if c["cohort"]=="legacy" else run.DATA)/"heldout"/(c["case"]+".npz")
        with np.load(d) as z:
            for key,e,y in (("held_tx","E","Y"),("held_rx","E_train_hemi","Y_train_opposite_rx")):
                if e not in z or y not in z or not np.any(np.isfinite(z[y])): continue
                model=run.la.PaperModel(dict(S=z["S"],E=z[e],centers=pub["centers"],k=pub["k"],radius=pub["radius"]),c,run.la.Clock(600))
                held.append((key,model,z[y]))
    out=[]; trajectories=[]
    for arm in (arms or expected(c)):
        path=prediction_path(c,arm) if directory is None else directory/f"{c['case']}__{arm}.npz"
        with np.load(path) as z:
            for row in json.loads(str(z["rows_json"])):
                k=f"m_{row['budget']}_{row['view']}"
                if k in z:
                    m=z[k]; row.update(shape_metrics(m,truth) if truth is not None else {})
                    for key,model,y in held:
                        row[key]=field_scores(model,m,y)
                counts=row.pop("counts")
                for key in ("LU","complex_RHS","Green_forward_FFT","Green_adjoint_FFT","PR_CG_steps","feedback_route_rank","strong_receiver_rank","feedback_setup_wall","reflection_cycles","DR_setup_wall","srap_cycles","fft_pcg_iterations","H_cholesky"):
                    row[key]=counts.get(key)
                row["counts_json"]=json.dumps(counts,sort_keys=True)
                out.append(row)
            if truth is not None:
                for t,phi,res,stage,step,m in zip(z["traj_wall"],z["traj_phi"],z["traj_res"],z["traj_stage"],z["traj_step"],z["traj_material"]):
                    trajectories.append(dict(case=c["case"],method=arm,seconds=t,phi=phi,residual=res,
                                             stage=stage,step=step,**shape_metrics(m,truth)))
    print("scored",c["case"],len(out),round(time.perf_counter()-start,2),flush=True)
    return out,trajectories


def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser();p.add_argument("--cohort",required=True);p.add_argument("--workers",type=int,default=2)
    a=p.parse_args();cs=[c for c in run.cases() if c["cohort"]==a.cohort]
    if not cs:raise SystemExit("No settings for "+a.cohort)
    missing=[(c["case"],arm) for c in cs for arm in expected(c) if not prediction_path(c,arm).is_file()]
    if missing:raise SystemExit(f"{len(missing)} predictions still outstanding: {missing[:4]}")
    start=time.perf_counter()
    with ProcessPoolExecutor(max_workers=a.workers) as ex: pairs=list(ex.map(score_case,cs))
    rows=[r for a1,_ in pairs for r in a1]; hist=[r for _,b in pairs for r in b]
    write_csv(ROOT/"results"/(a.cohort+"_metrics.csv"),rows)
    write_csv(ROOT/"results"/(a.cohort+"_trajectories.csv"),hist)
    print("offline scoring elapsed",round(time.perf_counter()-start,2),flush=True)

if __name__=="__main__":main()
