"""Paired fine-grid and laboratory scoring, with independent object as the unit of comparison."""
import csv
from collections import Counter,defaultdict
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE=Path(__file__).resolve().parent
OBJ=[f"b3m_{j:04d}" for j in range(1,13)]
NUMERIC=("L2","inside","air","IoU","train","held","wall","cpu","native_calls","coarse_calls","forward_rhs")
PAIRS=[("full","dual"),("full","diag"),("full","reset"),("pixel","dual"),("pixel","reset"),
       ("dual_then_pixel","dual_then_restart"),("pruned","matched_boxes"),("pruned","multistart"),
       ("pruned_reuse","matched_boxes"),("pruned_reuse","multistart"),("pruned_reuse","pruned")]
GOOD={"IoU":1}
def read(filename):
    return json.loads(filename.read_text())
def write(path,fields):
    path.write_text(json.dumps(fields,indent=2,allow_nan=False,ensure_ascii=False)+"\n")
def csvrows(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as out:
        wr=csv.DictWriter(out,fieldnames=keys);wr.writeheader();wr.writerows(rows)
def median(a):
    v=[x for x in a if x is not None]
    return float(np.median(v)) if v else None
def metrics(t):
    return {k:median([r[k] for r in t]) for k in NUMERIC}
def valmean(v):
    return float(np.mean(v)) if v and all(x is not None for x in v) else None

def main():
    results=read(BASE/"synthetic_scores.json")["rows"]
    rows=[]; lookup={}; cv=defaultdict(list); grouped=defaultdict(list)
    for r in results:
        name,arm=r["source"]["case"],r["source"]["arm"]
        hist=read(BASE/"predictions"/f"{name}__{arm}.json")
        oid,aperture,noise=name.split("__"); status=r["status"]
        m=r.get("material96",{})
        t=(r.get("native_fitting") or {}).get(name,{})
        held=(r.get("held_transmit_and_receive") or {})
        nr=hist["native_counts"];cr=hist["coarse_counts"]
        row=dict(case=name,object=oid,aperture=aperture,noise=noise,arm=arm,status=status,
                 L2=m.get("rel_L2"),inside=m.get("inside_L2"),air=m.get("background_L2"),
                 IoU=m.get("half_IoU"),train=t.get("relative_to_measured"),held=held.get("relative_to_measured"),
                 wall=hist["wall_seconds"],cpu=hist["cpu_seconds"],
                 coarse_calls=len(hist["coarse_audits"])-1,native_calls=nr.get("native_full_wave_calls",0),
                 forward_rhs=cr.get("complex_RHS",0)+nr.get("state_rhs",0),
                 steps=cr.get("PR_CG_steps",0),proof_wall=hist["proof_cost"]["wall"],
                 winner=hist["selected"],in_nominal_ball=hist["radius_status"],
                 other_stops=" ".join(rr["status"] for rr in hist["search_runs"]))
        rows.append(row);lookup[name,arm]=row;grouped[aperture,noise,arm].append(row)
        cv[oid,aperture,arm].append(row)
    csvrows(BASE/"per_observation.csv",rows)
    sums=[]
    for (ap,noise,arm),items in sorted(grouped.items()):
        sums.append(dict(aperture=ap,noise=noise,arm=arm,n=len(items),**metrics(items)))
    csvrows(BASE/"by_aperture_noise.csv",sums)
    objects={}
    for (oid,ap,arm),obs in cv.items():
        if len(obs)!=2:raise ValueError(f"Missing paired noise for {oid} {ap} {arm}.")
        objects[oid,ap,arm]={k:valmean([r[k] for r in obs]) for k in NUMERIC}
    apnames=sorted({r["aperture"] for r in rows})
    armlist=list(dict.fromkeys(r["arm"] for r in rows))
    combined=[]
    for ap in apnames:
        for arm in armlist:
            v=[objects[o,ap,arm] for o in OBJ]
            combined.append(dict(aperture=ap,arm=arm,n=len(v),**metrics(v)))
    csvrows(BASE/"independent_materials.csv",combined)
    comparisons=[];detailed=[]
    for ap in apnames:
        for arm,ref in PAIRS:
            dd=[]
            for oid in OBJ:
                a,b=objects[oid,ap,arm],objects[oid,ap,ref]
                d={k:a[k]-b[k] if a[k] is not None and b[k] is not None else None for k in NUMERIC}
                detailed.append(dict(object=oid,aperture=ap,arm=arm,reference=ref,**d));dd.append(d)
            summary=dict(aperture=ap,arm=arm,reference=ref,n_objects=len(dd),
                         both_inside_air_better=sum(v["inside"] < -1e-7 and v["air"] < -1e-7
                         for v in dd if v["inside"] is not None and v["air"] is not None))
            for k in NUMERIC:
                values=[w[k] for w in dd if w[k] is not None]
                sign=GOOD.get(k,-1)
                summary[k]=median(values)
                summary[k+"_better"]=sum(sign*x>1e-7 for x in values)
                summary[k+"_tie"]=sum(abs(x)<=1e-7 for x in values)
                summary[k+"_worse"]=sum(sign*x < -1e-7 for x in values)
            comparisons.append(summary)
    csvrows(BASE/"paired_differences.csv",comparisons);csvrows(BASE/"object_pairs.csv",detailed)
    counts={}
    for arm in armlist:
        items=[r for r in rows if r["arm"]==arm]
        count=Counter()
        for r in items:
            h=read(BASE/"predictions"/f"{r['case']}__{arm}.json")
            for cert in h["certificates"]:
                count[f"{cert['branch']}/{cert['decision']}"]+=1
            count["air_winner"]+=h["selected"]=="air"
        counts[arm]=dict(objects=12,conditions=len(items),native_in_ball=sum(
            r["in_nominal_ball"]=="numerical_in_ball" for r in items),sum_cpu=sum(r["cpu"] for r in items),
            sum_individual_wall=sum(r["wall"] for r in items),
            search_plus_physical_RHS=sum(r["forward_rhs"] for r in items),
            native_attempts=sum(r["native_calls"] for r in items),coarse_completed_calls=sum(r["coarse_calls"] for r in items),
            branch_counts=dict(count))
    write(BASE/"cohort_cost_and_exclusions.json",counts)
    labrows=[]
    for r in read(BASE/"chamber_scores.json")["rows"]:
        name,arm=r["source"]["case"],r["source"]["arm"]
        h=read(BASE/"predictions"/f"{name}__{arm}.json")
        t=r.get("train") or {};held=r.get("held_transmit") or {};nom=r.get("nominal_reference") or {}
        labrows.append(dict(observation=name,arm=arm,status=r["status"],
            held=held.get("relative_to_measured"),train=t.get("relative_to_measured"),
            wall=h["wall_seconds"],cpu=h["cpu_seconds"],**{k:v for k,v in nom.items() if isinstance(v,(float,int))}))
    csvrows(BASE/"laboratory_all.csv",labrows)
    visualize(detailed)
    print("Scored",len(rows),"synthetic predictions and",len(labrows),"chamber predictions.")

def visualize(differences):
    keep=[("full","dual","Cross / Dual"),("full","reset","Cross / reset"),
          ("dual_then_pixel","dual_then_restart","Pixel portfolio / restart"),
          ("pruned","multistart","Cold pruning / multistart"),
          ("pruned_reuse","multistart","Shared geometry / multistart")]
    fig,axs=plt.subplots(2,3,figsize=(11.5,5.8),layout="constrained",sharey=True)
    for j,ap in enumerate(("pw180_10x20","pw360_10x20")):
        for i,(key,title) in enumerate((("inside","Interior error"),("air","Air-background error"),("held","Unseen port error"))):
            ax=axs[j,i];ax.axvline(0,color="#565c61",lw=1)
            for a,(arm,ref,label) in enumerate(keep):
                v=[d[key] for d in differences if (d["aperture"],d["arm"],d["reference"])==(ap,arm,ref) and d[key] is not None]
                if not v:continue
                color="#145462" if a<3 else "#a34f2d"
                ax.scatter(v,a+np.linspace(-.14,.14,len(v)),s=10,alpha=.54,color=color,lw=0)
                ax.scatter(np.median(v),a,marker="D",s=26,color=color)
            ax.set_yticks(range(len(keep)),[f[2] for f in keep]);ax.grid(axis="x",alpha=.15)
            if j==0:ax.set_title(title)
            if i==0:ax.set_ylabel("180° aperture" if j==0 else "360° aperture")
            if j==1:ax.set_xlabel("Paired change in error")
    axs[0,0].invert_yaxis()
    for suffix in ("png","pdf"):
        fig.savefig(BASE/f"paired_materials.{suffix}",dpi=190,facecolor="white")
    plt.close(fig)
if __name__=="__main__":main()
