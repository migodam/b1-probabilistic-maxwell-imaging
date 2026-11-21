"""Per-array, per-material paired analysis with fixed-object bootstrap."""
import csv
import itertools
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
def load(name):
    with (ROOT/"results"/name).open() as f:return list(csv.DictReader(f))
def paired(rows,base="dual",new="srap",metric="rel_L2",items="case"):
    d={(r["method"],r[items]):r for r in rows}
    ids=sorted({k[1] for k in d})
    pairs=[(k,d.get((base,k)),d.get((new,k))) for k in ids]
    eligible=[(k,a,b) for k,a,b in pairs if a and b and a.get(metric) and b.get(metric)]
    if not eligible:return dict(n=len(ids),eligible=0)
    x=np.array([float(a[metric]) for _,a,_ in eligible]);y=np.array([float(b[metric]) for _,_,b in eligible])
    gain=(x-y)/x
    # For relative L2, log-ratio or gain medians are equivalent for odd N.
    rng=np.random.default_rng(20261121)
    boot=np.median(gain[rng.integers(len(gain),size=(10000,len(gain)))],axis=1)
    ci=np.quantile(boot,[.025,.975])
    tol=1e-10
    def med(key,col):
        return float(np.median([float((a,b)[col][key]) for _,a,b in eligible if (a,b)[col].get(key)]))
    o=dict(n=len(ids),eligible=len(x),wins=int(sum(gain>tol)),ties=int(sum(abs(gain)<=tol)),losses=int(sum(gain< -tol)),median_gain=float(np.median(gain)),ci95=ci.tolist(),baseline_median=float(np.median(x)),srap_median=float(np.median(y)))
    if metric=="rel_L2":
        o.update(dual_better_air=int(sum(x<1-tol)),srap_better_air=int(sum(y<1-tol)),
                 iou_dual=med("half_IoU",0),iou_srap=med("half_IoU",1),
                 background_dual=med("background_L2",0),background_srap=med("background_L2",1),
                 inside_dual=med("inside_L2",0),inside_srap=med("inside_L2",1),
                 rescued_past_air=[k for (k,_,_),a,b in zip(eligible,x,y) if b<1-tol and a>=1-tol],
                 lost_past_air=[k for (k,_,_),a,b in zip(eligible,x,y) if a<1-tol and b>=1-tol])
    return o

def main():
    summaries={}
    for cohort in ("confirmation","followup","development","legacy","measured_extra"):
        path=ROOT/"results"/(cohort+"_metrics.csv")
        if not path.exists(): continue
        records=load(path.name);store={}
        for (array,T,view), rows in itertools.groupby(sorted(records,key=lambda r:(r["array"],r["budget"],r["view"])),key=lambda r:(r["array"],r["budget"],r["view"])):
            rows=list(rows); metric="rel_L2" if rows[0]["suite"]=="LA_SYNTH" else "held_tx"
            entry=paired(rows,metric=metric)
            store[f"{array} {T} {view}"]=entry
        summaries[cohort]=store
    (ROOT/"results"/"paired_summary.json").write_text(json.dumps(summaries,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
