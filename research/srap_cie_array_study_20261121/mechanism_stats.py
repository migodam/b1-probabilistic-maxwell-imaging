"""Offline test of a priori first-cycle measurements against paired material gains."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from score_study import write_csv

ROOT=Path(__file__).resolve().parent
keys=("reinjection_change_error","commutator_vs_field","aperture_leak_vs_data","x_route_fraction")
def main():
    with (ROOT/"results"/"confirmation_metrics.csv").open() as f:
        rows={(r["case"],r["method"]):r for r in csv.DictReader(f) if r["budget"]=="30" and r["view"]=="raw"}
    pairs=[]
    for (case,arm),r in rows.items():
        if arm!="srap":continue
        ordinary=rows[case,"dual"]
        with np.load(ROOT/"results"/"predictions"/(case+"__srap.npz")) as z:
            trace=json.loads(str(z["route_json"]))
            ap=next((e for e in trace if e["event"]=="AP" and e["wall"]<=30),None)
            setup=next((e for e in trace if e["event"]=="setup"),{})
        pairs.append(dict(case=case,object_id=r["object_id"],array=r["array"],host=r["host"],had_cycle=ap is not None,
                          gain=(1-float(r["rel_L2"])/float(ordinary["rel_L2"])) if r["rel_L2"] and ordinary["rel_L2"] else None,
                          rank=setup.get("weak_rank"),first_wall=None if ap is None else ap["wall"],
                          **{k:None if ap is None else ap[k] for k in keys}))
    write_csv(ROOT/"results"/"first_cycle_vs_recovery.csv",pairs)
    out=[]
    for geom in sorted({r["array"] for r in pairs}):
        group=[r for r in pairs if r["array"]==geom]
        for host in ["all"]+sorted({r["host"] for r in group}):
            q=[r for r in group if r["had_cycle"] and r["gain"] is not None and (host=="all" or host==r["host"])]
            for k in keys:
                if len(q)<3: continue
                x=[r[k] for r in q];y=[r["gain"] for r in q]
                rho=spearmanr(x,y).statistic if len(set(x))>1 and len(set(y))>1 else None
                out.append(dict(array=geom,host=host,diagnostic=k,n=len(q),spearman=rho,
                                all_objects=sum(host=="all" or a["host"]==host for a in group)))
    write_csv(ROOT/"results"/"first_cycle_associations.csv",out)
if __name__=="__main__":main()
