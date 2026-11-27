"""Object-clustered uncertainty for paired, previously specified change-detection arms."""
from pathlib import Path
from itertools import product
import csv,json
import numpy as np

LAB=Path(__file__).resolve().parents[1];RES=LAB/"results"


def load(path):
    with path.open() as f:return list(csv.DictReader(f))


def compare(groups,axis,arm1,arm2,condition,seed):
    rows=[]
    for parent in sorted(set(x["parent"] for x in groups)):
        sub=[x for x in groups if x["parent"]==parent and x[axis]==condition]
        a={x["case"]:x for x in sub if (x.get("method") or x.get("mode"))==arm1}
        b={x["case"]:x for x in sub if (x.get("method") or x.get("mode"))==arm2}
        if set(a)!=set(b):raise ValueError(parent)
        win=[];err=[]
        for case in a:
            one,other=a[case],b[case]
            win.append(int(one["exact_support"])-int(other["exact_support"]))
            err.append(float(one.get("change_error") or one["change_relative_l2"])-
                       float(other.get("change_error") or other["change_relative_l2"]))
        rows.append(dict(parent=parent,pairs=len(win),difference_exact=np.mean(win),difference_L2=np.mean(err)))
    values=np.array([[x["difference_exact"],x["difference_L2"]] for x in rows]);n=len(rows)
    rng=np.random.default_rng(seed);index=rng.integers(n,size=(10000,n));means=values[index].mean(axis=1)
    return dict(first=arm1,second=arm2,condition=condition,parents=n,mean=values.mean(axis=0).tolist(),
                cluster_percentiles_95=np.percentile(means,[2.5,97.5],axis=0).tolist(),by_parent=rows)


if __name__=="__main__":
    main=load(RES/"interaction/offline.csv")
    main += [x for x in load(RES/"support/monitor_offline.csv") if x["mode"]=="conditional_readout"]
    other=load(RES/"cross_changes/offline.csv")
    result=[]
    for cond in ["clean","ten_percent_change","ten_percent_total"]:
        for first,second in [("interaction","diagonal"),("conditional_readout","interaction"),("interaction","music")]:
            result.append(compare(main,"noise",first,second,cond,81278))
    for cond in ["matched_clean","matched_noisy","cross_clean","cross_noisy"]:
        result.append(compare(other,"condition","interaction_closed","music_closed",cond,81278))
    (RES/"interaction/paired_object_statistics.json").write_text(json.dumps(dict(
        endpoints=["difference in exact-support fraction","difference in relative change L2"],seed=81278,
        interpretation="retrospective conditional cohort; paired object bootstrap, no new blind whole-image claim",
        comparisons=result),indent=2)+"\n")
