"""Longer common-budget physical polishing under separate per-layout balls."""
from argparse import ArgumentParser
import json

from fit import Attack, Observation
from loop import CFG, OUT, fit_family, old_inputs, row_fit
from ports import dump, elapsed, load, now, realify, save
import numpy as np
from threadpoolctl import threadpool_limits

CASE="fresh_bg_b"
ARMS=("B_interleaved","compensation")
CAP=220


class NormConsistentAttack(Attack):
    def state(self,x):
        st=self.engine.state(x)
        blocks=[o.residual(st["delta"]) for o in self.observations]
        st["ratios"]=[float(np.linalg.norm(r)) for r in blocks]
        norms=np.asarray(st["ratios"])
        p=2.8
        # IRLS weights produce the gradient of the same p-norm merit used
        # for acceptance. Division by the common norm does not change descent.
        st["weights"]=(np.maximum(norms,1e-12)/max(norms.max(),1e-12))**((p-2)/2)
        st["objective"]=float(np.linalg.norm(norms,ord=p)**2)
        st["residual"]=np.concatenate([w*realify(r) for r,w in zip(blocks,st["weights"])])
        return st


def run(arm):
    start=now()
    engine,old,roi,patch=old_inputs(CASE)
    p=OUT/CASE
    d=load(p/"new_observed.npz")
    design=json.loads((p/"design_frozen.json").read_text())["actions"][arm]
    new=Observation(design["tx"],design["rx"],d[arm+"_Y"],float(d[arm+"_r"]))
    baseline=load(p/"attacked_materials.npz",tuple(arm+"_"+f for f in ("unrestricted","zero","mean_low","mean_high","patch")))
    old_material=load(p/"design_candidates.npz",("initial","unrestricted_best","zero_best","patch_best","mean_low_best","mean_high_best"))
    out={}
    fits={}
    pre=NormConsistentAttack(engine,[old,new],CFG["fit_discrepancy_target"])
    free=None
    for family in ("unrestricted","zero","mean_low","mean_high","patch"):
        initial=[("first_PR_refit_own_data",baseline[arm+"_"+family])]
        if family=="unrestricted":
            other=("zero","mean_low","mean_high","patch")
            options=[(f,pre.state(baseline[arm+"_"+f])) for f in other]
            f,bs=min(options,key=lambda pair:max(pair[1]["ratios"]))
            initial.append(("best_constrained_own_data_"+f,bs["x"]))
        else:
            initial.append(("new_long_GN_unrestricted",free["x"]))
        st,logs=[],[]
        for name,x in initial:
            q,log=pre.fit(x,fit_family(family,roi,patch),CAP)
            st.append(q)
            logs.append(dict(seed=name,**log))
            print(CASE,arm,family,name,q["ratios"],log["stop"],flush=True)
            if max(q["ratios"])<.82:
                break
        i=min(range(len(st)),key=lambda z:max(st[z]["ratios"]))
        q=st[i]
        if family=="unrestricted":
            free=q
        fits[family]=row_fit(q["x"],q,dict(final=q["ratios"],runs=logs),
                             roi,patch,initial[i][1])
        out[family]=q["x"]
        out[family+"_j"]=q["j"]
        out[family+"_F"]=q["delta"]
        save(p/f"long_gn_{arm}.npz",**out)
        dump(p/f"long_gn_{arm}.json",dict(case=CASE,arm=arm,shots=0,fit_families=fits,
             optimized_target="weighted l_p aggregate line search with norm-consistent IRLS",
             cap_per_start=CAP,counts=dict(engine.counters),**elapsed(start)))


if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("arm",choices=ARMS)
    with threadpool_limits(limits=1,user_api="blas"):
        run(p.parse_args().arm)
