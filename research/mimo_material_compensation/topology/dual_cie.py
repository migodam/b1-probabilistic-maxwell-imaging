"""Source-block extension of the released Dual-CIE, followed by bounded physical polishing."""
from pathlib import Path
from argparse import ArgumentParser
from collections import Counter
import json
import sys

from ports import PRIOR, WATT, dump, elapsed, load, now, save
from loop import OUT, CFG, old_inputs, fit_family, row_fit
from fit import Attack, Observation

sys.path.insert(0,str(PRIOR/"imaging"))
from shared import FREEZE, NativeFourierAPCS, Regularizer, _readout

import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

BASE=FREEZE["baseline"]
ARMS=("B_interleaved","compensation")
CASE="fresh_bg_b"
CAP=2500
CAP_PHYSICAL=160
AUDIT_STRIDE=75


def energy(x):
    return np.sum(abs(x)**2,axis=0)


class CurrentGroup:
    def __init__(self,engine,obs):
        self.s=engine.c[obs.rx]
        self.sh=self.s.conj().T
        self.b=engine.b[:,obs.tx]*np.sqrt(WATT)
        self.y=obs.y
        u,v,vh=svd(self.s,full_matrices=False,check_finite=False)
        rank=int(np.count_nonzero(v>BASE["dpcs_relative_singular_cutoff"]*v[0]))
        self.vd=vh[:rank].conj().T.copy()
        self.jd=self.vd@((u[:,:rank].conj().T@self.y)/v[:rank,None])
        self.wd=1./np.maximum(energy(self.y),1e-30)
        self.ws=1./np.maximum(energy(self.b),1e-30)
        self.freq=None
        self.coeff=None
        self.rank=rank

    def stage(self,width,shape,counts,green):
        self.u=NativeFourierAPCS(shape,width,self.vd,counts)
        if self.freq is None:
            self.coeff=np.zeros((self.u.k,self.b.shape[1]),complex)
        else:
            self.coeff=self.u.transport(self.coeff,self.freq)
        self.freq=self.u.frequencies.copy()
        self.j=self.jd+self.u.apply(self.coeff)
        self.q=self.b+green.apply(self.j)
        self.d=self.s@self.j-self.y
        self.gradient=None
        self.direction=None

    def update(self,x,bj,green,counts):
        p=(bj/(1+bj*x))[:,None]
        r=p*(self.j-x[:,None]*self.q)
        v=p*r*self.ws
        g=self.u.adjoint(self.sh@(self.d*self.wd)+v-green.adjoint(x[:,None]*v))
        direction=-g
        if self.gradient is not None:
            pr=np.maximum(0., np.sum((g.conj()*(g-self.gradient)).real,axis=0)/np.maximum(energy(self.gradient),1e-30))
            direction+=self.direction*pr
            restart=(g.conj()*direction).real.sum(axis=0)>=0
            direction[:,restart]=-g[:,restart]
        jj=self.u.apply(direction)
        qq=green.apply(jj)
        dd=self.s@jj
        rr=p*(jj-x[:,None]*qq)
        denom=self.wd*energy(dd)+self.ws*energy(rr)
        a=-np.sum((g.conj()*direction).real,axis=0)/np.maximum(denom,1e-30)
        da=direction*a
        coeffs=energy(da)/np.maximum(energy(self.coeff),1e-30)
        self.coeff+=da
        self.j+=jj*a
        self.q+=qq*a
        self.d+=dd*a
        self.gradient,self.direction=g,direction
        counts["auxiliary_current_updates"]+=1
        counts["green_applies"]+=3
        return coeffs


def cie_seed(engine,observations,stage_callback):
    started=now()
    counts=Counter()
    groups=[CurrentGroup(engine,o) for o in observations]
    ws=np.concatenate([o.ws for o in groups])
    x=np.full(engine.n,BASE["initial_contrast"])
    shape=(engine.side,engine.side)
    regularizer=Regularizer(shape,smooth=1e-4,tv=1e-4,epsilon=.01)
    bounds=(np.zeros(engine.n),np.full(engine.n,4.5))
    model=Attack(engine,observations,CFG["fit_discrepancy_target"])
    best=model.state(x)
    audits=[]
    stages=[]
    for i,width in enumerate(BASE["fourier_corner_width"]):
        bj,bm=BASE["beta_current"][i],BASE["beta_material"][i]
        for group in groups:
            group.stage(width,shape,counts,engine.green)
        anchor=now()
        stop="update_budget"
        for it in range(1,CAP+1):
            coeffs=np.concatenate([o.update(x,bj,engine.green,counts) for o in groups])
            j=np.concatenate([o.j for o in groups],axis=1)
            q=np.concatenate([o.q for o in groups],axis=1)
            x,info=_readout(j,q,x,bm,ws,1.,bounds,
                     lambda v: regularizer.evaluate(v)[:2],True,
                     dict(material_inner_maxiter=BASE["material_inner_maxiter"]),counts)
            rms=float(np.sqrt(coeffs.mean()))
            stationarity=it>=BASE["minimum_updates_per_stage"] and rms<BASE["relative_coefficient_rms_stop_strictly_below"]
            if it%AUDIT_STRIDE==0 or it==CAP or stationarity:
                st=model.state(x)
                audits.append(dict(stage=i,iteration=it,aux_current_rms=rms,physical_ratios=st["ratios"],**elapsed(started)))
                if max(st["ratios"])<max(best["ratios"]):
                    best=st
                if max(st["ratios"])<.82:
                    stop="physical_discrepancy"
                    break
                if stationarity:
                    stop="published_current_stationarity"
                    break
        stages.append(dict(stage=i,width=width,beta_current=bj,beta_readout=bm,steps=it,
                           stop=stop,current_rms=rms,**elapsed(anchor)))
        stage_callback(i,it,max(best["ratios"]),stop)
        if stop=="physical_discrepancy":
            break
    return best,dict(stages=stages,full_wave_audits=audits,costs=dict(counts),
                     source_group_ranks=[o.rank for o in groups],**elapsed(started))


def run(arm):
    started=now()
    path=OUT/CASE
    data=load(path/"new_observed.npz")
    candidate=load(path/"design_candidates.npz")
    design=json.loads((path/"design_frozen.json").read_text())["actions"][arm]
    engine,old,roi,patch=old_inputs(CASE)
    new=Observation(design["tx"],design["rx"],data[arm+"_Y"],float(data[arm+"_r"]))
    model=Attack(engine,[old,new],CFG["fit_discrepancy_target"])
    best,report=cie_seed(engine,[old,new],lambda *args:print(CASE,arm,"cie",args,flush=True))
    records={}
    matrices={"cie_best":best["x"],"cie_best_j":best["j"],"cie_best_F":best["delta"]}
    first={}
    candidate_arrays=load(path/"attacked_materials.npz",tuple(arm+"_"+f for f in
                                ("unrestricted","zero","mean_low","mean_high","patch")))
    for family in ("unrestricted","zero","mean_low","mean_high","patch"):
        branches=[("CIE",best["x"]),("prior_own_data_GN",candidate_arrays[arm+"_"+family])]
        if family!="unrestricted":
            branches.append(("polished_free",first["unrestricted"]["x"]))
        fit_logs=[]
        final=[]
        for label,x in branches:
            st,fit=model.fit(x,fit_family(family,roi,patch),CAP_PHYSICAL)
            fit_logs.append(dict(seed=label,**fit))
            final.append(st)
            if max(st["ratios"])<.82:
                break
        i=min(range(len(final)),key=lambda i:max(final[i]["ratios"]))
        st=final[i]
        first[family]=st
        records[family]=row_fit(st["x"],st,dict(chosen_start=fit_logs[i]["seed"],final=st["ratios"],runs=fit_logs),
                                roi,patch,branches[i][1])
        matrices[family]=st["x"]
        matrices[family+"_j"]=st["j"]
        matrices[family+"_F"]=st["delta"]
        print(CASE,arm,family,st["ratios"],flush=True)
        save(path/f"dual_cie_{arm}.npz",**matrices)
        dump(path/f"dual_cie_{arm}.json",dict(case=CASE,arm=arm,scope="Same acquired data, one 400 MHz band, no truth access. Masked source groups reuse published R8 CIE readout and Fourier basis.",
             cie=report,cie_best_ratios=best["ratios"],families=records,
             cap_GN_each_start=CAP_PHYSICAL,total_counts=dict(engine.counters),**elapsed(started)))


if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("arm",choices=ARMS)
    args=p.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(args.arm)
