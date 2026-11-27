"""Try additional free exterior coordinates outside the public reference footprint."""
from attack import (HERE, ENV, BANDS, Circuit, Fit, R, RX, TX, arb, elapsed, lu_factor, lu_solve,
                    np, read, real, strict_replay, timer, write_json)
from probe import PRIOR
from scipy.fft import fft2, ifft2
from scipy.optimize import lsq_linear
import argparse
import json


def op_data(side,band):
    key=f"n{side}_f{band}_a08"+("" if band==400 else "_c048")
    with np.load(PRIOR/"ports/public/operators"/(key+".npz"),allow_pickle=False) as f:
        return {k:f[k] for k in ("spectrum","E0_bare","S_bare","Z_AA")}


def convolution(spec,side,v):
    f=fft2(v.reshape(side,side,-1),s=spec.shape,axes=(0,1),workers=1)
    return ifft2(spec[...,None]*f,axes=(0,1),workers=1)[:side,:side].reshape(v.shape)


class AllColumnJacobian:
    def __init__(self,side,band,ep):
        self.side=side
        self.e=(ep/8)**.5
        self.op=op_data(side,band)
        self.spec=self.op["spectrum"]
        self.lag=ifft2(self.spec)
        p,q=self.spec.shape
        self.reversed=self.spec[np.ix_((-np.arange(p))%p,(-np.arange(q))%q)]
        self.cost=dict(material_lu=0,forward_right_hand_sides=0,transpose_right_hand_sides=0,fft_count=0)

    def evaluate(self,x):
        ids=np.flatnonzero(x)
        n=len(ids)
        row,col=divmod(ids,self.side)
        p,q=self.spec.shape
        g=self.lag[(row[:,None]-row[None,:])%p,(col[:,None]-col[None,:])%q]
        a=np.eye(n)-x[ids,None]*g
        factors=lu_factor(a,check_finite=False)
        op=self.op
        b,c=op["E0_bare"],op["S_bare"]
        j=np.zeros_like(b)
        j[ids]=lu_solve(factors,x[ids,None]*b[ids],check_finite=False)
        nport=np.linalg.solve(op["Z_AA"]+R*np.eye(24)-1j*(c[:,ids]@j[ids]),np.eye(24))
        f=-2*R*self.e*nport[np.ix_(RX,TX)]
        right=(b+convolution(self.spec,self.side,j))@nport[:,TX]
        lower=nport[RX]@c
        u=np.zeros((len(x),8),complex)
        # Reuse the exact same bare LU for the transpose, rebuilding only
        # the receive rows of the dressed field.
        u[ids]=x[ids,None]*lu_solve(factors,lower[:,ids].T,trans=1,check_finite=False)
        left=lower+convolution(self.reversed,self.side,u).T
        d=-2j*R*self.e*(left[:,None,:]*right.T[None,:,:]).reshape(64,-1)
        self.cost["material_lu"]+=1
        self.cost["forward_right_hand_sides"]+=24
        self.cost["transpose_right_hand_sides"]+=8
        self.cost["fft_count"]+=64
        return f,np.concatenate((d.real,d.imag),axis=0)


def additional_candidates(side,x,columns,defects,valid,extent):
    # Rank exterior sites by first-order reduction of the two-ball objective.
    n=len(x)
    q=np.sum(np.stack([j.T@r for j,r in zip(columns,defects)]),axis=0)
    norm=np.sqrt(sum(np.sum(j*j,axis=0) for j in columns))
    direction=-q/np.maximum(norm,1e-13)
    active=np.ones(n,bool)
    active[valid]=False
    where=np.flatnonzero(active)
    chosen=[]
    for k in where[np.argsort(-direction[where])]:
        if direction[k]<=0.:
            break
        y,t=divmod(k,side)
        if any(abs(y-a)<=1 and abs(t-b)<=1 for a,b in chosen):
            continue
        chosen.append((y,t))
        if len(chosen)==extent:
            break
    return np.array([a*side+b for a,b in chosen],dtype=int)


def run(side,tag,family):
    prior=HERE/f"candidate_{side}_{tag}_{int(family)}.npz"
    if not prior.is_file():
        raise FileNotFoundError(prior)
    start=timer()
    _,core,_,data,labels,r,energy,paths=read(side)
    with np.load(prior,allow_pickle=False) as q:
        x=q["material"]
        support=q["original_support"]
    calc=[AllColumnJacobian(side,f,e) for f,e in zip(BANDS,energy)]
    trace=[]
    best=float("inf")
    best_x=x.copy()
    expansion_d=64 if side==96 else 160
    num_attempts=4
    for step in range(num_attempts):
        work=timer()
        outputs=[m.evaluate(x) for m in calc]
        j=[z[1]/r for z in outputs]
        residual=[real(z[0]-y)/r for z,y in zip(outputs,data[tag])]
        base=[float(np.linalg.norm(v)) for v in residual]
        if max(base)<best:
            best=max(base)
            best_x=x.copy()
        new=additional_candidates(side,x,j,residual,support,expansion_d)
        chosen=np.unique(np.r_[support,new])
        ext=chosen[~np.isin(chosen,core)]
        core_exact=x[core]
        assert np.all(core_exact==family)
        z0=x[ext]
        numerical_profile=Fit(side,chosen,core,family,data[tag],r,energy)
        init=numerical_profile.start(x,0)
        candidates=[]
        for cap in (.40,1.0):
            upper=np.minimum(4.5-z0,cap)
            lower=np.maximum(-z0,-cap)
            for bias in (1.,1.5,0.7):
                linear=np.vstack((j[0][:,ext]*bias,j[1][:,ext]))
                goal=-np.r_[residual[0]*bias,residual[1]]
                t=lsq_linear(linear,goal,bounds=(lower,upper),method="trf",
                             lsq_solver="lsmr",lsmr_tol=5e-5,lsmr_maxiter=150,
                             max_iter=14,tol=3e-4)
                proposal=np.clip(init+t.x,0,4.5)
                for rate in (1.,.5):
                    f=numerical_profile.evaluate(init+rate*(proposal-init),jac=False)
                    candidates.append((max(f["ratios"]),f))
            if min(z[0] for z in candidates)<.95:
                break
        value,current=min(candidates,key=lambda z:z[0])
        hist=dict(iteration=step,baseline=base,added_sites=new.tolist(),
                  trial_count=len(candidates),physics_cost_accepted=numerical_profile.ledger(),
                  best_numerical_residual=current["ratios"].tolist(),**elapsed(work))
        trace.append(hist)
        print(side,tag,family,step,base,current["ratios"],flush=True)
        if value>=max(base):
            break
        x=np.zeros(side**2)
        x[chosen]=current["x"]
        support=chosen
        if value<best:
            best=value
            best_x=x.copy()
        if value<.98:
            break
    audits,fields=strict_replay(side,best_x,data[tag],r,energy,core)
    is_feasible=max(v["strict_residual_over_r_upper"] for v in audits.values())<1.
    print("STRICT",side,tag,family,is_feasible,{k:v["strict_residual_over_r_upper"] for k,v in audits.items()},flush=True)
    name=f"expanded_{side}_{tag}_{int(family)}"
    write_json(name+".json",dict(label=labels[tag],sample_id=tag,side=side,core_fixed=family,
               is_opposite=(family==0.)==(labels[tag]=="high"),scope="Actual two balls, exact shared material. Additional exterior sites chosen only by data and frozen physical operators.",
               actual_sources=paths,strict_status="both balls verified" if is_feasible else "unresolved",
               independent_audit=audits,steps=trace,whole_grid_jacobian_cost={str(b):m.cost for b,m in zip(BANDS,calc)},
               **elapsed(start)))
    with (HERE/(name+".npz")).open("wb") as f:
        np.savez_compressed(f,material=best_x,roi_indices=core,
                            observed_400=data[tag][0],observed_460=data[tag][1],**fields)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("side",type=int,choices=(32,96))
    p.add_argument("label")
    p.add_argument("family",type=float,choices=(0.,2.))
    args=p.parse_args()
    run(args.side,args.label,args.family)
