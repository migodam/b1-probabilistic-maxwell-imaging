"""Necessary support size from full multistatic data; known-support readout follow-up."""
import csv
import json
from itertools import combinations
from pathlib import Path
from time import perf_counter,process_time

from interaction import (np, LAB, OUT, CFG, public_background, historical, prepare, forward,
                         predict, shortlist, atoms, NORM, threadpool_limits)
SUP = LAB / "results/support"
DESIGN = json.loads((LAB/"configs/support.json").read_text())


def save(name,rows):
    with (SUP/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)


def one_site(bg,r,lo,hi):
    """All cells and all stationary points of the exact real Möbius arc."""
    u,o=bg.u,bg.o
    v=np.sum((o.conj().T@r)*u.conj(),axis=1)
    b=np.sum(abs(o)**2,axis=0)*np.sum(abs(u)**2,axis=1)
    k=np.diag(bg.k)
    a=v.real
    h=b+2*(k*v).real
    qa=a*abs(k)**2-h*k.real
    disc=np.maximum(h*h+4*qa*a,0)
    q=-.5*(h+np.where(h>=0,1,-1)*np.sqrt(disc))
    candidates=[np.zeros_like(a),lo-bg.m,hi-bg.m]
    with np.errstate(divide="ignore",invalid="ignore"):
        candidates += [q/qa, -a/q]
    d=np.stack(candidates,axis=1)
    legal=(d>=lo-bg.m[:,None]) & (d<=hi-bg.m[:,None]) & np.isfinite(d)
    d=np.where(legal,d,0.)
    z=d/(1-k[:,None]*d)
    # Ranking by change in objective avoids subtracting two nearly equal large residuals.
    gain=b[:,None]*abs(z)**2-2*(z.conj()*v[:,None]).real
    i,p=np.unravel_index(np.argmin(gain),gain.shape)
    p1=z[i,p]*np.outer(o[:,i],u[i])
    return int(i),float(d[i,p]),float(NORM(r-p1))


def exact_reciprocal_span(bg,r):
    """Enumerate the union of unconstrained complex-symmetric two-site port images."""
    o,u=bg.o,bg.u
    h=o.conj().T@o
    v=u.conj()@u.T
    z=o.conj().T@r@u.conj().T
    hi=np.diag(h).real
    vi=np.diag(v).real
    n=len(vi)
    source, dest=np.triu_indices(n,1)
    best=-1.
    min_ratio=1.
    all_ok=True
    for start in range(0,len(source),2048):
        i=source[start:start+2048];j=dest[start:start+2048]
        g=np.empty((len(i),3,3),dtype=complex)
        g[:,0,0]=hi[i]*vi[i];g[:,1,1]=hi[j]*vi[j];g[:,0,1]=h[i,j]*v[i,j]
        g[:,0,2]=hi[i]*v[i,j]+h[i,j]*vi[i]
        g[:,1,2]=h[j,i]*vi[j]+hi[j]*v[j,i]
        g[:,2,2]=hi[i]*vi[j]+hi[j]*vi[i]+2*(h[i,j]*v[j,i]).real
        g[:,1,0]=g[:,0,1].conj();g[:,2,0]=g[:,0,2].conj();g[:,2,1]=g[:,1,2].conj()
        b=np.stack([z[i,i],z[j,j],z[i,j]+z[j,i]],axis=1)
        scale=np.sqrt(np.diagonal(g,axis1=1,axis2=2).real)
        gn=g/(scale[:,:,None]*scale[:,None,:]);bn=b/scale
        lam=np.linalg.eigvalsh(gn)
        min_ratio=min(min_ratio,float(np.min(lam[:,0]/lam[:,-1])))
        if np.min(lam[:,0]/lam[:,-1])<1e-11:
            all_ok=False
        coef=np.linalg.pinv(gn,hermitian=True,rcond=1e-13)@bn[:,:,None]
        energy=np.sum(bn.conj()*coef[:,:,0],axis=1).real
        idx=np.argmax(energy)
        if energy[idx]>best:
            best=energy[idx];win=(i[idx],j[idx])
    a=atoms(bg,win)
    diff=NORM(r.ravel()-a@np.linalg.lstsq(a,r.ravel(),rcond=None)[0])
    return dict(pair_i=int(win[0]),pair_j=int(win[1]),pairs=len(source),span_distance=float(diff),
                recovered_gram_distance=float(np.sqrt(max(0.,NORM(r)**2-best))),
                min_normalized_gram_ratio=min_ratio,well_resolved=all_ok)


def archived():
    rows=[];exact=[]
    for ix in range(48):
        data,m0=public_background(ix)
        name=f"r13_fresh_{ix:04d}"
        for arm in DESIGN["arms"]:
            folder=historical()/("supplement_aplus48/stack_run/RESULTS/aplus_stack112" if arm=="a_plus_arc_r7" else "RESULTS/confirmation_core48")
            with np.load(folder/(name+"__"+arm)/"primary.npz") as source:
                y=source["Y_prediction"]; m=source["m"].reshape(-1)
            r=data["Y"]-y
            sv=np.linalg.svd(r,compute_uv=False)
            tails=np.r_[np.sqrt(np.cumsum(sv[::-1]**2))[::-1],0.]
            delta=float(data["delta"])
            min_sites=int(np.flatnonzero(tails<=delta)[0])
            rows.append(dict(parent=ix,arm=arm,delta=delta,residual_delta=NORM(r)/delta,
                             rank_lower=min_sites,rank2_tail_delta=tails[2]/delta))
            if arm=="dual_cie" and ix in DESIGN["exact_screen_dual"]:
                wall,cpu=perf_counter(),process_time()
                bg=prepare(data,m)
                q=exact_reciprocal_span(bg,data["Y"]-bg.y)
                exact.append(dict(parent=ix,arm=arm,**q,delta=delta,
                   span_by_delta=q["span_distance"]/delta,rank2_by_delta=tails[2]/delta,
                   reference_agreement=NORM(y-bg.y)/NORM(y),cpu_s=process_time()-cpu,wall_s=perf_counter()-wall,
                   lu=1,rhs=m.size+y.shape[0]+y.shape[1]))
    save("archived_minimum.csv",rows);save("all_pairs.csv",exact)


def conditional(bg,y,lo,hi):
    sites=shortlist(bg,y)
    trial=[]
    for pair in combinations(sites,2):
        pair=tuple(sorted(pair))
        aa=atoms(bg,pair)
        q=np.linalg.qr(aa,mode="reduced")[0]
        trial.append((-NORM(q.conj().T@y.ravel()),pair))
    best=None
    for _,ix in sorted(trial)[:CFG["final_pairs"]]:
        ix=np.asarray(ix);o=bg.o[:,ix];u=bg.u[ix];k=bg.k[np.ix_(ix,ix)]
        v=np.linalg.lstsq(o,y,rcond=None)[0]
        e=u+k@v
        numerator=np.sum(e.conj()*v,axis=1).real
        den=np.sum(abs(e)**2,axis=1)
        d=np.clip(np.divide(numerator,den,out=np.zeros_like(den),where=den>0),lo-bg.m[ix],hi-bg.m[ix])
        t=np.linalg.solve(np.eye(2)-d[:,None]*k,np.diag(d))
        res=NORM(o@t@u-y)
        if best is None or res<best[0]:
            best=(res,ix,d)
    return best


def monitor():
    rows=[]
    for ix in CFG["parents"]:
        data,m=public_background(ix);bg=prepare(data,m)
        for distance in CFG["distances"]:
            case=f"r13_{ix:04d}_{distance}"
            with np.load(OUT/"public"/(case+".npz")) as raw:
                for noise in CFG["noise"]:
                    y=raw["data_"+noise];radius=float(raw["radius_"+noise])
                    for mode in ["one_site","conditional_readout"]:
                        wall,cpu=perf_counter(),process_time()
                        if mode=="one_site":
                            i,d,res=one_site(bg,y,float(data["lower"]),float(data["upper"]))
                            idx=[i];ds=[d]
                        else:
                            res,idx,ds=conditional(bg,y,float(data["lower"]),float(data["upper"]))
                        m1=m.copy();m1[idx]+=ds
                        audited=NORM(forward(data,m1)-bg.y-y)
                        rows.append(dict(parent=ix,case=case,noise=noise,distance=distance,mode=mode,
                          i=int(idx[0]),j=int(idx[-1]) if len(idx)>1 else -1,di=float(ds[0]),
                          dj=float(ds[-1]) if len(ds)>1 else 0,delta=radius,port_error=audited,
                          relative_data=audited/NORM(y),audit_gap=abs(audited-res),
                          in_ball=audited<=radius+1e-12*NORM(bg.y),cpu_s=process_time()-cpu,
                          wall_s=perf_counter()-wall,lu=1,rhs=y.shape[1]))
    save("monitor_predictions.csv",rows)
    (SUP/"monitor_seal.json").write_text(json.dumps({"count":len(rows),"truth_access":False})+"\n")
    scored=[]
    for row in rows:
        with np.load(OUT/"scoring"/(row["case"]+".npz")) as truth:
            pair=truth["pair"]; step=truth["step"]
        target=np.zeros(576);target[pair]=step
        got=np.zeros_like(target);got[row["i"]]=row["di"]
        if row["j"]>=0:got[row["j"]]=row["dj"]
        scored.append({**row,"exact_support":int(set([row["i"],row["j"]])==set(pair)),
                       "change_error":NORM(got-target)/NORM(target)})
    save("monitor_offline.csv",scored)


if __name__=="__main__":
    SUP.mkdir(exist_ok=True)
    with threadpool_limits(1):
        archived()
        monitor()
