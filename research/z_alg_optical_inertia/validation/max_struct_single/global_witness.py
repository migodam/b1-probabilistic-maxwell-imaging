"""Test simultaneous global coherent-source feasibility within the inherited pooled radius."""
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh, eigvalsh, qr, svd
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
BETAS=(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,3e-5,3e-4)
CAP,T=5,10000


def h(a):
    return (a+a.conj().T)/2


def proposals(G,S,B,Y,radii):
    p=(G-G.conj().T)/(2j)
    u,sv,v=svd(S,full_matrices=False)
    null=qr(S.conj().T,mode='full')[0][:,len(S):]
    hr=h(null.conj().T@p@null)
    e,rot=eigh(hr); null=null@rot
    dark_mat=h(null.conj().T @ (np.eye(len(G))-CAP*h(G)+T*p) @ null)
    value,vec=eigh(dark_mat,subset_by_index=(0,0))
    # Try the most sector-negative mode together with a small amount of the
    # most radiating receiver-invisible mode to obtain a robust optical sign.
    dark_coefs=vec[:,0]
    dark_candidates=[]
    for a in [0.,.05,.1,.2,.3,.4]:
        d=dark_coefs.copy();d[-1]+=a
        dark_candidates.append(null@d/np.linalg.norm(d))
    o_b=null.conj().T@(-.5j*B); o_g=null.conj().T@p
    outputs=[]
    for y,radius in zip(Y,radii):
        rhs0=u.conj().T@y
        best=None
        for rank in range(1,len(sv)+1):
            x0=v[:rank].conj().T@(rhs0[:rank]/sv[:rank,None])
            resid=float(np.linalg.norm(S@x0-y,2))
            if resid >= radius:continue
            rhs=o_g@x0+o_b
            for beta in BETAS:
                x=x0-null@(rhs/(e[:,None]+beta))
                je=x.conj().T@(B+G@x)
                optical=-h(je/1j); combined=T*optical-h(x.conj().T@x-CAP*je)
                margins=[eigvalsh(k,check_finite=False) for k in (optical,combined)]
                score=min((radius-resid)/radius,*[float(m[0]/max(1,abs(m[-1]))) for m in margins])
                if best is None or score>best[0]:
                    best=(score,x,rank,beta,resid,margins[0][0],margins[1][0])
        outputs.append(best)
    return dark_candidates,outputs


if __name__=='__main__':
    start=perf_counter(),process_time()
    with np.load(HERE.parent/'max_struct/public_inputs.npz',allow_pickle=False) as z:
        G,S,B,Y,radii,cases=[z[k] for k in ('G','S','E0','Y','delta','cases')]
    with threadpool_limits(1):
        ds,props=proposals(G,S,B,Y,radii)
    P=(G-G.conj().T)/(2j); H=np.eye(len(G))-CAP*h(G)
    meta=dict(scope='Floating proposals only, any strict result requires an independent interval check',
              cap=CAP,ratio=T, dark=[dict(p=float(np.vdot(d,P@d).real),kTp=float(np.vdot(d,(H+T*P)@d).real)) for d in ds],rows=[])
    for case,radius,record in zip(cases,radii,props):
        score,x,rank,beta,resid,o,c=record
        meta['rows'].append(dict(case=str(case),radius=float(radius),score=float(score),rank=int(rank),
                                 beta=float(beta),residual_2=float(resid),optical=float(o),combined=float(c)))
    meta['wall_seconds']=perf_counter()-start[0];meta['cpu_seconds']=process_time()-start[1]
    (HERE/'proposal.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(meta,flush=True)
