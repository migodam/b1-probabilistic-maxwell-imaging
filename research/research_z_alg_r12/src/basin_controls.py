"""Fixed random, feed-limited and cold-start controls for the two exposed stationary endpoints."""
import json
from time import perf_counter
import numpy as np
from scipy.linalg import eigh

from basin_sources import Wave, OUT, source_directions, load_endpoint, hermitian


def capped_beam(w,M,energy_cost,max_feed):
    I=np.eye(w.s)
    def beam(penalty):
        _,u=eigh(w.white@(M+penalty*I)@w.white,subset_by_index=[0,0])
        return np.sqrt(energy_cost)*w.white@u[:,0]
    left=0.;right=max(np.linalg.norm(M,2),1e-18)
    a=beam(0.)
    if np.linalg.norm(a)<=max_feed:return a
    while np.linalg.norm(beam(right))>max_feed:right*=2
    for _ in range(48):
        mid=(left+right)/2
        if np.linalg.norm(beam(mid))>max_feed:left=mid
        else:right=mid
    return beam(right)


def run(label):
    _,pub,m=load_endpoint(label)
    w=Wave(pub);I=np.eye(w.s)
    base=w.score(m)
    tic=perf_counter();lu=w.lu;rhs=w.rhs
    geo,weights,intern=source_directions(w,m)
    floor=geo['index']
    selected=np.flatnonzero(intern[6]==floor)[0]
    M=intern[8][selected]
    weights['profile_feed_l2_le_2_per_unit_power']=capped_beam(w,M,w.s,2*np.sqrt(w.s))
    r=intern[0];ix=int(np.argmax(np.linalg.norm(r,axis=0)))
    weights['single_largest_residual_shot']=np.sqrt(w.s)*I[:,ix]
    setup=dict(seconds=perf_counter()-tic,LU=w.lu-lu,RHS=w.rhs-rhs)
    for seed in range(1,13):
        rng=np.random.default_rng(78000+seed)
        z=rng.normal(size=w.s)+1j*rng.normal(size=w.s)
        weights[f'random_{seed:02d}']=np.sqrt(w.s)*w.white@(z/np.linalg.norm(z))
    records={};mats={'initial':m}
    for name,a in weights.items():
        start=perf_counter();lu=w.lu;rhs=w.rhs
        v,i1=w.fit(m,.04*I+.96*np.outer(a,a.conj()),48)
        x,i2=w.fit(v,I,96)
        sc=w.score(x)
        records[name]=dict(phase1=i1,phase2=i2,**sc,wall_s=perf_counter()-start,LU=w.lu-lu,RHS=w.rhs-rhs,
              unit_power_feed_L2=float(np.linalg.norm(a)/np.sqrt(w.s)),kept=bool(sc['phi']<base['phi']))
        mats[name]=x
    for name,v in [('continue',m),('cold_background',np.full(w.n,w.lo)),('cold_midbox',np.full(w.n,(w.lo+w.hi)/2))]:
        start=perf_counter();lu=w.lu;rhs=w.rhs
        x,info=w.fit(v,I,144);sc=w.score(x)
        records[name]=dict(phase2=info,**sc,wall_s=perf_counter()-start,LU=w.lu-lu,RHS=w.rhs-rhs,kept=bool(sc['phi']<base['phi']))
        mats[name]=x
    # Choose by the measured original objective, then compare the same further corrector on all alternatives.
    best_random=min((k for k in records if k.startswith('random')),key=lambda k:records[k]['phi'])
    for k in ['profile',best_random,'cold_background','cold_midbox','continue']:
        start=perf_counter();lu=w.lu;rhs=w.rhs
        x,info=w.fit(mats[k],I,500);sc=w.score(x)
        records[k]['long']=dict(fit=info,**sc,wall_s=perf_counter()-start,LU=w.lu-lu,RHS=w.rhs-rhs)
        mats[k+'_long']=x
    (OUT/f'basin_controls_{label}.json').write_text(json.dumps(dict(label=label,initial=base,geometry_cost=setup,
        best_random_by_phi=best_random,results=records),indent=2)+'\n')
    np.savez_compressed(OUT/f'basin_controls_{label}.npz',**mats)
    print(label,best_random,{k:round(records[k]['long']['phi'],7) for k in records if 'long' in records[k]},flush=True)


if __name__=='__main__':
    for label in ['shell','ring']:run(label)
