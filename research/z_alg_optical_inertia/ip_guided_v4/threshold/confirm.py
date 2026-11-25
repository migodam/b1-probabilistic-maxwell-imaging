"""Prospectively frozen data, online selection, and separate material scoring."""
from pathlib import Path
import json,sys,argparse
import numpy as np
from flint import arb,ctx
from threadpoolctl import threadpool_limits
from run import HERE,POOL,Border,clock,elapsed
from pool_common import FullWave,write_arrays,write_json
from structured_pool import RX,TX,N_STARTS,MAX_CALLS,proposals
# Avoid the other replay module on research paths.
import importlib.util
spec=importlib.util.spec_from_file_location('pool_timed_replay',POOL/'replay.py')
pool_replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(pool_replay)
from box_certificates import Bounds
OUT=HERE/'confirmation'
WEIGHTS=(1/32,1/128)
def geometry():
    with np.load(POOL/'structured_results/public/geometry.npz') as f:
        return f['green'],f['sensor'],f['incident']
def generate():
    p,q=OUT/'public',OUT/'scoring'
    p.mkdir(parents=True,exist_ok=True);q.mkdir(exist_ok=True)
    if (p/'cases.json').exists(): raise FileExistsError('Independent measurements already exist.')
    g,s,b=geometry(); rows=[]
    for k in range(8):
        t=clock();seed=20261125101+k;rng=np.random.default_rng(seed)
        name=f'confirm49_{seed}__cw4_rx16_n0125'
        yy,xx=np.meshgrid(np.linspace(-1,1,7),np.linspace(-1,1,7),indexing='ij')
        theta=rng.uniform(-1,1);x=np.cos(theta)*xx+np.sin(theta)*yy-rng.uniform(-.2,.2);y=-np.sin(theta)*xx+np.cos(theta)*yy-rng.uniform(-.2,.2)
        g1=np.exp(-((x/.65)**2+(y/.55)**2)); rr=np.hypot(x,y*1.15)
        base=(g1,np.exp(-((rr-.67)/.24)**2),np.exp(-((x-.48)**2+(y-.4)**2)/.26)+.72*np.exp(-((x+.48)**2+(y+.2)**2)/.18),1+np.tanh(3*x+y),
              np.exp(-((rr-.68)/.27)**2)*(.55+.45*np.cos(np.arctan2(y,x))),np.exp(-((x-.2*y)/.25)**2)+.75*np.exp(-((y+.1*x)/.28)**2),
              (1+np.tanh((.59-abs(x))/.19))*(1+np.tanh((.47-abs(y))/.21))*(1+.15*x),.7*np.exp(-((x-.58)**2+(y+.19)**2)/.35)+.45*np.exp(-((x+.23)**2+(y-.43)**2)/.23)+.35*g1)[k]
        l=rng.uniform(1.65,2.75);peak=rng.uniform(3.7,4.9);m=(l+(peak-l)*base/base.max()).ravel()
        model=FullWave(g,s,b,np.ones((32,12)));field=model.evaluate(m)[0];fit=field[np.ix_(RX,TX)]
        r=.0125*np.linalg.norm(fit); n=rng.standard_normal(fit.shape)+1j*rng.standard_normal(fit.shape)
        Y=fit+r*n/np.linalg.norm(n);radius=float(r*(1+1e-7))
        write_arrays(p/(name+'.npz'),Y=Y,radius=radius)
        write_arrays(q/(name+'.npz'),truth=m,ports=field)
        rows.append(dict(id=name,seed=seed,setup=elapsed(t)))
    write_json(p/'cases.json',rows)
def proposed(g,s,b,y,r,seed,checker):
    total=clock();m=FullWave(g,s,b,y)
    samples,lo,hi,prep=proposals(g,s,b,y,m,seed)
    t=clock();checker.measurement(y)
    decisions=[checker.check(l,h,y,r,weights=WEIGHTS,jvalues=(1,)) for l,h in zip(lo,hi)]
    time=elapsed(t);surv=[i for i,a in enumerate(decisions) if not a['rejected']]
    auditor=Bounds(g,s,b,y);winner=None;resid=np.inf
    mats={};physical=[]
    for i in surv:
        mat,pred,fit=m.refine(samples[N_STARTS*i],lo[i],hi[i],MAX_CALLS)
        physical.append(fit);mats[i]=mat
        if fit['residual']<resid: winner,resid=i,fit['residual']
    v=auditor.point_upper(mats[winner]) if winner is not None else None
    row=dict(winner=winner,visited=surv,preparation=prep,screening=time,
             decisions=decisions,verified_residual_upper=v,estimated_residual=None if winner is None else resid,
             dense_optimization_counts=dict(m.counts),**elapsed(total))
    return row,mats[winner] if winner is not None else None,lo,hi
def online():
    ctx.prec=112
    g,fulls,fullb=geometry();s,b=fulls[RX],fullb[:,TX]
    start=clock();checker=Border(g,s,b);setup=elapsed(start)
    records={}; seal={}; refs=json.loads((OUT/'public/cases.json').read_text())
    names=['all_cluster','top_early','fast_boxes']
    for k,entry in enumerate(refs):
        name=entry['id'];target=entry['seed']
        with np.load(OUT/'public'/(name+'.npz')) as f:y,r=f['Y'],float(f['radius'])
        rec={}; order=list(np.random.default_rng(19361+k).permutation(names))
        for method in order:
            if method=='fast_boxes':
                row,mat,lo,hi=proposed(g,s,b,y,r,target+1007,checker)
                seal[name+'__lower']=lo;seal[name+'__upper']=hi
            else:
                row,chosen=pool_replay.run_one(method,g,s,b,y,r,target+1007)
                mat=chosen[0] if chosen else None
            rec[method]=row
            if mat is not None: seal[name+'__'+method]=mat
        records[name]=dict(radius=r,order=order,methods=rec)
    write_arrays(OUT/'predictions.npz',**seal)
    write_json(OUT/'online.json',dict(setup=setup,cases=records))
def score():
    ctx.prec=112
    g,fs,fb=geometry();s,b=fs[RX],fb[:,TX]; checker=Border(g,s,b)
    results={};raw=json.loads((OUT/'online.json').read_text())
    with np.load(OUT/'predictions.npz') as seal:
        for name,test in raw['cases'].items():
            with np.load(OUT/'public'/(name+'.npz')) as p, np.load(OUT/'scoring'/(name+'.npz')) as f:
                y,r=p['Y'],float(p['radius']);m,ports=f['truth'],f['ports']
            lo,hi=seal[name+'__lower'],seal[name+'__upper']
            contained=np.all((lo<=m)&(m<=hi),axis=1);truth_upper=Bounds(g,s,b,y).point_upper(m)
            pos=test['methods']['fast_boxes']['decisions']
            excluded=np.array([d['rejected'] for d in pos]);t=clock();checker.measurement(y)
            independent=[checker.check(l,h,y,r,verify=True,weights=WEIGHTS,jvalues=(1,)) for l,h,e in zip(lo,hi,excluded) if e]
            replay=elapsed(t)
            assert truth_upper<=r and not np.any(contained&excluded)
            assert all(a.get('arb_congruence_lower',0)>0 for a in independent)
            scores={}
            for method,row in test['methods'].items():
                est=seal[name+'__'+method];model=FullWave(g,fs,fb,ports);ph=model.evaluate(est)[0]
                rx=np.ix_(np.setdiff1d(np.arange(fs.shape[0]),RX),TX)
                tx=np.ix_(RX,np.setdiff1d(np.arange(fb.shape[1]),TX))
                scores[method]=dict(relative_material=float(np.linalg.norm(est-m)/np.linalg.norm(m)),background_rms=0,
                   heldout_receivers=float(np.linalg.norm((ph-ports)[rx])/np.linalg.norm(ports[rx])),
                   heldout_sources=float(np.linalg.norm((ph-ports)[tx])/np.linalg.norm(ports[tx])),
                   inside_noise=row['verified_residual_upper']<=r,runtime=row['wall_seconds'],
                   lu=row['dense_optimization_counts']['lu'],winner=int(row['winner']))
            results[name]=dict(covered=np.where(contained)[0].tolist(),excluded=np.where(excluded)[0].tolist(),
                               physical_reference_upper=truth_upper,radius=r,models=scores,
                               independent_ball_recheck=replay,checked=independent)
    write_json(OUT/'offline.json',results)
if __name__=='__main__':
    with threadpool_limits(1):
        cmd=argparse.ArgumentParser();cmd.add_argument('stage',choices=['generate','online','score'])
        {'generate':generate,'online':online,'score':score}[cmd.parse_args().stage]()
