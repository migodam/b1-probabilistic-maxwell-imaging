"""Attack two-frequency class-wide claims with shared nondispersive shell materials."""
import argparse,json
from time import perf_counter,process_time
import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from probe import HERE,PRIOR,Circuit,material_task,ACTIVE_PORTS,R,summarize_pair

def search(side,t,mode,nsteps):
    clock=perf_counter(),process_time();ix,group=material_task(side)
    groups=np.unique(group[group>0]); n=len(groups)
    masks=[group==g for g in groups]
    models=[Circuit(side,f,ix) for f in (400,460)]
    r=float(np.load(PRIOR/'ports/public/cases'/f'n{side}_c01.npz')['radius'])
    local=HERE/f'phase_{side}_two_shells.npz'
    s0=np.load(local)[f'v_{t}']
    prev=HERE/f'common_32_{t}_{mode}.npz'
    if side>32 and prev.exists():
        q=np.load(prev);s0=q['v']
    last={}
    def calc(y):
        if 'y' in last and np.array_equal(last['y'],y):return
        xp=np.full(len(ix),t);xm=np.zeros(len(ix))
        for j,mask in enumerate(masks):xp[mask]=y[j];xm[mask]=y[n+j]
        vec=[];allj=[];pall=[];mall=[]
        for w in models:
            p,lp,ep=w.solve(xp,True);m,lm,em=w.solve(xm,True)
            # The fitting norm is on the entire common augmented interface,
            # including the extra test terminal that need not be observed.
            pp=np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)
            dif=(p-m)[pp]
            deriv=[-2j*R*(lp[:,h]@ep[h])[pp] for h in masks]+[2j*R*(lm[:,h]@em[h])[pp] for h in masks]
            jac=np.stack(deriv,axis=-1).reshape(dif.size,2*n)
            allj.append(np.vstack((jac.real,jac.imag)));vec.append(np.r_[dif.real.ravel(),dif.imag.ravel()])
            pall.append(p);mall.append(m)
        last.update(y=y.copy(),jac=np.vstack(allj),res=np.concatenate(vec),p=pall,m=mall,plus=xp,minus=xm)
    def f(y):calc(y);return last['res']/(2*r)
    def j(y):calc(y);return last['jac']/(2*r)
    pre=float(np.linalg.norm(f(s0)))
    y=least_squares(f,s0,j,bounds=(0,4.5),ftol=1e-10,xtol=1e-10,gtol=1e-8,max_nfev=nsteps);calc(y.x)
    prefix=HERE/f'common_{side}_{t}_{mode}'
    data={str(k):summarize_pair(a,b,r) for k,a,b in zip([400,460],last['p'],last['m'])}
    np.savez_compressed(prefix.with_suffix('.npz') if '.' not in prefix.name else str(prefix)+'.npz',ix=ix,group=group,groups=groups,v=y.x,plus=last['plus'],minus=last['minus'],r=r)
    out=dict(nfev=y.nfev,stop=y.message,initial_combined=pre,fit=float(np.linalg.norm(y.fun)),diagnostics=data,
             costs=[w.count for w in models],wall_s=perf_counter()-clock[0],cpu_s=process_time()-clock[1])
    (HERE/(prefix.name+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(side,t,y.nfev,[(k,v['open']['max_cut'],v['all_loads_grid_sup_transfer_over_2r']) for k,v in data.items()],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',type=int);p.add_argument('t',type=float);p.add_argument('--steps',type=int,default=100);a=p.parse_args()
    with threadpool_limits(1):search(a.side,a.t,'shared',a.steps)
