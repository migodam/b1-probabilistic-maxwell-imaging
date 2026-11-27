"""Re-investigate the failed shell construction with independent material endpoints."""
from probe import Circuit, material_task, summarize_pair, ACTIVE_PORTS, R, HERE, PRIOR
from scipy.optimize import least_squares
from time import perf_counter,process_time
from threadpoolctl import threadpool_limits
import argparse,json
import numpy as np


def compute(side,mode):
    ix,group=material_task(side)
    wave=Circuit(side,400,ix)
    # Preserve the same physical annuli. The sector version relaxes radial symmetry
    # in the exterior and can compensate departures of the sampled core from a disk.
    if mode.endswith('sectors'):
        xy=wave.disk.centers[ix]-np.array((.125,.0625))
        sector=(np.arctan2(xy[:,1],xy[:,0])%(2*np.pi)//(np.pi/2)).astype(int)
        group=np.where(group,1+(group-1)*4+sector,0)
    groups=np.unique(group[group>0])
    masks=[group==i for i in groups]
    q=len(groups)
    pair=not mode.startswith('one')
    pars=(2 if pair else 1)*q
    y=np.full(pars,2.25)
    r=float(np.load(PRIOR/'ports/public/cases'/f'n{side}_c01.npz')['radius'])
    result={}; data={}; clock=perf_counter(),process_time()
    mids=(.5,1.,1.5,2.,2.5,3.,3.5)
    for t in mids:
        saved={}
        def val(v):
            if 'v' in saved and np.array_equal(v,saved['v']):return
            plus=np.full(len(ix),float(t));minus=np.zeros(len(ix))
            for i,mask in enumerate(masks):
                plus[mask]=v[i];minus[mask]=v[i+q] if pair else 2.25
            p,l1,e1=wave.solve(plus,True);m,l2,e2=wave.solve(minus,True)
            d=(p-m)[np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)]
            rows=[-2j*R*(l1[:,mask]@e1[mask])[np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)] for mask in masks]
            if pair: rows += [2j*R*(l2[:,mask]@e2[mask])[np.ix_(ACTIVE_PORTS,ACTIVE_PORTS)] for mask in masks]
            j=np.stack(rows,axis=-1).reshape(d.size,pars)
            saved.update(v=v.copy(),plus=plus,minus=minus,p=p,m=m,res=np.r_[d.real.ravel(),d.imag.ravel()],jac=np.vstack((j.real,j.imag)))
        def fun(v):val(v);return saved['res']
        def jac(v):val(v);return saved['jac']
        seed_file=HERE/f'phase_32_{mode}.npz'
        if side>32 and seed_file.exists() and t==mids[0]:
            with np.load(seed_file) as f:
                old_groups=f['groups']
                for i,g in enumerate(groups):
                    oldidx=np.flatnonzero(old_groups==g)
                    if len(oldidx):
                        y[i]=f[f'v_{t}'][oldidx[0]]
                        if pair:y[q+i]=f[f'v_{t}'][len(old_groups)+oldidx[0]]
        tol=2e-9
        opt=least_squares(fun,y,jac=jac,bounds=(0.,4.5),ftol=tol,xtol=tol,gtol=2e-10,max_nfev=55 if side==32 else 25)
        y=opt.x;val(y)
        result[str(t)]=dict(nfev=opt.nfev,stop=opt.message,optimality=opt.optimality,
                      distances=summarize_pair(saved['p'],saved['m'],r),
                      unobserved_coherent_norm=float(np.linalg.norm(saved['p']-saved['m'],2))/(2*r))
        for name in ('v','plus','minus','p','m'):data[f'{name}_{t}']=saved[name].copy()
        result['cost']=dict(wall_s=perf_counter()-clock[0],cpu_s=process_time()-clock[1],**wave.count)
        dat=dict(ix=ix,group=group,groups=groups,side=side,r=r,**data)
        np.savez_compressed(HERE/f'phase_{side}_{mode}.npz',**dat)
        (HERE/f'phase_{side}_{mode}.json').write_text(json.dumps(result,indent=2)+'\n')
        print(side,mode,t,opt.nfev,round(result[str(t)]['distances']['open']['max_cut'],5),
              round(result[str(t)]['distances']['all_loads_grid_sup_transfer_over_2r'],5),
              round(result['cost']['wall_s'],2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('side',type=int);ap.add_argument('modes',nargs='+')
    a=ap.parse_args()
    with threadpool_limits(1):
        for mode in a.modes:compute(a.side,mode)
