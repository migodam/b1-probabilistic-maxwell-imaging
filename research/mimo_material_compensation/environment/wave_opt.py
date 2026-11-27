"""Refine minimax full-aperture waveforms at a fixed physical core contrast."""
import argparse,json
from time import perf_counter,process_time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from probe import Circuit, material_task, summarize_pair, at_load, HERE, PRIOR, R


def minimize_wave(side,t,sectors,startfile,iterations,order):
    timer=perf_counter(),process_time()
    ix,ann=material_task(side)
    wave=Circuit(side,400,ix)
    rho=wave.disk.centers[ix]-[.125,.0625]
    theta=np.arctan2(rho[:,1],rho[:,0])%(2*np.pi)
    group=np.where(ann>0,1+(ann-1)*sectors+np.minimum(sectors-1,(theta*sectors/(2*np.pi)).astype(int)),0)
    groups=np.unique(group[group>0]); masks=[group==g for g in groups]; n=len(groups)
    r=float(np.load(PRIOR/'ports/public/cases'/f'n{side}_c01.npz')['radius'])
    def read_start():
        with np.load(startfile) as p:
            a=p['groups']; l=len(a); b=p[f'v_{t}'] if f'v_{t}' in p else p['v']
            q=int(p['sectors']) if 'sectors' in p else 1
            def value(g,j):
                ann0=(g-1)//sectors; sec0=(g-1)%sectors
                old=1+ann0*q+min(q-1,sec0*q//sectors)
                hit=np.flatnonzero(a==old)
                return float(b[hit[0]+j*l]) if len(hit) else 2.25
            return np.array([value(g,j) for j in range(2) for g in groups])
    initial=read_start()
    loads=np.exp(2j*np.pi*np.arange(4)/4)
    def objective(v,grad=True):
        x0=np.zeros(len(ix));x1=np.zeros(len(ix));x1[group==0]=t
        for j,m in enumerate(masks):x1[m]=v[j];x0[m]=v[n+j]
        mats=[]
        for x in (x1,x0):
            s,l,e=wave.solve(x,True)
            direct=np.asarray([-2j*R*(l[:,m]@e[m,:]) for m in masks])
            mats.append((s,direct))
        score=0.; derivative=np.zeros(2*n)
        for z in loads:
            (s1,h1),(s0,h0)=mats
            delta=(at_load(s1,z)-at_load(s0,z))/(2*r)
            u,s,vh=np.linalg.svd(delta,full_matrices=False)
            score+=sum(s**order)
            g=(u*(s**(order-1)))@vh
            for j,(p,ds) in enumerate(mats):
                f=z/(1-z*p[24,24]); left=np.zeros((16,25),complex);left[:,:16]=np.eye(16)
                right=left.T.copy();left[:,24]=f*p[:16,24];right[24,:]=f*p[24,:16]
                der=(left@ds@right)/(2*r)
                derivative[j*n:(j+1)*n]+=(1-2*j)*order*np.einsum('ij,kij->k',g.conj(),der).real
        pnorm=(score/len(loads))**(1/order)
        return pnorm,derivative*(pnorm/(order*score))
    result=minimize(objective,initial,jac=True,bounds=[(0,4.5)]*(2*n),method='L-BFGS-B',
        options=dict(maxiter=iterations,ftol=1e-14,gtol=1e-8,maxcor=15))
    xp=np.full(len(ix),t);xm=np.zeros(len(ix))
    for i,m in enumerate(masks):xp[m]=result.x[i];xm[m]=result.x[n+i]
    sp,sm=wave.solve(xp),wave.solve(xm)
    name=f'wave_{side}_{t}_{sectors}'
    pairs=summarize_pair(sp,sm,r)
    xpares=dict(ix=ix,group=group,groups=groups,sectors=sectors,v=result.x,plus=xp,minus=xm,p=sp,m=sm,r=r)
    np.savez_compressed(HERE/(name+'.npz'),**xpares)
    out=dict(side=side,core=t,sectors=sectors,start=Path(startfile).name,nit=int(result.nit),nfev=int(result.nfev),
             reason=str(result.message),objective=float(result.fun),distances=pairs,ops=wave.count,wall_s=perf_counter()-timer[0],cpu_s=process_time()-timer[1])
    (HERE/(name+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(side,t,sectors,round(pairs['all_loads_grid_sup_transfer_over_2r'],6),result.nit,out['wall_s'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',type=int);p.add_argument('t',type=float);p.add_argument('sectors',type=int);p.add_argument('start');p.add_argument('--iterations',type=int,default=200);p.add_argument('--order',type=int,default=8);a=p.parse_args()
    with threadpool_limits(1):minimize_wave(a.side,a.t,a.sectors,HERE/a.start,a.iterations,a.order)
