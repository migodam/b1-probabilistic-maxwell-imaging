"""Cache only physical port fields and extract the needed host inverse columns after screening."""
import argparse
import csv
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.stats import chi2
from threadpoolctl import threadpool_limits
from ndt import Air, BACKGROUND, BOUNDS, Host, debye_epsilon, dielectric_factor, grid
from designs import SIGMA_FRACTION, read_arm
from sparse_features import OUT, ROOT, estimate, prepare, real_stack, remove_shared_gain, schedule
sys.path.insert(0,str(ROOT))
from src.inversion import _bound_step, fullwave
from src.roi_refine import region, RestrictedPenalty
from run_inversion import pair_centers, score

SELECTED=['f1_t8_r8','f4_t2_r8','loc_nominal']


class PortFields:
    def __init__(self,n,f,debye):
        self.n=n; self.factor=dielectric_factor(f) if debye else 1.
        self.bg=debye_epsilon(f)-1 if debye else BACKGROUND
        air=Air(n,f); self.lu=lu_factor(np.eye(n*n)-self.bg*air.dense(),check_finite=False)
        self.e=lu_solve(self.lu,air.E0,check_finite=False)
        self.s=lu_solve(self.lu,air.S.T,trans=1,check_finite=False).T
        self.xy=air.xy

    def tangent(self,d):
        stride=self.n//24
        j=np.einsum('ri,it->rti',self.s[list(d.rx)],self.e[:,list(d.tx)])*self.factor
        return j.reshape(len(d.rx),len(d.tx),24,stride,24,stride).sum(axis=(3,5)).reshape(-1,576)

    def extract(self,d,ids):
        # (I-cG)^-1 G = ((I-cG)^-1-I)/c for a uniform scalar background c.
        n=self.n; a=n//24
        y,x=np.unravel_index(ids,(24,24))
        u,v=np.indices((a,a)).reshape(2,-1)
        sub=((a*y[:,None]+u)*n+a*x[:,None]+v).ravel()
        rhs=np.zeros((n*n,len(sub)),dtype=complex,order='F')
        rhs[sub,np.arange(len(sub))]=1.
        z=lu_solve(self.lu,rhs,overwrite_b=True,check_finite=False)[sub]
        z.flat[::len(sub)+1]-=1.
        return Host(self.e[np.ix_(sub,d.tx)],self.s[np.ix_(d.rx,sub)],z/self.bg,self.xy[sub])


def refine_with_discrepancy(blocks,d,m0,ids,obs,blank,ref,sd,ft,n,profile,stop_at_noise=True):
    a=n//24; qt=m0.ravel()[ids].copy();reg=RestrictedPenalty(ids);weight=1/(ref*np.sqrt(obs.size));calls=0
    # If hypothesis is truly compatible with known noise, further fitting is optional.
    # This threshold requires no independence between the fitted estimate and residual.
    limit=sd*np.sqrt(chi2.ppf(.95,2*(obs.size-int(profile)))/2)
    def evaluate(x,der=False):
        nonlocal calls
        calls+=1
        out=fullwave(blocks,np.repeat(x,a*a),jacobian=der,material_factors=ft)
        field,j=out if der else (out,None)
        rz=field.ravel()-obs
        if profile:
            rz=remove_shared_gain(rz[:,None],blank).ravel()
        if der:
            j=j.reshape(-1,len(qt),a*a).sum(axis=2)
            if profile:j=remove_shared_gain(j,blank)
        return rz,j
    r,j=evaluate(qt,True);initial=float(np.linalg.norm(r));steps=0;stop='six_steps'
    for i in range(6):
        if stop_at_noise and np.linalg.norm(r)<=limit:
            stop='residual_within_noise';break
        A=real_stack(j)*weight;v=real_stack(r)*weight;H=A.T@A;g=A.T@v
        regval,der=reg.evaluate(qt,H);g+=der;obj=.5*np.dot(v,v)+regval
        step=_bound_step(H,g,qt);slope=g@step
        for e in range(9):
            alpha=2.**-e;x=np.clip(qt+alpha*step,*BOUNDS);rr,_=evaluate(x)
            cost=.5*np.linalg.norm(rr*weight)**2+reg.evaluate(x)[0]
            if cost<=obj+1.e-4*alpha*slope:break
        else:
            stop='line_search';break
        qt,r=x,rr;steps=i+1
        if steps<6 and (not stop_at_noise or np.linalg.norm(r)>limit):r,j=evaluate(qt,True)
    image=np.zeros((24,24));image.ravel()[ids]=qt
    return image,dict(steps=steps,reason=stop,initial_norm_over_noise_rms=initial/(np.sqrt(len(obs))*sd),
             final_norm_over_noise_rms=float(np.linalg.norm(r)/(np.sqrt(len(obs))*sd)),tolerance_ratio=float(limit/(np.sqrt(len(obs))*sd)),
             model_calls=calls)


def run(args):
    with np.load(OUT/(args.data+'_observations.npz')) as z:
        keys=('frequencies_ghz','case_ids','splits','reference_amplitude','clean','clean_same','noise','background')
        data={k:z[k] for k in keys}
    ref=float(data['reference_amplitude']);ds=schedule();ports={};costs={};numbers=[];all_images=[]
    for name in SELECTED:
        d=ds[name];c=process_time()
        for f in d.freqs:
            if f not in ports:
                s=perf_counter();ports[f]=PortFields(48,f,args.debye);costs[f]=perf_counter()-s
        ps=[ports[f] for f in d.freqs]; f_index=[int(np.flatnonzero(data['frequencies_ghz']==f)[0]) for f in d.freqs]
        blank=data['background'][np.ix_(f_index,d.rx,d.tx)].ravel();sd=ref*SIGMA_FRACTION/np.sqrt(d.shot_energy)
        s=perf_counter();tangent=np.vstack([p.tangent(d) for p in ps]);op=prepare(tangent,blank,sd,True);prep=perf_counter()-s
        fits=[];ft=[p.factor for p in ps]
        for i,cid in enumerate(data['case_ids']):
            tic=perf_counter();y,_=read_arm(data,i,d);obs=y.ravel();Y=real_stack(remove_shared_gain(obs[:,None],blank))*np.sqrt(2.)/sd
            seed,ev,crit=estimate(*op,Y);seed=seed.reshape(24,24);ids,_=region(seed);t1=perf_counter();k=len(ids)*4
            if k:
                blocks=[p.extract(d,ids) for p in ps];t2=perf_counter()
                result,record=refine_with_discrepancy(blocks,d,seed,ids,obs,blank,ref,sd,ft,48,
                                                      not args.unprofile_refinement,not args.fixed_steps)
            else:
                result=seed;t2=t1;record=dict(steps=0,reason='empty_region',model_calls=0)
            t3=perf_counter();fits.append([seed,result]);row=dict(case=str(cid),design=name,n_pixels=len(ids),fine_states=k,
                screen_seconds=t1-tic,column_solve_seconds=t2-t1,nonlinear_seconds=t3-t2,
                warm_seconds=t3-tic,cold_seconds=t3-tic+sum(costs[f] for f in d.freqs)+prep,
                preparation_s=prep,field_rhs=len(d.freqs)*24,inverse_column_rhs=k*len(d.freqs),matched_score=float(ev[0]),threshold=float(crit),**record)
            numbers.append(row)
        all_images.append(fits);print(name,'CPU',round(process_time()-c,2),'shared cold',round(sum(costs[f] for f in d.freqs)+prep,2),flush=True)
    prefix=OUT/('port_online_'+args.data+args.suffix)
    np.savez_compressed(str(prefix)+'.npz',names=SELECTED,case_ids=data['case_ids'],images=all_images)
    (Path(str(prefix)+'.json')).write_text(json.dumps(dict(methods=['sparse','fine_discrepancy'],host_costs_s=costs,cases=numbers),indent=2)+'\n')


def score_fits(data,suffix):
    with np.load(OUT/(data+'.npz')) as z:gt=z['truth']
    pairs=pair_centers(OUT/(data+'.json'))
    rows=[]
    with np.load(OUT/('port_online_'+data+suffix+'.npz')) as z:
        for a,(arm,im) in enumerate(zip(z['names'],z['images'])):
            for i,(cid,both) in enumerate(zip(z['case_ids'],im)):
                for method,pic in zip(('sparse','fine_discrepancy'),both):
                    rows.append(dict(arm=str(arm),case=str(cid),method=method,**score(pic.ravel(),gt[i].ravel(),grid(24)*1000,pairs.get(str(cid)))))
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/('port_scores_'+data+suffix+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',default='holdout');p.add_argument('--debye',action='store_true');p.add_argument('--score',action='store_true')
    p.add_argument('--fixed-steps',action='store_true');p.add_argument('--unprofile-refinement',action='store_true')
    p.add_argument('--suffix',default='');args=p.parse_args()
    with threadpool_limits(1):
        if args.score:score_fits(args.data,args.suffix)
        else:run(args)
