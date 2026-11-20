"""Compare frequency and spatial degrees of freedom at fixed acquisition cost."""
import csv
import json
from itertools import combinations
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits
from designs import Design, all_designs, spaced, SIGMA_FRACTION

OUT=Path(__file__).resolve().parents[1]/'results'


def tangent_grams(a, f, sigma):
    # Columns dx/mm, dy/mm, permittivity amplitude, and a real shared grid-discrepancy direction.
    j=a['jac'+str(f)]
    bias=a['y48']-a['y24']
    v=np.concatenate((j,bias[...,None]),axis=-1)
    # One complex receiver measurement supplies two real Gaussian components.
    return (2/sigma**2)*(v.conj()[...,None]*v[...,None,:]).real


def schur_info(grams, nuisance):
    sub=grams[:,:2,2:2+nuisance];n=grams[:,2:2+nuisance,2:2+nuisance]
    inv=np.linalg.pinv(n,hermitian=True,rcond=1.e-12)
    return grams[:,:2,:2]-sub@inv@sub.swapaxes(-1,-2)


def loc_bound(g, count):
    info=schur_info(g,count)
    eigen=np.linalg.eigvalsh(info)
    if np.any(eigen[:,0]<=1.e-11*np.maximum(eigen[:,1],1.e-20)):
        return float('inf')
    return float(1./np.sqrt(eigen[:,0].min()))


def at(grams,d,fi):
    # Order f,case,rx,tx,par,par.
    idx=[fi[f] for f in d.freqs]
    slab=grams[idx][:,:,d.rx][:,:,:,d.tx]
    return slab.sum(axis=(0,2,3))*d.shot_energy


def image_information(j,pixel_sd,sigma,energy):
    v=j.reshape(-1,j.shape[-1]);a=np.vstack((v.real,v.imag))*(np.sqrt(2*energy)*pixel_sd/sigma)
    m=a@a.T if len(a)<=a.shape[1] else a.T@a
    e=np.maximum(eigvalsh(m,check_finite=False),0.)
    return dict(effective=float(np.sum(e/(1.+e))),above_one=int(np.sum(e>1.)),
                information_nats=float(np.log1p(e).sum()/2),largest=float(np.sqrt(e[-1])))


def measurement_row(d, bank, fi, sigma, grams):
    x=at(grams['24'],d,fi);z=at(grams['48'],d,fi)
    idx=[fi[f] for f in d.freqs]
    j=bank['pixels48'][idx][:,d.rx][:,:,d.tx];jc=bank['pixels24'][idx][:,d.rx][:,:,d.tx]
    a=image_information(j,.30,sigma,d.shot_energy);b=image_information(jc,.30,sigma,d.shot_energy)
    row=dict(design=d.name,frequencies=len(d.freqs),transmitters=len(d.tx),receivers=len(d.rx),
             shots=d.shots,complex_samples=d.samples,energy=1.,min_ghz=min(d.freqs),max_ghz=max(d.freqs),
             noise_sd=sigma,prior_pixel_sd=.3,eff_48=a['effective'],eff_24=b['effective'],
             modes_above_noise=a['above_one'],info_nats=a['information_nats'],largest_singular=a['largest'],
             location_sd_known_model_mm=loc_bound(z,1),location_sd_coarse_mm=loc_bound(x,1),
             location_sd_profiled_model_mm=loc_bound(z,2))
    # Reweight the same spectrum for two prior/noise ratio sensitivities on the core trade.
    if d.shots==8 and len(d.rx)==8:
        for scale in (.1,1.):
            row[f'eff_sd_{scale}']=image_information(j,scale,sigma,d.shot_energy)['effective']
    return row


def scheduled_grid():
    opts=(2.,6.,10.,12.,16.,20.,24.,28.,30.)
    for nf in (1,2,4,8):
        nt=8//nf
        for fr in combinations(opts,nf):
            for rot in range(8//nt):
                tx=tuple((i+rot)%8 for i in spaced(8,nt))
                yield Design('',fr,tx,spaced(16,8))


def design_record(d,name):
    return dict(name=name,freqs=d.freqs,tx=d.tx,rx=d.rx,energy=1.,shots=d.shots,complex_samples=d.samples)


def main():
    t,c=perf_counter(),process_time()
    with np.load(OUT/'task_fields.npz') as z: bank=dict(z)
    with np.load(OUT/'acquisition.npz') as z: sigma=float(z['reference_amplitude'])*SIGMA_FRACTION
    fi={f:i for i,f in enumerate(bank['frequencies_ghz'])}; gg={k:tangent_grams(bank,int(k),sigma) for k in ('24','48')}
    # Pre-aggregate identical Rx sums and all physical Tx subsets across candidate schedules.
    gr={k:g[:,:,spaced(16,8)].sum(axis=2) for k,g in gg.items()}
    tc=perf_counter();best=[(np.inf,None),(np.inf,None)];evaluated=0
    for d in scheduled_grid():
        i=[fi[f] for f in d.freqs]; g=[v[i][:,:,d.tx].sum(axis=(0,2))/8 for v in gr.values()]
        scores=(loc_bound(g[0],1),max(loc_bound(g[0],2),loc_bound(g[1],2)))
        for k,s in enumerate(scores):
            if s<best[k][0]:best[k]=(s,d)
        evaluated+=1
    chosen=[design_record(d,key) for (_,d),key in zip(best,('loc_nominal','loc_grid_profiled'))]
    meta=dict(designs=chosen,criterion_values_mm=[v[0] for v in best],virtual_centers_mm=bank['centers_mm'].tolist(),
              candidates=evaluated,search_wall_s=perf_counter()-tc,available_frequencies_ghz=[2,6,10,12,16,20,24,28,30],
              inputs='Known host and public nine virtual defects. No target observations or target material arrays used.')
    (OUT/'task_design.json').write_text(json.dumps(meta,indent=2)+'\n');print('selected:',meta,flush=True)
    designs=all_designs()
    for item in chosen:
        d=Design(item['name'],tuple(item['freqs']),tuple(item['tx']),tuple(item['rx']));designs[d.name]=d
    rows=[measurement_row(d,bank,fi,sigma,gg) for d in designs.values()]
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/'dof_trade.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=cols,lineterminator='\n');w.writeheader();w.writerows(rows)
    (OUT/'dof_cost.json').write_text(json.dumps(dict(cpu_s=process_time()-c,wall_s=perf_counter()-t,designs=len(rows),design_time_s=meta['search_wall_s']),indent=2)+'\n')
    for r in rows:
        if r['shots']==8 and r['receivers']==8:print(r['design'],round(r['eff_48'],2),r['location_sd_known_model_mm'],flush=True)

if __name__=='__main__':
    with threadpool_limits(1):main()
