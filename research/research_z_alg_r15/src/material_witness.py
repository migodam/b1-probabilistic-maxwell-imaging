"""Seek opposite preferences between two equal-cost MIMO schedules on finite material pairs."""
from pathlib import Path
from time import perf_counter, process_time
import json

import numpy as np
from scipy.linalg import eigh,lu_factor,lu_solve
from scipy.special import ndtr
from threadpoolctl import threadpool_limits
from designs import all_designs, SIGMA_FRACTION
from ndt import Air,BACKGROUND,grid

OUT=Path(__file__).resolve().parents[1]/'results'
PAIRS=(('f2_t4_r8','f4_t2_r8'),('f1_t8_r8','f8_t1_r8'),('f2_t4_r8','oneside_f2_t4_r8'),('high_t8_r8','f4_t2_r8'))
ARMS=all_designs();FREQUENCIES=sorted(set(f for p in PAIRS for name in p for f in ARMS[name].freqs))


def defect_patch():
    xy=1000*grid(24);r=np.linalg.norm(xy-np.array([1.,-1.]),axis=1)
    t=np.clip((5.5-r)/1.5,0,1);w=t*t*(3-2*t);ids=np.flatnonzero(w>0);w=w[ids]
    dist=np.linalg.norm(xy[ids,None]-xy[None,ids],axis=2)
    edge=(dist>.9)&(dist<1.1);L=np.diag(edge.sum(axis=1))-edge
    return ids,w,np.eye(len(ids))+2*L


def physical_patch(n, f, ids, w):
    ratio=n//24
    coarse=(np.arange(n*n).reshape(n,n)//n//ratio*24+(np.arange(n*n).reshape(n,n)%n)//ratio).ravel()
    select=np.flatnonzero(np.isin(coarse,ids));order=np.argsort(coarse[select],kind='stable');select=select[order]
    air=Air(n,f);G=air.dense();fac=lu_factor(np.eye(n*n)-BACKGROUND*G,check_finite=False)
    E=lu_solve(fac,air.E0,check_finite=False)[select]
    S=lu_solve(fac,air.S.T,trans=1,check_finite=False).T[:,select]
    K=lu_solve(fac,G[:,select],check_finite=False)[select]
    q=np.repeat(-.625*w,ratio**2)
    local=lu_factor(np.eye(len(q))-q[:,None]*K,check_finite=False)
    J=lu_solve(local,q[:,None]*E,check_finite=False)
    et=E+K@J; st=lu_solve(local,S.T,trans=1,check_finite=False).T
    # Host eliminated again about a legal, partially depleted material region.
    gt=lu_solve(local,K.T,trans=1,check_finite=False).T  # K (I-Dq K)^-1
    tangent=(st[:,:,None]*et[None,:,:]).reshape(16,len(ids),ratio**2,8).sum(axis=2).transpose(0,2,1)*w
    return dict(E=E,S=S,K=K,et=et,st=st,gt=gt,base=q,ratio=ratio,tangent=tangent)


def gram(terms,design,sd):
    grad=np.concatenate([terms[f]['tangent'][np.ix_(design.rx,design.tx)].reshape(-1,terms[f]['tangent'].shape[-1]) for f in design.freqs])*np.sqrt(design.shot_energy)/sd
    return (grad.conj().T@grad).real


def gap(terms,direct,z,w,t,sd):
    values={};identity=0.
    for f,x in terms.items():
        k=np.repeat(z*w,x['ratio']**2);q=x['base']; eye=np.eye(len(q));a=t*k
        j1=lu_solve(lu_factor(eye-(q+a)[:,None]*x['K'],check_finite=False),(q+a)[:,None]*x['E'],check_finite=False)
        j2=lu_solve(lu_factor(eye-(q-a)[:,None]*x['K'],check_finite=False),(q-a)[:,None]*x['E'],check_finite=False)
        diff=x['S']@(j1-j2);values[f]=diff
        if direct:
            op=k[:,None]*x['gt'];mid=2*t*x['st']@lu_solve(lu_factor(eye-t*t*(op@op),check_finite=False),k[:,None]*x['et'],check_finite=False)
            identity=max(identity,float(np.linalg.norm(mid-diff)/max(1.e-20,np.linalg.norm(diff))))
    def outcome(d):
        distance=np.sqrt(d.shot_energy*sum(np.linalg.norm(values[f][np.ix_(d.rx,d.tx)])**2 for f in d.freqs))/sd
        return dict(gap_over_sigma=float(distance),equal_prior_error=float(ndtr(-distance/np.sqrt(2.))))
    return outcome,identity


def main():
    t,c=perf_counter(),process_time();ids,w,pen=defect_patch()
    with np.load(OUT/'acquisition.npz') as a:sd=float(a['reference_amplitude'])*SIGMA_FRACTION
    all_terms={n:{f:physical_patch(n,f,ids,w) for f in FREQUENCIES} for n in (24,48)}
    rows=[];signs=[];examples=[]
    for left,right in PAIRS:
        gs={n:[gram(terms,ARMS[name],sd) for name in (left,right)] for n,terms in all_terms.items()}
        signs.append(dict(first=left,second=right,eigen_extremes={str(n):np.linalg.eigvalsh(g[1]-g[0])[[0,-1]].tolist() for n,g in gs.items()}))
        for small,large,index in ((left,right,0),(right,left,1)):
            g=gs[24];_,vec=eigh(g[1-index],g[index]+pen,subset_by_index=(len(ids)-1,len(ids)-1));z=vec[:,-1];z=z/np.max(abs(z))
            # Choose a finite perturbation on the coarse model, before finer prediction.
            # A gap of .8 sigma corresponds to equal-prior error about 28.6 percent.
            tmax=.625; amplitude=min(tmax,.8/max(2*np.sqrt(max(z@g[index]@z,0)),1.e-15))
            for _ in range(3):
                opt,_=gap({f:all_terms[24][f] for f in ARMS[small].freqs},False,z,w,amplitude,sd)
                sg=opt(ARMS[small])['gap_over_sigma'];amplitude=min(tmax,amplitude*.8/max(sg,1.e-15))
            row=dict(hard=small,easy=large,amplitude=float(amplitude),coarse_indices=ids.tolist(),
                    lower_material=(-.625*w-amplitude*w*z).tolist(),upper_material=(-.625*w+amplitude*w*z).tolist(),
                    patch_rms=float(np.linalg.norm(2*amplitude*w*z)/np.sqrt(len(ids))),global_rms=float(np.linalg.norm(2*amplitude*w*z)/24),
                    max_cell_difference=float(max(abs(2*amplitude*w*z))),frequency_ghz=FREQUENCIES,grid={})
            for n,terms in all_terms.items():
                used={f:terms[f] for f in sorted(set(ARMS[small].freqs+ARMS[large].freqs))}
                opt,err=gap(used,True,z,w,amplitude,sd)
                row['grid'][str(n)]=dict(hard=opt(ARMS[small]),easy=opt(ARMS[large]),identity_relative_error=err)
            rows.append(row);print(small,'vs',large,round(amplitude,3),row['grid']['48'],flush=True)
    result=dict(noise_sigma=sd,patch_cells=len(ids),full_material_cells=576,fine_em_states=2304,
                tangent_crossings=signs,finite_cases=rows,cpu_s=process_time()-c,wall_s=perf_counter()-t,
                hypotheses='Numerically validated paired partially depleted material maps. Same-map 48x48 verification uses four physical currents for each real material pixel.')
    (OUT/'material_witnesses.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    with threadpool_limits(1):main()
