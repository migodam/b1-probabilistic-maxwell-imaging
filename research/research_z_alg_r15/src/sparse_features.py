"""Unknown-support detection and fast signed material imaging in the intact-plate tangent."""
import argparse
import csv
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigvalsh
from scipy.stats import norm
from threadpoolctl import threadpool_limits
from designs import Design, RECOVERY_ARMS, all_designs, read_arm, SIGMA_FRACTION
from ndt import BOUNDS, grid, read_data

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'
ALPHA=.05; SMOOTH=.025; ITERATIONS=160


def real_stack(a):
    return np.concatenate((a.real,a.imag),axis=0)


def remove_shared_gain(values,b):
    return values-b[:,None]*((b.conj()@values)/np.vdot(b,b))[None,:]


def laplacian(x):
    a=x.reshape(24,24,-1);g=np.zeros_like(a)
    dy=a[1:]-a[:-1];dx=a[:,1:]-a[:,:-1]
    g[1:]+=dy;g[:-1]-=dy;g[:,1:]+=dx;g[:,:-1]-=dx
    return g.reshape(576,-1)


def schedule():
    arms=all_designs()
    names=list(RECOVERY_ARMS)+['low_t8_r8','f8_t8_r16']
    for v in json.loads((OUT/'task_design.json').read_text())['designs'][:1]:
        d=Design(v['name'],tuple(v['freqs']),tuple(v['tx']),tuple(v['rx']));arms[d.name]=d;names.append(d.name)
    return {k:arms[k] for k in names}


def prepare(j,blank,sd,project):
    if project:j=remove_shared_gain(j,blank)
    A=np.ascontiguousarray(real_stack(j)*(np.sqrt(2.)/sd));w=np.maximum(np.linalg.norm(A,axis=0),1.e-20)
    A/=w
    gram=A@A.T if len(A)<=A.shape[1] else A.T@A
    top=eigvalsh(gram,subset_by_index=(len(gram)-1,len(gram)-1),check_finite=False)[0]
    sm=np.median(w)/w
    L=top+8.*SMOOTH*np.max(sm)**2
    return A,w,sm,L


def estimate(A,w,sm,L,Y):
    # In u_i=||A_i|| q_i units, each real matched filter has null variance one.
    k=norm.isf(ALPHA/(2.*len(w)))
    evidence=A.T@Y
    threshold=np.max(abs(evidence),axis=0)
    active=np.flatnonzero(threshold>k)
    u=np.zeros((len(w),len(active)));z=u.copy();moment=1.
    lo,hi=BOUNDS[0]*w[:,None],BOUNDS[1]*w[:,None]
    for _ in range(ITERATIONS):
        gradient=A.T@(A@z-Y[:,active])+SMOOTH*sm[:,None]*laplacian(sm[:,None]*z)
        p=z-gradient/L
        nxt=np.clip(np.sign(p)*np.maximum(abs(p)-k/L,0.),lo,hi)
        new=(1.+np.sqrt(1.+4*moment**2))/2
        z=nxt+(moment-1)/new*(nxt-u);u=nxt;moment=new
    image=np.zeros((len(w),Y.shape[1]));image[:,active]=u/w[:,None]
    return image.T,threshold,k


def null_trials(j,blank,sigma,seed):
    rng=np.random.default_rng(seed);trials=2000;critical=norm.isf(ALPHA/(2*j.shape[1]))
    cols=[j,remove_shared_gain(j,blank)];units=[v/np.maximum(np.linalg.norm(v,axis=0),1.e-30) for v in cols]
    counts=np.zeros((2,3),int)
    for _ in range(trials//100):
        noise=(rng.normal(size=(len(blank),100))+1j*rng.normal(size=(len(blank),100)))*(sigma/np.sqrt(2.))
        phases=rng.uniform(0,2.*np.pi,100);gain=.01*np.exp(1j*phases)
        # A separate stress arm uses independent per-channel phase/amplitude drift.
        channel=.01*np.exp(1j*rng.uniform(0.,2*np.pi,(len(blank),100)))
        for i,dic in enumerate(units):
            for k,obs in enumerate((noise,noise+blank[:,None]*gain,noise+blank[:,None]*channel)):
                tops=np.max(abs((dic.conj().T@obs).real),axis=0)*(np.sqrt(2.)/sigma)
                counts[i,k]+=np.sum(tops>critical)
    return dict(trials=trials,critical=critical,noise_only=counts[:,0].tolist(),shared_gain_1pct=counts[:,1].tolist(),independent_channel_gain_1pct=counts[:,2].tolist(),methods=['raw','profile_shared_gain'])


def run(data_path,tag):
    with np.load(data_path or OUT/'acquisition.npz') as source:
        data={key:source[key] for key in ('frequencies_ghz','case_ids','reference_amplitude',
                                          'clean','clean_same','noise','background')}
    with np.load(OUT/'task_fields.npz') as z:
        freq=z['frequencies_ghz'];banks={n:z['pixels'+str(n)] for n in (24,48)}
    arms=schedule();names=list(arms);cases=data['case_ids'];shape=(len(arms),len(cases),24,24)
    images={24:np.zeros(shape),48:np.zeros(shape)};states=[];null={}
    for idx,(name,d) in enumerate(arms.items()):
        fi=[int(np.flatnonzero(freq==f)[0]) for f in d.freqs];sl=np.ix_(fi,d.rx,d.tx)
        blank=data['background'][sl].ravel();rows=[read_arm(data,i,d)[0].ravel() for i in range(len(cases))]
        sd=float(data['reference_amplitude'])*SIGMA_FRACTION/np.sqrt(d.shot_energy)
        Y=real_stack(remove_shared_gain(np.array(rows).T,blank))*(np.sqrt(2.)/sd)
        for n in (24,48):
            t,c=perf_counter(),process_time();j=banks[n][sl].reshape(-1,576)
            A,w,sm,L=prepare(j,blank,sd,True);pre=perf_counter()-t
            s=perf_counter();im,test,crit=estimate(A,w,sm,L,Y);joint=perf_counter()-s
            images[n][idx]=im.reshape(-1,24,24)
            states.append(dict(arm=name,grid=n,preparation_wall_s=pre,batched_fit_wall_s=joint,
                fit_wall_s_per_material=joint/len(cases),cpu_s=process_time()-c,critical=crit,
                max_score=test.tolist(),per_material_estimation_includes=ITERATIONS))
        if tag=='' and name in ('f1_t8_r8','f2_t4_r8','f4_t2_r8','f8_t1_r8','loc_nominal','oneside_f2_t4_r8'):
            null[name]=null_trials(j,blank,sd,145316+idx)
    for n in (24,48):
        out=OUT/('sparse_materials'+('' if n==48 else '_coarse')+tag+'.npz')
        np.savez_compressed(out,names=names,case_ids=cases,images=images[n])
    (OUT/('sparse_cost'+tag+'.json')).write_text(json.dumps(states,indent=2)+'\n')
    if null:(OUT/'null_false_alarms.json').write_text(json.dumps(null,indent=2)+'\n')


def score_results(data_path,tag):
    sys.path.insert(0,str(ROOT));from run_inversion import score, pair_centers
    data=read_data(data_path);xy=grid(24)*1000;pc=pair_centers() if not tag else {};output=[]
    for n in (24,48):
        with np.load(OUT/('sparse_materials'+('' if n==48 else '_coarse')+tag+'.npz')) as z:
            for label,stack in zip(z['names'],z['images']):
                for i,(cid,image) in enumerate(zip(z['case_ids'],stack)):
                    met=score(image.ravel(),data['truth'][i].ravel(),xy,pc.get(str(cid)))
                    output.append(dict(arm=str(label),case=str(cid),grid=n,**met))
    keys=list(dict.fromkeys(k for row in output for k in row))
    with (OUT/('sparse_scores'+tag+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys,lineterminator='\n');w.writeheader();w.writerows(output)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path);p.add_argument('--tag',default='');p.add_argument('--score',action='store_true');a=p.parse_args()
    with threadpool_limits(1):
        if a.score:score_results(a.data,a.tag)
        else:run(a.data,a.tag)
