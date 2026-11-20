"""Spend a common first acquisition on a frequency, transmitter, or receiver extension."""
import argparse
import csv
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigvalsh
from threadpoolctl import threadpool_limits

from designs import Design, SIGMA_FRACTION
from ndt import grid
from port_online import PortFields
from sparse_features import OUT, prepare, estimate, real_stack, remove_shared_gain
from run_inversion import pair_centers, score

BASE_TX = (0, 2, 4, 6)
OTHER_TX = (1, 3, 5, 7)
BASE_RX = tuple(range(0, 16, 2))
OTHER_RX = tuple(range(1, 16, 2))
BASE = (16., BASE_TX, BASE_RX)


def gain(a):
    gram = a @ a.T if len(a) <= a.shape[1] else a.T @ a
    eigen = np.maximum(eigvalsh(gram, check_finite=False), 0.)
    return np.log1p(eigen).sum()/2, (eigen/(1+eigen)).sum()


def fitted_inputs(dataset):
    with np.load(OUT/(dataset+'_observations.npz')) as z:
        return {k: z[k] for k in ('frequencies_ghz','clean','noise','background','case_ids','reference_amplitude')}


def run(dataset):
    d=fitted_inputs(dataset)
    sigma=float(d['reference_amplitude'])*SIGMA_FRACTION
    ix={f:i for i,f in enumerate(d['frequencies_ghz'])}
    append=[('pilot_alone',None,4),('add_tx',(16.,OTHER_TX,BASE_RX),8),
            ('switch_rx',(16.,BASE_TX,OTHER_RX),8),('parallel_rx',(16.,BASE_TX,OTHER_RX),4)]
    append += [('add_'+str(f)+'GHz',(float(f),BASE_TX,BASE_RX),8) for f in (6,12,20,24,30)]
    cache={};host_times={};records=[];maps=[]
    for name,extension,shots in append:
        legs=[BASE] if extension is None else [BASE,extension]
        all_j=[]; blanks=[]; ys=[]
        sd=sigma*np.sqrt(shots)
        for f,tx,rx in legs:
            if f not in cache:
                t=perf_counter();cache[f]=PortFields(48,f,dataset=='holdout_debye');host_times[f]=perf_counter()-t
            a=cache[f].tangent(Design('',(f,),tx,rx))
            all_j.append(a); k=ix[f]
            blanks.append(d['background'][k][np.ix_(rx,tx)].ravel())
            v=(d['clean'][:,k][:,rx][:,:,tx]+sd*d['noise'][:,k][:,rx][:,:,tx])
            ys.append(v.reshape(len(v),-1).T)
        j=np.vstack(all_j);blank=np.concatenate(blanks);y=np.vstack(ys)
        before=perf_counter()
        projected=real_stack(remove_shared_gain(y,blank))*np.sqrt(2.)/sd
        imgs,_,_=estimate(*prepare(j,blank,sd,True),projected)
        fit=perf_counter()-before
        a_pilot=real_stack(all_j[0])*np.sqrt(2./4)*.30/sigma
        I0,_=gain(a_pilot)
        Ihalf,_=gain(a_pilot*np.sqrt(4./shots))
        I,d_eff=gain(real_stack(j)*np.sqrt(2./shots)*.30/sigma)
        records.append(dict(scheme=name,shots=shots,samples=len(j),power=1.,effective=float(d_eff),
                            info_nats=float(I),pilot_info_nats=float(I0),
                            dilution_nats=float(Ihalf-I0),conditional_new_nats=float(I-Ihalf),
                            net_nats=float(I-I0),batch_fit_s=fit))
        maps.append(imgs)
    np.savez_compressed(OUT/('nested_'+dataset+'.npz'),schemes=[r['scheme'] for r in records],
                       case_ids=d['case_ids'],maps=maps)
    (OUT/('nested_'+dataset+'.json')).write_text(json.dumps(dict(rows=records,port_cold_s=host_times,
        conventions='All schedules have energy one. Parallel Rx collects both receiver groups in four shots. Others with a second group use eight shots and half the per-shot energy of the pilot-only control. Known-gain local information and gain-profiled image scores are separate.'),indent=2)+'\n')


def evaluate(dataset):
    with np.load(OUT/(dataset+'.npz')) as z:
        truth=z['truth']
    pairs=pair_centers(OUT/(dataset+'.json'));xy=1000*grid(24)
    with np.load(OUT/('nested_'+dataset+'.npz')) as z:
        rows=[dict(scheme=str(name),case=str(key),**score(image,truth[i].ravel(),xy,pairs.get(str(key))))
              for name,batch in zip(z['schemes'],z['maps']) for i,(key,image) in enumerate(zip(z['case_ids'],batch))]
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/('nested_scores_'+dataset+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--data',default='holdout72');parser.add_argument('--score',action='store_true')
    args=parser.parse_args()
    with threadpool_limits(1):
        if args.score: evaluate(args.data)
        else: run(args.data)
