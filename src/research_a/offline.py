"""Truth process only. Never imported by search.py or target.py."""
import json
import numpy as np
from scipy.stats import truncnorm
from b1_joint.material import prior_sample
from b1_joint.accounting import dump
from .target import OUT,to_p
from .material_metrics import metrics

def generate(target,n=32):
    rng=np.random.default_rng(261015102);noise=np.random.default_rng(261015103)
    x=prior_sample(rng,n);gain=truncnorm.rvs(-2,2,loc=1,scale=.15,size=n,random_state=rng)
    online=[];truth=[]
    for i,(xx,g) in enumerate(zip(x,gain)):
        p=to_p(xx);target.budget.context=dict(parent=i,role='fresh_truth')
        f=target.forward(p,kind='FRESH_TRUTH')[0];y=g*f+target.sigma[:32]*noise.normal(size=32)
        online.append(dict(id=f'confirm_{i:02d}',y=y,sigma=target.sigma[:32],target=target.physics_identity()))
        truth.append(dict(id=f'confirm_{i:02d}',p=p,latent=xx,gain=g,measurement=f))
    dump(OUT/'results/CONFIRM_OBSERVATIONS.json',online)
    dump(OUT/'results/offline_truth/CONFIRM_TRUTH.json',truth)
    return online

def score(campaign,truth,threshold=None):
    truths={c['id']:c for c in truth};rows=[]
    for item in campaign:
        p=truths[item['id']]['p']
        for result in item['methods']:
            archive=result['archive'];idx=result['selected_indices'];e=[metrics(c['p'],p) for c in archive]
            row=dict(parent=item['id'],arm=result['arm'],valid_k=len(idx),distinct_archive=len(archive),base=metrics(item['base']['p'],p),selected=e[result['selected_online_index']] if e else None,by_k={},calls=result['endpoint_calls'],wall=result['wall_seconds'],proxy_evaluations=result['proxy_evaluations'])
            for k in (1,2,4,8):
                ids=idx[:k]
                if not ids:row['by_k'][str(k)]=None;continue
                best=min(ids,key=lambda j:e[j]['normalized_mae']);entry=dict(best=e[best],candidate_id=best,actual_k=len(ids))
                if threshold is not None:entry['coverage']=bool(e[best]['normalized_mae']<=threshold and e[best]['central_radius_error']<=.05 and e[best]['ring_halfwidth_error']<=.025)
                row['by_k'][str(k)]=entry
            rows.append(row)
    return rows
