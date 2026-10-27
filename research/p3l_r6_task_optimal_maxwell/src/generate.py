"""Dedicated prior-predictive generation; never imported by online inference."""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.special import ndtr
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from shared.physics import OUT,Book,Target,SIGMA,write,save

def generate(group,sigma_scale=1.):
    manifest=json.loads((OUT/f'configs/parents_{group}.json').read_text())
    b=Book('C' if group=='CDEV' else 'D',f'generate_{group}',role='DATA_GENERATION',max_wall=3600)
    status='COMPLETE';error=None
    try:
        target=Target(b)
        for i,row in enumerate(manifest['parents']):
            parent=row['parent'];dst=OUT/'raw/online'/f'{parent}.npz'
            if dst.exists():raise FileExistsError('generation never replaces an existing observation')
            rng=np.random.default_rng(row['prior_seed']);x=rng.normal(size=66);q=x.copy();q[64:]=ndtr(q[64:])
            p=target.evaluate(q,purpose='new_prior_predictive_truth')
            noise_rng=np.random.default_rng(row['noise_seed'])
            innovation=(noise_rng.normal(size=(32,8))+1j*noise_rng.normal(size=(32,8)))/np.sqrt(2)
            sigma=SIGMA*float(sigma_scale);y=p.state.data+sigma*innovation
            save(OUT/'raw/offline_truth'/f'{parent}.npz',x=x,q=q,chi=p.state.chi,mean=p.state.data,noise_innovation=innovation,
                 parent=np.array(parent),prior_seed=np.array(row['prior_seed']),noise_seed=np.array(row['noise_seed']))
            save(dst,observation=y,sigma=np.array(sigma),parent=np.array(parent),target_id=np.array(target.cfg['target_id']))
            if i%16==0:print(json.dumps(dict(event='GENERATED',group=group,n=i+1,RHS=Book.totals()['RHS'])),flush=True)
        write(OUT/f'raw/GENERATION_{group}.json',dict(status='COMPLETE',n=len(manifest['parents']),sigma_scale=sigma_scale))
    except BaseException as exc:status='FAILED';error=repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('group');p.add_argument('--sigma-scale',type=float,default=1.);a=p.parse_args();generate(a.group,a.sigma_scale)
