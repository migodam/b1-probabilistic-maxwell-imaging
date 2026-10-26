"""Dedicated prior/truth generation process. No model evaluation or scoring."""
from pathlib import Path
import argparse,json,gc
import numpy as np
from physics_confirm import OUT,Book,Target,save,write,plain

def main():
    ap=argparse.ArgumentParser();ap.add_argument('group',choices=['probe','A','B']);args=ap.parse_args()
    manifest=json.loads((OUT/f'configs/parents_{args.group}.json').read_text())
    for index,row in enumerate(manifest['parents']):
        parent=row['parent'];dest=OUT/'raw/online'/f'{parent}.npz';truth=OUT/'raw/offline_only'/f'{parent}.npz'
        if dest.exists() and truth.exists():continue
        if dest.exists() or truth.exists():raise RuntimeError('partial generator endpoint: reconcile, no implicit overwrite')
        b=Book('generate_'+parent,parent,'TRUTH_GENERATION');status='COMPLETE';error=None
        try:
            t=Target(b)
            rng=np.random.default_rng(row['latent_seed']);q=t.prior.draw(rng)[0]
            nr=np.random.default_rng(row['noise_seed'])
            noise=t.sigma/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
            p=t.evaluate(q,'TRUTH_GENERATION_'+parent)
            save(truth,q=q,chi=p.state.chi,prediction=p.state.data,noise=noise,parent=np.array(parent))
            save(dest,observation=p.state.data+noise,sigma=np.array(t.sigma),target_id=np.array(t.cfg['target_id']))
            write(dest.with_suffix('.json'),dict(parent=parent,group=args.group,latent_seed=row['latent_seed'],noise_seed=row['noise_seed'],
                  prior_draw='original MaterialPrior.draw(default_rng(latent_seed))',RHS=8,truth_access_role='generation only'))
            del p,t;gc.collect()
        except BaseException as exc:status='FAILED';error=repr(exc);raise
        finally:b.close(status,error)
        if (index+1)%32==0 or args.group=='probe':print(json.dumps(dict(generated=index+1,group=args.group,cumulative_RHS=Book.totals()['RHS'])),flush=True)
    write(OUT/f'raw/GENERATED_{args.group}.json',dict(status='GENERATED',parents=len(manifest['parents']),group=args.group,scoring=False))
if __name__=='__main__':main()
