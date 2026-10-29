"""Separate paid data generator. Online runner never imports this module."""
from pathlib import Path
import argparse,json,time
import numpy as np
import possibility_volume as pv

def generate(c):
    root=Path(c['output']);private=Path(c['private_output'])
    private.mkdir(parents=True,exist_ok=True)
    (private/'.gitignore').write_text('*\n!.gitignore\n')
    if (root/'GENERATION_SEALED.json').exists():raise RuntimeError('no generation replay')
    with pv.Book(c,'GENERATE_FRESH32','DATA_GENERATION') as book:
        # All attempted truth states and solver arrays remain private too.
        book.jobdir=private/'physics_jobs'/book.job
        book.jobdir.mkdir(parents=True,exist_ok=False)
        ph=pv.ExactROI(book)
        try:
            for row in c['parents']:
                parent=row['parent'];online=root/'inputs/online'/f'{parent}.npz'
                truth=private/'truth'/f'{parent}.npz'
                if online.exists() or truth.exists():raise RuntimeError('existing generated parent no replay')
                x=np.random.default_rng(np.random.SeedSequence(row['prior_seed'])).normal(size=66)
                p=ph.evaluate(x,parent+'_private_generation')
                nr=np.random.default_rng(np.random.SeedSequence(row['noise_seed']))
                noise=pv.SIGMA/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
                observation=p.prediction+noise
                book.save(truth,x=x,q=p.q,material=p.material,prediction=p.prediction,
                          observation=observation,noise=noise,sigma=np.array(pv.SIGMA),parent=np.array(parent))
                book.save(online,observation=observation,sigma=np.array(pv.SIGMA),parent=np.array(parent))
                print('VOLUME_GENERATED',parent,book.total['reduced_RHS256'],flush=True)
            pv.write(root/'GENERATION_SEALED.json',dict(protocol_id=c['id'],utc=time.time(),
                parents=[r['parent'] for r in c['parents']],online_fields=['observation','sigma','parent'],
                truth_private=True,truth_initialized_online=False,ledger=book.total))
        finally:ph.close()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',default=str(pv.CONFIG));a=p.parse_args()
    generate(pv.read_protocol(a.config))
