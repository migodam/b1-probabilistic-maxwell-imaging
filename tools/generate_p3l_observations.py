from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save

def main():
    group=sys.argv[1] if len(sys.argv)>1 else 'development'
    if group=='locked' and not (OUT/'METHODS_FROZEN_FOR_LOCKED.json').exists():raise RuntimeError('methods not frozen; locked generation prohibited')
    stage='A' if group=='development' else 'D'
    book=Book(stage,'dataset_'+group,max_wall=1200);status='COMPLETE';error=None
    try:
        target=Target(book,keep_fields=True);rows=json.loads((OUT/'configs/PARENTS.json').read_text())
        for row in [v for v in rows if v['group']==group]:
            parent=row['parent'];rng=np.random.default_rng(np.random.SeedSequence(row['seed']+[0]));q=target.prior.draw(rng)[0]
            noise_rng=np.random.default_rng(np.random.SeedSequence(row['seed']+[1]));noise=target.sigma/np.sqrt(2)*(noise_rng.normal(size=(32,8))+1j*noise_rng.normal(size=(32,8)))
            p=target.evaluate(q,'TRUTH_GENERATION_'+parent)
            save(OUT/f'inputs/offline_only/{parent}.npz',q=q,chi=p.state.chi,prediction=p.state.data,noise=noise,endpoint_key=np.array(p.key))
            save(OUT/f'inputs/online/{parent}.npz',observation=p.state.data+noise,target_id=np.array(target.cfg['target_id']),sigma=np.array(target.sigma))
            print(json.dumps(dict(parent=parent,generated=True,RHS=8)),flush=True)
        write(OUT/f'AUDIT/{group}_DATA_GENERATOR_EXIT_REQUIRED.json',dict(group=group,generated=sum(v['group']==group for v in rows),GT_inference=False,process='dedicated truth generator',no_inference_in_this_process=True))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:book.close(status,error)

if __name__=='__main__':main()
