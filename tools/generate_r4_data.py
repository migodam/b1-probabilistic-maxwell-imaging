"""Dedicated prior simulation; online inversion never imports this module."""
from pathlib import Path
import argparse,sys,json,gc,os,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,Target,save,write

def main():
 ap=argparse.ArgumentParser();ap.add_argument('group',choices=['dev','locked','confirmation']);ap.add_argument('--count',type=int,required=True);ap.add_argument('--start',type=int,default=0);a=ap.parse_args()
 plan=json.loads((OUT/'configs/PROTOCOL.json').read_text())
 key={'dev':'new_development','locked':'new_locked','confirmation':'independent_confirmation'}[a.group];cfg=plan[key]
 if a.start<0 or a.count<1 or a.start+a.count>cfg['n']:raise ValueError('frozen parent count')
 if a.group!='dev' and not (OUT/f'configs/FROZEN_{a.group.upper()}.json').exists():raise ValueError('freeze methods before new test observations')
 stage='I' if a.group=='dev' else 'IV'
 b=Book(stage,f'generate_{a.group}_{a.start}_{a.count}',max_wall=900);status='COMPLETE';error=None
 try:
  t=Target(b,keep_fields=True)
  for k in range(a.start,a.start+a.count):
   parent=f"{cfg['prefix']}_{k:05d}";dest=OUT/f'inputs/online/{parent}.npz'
   if dest.exists() or (OUT/f'inputs/offline_only/{parent}.npz').exists():raise FileExistsError('never overwrite a parent')
   rng=np.random.default_rng(np.random.SeedSequence([cfg['seed'],k,0]));q=t.prior.draw(rng)[0]
   nr=np.random.default_rng(np.random.SeedSequence([cfg['seed'],k,1]));noise=t.sigma/np.sqrt(2)*(nr.normal(size=(32,8))+1j*nr.normal(size=(32,8)))
   b.context=dict(parent=parent,split=a.group,evidence_role='TRUTH_GENERATION')
   point=t.evaluate(q,'TRUTH_GENERATION_'+parent);y=point.state.data+noise
   save(OUT/f'inputs/offline_only/{parent}.npz',q=q,chi=point.state.chi,prediction=point.state.data,noise=noise)
   save(dest,observation=y,sigma=np.array(t.sigma),target_id=np.array(t.cfg['target_id']))
   write(dest.with_suffix('.json'),dict(parent=parent,split=a.group,seed=[cfg['seed'],k],new_parent=True,GT_scope='generation/training/offline only'))
   del point
   if (k-a.start+1)%32==0:t.backend()._transfers.clear();gc.collect();print(parent,flush=True)
  write(b.jobdir/'DATA_COMPLETE.json',dict(group=a.group,start=a.start,count=a.count,target_id=t.cfg['target_id'],new_parents=True))
 except BaseException as e:status='FAILED';error=repr(e);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
