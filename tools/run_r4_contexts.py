"""Truth-free strong native/probit contexts; raw endpoint persistence is immediate."""
from pathlib import Path
import sys,argparse,json,gc,subprocess,os,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,Target,save,write
from p3r4.inversion import reconstruct
from p3r4.anchors import probit_optimize

def main():
 ap=argparse.ArgumentParser();ap.add_argument('group',choices=['dev','locked','confirmation','train','validation']);ap.add_argument('--count',required=True,type=int);ap.add_argument('--start',default=0,type=int)
 ap.add_argument('--native',default=96,type=int);ap.add_argument('--probit',default=64,type=int);a=ap.parse_args()
 if (a.native,a.probit)!=(96,64):raise ValueError('frozen common context budget96+64')
 stage='IV' if a.group in ('locked','confirmation') else 'I';cfg=json.loads((OUT/'configs/PROTOCOL.json').read_text())
 if a.group in ('train','validation'):
  prefix='R2_'+a.group;old=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
  cap=cfg['expanded_anchor_context_max'][0 if a.group=='train' else 1]
 else:
  d=cfg[{'dev':'new_development','locked':'new_locked','confirmation':'independent_confirmation'}[a.group]];prefix=d['prefix'];cap=d['n']
 if a.start<0 or a.count<1 or a.start+a.count>cap:raise ValueError('frozen groupcap')
 if stage=='IV' and not (OUT/f'configs/FROZEN_{a.group.upper()}.json').exists():raise ValueError('testfreeze')
 b=Book(stage,f'contexts_{a.group}_{a.start}_{a.count}',max_wall=7200);status='COMPLETE';error=None
 try:
  t=Target(b,keep_fields=True)
  for k in range(a.start,a.start+a.count):
   parent=f'{prefix}_{k:05d}';dest=OUT/f'results/contexts/{a.group}/{parent}.npz'
   if dest.exists():raise FileExistsError('context overwrite prohibited')
   if a.group in ('train','validation'):
    # Copy ONLY online keys; old offline/anchor truths never opened.
    path=old/f'inputs/online/{parent}.npz'
    with np.load(path,allow_pickle=False) as z:
     if set(z.files)!={'observation','target_id','sigma'}:raise ValueError('online source allowlist')
     if not (OUT/f'inputs/online/{parent}.npz').exists():save(OUT/f'inputs/online/{parent}.npz',**{key:z[key] for key in z.files})
   y=t.observation(parent);folder=OUT/f'results/contexts/details/{parent}';folder.mkdir(parents=True,exist_ok=True)
   b.context=dict(parent=parent,split=a.group,evidence_role='ONLINE_INFERENCE')
   before=json.loads((OUT/'RESOURCE_LEDGER.json').read_text())['stages'][stage];tick=time.perf_counter()
   native,nstat=reconstruct(t,y,parent,method='full_TRF',max_evaluations=a.native)
   save(folder/'NATIVE.npz',q=native.q,chi=native.state.chi,prediction=native.state.data)
   point,x,H,pstat=probit_optimize(t,y,parent,folder/'probit',native,cap=a.probit)
   J=t.jacobian(point)
   save(dest,observation=y,q0=point.q,chi0=point.state.chi,J=J,prediction=point.state.data,residual=y-point.state.data,
        phi=np.array(t.phi(point,y)),x=x,H=H,parent=np.array(parent))
   t.persist_cache(point,OUT/f'results/contexts/cache/{parent}.npz')
   after=json.loads((OUT/'RESOURCE_LEDGER.json').read_text())['stages'][stage]
   write(dest.with_suffix('.json'),dict(parent=parent,GT_read=False,native_cap=96,probit_cap=64,Phi=t.phi(point,y),
       native_stationarity=nstat['kkt'],probit_stationarity=pstat,
       new_RHS=after['RHS']-before['RHS'],new_states=after['states']-before['states'],new_LU=after['LU']-before['LU'],wall_seconds=time.perf_counter()-tick,
       target_id=t.cfg['target_id'],context_shift='new common96+64; old512/128 were48+48'))
   del point,native,J;t.backend()._transfers.clear();gc.collect()
   rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(os.getpid())],text=True))
   if rss>3*1024*1024:raise MemoryError('3GB physical worker guard')
   print(json.dumps(dict(parent=parent,completed=True,RSS_MiB=rss/1024)),flush=True)
 except BaseException as e:status='FAILED';error=repr(e);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
