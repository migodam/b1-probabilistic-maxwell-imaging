"""Post-score timing control; fixed metadata subset, no truth/model/risk selection.

Replays three unchanged physical algorithms in cyclic order in ONE window.
No original prediction/score is replaced. All replay solves are genuinely paid.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import sys,time,json,gc,hashlib
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P));R=P/'track_c/common_window_timing'
from shared.physics import Book,Target,write
from track_c.anchors import run as staged_run,CONFIGS,REVISION1_CONFIGS
from track_c.probit_prior_baseline import run as direct_run
from track_d.common import identity
N=16;CONFIG=('full96p64','low24p64','priorp128');MAX_RHS=(6400+3520+5160)*N*2

def main():
 if R.exists():raise FileExistsError('never replay/overwrite existing timing control')
 assert (P/'track_d/scores/SCORING_COMPLETE.json').exists()
 assert (P/'track_c/supplementary_prior_baseline/scores/SCORING_COMPLETE.json').exists()
 t=Book.totals();auth=json.loads((P/'configs/AUTHORIZATION.json').read_text())
 assert not t['active_jobs'] and not t['active_actions']
 assert t['RHS']+MAX_RHS<=auth['total_complex_RHS']
 assert t['stages']['D']['RHS']+MAX_RHS<=auth['stage_caps']['D']
 ids={g:[r['parent'] for r in json.loads((P/f'configs/parents_{g}.json').read_text())['parents']][:N] for g in ('DLOW','DHIGH')}
 write(R/'PROTOCOL_PRE_TIMING.json',dict(status='POST_SCORE_TIMING_REPLAY_FROZEN',new_risk_test=False,new_training=False,new_method=False,
  parents=ids,selection='first16 metadata-index parents per condition; never risk-selected',configs=CONFIG,
  order='cyclic permutations per index and condition; independent target/prior start each path',max_RHS=MAX_RHS,
  RHS_before=t['RHS'],scope='common-window physical timing; does not change global8 primary risk/cost outcomes',
  no_GT_reader=True,source_identity=identity(Path(__file__)),source_dependencies=[identity(P/'track_c/anchors.py'),identity(P/'track_c/probit_prior_baseline.py'),identity(P/'shared/physics.py')],
  statistical_qualification='20000 parent-paired bootstrap, descriptive two-sided95%; n16 per noise condition, not new independent risk confirmation',utc=time.time()))
 rows=[]
 for gi,(group,parents) in enumerate(ids.items()):
  for i,parent in enumerate(parents):
   with np.load(P/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:y,sigma=z['observation'].copy(),float(z['sigma'])
   order=list(CONFIG[(i+gi)%3:]+CONFIG[:(i+gi)%3])
   for name in order:
    b=Book('D',f'TIMING_{parent}_{name}',parent,'POST_SCORE_NO_GT_TIMING_CONTROL',max_wall=120);status,error='COMPLETE',None
    before=Book.totals();tick=time.perf_counter();destination=R/'contexts'/group/parent/name
    try:
     target=Target(b,sigma=sigma)
     if name=='priorp128':context,report=direct_run(target,y,parent,destination,cap=128)
     else:
      cfg=CONFIGS['native96_probit64'] if name=='full96p64' else REVISION1_CONFIGS['native24_probit64']
      context,report=staged_run(target,y,parent,destination,cfg)
     target.backend()._transfers.clear();del context,target;gc.collect();wall=time.perf_counter()-tick
     after=Book.totals()
     original=(P/'track_c/supplementary_prior_baseline/contexts'/group/parent/'prior_x0_probit128/CONTEXT.npz') if name=='priorp128' else P/'track_d/contexts'/group/parent/({'full96p64':'native96_probit64','low24p64':'native24_probit64'}[name])/'CONTEXT.npz'
     with np.load(original,allow_pickle=False) as a,np.load(destination/'CONTEXT.npz',allow_pickle=False) as c:
      differences={k:float(np.max(abs(a[k]-c[k]))) for k in ('q0','chi0','prediction')}
     row=dict(group=group,parent=parent,config=name,order=order,wall_seconds=wall,
       RHS=after['RHS']-before['RHS'],states=after['states']-before['states'],LU=after['LU']-before['LU'],
       physics_wall_seconds=after['physics_wall_seconds']-before['physics_wall_seconds'],original_context_identity=identity(original),endpoint_difference=differences,GT_read=False)
     rows.append(row);write(destination/'TIMING_REPORT.json',row)
    except BaseException as exc:status,error='FAILED',repr(exc);raise
    finally:b.close(status,error)
   if (i+1)%8==0:print(json.dumps(dict(group=group,parents_complete=i+1,RHS=Book.totals()['RHS'])),flush=True)
 result={};rng=np.random.default_rng(202610270802)
 for group,parents in ids.items():
  wall=np.array([[next(r['wall_seconds'] for r in rows if r['group']==group and r['parent']==p and r['config']==c) for c in CONFIG] for p in parents])
  ix=rng.integers(N,size=(20000,N));draw=wall[ix].mean(1);comparisons={}
  for a,bb in ((1,0),(2,0),(2,1)):
   saving=1-draw[:,a]/draw[:,bb]
   comparisons[f'{CONFIG[a]}_vs_{CONFIG[bb]}']=dict(saving=1-float(wall[:,a].mean()/wall[:,bb].mean()),descriptive95_interval=np.quantile(saving,[.025,.975]).tolist())
  result[group]=dict(n=N,anchor_wall_mean=dict(zip(CONFIG,wall.mean(0).tolist())),comparisons=comparisons,
   mean_RHS={c:float(np.mean([r['RHS'] for r in rows if r['group']==group and r['config']==c])) for c in CONFIG})
  np.savez_compressed(R/f'{group}_TIMING_DRAWS.npz',parent=np.array(parents),configs=np.array(CONFIG),wall_seconds=wall,bootstrap_indices=ix,bootstrap_wall_means=draw)
 write(R/'RESULTS.json',dict(status='COMPLETE_POST_SCORE_TIMING_CONTROL',conditions=result,records=rows,
  actual_new_RHS=Book.totals()['RHS']-t['RHS'],worst_pre_reserved_RHS=MAX_RHS,no_GT_reader=True,no_new_model_or_risk_claim=True,
  qualification='same physical algorithms/canonical environment; replay outputs compared to original sealed contexts; measured anchor wall includes context persistence+cleanup; no head re-fit or change to original global8 results'))

if __name__=='__main__':main()
