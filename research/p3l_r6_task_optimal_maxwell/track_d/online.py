"""Root-gated D generation/anchor/prediction process. Never opens test truth.

The generation subcommand delegates to the separate dedicated generator process.
This online module has no offline truth reader and no scoring entry point.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
from pathlib import Path
import argparse,json,sys,time,subprocess,gc
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_d.common import freeze,parents,context_path,identity,write,GROUPS,SEEDS,ROOT

def generate(group,d):
 # A separate generation process owns original prior-predictive truth generation.
 subprocess.run([sys.executable,'-B',str(OUT/'src/generate.py'),group,'--sigma-scale',str(GROUPS[group])],check=True,cwd=OUT)

def anchors(group,d):
 import shared.physics as physics
 from track_c.anchors import run,CONFIGS,REVISION1_CONFIGS
 registry={**CONFIGS,**REVISION1_CONFIGS}
 if d['selected_low_config'] not in registry:raise ValueError('selected low config absent from Revision1 registry')
 ids=parents(group)
 for index,parent in enumerate(ids):
  with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:
   y=z['observation'].copy();sigma=float(z['sigma'])
   if str(z['parent'])!=parent:raise ValueError('online parent mismatch')
  if abs(sigma/physics.SIGMA-GROUPS[group])>1e-13:raise ValueError('noise condition mismatch')
  pair=[('full','native96_probit64'),('low',d['selected_low_config'])]
  if index%2:pair.reverse()
  pair_order=[label for label,config in pair]
  for label,config in pair:
   dest=ROOT/'contexts'/group/parent/config;job=f'D_anchor_{parent}_{label}'
   reportpath=dest/'REPORT.json';closed=OUT/'raw/jobs'/job/'CLOSED.json'
   if reportpath.exists():
    report=json.loads(reportpath.read_text())
    complete=closed.exists() and json.loads(closed.read_text())['status']=='COMPLETE'
    recovery_path=dest/'ADMINISTRATIVE_RECOVERY.json'
    if recovery_path.exists():
     recovery=json.loads(recovery_path.read_text())
     complete=(recovery.get('status')=='RECOVERED_REPORT_ONLY_NO_NEW_PHYSICS' and recovery.get('original_closed_identity')==identity(closed)
       and recovery.get('context_identity')==identity(dest/'CONTEXT.npz') and recovery.get('recovered_report_identity')==identity(reportpath)
       and recovery.get('QA_identity')==identity(recovery['QA_identity']['path']))
    if not (complete and report.get('GT_read') is False and report.get('path_complete') is True and report.get('freeze_identity')==identity(ROOT/'D_FREEZE.json') and report.get('context_identity')==identity(dest/'CONTEXT.npz')):raise RuntimeError('completed context cannot be safely resumed; root audit required')
    continue
   if dest.exists():raise RuntimeError('partial context: no implicit replay; root audit required')
   b=physics.Book('D',job,parent,'ONLINE_INFERENCE',max_wall=d['anchor_max_wall_seconds']);status,error='COMPLETE',None;start=time.perf_counter()
   try:
    target=physics.Target(b,sigma=sigma)
    _,report=run(target,y,parent,dest,registry[config])
    cleanup_start=time.perf_counter();target.backend()._transfers.clear();del target
    collected=gc.collect();cleanup_seconds=time.perf_counter()-cleanup_start
    write(reportpath,dict(report,sigma=sigma,GT_read=False,context_identity=identity(dest/'CONTEXT.npz'),freeze_identity=identity(ROOT/'D_FREEZE.json'),cleanup_seconds=cleanup_seconds,GC_collected_objects=collected,pair_order=pair_order,frozen_parent_index=index,independent_prior_initialization=True,end_to_end_anchor_wall_seconds=time.perf_counter()-start,wall_scope='target setup + physical path + context persistence + cleanup; excludes Book close/administrative receipt write'))
   except BaseException as exc:status,error='FAILED',repr(exc);raise
   finally:b.close(status,error)
   progress=dict(event='D_ANCHOR_CONFIG_COMPLETE',group=group,parent=parent,config=config,completed_parents_before=index,pair_order=pair_order,RHS=physics.Book.totals()['RHS'],utc=time.time())
   with (ROOT/'ONLINE_PROGRESS.jsonl').open('a') as f:f.write(json.dumps(progress)+'\n')
   print(json.dumps(progress),flush=True)
 write(ROOT/f'ANCHORS_COMPLETE_{group}.json',dict(status='COMPLETE',parents=ids,new_probability_claim=False,freeze_identity=identity(ROOT/'D_FREEZE.json')))

def predict(group,d):
 # No estimator prediction until all512 independent Full/Low contexts complete.
 for condition in GROUPS:
  marker=ROOT/f'ANCHORS_COMPLETE_{condition}.json'
  if not marker.exists():raise RuntimeError('all512 contexts must complete before any D prediction')
  ready=json.loads(marker.read_text())
  if ready.get('status')!='COMPLETE' or ready.get('parents')!=parents(condition) or ready.get('freeze_identity')!=identity(ROOT/'D_FREEZE.json'):raise ValueError('all512 context completion identity')
 from track_c.frozen_heads import load_action,predict_contexts,ONLINE_KEYS
 ids=parents(group);destination=ROOT/'predictions'/group
 if destination.exists():
  sealpath=destination/'PREDICTIONS_SEALED.json'
  if not sealpath.exists():raise RuntimeError('partial prediction job; root audit required')
  existing=json.loads(sealpath.read_text())
  if existing['freeze_identity']!=identity(ROOT/'D_FREEZE.json') or existing['parents']!=ids:raise ValueError('existing seal freeze/parents differ')
  for record in existing['records']:
   if identity(record['identity']['path'])!=record['identity']:raise ValueError('existing prediction seal changed')
  if identity(existing['phases_identity']['path'])!=existing['phases_identity']:raise ValueError('existing phase seal changed')
  seal_global(d)
  return
 destination.mkdir(parents=True);records=[];context_artifacts=[]
 for kind in ('full','low'):
  for parent in ids:
   path=context_path(group,parent,kind,d);reportpath=path.with_name('REPORT.json');report=json.loads(reportpath.read_text())
   if report.get('context_identity')!=identity(path) or report.get('freeze_identity')!=identity(ROOT/'D_FREEZE.json') or report.get('GT_read') is not False or report.get('path_complete') is not True:raise ValueError('context report identity/freeze/isolation failed before prediction')
   context_artifacts.append(dict(parent=parent,context=kind,context_identity=identity(path),report_identity=identity(reportpath)))
 systems=d['estimator_arms']
 for name,arm in systems.items():
  contexts=[];anchor_cost=[]
  for parent in ids:
   path=context_path(group,parent,arm['context'],d)
   with np.load(path,allow_pickle=False) as z:
    if str(z['parent'])!=parent:raise ValueError('context order')
    contexts.append({k:z[k].copy() for k in (*ONLINE_KEYS,'parent')})
   report=json.loads(path.with_name('REPORT.json').read_text())
   anchor_cost.append(report['end_to_end_anchor_wall_seconds'])
  for seed in SEEDS:
   start=time.perf_counter()
   if arm.get('feature_mode')=='kkt_q':
    from track_c.kkt_closure.closure_action import Action as KKTAction
    action=KKTAction(method=arm['original'],seed=seed,closure_mode='q',noise_ratio=GROUPS[group])
   else:action=load_action(arm['original'],seed)
   load=time.perf_counter()-start
   start=time.perf_counter();prediction=predict_contexts(contexts,action,ids,batch_size=64);duration=time.perf_counter()-start
   if prediction.shape!=(256,576) or not np.isfinite(prediction).all():raise ValueError('finite full-grid material prediction required; no GT fallback')
   path=destination/f'{name}_s{seed}.npz';np.savez_compressed(path,parent=np.array(ids),prediction=prediction)
   rec=dict(method=name,context=arm['context'],feature_mode=arm.get('feature_mode','actual'),noise_ratio=GROUPS[group],seed=seed,identity=identity(path),load_wall_seconds=load,
      prediction_wall_seconds=duration,per_parent_warm_learning_seconds=duration/256,
      per_parent_cold_anchor_seconds=anchor_cost,GT_read=False,new_physics_RHS=0,
      timing_qualification='full/low each separately paid actual path plus measured batch64 original frozen action inference')
   records.append(rec);del action
 # Phases use full-anchor fraction, never true material phase.
 from track_c.evaluate_dev import make_prior
 prior=make_prior();phases=[];codes=[]
 for parent in ids:
  with np.load(context_path(group,parent,'full',d),allow_pickle=False) as z:q=z['q0'].copy()
  field=prior.config['field_mean']+prior.A@q[:64];s=field[prior.mask]/prior.temperature
  phases.append([float(np.mean(s<=-1)),float(np.mean((s>-1)&(s<1))),float(np.mean(s>=1))]);code=np.full(576,-1,dtype=np.int8);code[prior.mask]=np.where(s<=-1,0,np.where(s>=1,2,1));codes.append(code)
 phasepath=destination/'FULL_ANCHOR_PHASES.npz';np.savez_compressed(phasepath,parent=np.array(ids),fractions=np.array(phases),phase_code=np.array(codes))
 write(destination/'PREDICTIONS_SEALED.json',dict(status='GROUP_PREDICTIONS_SEALED',parents=ids,records=records,
   phases_identity=identity(phasepath),context_artifacts=context_artifacts,GT_read=False,freeze_identity=identity(ROOT/'D_FREEZE.json')))
 seal_global(d)

def seal_global(d):
 if all((ROOT/'predictions'/g/'PREDICTIONS_SEALED.json').exists() for g in GROUPS):
  seals=[]
  for g in GROUPS:
   sealpath=ROOT/'predictions'/g/'PREDICTIONS_SEALED.json'
   groupseal=json.loads(sealpath.read_text())
   if groupseal['freeze_identity']!=identity(ROOT/'D_FREEZE.json') or groupseal['parents']!=parents(g):raise ValueError('both noise conditions must use identical root freeze/parents')
   expected={(name,seed) for name in d['estimator_arms'] for seed in SEEDS}
   actual=[(r['method'],r['seed']) for r in groupseal['records']]
   if len(actual)!=len(expected) or set(actual)!=expected:raise ValueError('complete frozen action/seed family required')
   for record in groupseal['records']:
    if identity(record['identity']['path'])!=record['identity']:raise ValueError('prediction changed before global seal')
   if identity(groupseal['phases_identity']['path'])!=groupseal['phases_identity']:raise ValueError('phase changed before global seal')
   for ctx in groupseal['context_artifacts']:
    for key in ('context_identity','report_identity'):
     if identity(ctx[key]['path'])!=ctx[key]:raise ValueError('context/report changed before global seal')
   seals.append(identity(sealpath))
  globalpath=ROOT/'ALL_512_PREDICTIONS_SEALED.json'
  if globalpath.exists():
   old=json.loads(globalpath.read_text())
   if old.get('group_seals')!=seals or old.get('freeze_identity')!=identity(ROOT/'D_FREEZE.json') or old.get('status')!='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT':raise ValueError('global seal already exists with different immutable inputs')
   return
  write(globalpath,dict(status='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT',group_seals=seals,freeze_identity=identity(ROOT/'D_FREEZE.json'),GT_read=False,utc=time.time()))

def main():
 p=argparse.ArgumentParser();p.add_argument('operation',choices=['generate','anchors','predict']);p.add_argument('group',choices=list(GROUPS));a=p.parse_args();d=freeze();
 if (ROOT/'GT_UNLOCK.json').exists() or (ROOT/'scores').exists():raise RuntimeError('online phase closed after GT unlock/scoring; no prediction or seal rewrites')
 {'generate':generate,'anchors':anchors,'predict':predict}[a.operation](a.group,d)
if __name__=='__main__':main()
