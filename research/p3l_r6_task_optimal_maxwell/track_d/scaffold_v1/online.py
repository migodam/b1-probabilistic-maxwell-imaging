"""Root-gated D generation/anchor/prediction process. Never opens test truth.

The generation subcommand delegates to the separate dedicated generator process.
This online module has no offline truth reader and no scoring entry point.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
from pathlib import Path
import argparse,json,sys,time,subprocess
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_d.common import freeze,parents,context_path,identity,write,GROUPS,SEEDS,ROOT

def generate(group,d):
 # A separate generation process owns original prior-predictive truth generation.
 subprocess.run([sys.executable,'-B',str(OUT/'src/generate.py'),group,'--sigma-scale',str(GROUPS[group])],check=True,cwd=OUT)

def anchors(group,d):
 # Imports occur only after root freeze has passed.
 import shared.physics as physics
 from track_c.anchors import run,CONFIGS,_ObservedTarget,_Trace,_accepted_iteration_hook,build_online_context
 if d['selected_low_config'] not in CONFIGS:raise ValueError('selected low config absent from frozen Track C implementation')
 for parent in parents(group):
  with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:
   y=z['observation'].copy();sigma=float(z['sigma'])
   if str(z['parent'])!=parent:raise ValueError('online parent mismatch')
  if abs(sigma/physics.SIGMA-GROUPS[group])>1e-13:raise ValueError('noise condition mismatch')
  full=ROOT/'contexts'/group/parent/'native96_probit64';native=ROOT/'contexts'/group/parent/'native96_only'
  b=physics.Book('D',f'D_anchor_{parent}_full',parent,'ONLINE_INFERENCE',max_wall=d['anchor_max_wall_seconds'])
  status,error='COMPLETE',None;start=time.perf_counter();original=physics.probit_optimize
  target=physics.Target(b,sigma=sigma)
  def capture_native(point,native_report):
   # Original Native optimizer has already paid its endpoint J/KKT. Never solve
   # a counterfactual future state to build the native-only action context.
   if point.J is None:raise RuntimeError('paid Native J missing; stop rather than unbilled solve')
   if native.exists():raise FileExistsError('never overwrite Native prefix')
   native.mkdir(parents=True)
   physics.save(native/'CONTEXT.npz',parent=np.array(parent),observation=y,residual=y-point.state.data,
      prediction=point.state.data,q0=point.q,chi0=point.state.chi,J=point.J)
   target.persist_cache(point,native/'PAID_NATIVE_STATE.npz')
   ledger=physics.Book.totals()
   write(native/'REPORT.json',dict(parent=parent,GT_read=False,context='native96_only',
      captured_before_Probit=True,cold_native_phase_wall_seconds=native_report['wall_seconds'],
      cold_native_counterfactual_wall_seconds=time.perf_counter()-start,
      cold_native_RHS=8*native_report['fresh_forward_states']+32*native_report['receiver_J32_calls'],
      actual_prefix_reserved_RHS=ledger['RHS']-ledger_start['RHS'],
      cold_native_cost_qualification='paid full-run Native phase plus context persistence; no separate standalone timing replay',
      context_capture_wall_seconds=time.perf_counter()-start,incremental_reuse_RHS=0))
  ledger_start=physics.Book.totals()
  observed=_ObservedTarget(target,full);trace=_Trace(observed,y,parent,full/'path',CONFIGS['native96_probit64'])
  try:
   with _accepted_iteration_hook(physics.reconstruct,trace,'native'):
    native_point,native_report=physics.reconstruct(observed,y,parent,method='full_TRF',max_evaluations=96)
   capture_native(native_point,native_report)
   native_q=native_point.q.copy();native_chi=native_point.state.chi.copy()
   with _accepted_iteration_hook(original,trace,'probit'):
    point,x,H,probit_report=original(observed,y,parent,full/'original_optimizer/probit',native_point,cap=64)
   context=build_online_context(observed,y,parent,point,x,H);context.update(native_q=native_q,native_chi=native_chi)
   physics.save(full/'CONTEXT.npz',**context)
   after=physics.Book.totals();report=dict(parent=parent,config='native96_probit64',GT_read=False,
      original_native_report=native_report,original_probit_report=probit_report,path_complete=True,
      output_context=str(full/'CONTEXT.npz'),captured_native_before_Probit=True,
      path_resource_cost={k:after[k]-ledger_start[k] for k in ('RHS','forward_RHS','derivative_RHS','states','LU','physics_wall_seconds')})
   write(full/'ANCHOR_REPORT.json',report)
   write(full/'REPORT.json',dict(report,sigma=sigma,end_to_end_anchor_wall_seconds=time.perf_counter()-start))
  except BaseException as exc:status,error='FAILED',repr(exc);raise
  finally:trace.close();b.close(status,error)
  # No second Native solve for Native-selected candidate; a genuinely different
  # low path is paid once separately and its actual cost stays in singleton Book.
  if any(arm['context']=='low' for arm in d['matched_arms'].values()):
   dest=ROOT/'contexts'/group/parent/d['selected_low_config'];b=physics.Book('D',f'D_anchor_{parent}_low',parent,'ONLINE_INFERENCE',max_wall=d['anchor_max_wall_seconds']);status,error='COMPLETE',None;start=time.perf_counter()
   try:
    target=physics.Target(b,sigma=sigma);_,report=run(target,y,parent,dest,d['selected_low_config']);write(dest/'REPORT.json',dict(report,sigma=sigma,GT_read=False,end_to_end_anchor_wall_seconds=time.perf_counter()-start))
   except BaseException as exc:status,error='FAILED',repr(exc);raise
   finally:b.close(status,error)
 write(ROOT/f'ANCHORS_COMPLETE_{group}.json',dict(status='COMPLETE',parents=parents(group),new_probability_claim=False))

def predict(group,d):
 from track_c.frozen_heads import load_action,predict_contexts,ONLINE_KEYS
 from track_c.matched_heads import MatchedAction
 ids=parents(group);destination=ROOT/'predictions'/group
 if destination.exists():raise FileExistsError('no prediction overwrite or post-score replay')
 destination.mkdir(parents=True);records=[]
 systems={name:dict(context='full',original=name) for name in ('strongCNN','distilled_w128_n32768','deterministic_matched')}
 systems.update(d['matched_arms'])
 for name,arm in systems.items():
  contexts=[];anchor_cost=[]
  for parent in ids:
   path=context_path(group,parent,arm['context'],d)
   with np.load(path,allow_pickle=False) as z:
    if str(z['parent'])!=parent:raise ValueError('context order')
    contexts.append({k:z[k].copy() for k in (*ONLINE_KEYS,'parent')})
   report=json.loads(path.with_name('REPORT.json').read_text())
   anchor_cost.append(report['cold_native_counterfactual_wall_seconds'] if arm['context']=='native' else report['end_to_end_anchor_wall_seconds'])
  for seed in SEEDS:
   start=time.perf_counter();action=load_action(arm['original'],seed) if 'original' in arm else MatchedAction(arm['checkpoints'][str(seed)]['path']);load=time.perf_counter()-start
   start=time.perf_counter();prediction=predict_contexts(contexts,action,ids,batch_size=64);duration=time.perf_counter()-start
   if prediction.shape!=(256,576) or not np.isfinite(prediction).all():raise ValueError('finite full-grid material prediction required; no GT fallback')
   path=destination/f'{name}_s{seed}.npz';np.savez_compressed(path,parent=np.array(ids),prediction=prediction)
   rec=dict(method=name,context=arm['context'],seed=seed,identity=identity(path),load_wall_seconds=load,
      prediction_wall_seconds=duration,per_parent_warm_learning_seconds=duration/256,
      per_parent_cold_anchor_seconds=anchor_cost,GT_read=False,new_physics_RHS=0,
      timing_qualification='full/low actual paid path; Native paid prefix counterfactual standalone cost, shared campaign charged once')
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
   phases_identity=identity(phasepath),GT_read=False,freeze_identity=identity(ROOT/'D_FREEZE.json')))
 if all((ROOT/'predictions'/g/'PREDICTIONS_SEALED.json').exists() for g in GROUPS):
  seals=[]
  for g in GROUPS:
   sealpath=ROOT/'predictions'/g/'PREDICTIONS_SEALED.json'
   if json.loads(sealpath.read_text())['freeze_identity']!=identity(ROOT/'D_FREEZE.json'):raise ValueError('both noise conditions must use identical root freeze')
   seals.append(identity(sealpath))
  write(ROOT/'ALL_512_PREDICTIONS_SEALED.json',dict(status='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT',group_seals=seals,freeze_identity=identity(ROOT/'D_FREEZE.json'),GT_read=False,utc=time.time()))

def main():
 p=argparse.ArgumentParser();p.add_argument('operation',choices=['generate','anchors','predict']);p.add_argument('group',choices=list(GROUPS));a=p.parse_args();d=freeze();{'generate':generate,'anchors':anchors,'predict':predict}[a.operation](a.group,d)
if __name__=='__main__':main()
