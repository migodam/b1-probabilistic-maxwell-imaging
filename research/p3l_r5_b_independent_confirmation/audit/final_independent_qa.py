"""Authorized offline independent mechanical QA; no solvers or model imports."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib,collections,subprocess,time
import numpy as np
S=Path(__file__).resolve().parents[1];ROOT=S.parents[1];AUD=S/'audit';start=time.perf_counter()
def rd(p):return json.loads(p.read_text())
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
snap=rd(S/'RESOURCE_LEDGER_SNAPSHOT.json');total=collections.Counter();actual=collections.Counter();active={};seq=0;jobends=0;gt_flags=[];jobs={}
with (S/'RESOURCE_LEDGER.jsonl').open() as f:
 for line in f:
  r=json.loads(line);seq+=1;assert r['sequence']==seq
  for k in ['RHS','state_attempts','LU_attempts']:total[k]+=r.get(k,0)
  total['forward_RHS' if r.get('kind')=='forward' else 'derivative_RHS']+=r.get('RHS',0)
  assert r['cumulative_RHS']==total['RHS']
  if r['event']=='PHYSICS_BEGIN':assert not active;active[r['action']]=r
  if r['event']=='PHYSICS_END':
   assert r['action'] in active;active.pop(r['action']);total['physics_wall_seconds']+=r['wall_seconds'];total['physics_CPU_seconds']+=r['CPU_seconds'];total['failures']+=r['error'] is not None
  if r['event']=='JOB_END':
   jobends+=1;actual.update(r['actual_counters']);total['job_wall_seconds']+=r['wall_seconds'];total['job_CPU_seconds']+=r['CPU_seconds'];assert r['unused_prepaid_forward_RHS']==0
  jobs.setdefault(r['job'],collections.Counter())['RHS']+=r.get('RHS',0)
ledger=dict(sequence=seq,job_end_count=jobends,active_actions=active,totals=dict(total),snapshot_differences={k:float(total[k])-float(v) for k,v in snap.items() if isinstance(v,(int,float)) and k!='sequence' and abs(float(total[k])-float(v))>1e-8},counter_differences={k:actual[k]-v for k,v in snap['actual_counters'].items() if abs(actual[k]-v)>1e-8})
ledger['actual_solver_checks']=dict(actual_full_forward_RHS=actual['actual_full_forward_RHS'],actual_derivative_RHS=actual['actual_derivative_RHS'],linear_solve_rhs=actual['linear_solve_rhs'],LU=actual['LU_factorizations'],full_forward_calls=actual['full_forward_calls'],reserved_RHS=total['RHS'])
assert not ledger['snapshot_differences'] and not ledger['counter_differences'] and not active
assert actual['linear_solve_rhs']==total['RHS']==1999504 and actual['actual_full_forward_RHS']==total['forward_RHS'] and actual['actual_derivative_RHS']==total['derivative_RHS']
parents=[]
for g in ['A','B']:parents += [r['parent'] for r in rd(S/f'configs/parents_{g}.json')['parents']]
sysroot=ROOT/'src';import sys;sys.path.insert(0,str(sysroot));from p1_material.prior import MaterialPrior
cfg=rd(ROOT/'research/p3l_r4_conditional_transport/configs/TARGET.json');geom=rd(ROOT/'research/p1_material_target_20261018/inputs/geometry.json');prior=MaterialPrior(cfg,geom['points'])
contexts=[];native=[];backmax=0.;decmax=0.;probit_success=0;native_failure=collections.Counter()
for p in parents:
 d=rd(S/f'raw/contexts/{p}.json');assert d['GT_read'] is False and d['native_cap']==96 and d['probit_cap']==64
 n=d['native_report'];native_failure[n['status']]+=1;assert n['fresh_forward_states']<=96
 assert d['probit_stationarity']['fresh_state_attempts']<=64
 probit_success+=d['probit_stationarity']['optimizer'].get('success') is True
 with np.load(S/f'raw/contexts/{p}.npz',allow_pickle=False) as z:
  q=z['q0'];chi=z['chi0'];nchi=z['native_chi'];nq=z['native_q'];assert prior.legal(q) and prior.legal(nq);decmax=max(decmax,float(np.max(np.abs(prior.decode(q)[0]-chi))),float(np.max(np.abs(prior.decode(nq)[0]-nchi))));native.append(nchi.copy());assert np.isfinite(z['J']).all() and z['J'].shape==(512,66)
 with np.load(S/f'raw/physics/anchors_{p}/details/probit/RAW_MODE_CACHE.npz',allow_pickle=False) as z:backmax=max(backmax,float(z['backward_residual']))
 assert d['new_RHS']==jobs['anchors_'+p]['RHS']
 contexts.append(d['new_RHS'])
assert decmax==0 and backmax<=1e-10
anchor=dict(parents=512,all_GT_read_false=True,all_native96_probit64_caps=True,decoder_max_absolute_error=decmax,backward_residual_max=backmax,native_status_counts=dict(native_failure),probit_SciPy_success=probit_success,context_RHS=sum(contexts),context_RHS_min=min(contexts),context_RHS_max=max(contexts))
# Algebra only on two fixed first parents, no factorization or solve.
algebra=[]
with np.load(ROOT.parent/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/inputs/background/weak.npz',allow_pickle=False) as z:G=z['G'];SS=z['S'];inc=z['incident']
for p in [parents[0],parents[256]]:
 with np.load(S/f'raw/physics/anchors_{p}/details/probit/RAW_MODE_CACHE.npz',allow_pickle=False) as z:
  chi=z['chi'];cur=z['current'];field=z['field'];pred=z['prediction'];L=np.eye(576)-chi[:,None]*G;rhs=chi[:,None]*inc
  rel=lambda x,y:float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-30))
  r=dict(parent=p,current_equation_relative=rel(L@cur,rhs),field_equation_relative=rel(inc+G@cur,field),prediction_equation_relative=rel(SS@cur,pred));assert max(v for k,v in r.items() if k!='parent')<=1e-10;algebra.append(r)
anchor['fixed_two_algebra']=algebra;anchor['algebra_new_RHS']=0;anchor['algebra_scope']='matrix applications only; no solve/factorization'
res=rd(S/'RESULTS_AND_STATISTICS.json')
with np.load(S/'raw/offline_scores/OFFLINE_SCOREPACK.npz',allow_pickle=False) as z:
 ids=z['parent_ids'].astype(str).tolist();methods=z['methods'].astype(str).tolist();truth=z['truth_chi'];raw=z['prediction_raw'];clipped=z['prediction_clipped'];assert ids==parents
 assert raw.shape==(512,3,len(methods),576)
 mse=np.mean(abs(raw-truth[:,None,None,:])**2,axis=-1);mae=np.mean(abs(raw-truth[:,None,None,:]),axis=-1);cmse=np.mean(abs(clipped-truth[:,None,None,:])**2,axis=-1)
 ensemble=np.mean(abs(raw.mean(axis=1)-truth[:,None,:])**2,axis=-1)
 metrics=dict(methods=methods,raw_seed_mean_MSE={m:float(mse[:,:,i].mean()) for i,m in enumerate(methods)},raw_seed_mean_MAE={m:float(mae[:,:,i].mean()) for i,m in enumerate(methods)},raw_ensemble_MSE={m:float(ensemble[:,i].mean()) for i,m in enumerate(methods)},max_difference_published_means=max(abs(float(mse[:,:,i].mean())-res['raw_primary']['method_means'][m]) for i,m in enumerate(methods[:5])))
 with np.load(S/'raw/offline_scores/PARENT_LOSSES.npz',allow_pickle=False) as p:metrics['max_parent_seed_raw_loss_difference']=float(np.max(abs(mse[:,:,:5]-p['raw_losses'])));metrics['max_parent_seed_clipped_loss_difference']=float(np.max(abs(cmse[:,:,:5]-p['clipped_losses'])))
 with np.load(S/'raw/offline_scores/PARENT_METRICS.npz',allow_pickle=False) as p:metrics['max_parent_seed_MAE_difference']=float(np.max(abs(mae-p['raw_mae_full'])))
 # Native baseline descriptive only, not a new primary contrast.
 native=np.asarray(native);nmse=np.mean(abs(native-truth)**2,axis=1);nmae=np.mean(abs(native-truth),axis=1)
 np.savez_compressed(AUD/'NATIVE_BASELINE_DESCRIPTIVE.npz',parent_ids=np.array(ids),mse_full=nmse,mae_full=nmae)
 metrics['native_baseline_descriptive']=dict(MSE=float(nmse.mean()),MAE=float(nmae.mean()),scope='additional descriptive baseline; no new primary hypothesis')
 # Re-read independently sealed predictions, check pack alignment and identity.
 predmax=0.;artifacts=0
 global_seal=rd(S/'raw/predictions/ALL_PREDICTIONS_SEALED.json')
 for a in global_seal['artifacts']:
  p=Path(a['path']);assert p.stat().st_size==a['bytes'];assert hashlib.blake2b(p.read_bytes()).hexdigest()==a['blake2b'];artifacts+=1
 for gi,g in enumerate(['A','B']):
  for mi,m in enumerate(methods):
   for si,seed in enumerate([2026102601,2026102602,2026102603]):
    sd=2026102601 if m=='physics_probit' else seed
    with np.load(S/f'raw/predictions/{g}/{m}_s{sd}.npz',allow_pickle=False) as p:predmax=max(predmax,float(np.max(abs(p['prediction_raw']-raw[gi*256:(gi+1)*256,si,mi]))))
 metrics['sealed_artifacts_checked']=artifacts;metrics['max_sealed_prediction_pack_difference']=predmax
means=metrics['raw_seed_mean_MSE'];gain=means['strong_cnn']-means['k32'];ret=(means['strong_cnn']-means['distilled_w128'])/gain;retrow=res['raw_primary']['comparisons']['retention_80pct'];metrics['retention']=dict(point_fraction=ret,margin=means['strong_cnn']-means['distilled_w128']-.8*gain,original_adjusted_one_sided_lower=retrow['lower_one_sided_98_75'],confidence_gate=retrow['confidence_gate'],status='HOLD_POINT_ABOVE_80_NOT_STATISTICAL_PASS')
assert ret>.8 and not retrow['confidence_gate'] and retrow['lower_one_sided_98_75']<0
sources=[]
for r in rd(S/'configs/FREEZE.json')['files']:
 p=Path(r['path'])
 if '/src/' in str(p):sources.append(dict(path=str(p),bytes_match=p.stat().st_size==r['bytes'],freeze_identity_match=hashlib.blake2b(p.read_bytes()).hexdigest()==r['blake2b']))
old=[]
for name in ['b1-p3l-r4-conditional-transport','b1-p3l-r5-parallel-ab']:
 p=ROOT.parent/name
 old.append(dict(path=str(p),HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=p,text=True).strip(),tracked_diff=subprocess.check_output(['git','diff','--stat'],cwd=p,text=True).strip(),status_porcelain=subprocess.check_output(['git','status','--porcelain'],cwd=p,text=True).strip(),scope='current snapshot; no pre-run git status baseline recovered by this worker'))
qa=dict(status='PASS_MECHANICAL_QA_WITH_RETAINED_SCIENTIFIC_HOLD',ledger=ledger,anchors=anchor,metrics=metrics,frozen_sources=sources,old_worktrees=old,new_physics_RHS=0,new_training=0,offline_truth_read_authorized=True,new_primary_hypotheses=0,wall_seconds=time.perf_counter()-start)
assert all(r['bytes_match'] and r['freeze_identity_match'] for r in sources)
dump(AUD/'FINAL_INDEPENDENT_QA.json',qa)
print(json.dumps(dict(status=qa['status'],anchor=anchor,metrics=metrics,old=old,source_match=len(sources),wall=qa['wall_seconds']),indent=2))
