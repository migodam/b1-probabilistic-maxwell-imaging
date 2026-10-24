#!/usr/bin/env python3
"""ZERO-NEW-PHYSICS paid operator teacher preparation; no NN optimization."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,time,traceback,gc,resource
from collections import Counter
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,save,pack
from p1_material.prior import MaterialPrior
from p3r4.operator_surrogate import (BASELINES,ACTION_METHODS,FEATURE_NAMES,ERROR_SCALE,
 coordinate_derivative,compose_from_transfer,one_secant,ordered_evolution,make_features,
 compare_actions,ranker_parameter_count)

JOB='operator_transport_pack_v1'
DEST=OUT/'results/II/operator_pack_v1'
AUDIT=OUT/'AUDIT/operator_transport'
SOURCES=('dynamic_geometry_first','dynamic_geometry_heun_r1')
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(32))


def split(parent):
 i=int(parent.split('_')[-1]);return 'train' if i<16 else 'validation' if i<24 else 'operator_test'


def append(path,row):
 with Path(path).open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n');f.flush()


def rss():
 a=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;return a if sys.platform=='darwin' else a*1024


def prepare():
 if DEST.exists() or AUDIT.exists():raise FileExistsError('operator outputs exist; no implicit rewrite')
 for job in SOURCES:
  if json.loads((OUT/'III/jobs'/job/'CLOSED.json').read_text())['status']!='COMPLETE':raise ValueError('paid teacher source not closed')
 DEST.mkdir(parents=True);AUDIT.mkdir(parents=True)
 b=Book('II',JOB,max_wall=900,analysis_only=True);started=time.perf_counter();status='COMPLETE';error=None
 inventory=Counter();missing=[];parent_profiles=[];quality=Counter();global_metrics=[]
 try:
  cfg=json.loads((OUT/'configs/TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
  prior=MaterialPrior(cfg,points);sigma=cfg['noise']['sigma']
  protocol=dict(status='FROZEN_BEFORE_LABEL_PREPARATION',parents=PARENTS,
    train_parents=PARENTS[:16],validation_parents=PARENTS[16:24],operator_test_parents=PARENTS[24:],
    exposure='ALL_EXPOSED_DEVELOPMENT_NOT_LOCKED',split_unit='original parent',
    pairing='same paid StageI context anchor to every saved paid-J state, not adjacent path pairs',
    coordinates='probit66 primary; fixed-property Gaussian64 separate pack',
    keep='all available legal teacher states, no Phi/GT/success selection; exact q+x dedup within parent and coordinate',
    current_residual='same exact already-paid candidate Full residual for every action baseline',
    material_action='D_M1*delta;H=I+J.TJ;delta=-solve(H,J.Tr+x1)',
    baselines=ACTION_METHODS,baseline_conditional_costs=dict(anchor_forward_RHS=8,anchor_J_RHS=32,
      candidate_Full_state_RHS=8,teacher_candidate_J_RHS=32),
    new_receiver_T1='MISSING unless explicitly saved; never reconstruct by solving',
    alternate_decomposition='ordered decoder/field/Full-teacher-inferred propagation residual; original T1 ordering MISSING_T1',
    feature_dimension=len(FEATURE_NAMES),feature_names=FEATURE_NAMES,error_scale=ERROR_SCALE,
    feature_no_J1_or_teacher_error=True,model_family='empirical action-error/refreshranker, not DINO replication',
    max_parameters=800000,proposed_parameters=ranker_parameter_count(),NN_optimization=False,
    source_geometry_paid_RHS=459712,source_geometry_paid_states=16060,
    source_J_array_count_expected=10351,source_J_state_forward_plus_J_RHS=414040,
    source_remaining_nonJ_state_RHS=45672,source_cost_owner='R4III dynamic_geometry_first and dynamic_geometry_heun_r1; common contexts R4I separately retained',
    pack_storage='all matched action/delta arrays and Broyden factors; exact FullJ/source pointers retained; deterministic dense-J example per parent/coordinate',
    no_new_complex_RHS=True,no_new_physical_LU=True,real_GN_solves='5 per unique pair; real <=66 parameter systems, not complex Maxwell/current RHS')
  write(DEST/'FROZEN_PROTOCOL.json',protocol)
  for parent in PARENTS:
   tick=time.perf_counter();cachepath=OUT/'results/contexts/cache'/f'{parent}.npz'
   if not cachepath.exists():missing.append(dict(parent=parent,role='base_cache',status='MISSING'));continue
   with np.load(cachepath,allow_pickle=False) as z:
    baseq=z['q'].copy();chi0=z['chi'].copy();T0=z['receiver_T'].copy();E0=z['field'].copy();F0=z['prediction'].copy()
   if T0.shape!=(32,576):missing.append(dict(parent=parent,role='base_receiver_T0',status='MISSING'));continue
   with np.load(OUT/'inputs/online'/f'{parent}.npz',allow_pickle=False) as z:
    if set(z.files)!={'observation','sigma','target_id'} or str(z['target_id'])!=cfg['target_id'] or float(z['sigma'])!=sigma:raise ValueError('observation contract')
    observation=z['observation'].copy()
   save(DEST/'anchors'/f'{parent}.npz',q0=baseq,chi0=chi0,T0=T0,E0=E0,F0=F0,observation=observation,sigma=np.array(sigma))
   for dimension,label in ((66,'probit66'),(64,'field64')):
    basepath=OUT/'results/III/dynamic_geometry_first'/parent/f'd{dimension}_r2/BASE_FRAME.npz'
    if not basepath.exists():missing.append(dict(parent=parent,dimension=dimension,role='base_frame',status='MISSING'));continue
    with np.load(basepath,allow_pickle=False) as z:x0=z['x'].copy();D0=z['Dchi'].copy();J0=z['J'].copy()
    recomposed=compose_from_transfer(T0,E0,D0,sigma)
    base_error=float(np.linalg.norm(recomposed-J0)/max(np.linalg.norm(J0),1e-30))
    if base_error>1e-10:raise ValueError('base T0E0D0/J0 mismatch')
    buffers={key:[] for key in ('features','log_error_target','material_action_error','relative_material_action_error',
      'relative_J_F_error','delta','material_action','broyden_u','broyden_v','term_norms','cross_inner_products',
      'squared_norm_identity_error','telescoping_error','source_J','source_endpoint','source_job','rank','x','q')}
    seen=set();source_count=0;duplicates=0;example_saved=False
    for job in SOURCES:
     source=OUT/'results/III'/job/parent
     for path in sorted(source.glob(f'd{dimension}_r*/derivatives/*.npz')):
      inventory['source_J_arrays']+=1;source_count+=1
      endpoint=OUT/'III/jobs'/job/'endpoints'/path.name
      if not endpoint.exists():missing.append(dict(parent=parent,source=str(path),role='Full_endpoint',status='MISSING'));continue
      with np.load(path,allow_pickle=False) as z:
       x1=z['x'].copy();q1=z['q'].copy();Jq=z['Jq'].copy();J1=z['Jcoordinate'].copy();has_T1='receiver_T' in z.files
      inventory['new_receiver_T1_PRESENT' if has_T1 else 'new_receiver_T1_MISSING']+=1
      key=(x1.tobytes(),q1.tobytes())
      if key in seen:
       duplicates+=1;inventory['exact_duplicate_pair']+=1
       append(DEST/'SOURCE_MANIFEST.jsonl',dict(parent=parent,dimension=dimension,source=str(path),status='EXACT_DUPLICATE_WITHIN_PARENT_COORDINATE',receiver_T1='PRESENT' if has_T1 else 'MISSING'))
       continue
      if x1.shape!=(dimension,) or J1.shape!=(512,dimension) or Jq.shape!=(512,66) or not np.isfinite(J1).all():raise ValueError('saved teacher layout')
      expected=Jq[:,:64] if dimension==64 else Jq*np.r_[np.ones(64),np.exp(-.5*x1[64:]**2)/np.sqrt(2*np.pi)][None,:]
      transform_error=float(np.linalg.norm(expected-J1)/max(np.linalg.norm(J1),1e-30))
      if transform_error>1e-12:raise ValueError('saved coordinate-J/probit chain mismatch')
      with np.load(endpoint,allow_pickle=False) as z:
       if not np.array_equal(z['q'],q1):raise ValueError('endpoint and J candidate mismatch')
       chi1=z['chi'].copy();E1=z['field'].copy();F1=z['prediction'].copy();backward=float(z['backward_residual'])
      if dimension==64 and not np.array_equal(q1[64:],baseq[64:]):raise ValueError('fixed-u64 property drift')
      if not prior.legal(q1):raise ValueError('illegal saved candidate')
      D1=coordinate_derivative(prior,q1,x1,dimension)
      Jdecoder=compose_from_transfer(T0,E0,D1,sigma)
      Jfield=compose_from_transfer(T0,E1,D1,sigma)
      h=x1-x0;df=pack(F1-F0)*np.sqrt(2)/sigma
      Jb,u,v=one_secant(J0,h,df);r=pack(F1-observation)*np.sqrt(2)/sigma
      comparisons=compare_actions(J0,Jdecoder,Jfield,Jb,J1,r,x1,D1)
      evolution=ordered_evolution(J0,Jdecoder,Jfield,J1)
      features=make_features(prior=prior,q0=baseq,q1=q1,x0=x0,x1=x1,dimension=dimension,
        D0=D0,D1=D1,E0=E0,E1=E1,J0=J0,Jdecoder=Jdecoder,Jfield=Jfield,Jb=Jb,
        F0=F0,F1=F1,observation=observation,sigma=sigma)
      seen.add(key);inventory['unique_pairs']+=1;quality['real_GN_solves']+=5
      for k in ('features',):buffers[k].append(features.astype(np.float32))
      for k in ('log_error_target','material_action_error','relative_material_action_error','relative_J_F_error'):buffers[k].append(comparisons[k])
      for k in ('delta','material_action'):buffers[k].append(comparisons[k].astype(np.float32))
      buffers['broyden_u'].append(u.astype(np.float32));buffers['broyden_v'].append(v.astype(np.float32))
      for k in ('term_norms','cross_inner_products','squared_norm_identity_error','telescoping_error'):buffers[k].append(evolution[k])
      buffers['source_J'].append(str(path));buffers['source_endpoint'].append(str(endpoint));buffers['source_job'].append(job)
      rank=int(path.parent.parent.name.split('_r')[-1]);buffers['rank'].append(rank);buffers['x'].append(x1);buffers['q'].append(q1)
      append(DEST/'SOURCE_MANIFEST.jsonl',dict(parent=parent,dimension=dimension,split=split(parent),rank=rank,
        source=str(path),endpoint=str(endpoint),status='USED',receiver_T1='PRESENT' if has_T1 else 'MISSING',
        coordinate_transform_relative_error=transform_error,base_recomposition_relative_error=base_error,backward_residual=backward,
        exact_D1='ANALYTIC_RECOMPUTED_ZERO_PHYSICS',original_T1_ordering='MISSING_T1',ordered_propagation='FULL_TEACHER_INFERRED_NOT_ONLINE_FEATURE'))
      if not example_saved:
       save(DEST/'matched_dense_examples'/f'{parent}_d{dimension}.npz',baseline_names=np.array(ACTION_METHODS),
        frozenJ=J0,decoder_exact=Jdecoder,field_updated=Jfield,one_secant_Broyden=Jb,FullJ=J1,
        D1=D1,x0=x0,x1=x1,residual=r,delta=comparisons['delta'],material_action=comparisons['material_action'],
        material_action_error=comparisons['material_action_error'],evolution_term_norms=evolution['term_norms'],evolution_cross_inner_products=evolution['cross_inner_products'])
       example_saved=True
      if rss()>3*1024**3:raise MemoryError('operator pack3GB memory cap')
    if buffers['features']:
     arrays={key:np.asarray(value) for key,value in buffers.items()}
     save(DEST/label/f'{parent}.npz',parent=np.array(parent),split=np.array(split(parent)),dimension=np.array(dimension),
       baseline_names=np.array(ACTION_METHODS),feature_names=np.array(FEATURE_NAMES),**arrays)
     errors=np.asarray(buffers['material_action_error']);relative=np.asarray(buffers['relative_J_F_error'])
     summary=dict(parent=parent,split=split(parent),dimension=dimension,source_records=source_count,unique_pairs=len(errors),exact_duplicates=duplicates,
       base_recomposition_relative_error=base_error,material_action_error_mean=errors.mean(axis=0).tolist(),material_action_error_median=np.median(errors,axis=0).tolist(),
       relative_J_F_error_mean=relative.mean(axis=0).tolist(),maximum_telescoping_error=max(buffers['telescoping_error']),
       maximum_squared_norm_identity_error=max(buffers['squared_norm_identity_error']),FullJ1='PAID_PERSISTED',T1='MISSING',
       D1='ANALYTIC_WITH_CORRECT_COORDINATE_CHAIN',new_complex_RHS=0,new_physical_LU=0)
     parent_profiles.append(summary);write(AUDIT/f'{parent}_d{dimension}.json',summary)
    del buffers;gc.collect()
   write(AUDIT/'PROGRESS.json',dict(completed_parent=parent,inventory=dict(inventory),profiles=len(parent_profiles),elapsed_seconds=time.perf_counter()-started))
   print(json.dumps(dict(parent=parent,inventory=dict(inventory),elapsed_seconds=time.perf_counter()-started)),flush=True)
  manifest=dict(status='READY_NO_NN_OPTIMIZATION',inventory=dict(inventory),missing=missing,profiles=parent_profiles,
    coordinates=dict(primary='probit66',separate_conditional='field64'),feature_dimension=len(FEATURE_NAMES),parameter_count=ranker_parameter_count(),
    train_parents=PARENTS[:16],validation_parents=PARENTS[16:24],operator_test_parents=PARENTS[24:],
    exposure='ALL_EXPOSED_DEVELOPMENT_NOT_LOCKED',feature_precision='float32; baseline error metrics float64; FullJ original float64 source retained',
    matched_actions='5 methods per unique pair; delta/action arrays float32 for training, exact float64 dense examples and source references',
    new_complex_RHS=0,new_Maxwell_states=0,new_physical_LU=0,real_GN_solves=quality['real_GN_solves'],
    original_geometry_acquisition_RHS=459712,original_J_source_acquisition_RHS=414040,
    original_StageI_context_cost='retained separately in R4I ledger, not zero',
    Full_receiver_T1='MISSING; no solve or pseudoinverse reconstruction',
    decomposition='ordered telescope with Full-teacher-inferred propagation residual; not independent T1-based audit',
    NN_optimization=False,scientific_scope='local Gaussianized GN material action fidelity, not material GT risk or posterior mixing',
    elapsed_seconds=time.perf_counter()-started,RSS_peak_bytes=rss())
  write(DEST/'MANIFEST.json',manifest);write(AUDIT/'RESULTS.json',manifest)
  b.counts['source_J_arrays']=inventory['source_J_arrays'];b.counts['unique_operator_pairs']=inventory['unique_pairs'];b.counts['real_parameter_GN_solves']=quality['real_GN_solves']
 except BaseException as exc:
  status='ERROR';error=repr(exc);write(AUDIT/'ERROR.json',dict(error=error,traceback=traceback.format_exc(),inventory=dict(inventory)));raise
 finally:
  b.close(status,error);write(AUDIT/'FINAL_JOB_COST.json',dict(status=status,error=error,elapsed_seconds=time.perf_counter()-started,
    new_complex_RHS=0,new_Maxwell_states=0,new_physical_LU=0,real_GN_solves=quality['real_GN_solves'],receipt=str(b.jobdir/'CLOSED.json')))

if __name__=='__main__':prepare()
