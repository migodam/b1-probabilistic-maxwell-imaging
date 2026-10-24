#!/usr/bin/env python3
"""ZERO-PHYSICS forward/JVP data preparation; no optimizer or neural model fit."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,time,traceback,subprocess,gc
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,save
from p3r4.forward_jvp_surrogate import (V_SEED,SEEDS,whiten,fixed_directions,input_tangent_matrix,parameter_count)

JOB='forward_jvp_prepare_v1'
SOURCE=OUT/'results/II/operator_pack_v1'
DEST=OUT/'results/II/forward_jvp_pack_v1'
AUDIT=OUT/'AUDIT/forward_jvp_surrogate'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(32))


def prepare():
 if DEST.exists() or AUDIT.exists():raise FileExistsError('forward/JVP preparation already exists')
 if json.loads((SOURCE/'MANIFEST.json').read_text())['status']!='READY_NO_NN_OPTIMIZATION':raise ValueError('paid operator source pack required')
 DEST.mkdir(parents=True);AUDIT.mkdir(parents=True)
 b=Book('II',JOB,max_wall=900,analysis_only=True);started=time.perf_counter();status='COMPLETE';error=None;profiles=[]
 try:
  protocol=dict(status='PREPARE_ONLY_NO_OPTIMIZER',target='whitened nonlinear response remainder R=f1-f0-J0*h',
    input='x0,h,f0_white,residual0_white,J0*h; no queried candidate Full quantity',
    coordinates='full66 probit primary; fixed-property64 separate',input_dim66=1668,input_dim64=1664,
    architecture='256LayerNormSiLU ->256LayerNormSiLU ->512; anchored net(h)-net(0)',parameters66=parameter_count(66),parameters64=parameter_count(64),cap_parameters=800000,
    derivative='candidate raw h only; x0,f0,residual0,J0 fixed context; input/output normalizer chain rule explicit',
    JVP_target='(true_paid_J1-J0)*V',JVP_directions=2,V_seed=V_SEED,V_distribution='fixed Rademacher +/-1/sqrt(d), no QR or re-draw',
    controls=dict(value_only=dict(value_weight=1,derivative_weight=0),derivative_informed=dict(value_weight=1,derivative_weight=1,JVP_directions=2)),
    seeds=SEEDS,initial_updates=1000,maximum_updates=3000,batch=64,optimizer='AdamW lr3e-4 weight_decay1e-5',
    normalizers='parent-balanced TRAIN only; input mean/std floor1e-6; output RMS about zero floor1e-4; validation/test frozen',
    folds=dict(train=PARENTS[:16],validation=PARENTS[16:24],operator_test=PARENTS[24:]),exposure='ALL_EXPOSED_DEVELOPMENT_NOT_LOCKED',
    validation_gates=dict(max_whitened_L2_response_error=.05,max_absolute_Phi_error=.05),
    gate_scope='empirical method screening, not pointwise certificates; failed gate implies Full fallback8RHS/candidate, never forward saving claim',
    observation_potential='.5*||fhat_white-y_white||², original sigma/packing; no density/posterior correctness inferred',
    source_teacher_geometry_RHS=459712,source_J_state_RHS=414040,source_context_owner='R4I paid96+64 contexts; retained separately',
    cold_per_pair_label_cost=dict(value_only_forward_at_given_anchor=8,FullJ_derivative_teacher_at_given_anchor=40),
    actual_data_reuse='no new RHS; shared original costs retained, not free cold labels',T1='MISSING_IRRELEVANT_FOR_SAVED_JVP_TARGET',
    inference_timing='warm single-candidate context assembly/normalization/anchored forward/reconstruction; cold load and Fullfallback separate',
    inference_wall_measured='NOT_YET_UNTRAINED',new_physics_RHS=0,NN_optimization=0,GT_read=False)
  write(DEST/'FROZEN_PROTOCOL.json',protocol)
  for d,label in ((66,'probit66'),(64,'field64')):
   V=fixed_directions(d);train_moments=[];train_second=[];train_output_second=[]
   for parent in PARENTS:
    with np.load(SOURCE/label/f'{parent}.npz',allow_pickle=False) as z:
     xs=z['x'].copy();sources_J=z['source_J'].astype(str);sources_ep=z['source_endpoint'].astype(str);split=str(z['split'])
    with np.load(SOURCE/'anchors'/f'{parent}.npz',allow_pickle=False) as z:
     F0=z['F0'].copy();y=z['observation'].copy();sigma=float(z['sigma'])
    with np.load(OUT/'results/III/dynamic_geometry_first'/parent/f'd{d}_r2/BASE_FRAME.npz',allow_pickle=False) as z:x0=z['x'].copy();J0=z['J'].copy()
    h=xs-x0;linear=h@J0.T;f0=whiten(F0,sigma);yw=whiten(y,sigma);r0=f0-yw
    responses=[];jvps=[];phi=[]
    for Jpath,Epath in zip(sources_J,sources_ep):
     with np.load(Jpath,allow_pickle=False) as z:J1=z['Jcoordinate'].copy()
     with np.load(Epath,allow_pickle=False) as z:prediction=z['prediction'].copy()
     f1=whiten(prediction,sigma);responses.append(f1);jvps.append((J1-J0)@V)
     diff=prediction-y;phi.append(float(np.vdot(diff,diff).real/sigma**2))
    responses=np.asarray(responses);jvps=np.asarray(jvps);target=responses-f0-linear
    save(DEST/label/f'{parent}.npz',parent=np.array(parent),split=np.array(split),dimension=np.array(d),
      x0=x0,h=h,f0_white=f0,residual0_white=r0,y_white=yw,J0=J0,V=V,
      nonlinear_remainder_target=target.astype(np.float32),response_target_white=responses,
      JVP_remainder_target=jvps.astype(np.float32),Phi_true=np.array(phi),source_J=sources_J,source_endpoint=sources_ep)
    if split=='train':
     n=len(h);raw=np.c_[np.broadcast_to(x0,(n,d)),h,np.broadcast_to(f0,(n,512)),np.broadcast_to(r0,(n,512)),linear]
     train_moments.append(raw.mean(axis=0));train_second.append((raw*raw).mean(axis=0));train_output_second.append((target*target).mean(axis=0))
    profiles.append(dict(parent=parent,dimension=d,split=split,pairs=len(h),queries_new=0,source_saved_JVP=True))
    del responses,jvps,target;gc.collect()
   mean=np.mean(train_moments,axis=0);second=np.mean(train_second,axis=0);scale=np.sqrt(np.maximum(second-mean*mean,1e-12));out_scale=np.sqrt(np.maximum(np.mean(train_output_second,axis=0),1e-8))
   save(DEST/label/'TRAIN_NORMALIZATION.npz',input_mean=mean,input_scale=scale,output_scale=out_scale,V=V,
      train_parents=np.array(PARENTS[:16]),V_seed=np.array(V_SEED))
   with np.load(DEST/label/f'{PARENTS[0]}.npz',allow_pickle=False) as z:J0=z['J0'].copy()
   save(AUDIT/f'CHAIN_RULE_MATRICES_d{d}.npz',input_raw_tangent_matrix=input_tangent_matrix(J0),
      input_normalized_tangent_matrix=input_tangent_matrix(J0)/scale[:,None],V=V,
      normalized_input_direction_matrix=(input_tangent_matrix(J0)@V)/scale[:,None],output_scale=out_scale)
   print(json.dumps(dict(coordinates=d,pairs=sum(p['pairs'] for p in profiles if p['dimension']==d),parameters=parameter_count(d))),flush=True)
  test=subprocess.run([sys.executable,'-m','unittest','discover','-s','tests/p3r4','-p','test_forward_jvp_algebra.py','-v'],
    cwd=ROOT,env={**os.environ,'PYTHONPATH':str(ROOT/'src')},capture_output=True,text=True)
  (AUDIT/'ALGEBRA_TEST_LOG.txt').write_text(test.stdout+test.stderr)
  if test.returncode:raise RuntimeError('packing/JVP algebra check failed')
  manifest=dict(status='READY_NO_NN_OPTIMIZATION',profiles=profiles,total_pairs=sum(z['pairs'] for z in profiles),
    physical_RHS_new=0,Maxwell_states_new=0,physical_LU_new=0,NN_optimizer_steps=0,
    teacher_geometry_owner_RHS=459712,teacher_J_state_owner_RHS=414040,
    normalization='TRAIN only, equally weighted parent moments',V_seed=V_SEED,heldout_folds_exposed=True,
    module='src/p3r4/forward_jvp_surrogate.py',trainer='tools/train_forward_jvp.py',tests=3,
    inference_wall='NOT_MEASURED_UNTRAINED',elapsed_seconds=time.perf_counter()-started)
  write(DEST/'MANIFEST.json',manifest);write(AUDIT/'PREPARATION.json',manifest)
 except BaseException as exc:status='ERROR';error=repr(exc);write(AUDIT/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:b.close(status,error)

if __name__=='__main__':prepare()
