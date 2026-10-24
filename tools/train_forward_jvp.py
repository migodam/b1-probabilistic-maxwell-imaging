#!/usr/bin/env python3
"""Bounded value-only vs derivative-informed local response controls.

Preparation and --help do not optimize. Training requires explicit --train,
mode/seed/job. Inference uses fixed anchor and candidate delta only. Evaluation
uses already-paid Full/JVP labels, reports white/Phi errors and single-query
compute wall. It never silently executes a physical Full fallback.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,argparse,time,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,BudgetStop
from p3r4.forward_jvp_surrogate import (SEEDS,build_local_forward_model,value_jvp_loss,
 residual_jvp,parameter_count)
PACK=OUT/'results/II/forward_jvp_pack_v1'


def names(split):return [f'R4DEV_{i:05d}' for i in (range(16) if split=='train' else range(16,24) if split=='validation' else range(24,32))]


def load_data(d,split):
 folder=PACK/('probit66' if d==66 else 'field64');rows=[]
 for parent in names(split):
  with np.load(folder/f'{parent}.npz',allow_pickle=False) as z:
   if str(z['split'])!=split or int(z['dimension'])!=d:raise ValueError('fixed fold/coordinate mismatch')
   data={key:z[key].copy() for key in ('x0','h','f0_white','residual0_white','y_white','J0',
     'nonlinear_remainder_target','response_target_white','JVP_remainder_target','Phi_true','V')}
  data['parent']=parent;rows.append(data)
 return rows


def context_tensors(data,torch,batch=1):
 return dict(x0=torch.from_numpy(data['x0'].astype(np.float32))[None,:].expand(batch,-1),
  f0=torch.from_numpy(data['f0_white'].astype(np.float32))[None,:].expand(batch,-1),
  r0=torch.from_numpy(data['residual0_white'].astype(np.float32))[None,:].expand(batch,-1),
  J0=torch.from_numpy(data['J0'].astype(np.float32))[None,:,:].expand(batch,-1,-1))


def evaluate_groups(model,data,torch,*,jvp_diagnostics=False):
 model.eval();parent_rows=[];predictions={};maximum_white=0.;maximum_phi=0.;all_times=[]
 for p in data:
  tick=time.perf_counter();context=context_tensors(p,torch);context_wall=time.perf_counter()-tick
  values=[];times=[];jerrors=[]
  for i,h_raw in enumerate(p['h']):
   tick=time.perf_counter();h=torch.from_numpy(h_raw.astype(np.float32))[None,:]
   with torch.no_grad():remainder=model(h,context['x0'],context['f0'],context['r0'],context['J0']).numpy()[0].astype(np.float64)
   # Mixed precision is explicit: exact cached anchor/linear reconstruction
   # in float64, learned remainder in float32; no new Full datum is read.
   prediction=p['f0_white']+p['J0']@h_raw+remainder
   times.append(time.perf_counter()-tick);values.append(prediction)
   if jvp_diagnostics:
    errors=[]
    for k in range(2):
     direction=torch.from_numpy(p['V'][:,k].astype(np.float32))[None,:]
     got=residual_jvp(model,h,context,direction,create_graph=False).detach().numpy()[0]
     errors.append(float(np.linalg.norm(got-p['JVP_remainder_target'][i,:,k])))
    jerrors.append(errors)
  values=np.asarray(values);errors=np.linalg.norm(values-p['response_target_white'],axis=1)
  predicted_phi=.5*np.sum((values-p['y_white'])**2,axis=1);phi_error=predicted_phi-p['Phi_true']
  frozen=p['f0_white']+p['h']@p['J0'].T
  frozen_error=np.linalg.norm(frozen-p['response_target_white'],axis=1)
  frozen_phi_error=.5*np.sum((frozen-p['y_white'])**2,axis=1)-p['Phi_true']
  maximum_white=max(maximum_white,float(errors.max()));maximum_phi=max(maximum_phi,float(abs(phi_error).max()));all_times+=times
  parent_rows.append(dict(parent=p['parent'],pairs=len(values),white_MSE=float(np.mean(errors**2)/512),white_L2_mean=float(errors.mean()),white_L2_max=float(errors.max()),
   absolute_Phi_error_mean=float(abs(phi_error).mean()),absolute_Phi_error_max=float(abs(phi_error).max()),
   single_query_compute_seconds_mean=float(np.mean(times)),single_query_compute_seconds_max=float(np.max(times)),warm_context_setup_seconds=context_wall,
   JVP_error_mean=float(np.mean(jerrors)) if jerrors else None,
   frozenJ_white_L2_mean=float(frozen_error.mean()),frozenJ_white_L2_max=float(frozen_error.max()),
   frozenJ_absolute_Phi_error_mean=float(abs(frozen_phi_error).mean()),frozenJ_absolute_Phi_error_max=float(abs(frozen_phi_error).max())))
  predictions[p['parent']]=dict(prediction_white=values,white_L2_error=errors,Phi_prediction=predicted_phi,Phi_error=phi_error,JVP_L2_error=np.asarray(jerrors),
   frozenJ_prediction_white=frozen,frozenJ_white_L2_error=frozen_error,frozenJ_Phi_error=frozen_phi_error)
 gate=bool(maximum_white<=.05 and maximum_phi<=.05)
 summary=dict(parent_balanced_white_MSE=float(np.mean([z['white_MSE'] for z in parent_rows])),max_white_L2_error=maximum_white,max_absolute_Phi_error=maximum_phi,
  empirical_validation_gate=gate,gate_scope='empirical screening, not a certified bound',
  single_query_compute_seconds_mean=float(np.mean(all_times)),single_query_compute_seconds_max=float(np.max(all_times)),
  timing='CPU single candidate, warm model/anchor; tensor delta assembly+normalized anchored net+float64 reconstruction; labels/error analysis and cold load separate',
  Full_fallback_required_if_gate_fails=not gate,Full_fallback_complex_RHS_per_candidate=8,
  counterfactual_all_Full_fallback_RHS=8*sum(z['pairs'] for z in parent_rows) if not gate else None,
  actual_evaluation_physics_RHS=0,acceleration_claim=False,posterior_correctness_inferred=False,parents=parent_rows)
 return summary,predictions


def source_snapshot(book):
 # Book's generic tool glob need not cover this differently named entrypoint.
 target=book.jobdir/'source_snapshot/tools'/Path(__file__).name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(Path(__file__).read_bytes())
 write(book.jobdir/'ADDITIONAL_SOURCE_SNAPSHOT.json',dict(entrypoint=str(Path(__file__).resolve()),snapshot=str(target),bytes=target.stat().st_size,phase='before any optimizer/model fitting',checksum_checks=False))


def train(args):
 dest=OUT/'results/II/forward_jvp_training'/args.job
 if dest.exists():raise FileExistsError('frozen training output exists')
 b=Book('II',args.job,max_wall=args.max_wall,analysis_only=True,kind='training_cpu');status='COMPLETE';error=None;started=time.perf_counter()
 try:
  source_snapshot(b);dest.mkdir(parents=True)
  import torch
  torch.set_num_threads(1);torch.manual_seed(args.seed)
  proto=json.loads((PACK/'FROZEN_PROTOCOL.json').read_text())
  training=load_data(args.coordinates,'train');validation=load_data(args.coordinates,'validation')
  normpath=PACK/('probit66' if args.coordinates==66 else 'field64')/'TRAIN_NORMALIZATION.npz'
  with np.load(normpath,allow_pickle=False) as z:normalization={key:z[key].copy() for key in ('input_mean','input_scale','output_scale','V')}
  model=build_local_forward_model(args.coordinates,normalization['input_mean'],normalization['input_scale'],normalization['output_scale'])
  optimizer=torch.optim.AdamW(model.parameters(),lr=.0003,weight_decay=.00001);rng=np.random.default_rng(args.seed)
  derivative_weight=1. if args.mode=='jvp' else 0.;best=np.inf;best_update=0
  write(dest/'RUN_PROTOCOL.json',dict(mode=args.mode,seed=args.seed,coordinates=args.coordinates,updates=args.updates,batch=64,
    derivative_weight=derivative_weight,parameters=sum(p.numel() for p in model.parameters()),source_protocol=proto,
    training_physics_RHS=0,teacher_geometry_owner_RHS=459712,operator_test_read=False,device='cpu'))
  def checkpoint(update):
   nonlocal best,best_update
   metrics,_=evaluate_groups(model,validation,torch)
   with (dest/'VALIDATION.jsonl').open('a') as f:f.write(json.dumps(dict(update=update,**metrics))+'\n')
   if metrics['parent_balanced_white_MSE']<best:
    best=metrics['parent_balanced_white_MSE'];best_update=update
    torch.save(dict(model=model.state_dict(),normalization=normalization,coordinates=args.coordinates,mode=args.mode,seed=args.seed,
      selected_update=update,validation=metrics,protocol=proto,operator_test_read=False),dest/'BEST.pt')
   model.train()
  checkpoint(0)
  for update in range(1,args.updates+1):
   choices=rng.integers(0,len(training),size=64);items=[]
   for choice in choices:
    p=training[int(choice)];i=int(rng.integers(len(p['h'])));items.append((p,i))
   h=torch.from_numpy(np.stack([p['h'][i] for p,i in items]).astype(np.float32))
   context=dict(x0=torch.from_numpy(np.stack([p['x0'] for p,i in items]).astype(np.float32)),
    f0=torch.from_numpy(np.stack([p['f0_white'] for p,i in items]).astype(np.float32)),
    r0=torch.from_numpy(np.stack([p['residual0_white'] for p,i in items]).astype(np.float32)),
    J0=torch.from_numpy(np.stack([p['J0'] for p,i in items]).astype(np.float32)))
   target=torch.from_numpy(np.stack([p['nonlinear_remainder_target'][i] for p,i in items]))
   jtarget=torch.from_numpy(np.stack([p['JVP_remainder_target'][i] for p,i in items]))
   V=torch.from_numpy(normalization['V'].astype(np.float32))
   optimizer.zero_grad(set_to_none=True)
   loss,parts=value_jvp_loss(model,h,context,target,jtarget,V,derivative_weight=derivative_weight)
   if not torch.isfinite(loss):raise FloatingPointError('nonfinite training loss; stop without tuning')
   loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
   if update%100==0 or update==args.updates:
    checkpoint(update);print(json.dumps(dict(update=update,loss=float(loss.detach()),value=float(parts['value'].detach()),derivative=float(parts['derivative'].detach()),best_validation_whiteMSE=best)),flush=True)
  write(dest/'RESULTS.json',dict(status='COMPLETE',mode=args.mode,seed=args.seed,selected_update=best_update,validation_white_MSE=best,
    training_physics_RHS=0,teacher_geometry_owner_RHS=459712,elapsed_seconds=time.perf_counter()-started,operator_test_read=False,acceleration_claim=False))
 except BaseException as exc:status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc);write(dest/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:b.close(status,error)


def evaluate(args):
 dest=OUT/'results/II/forward_jvp_evaluation'/args.job
 if dest.exists():raise FileExistsError('evaluation output exists')
 b=Book('II',args.job,max_wall=args.max_wall,analysis_only=True);status='COMPLETE';error=None
 try:
  source_snapshot(b);dest.mkdir(parents=True);tick=time.perf_counter()
  import torch
  torch.set_num_threads(1);state=torch.load(args.checkpoint,map_location='cpu',weights_only=False)
  checkpoint_load=time.perf_counter()-tick;tick=time.perf_counter();n=state['normalization']
  model=build_local_forward_model(state['coordinates'],n['input_mean'],n['input_scale'],n['output_scale']);model.load_state_dict(state['model']);cold_load=time.perf_counter()-tick
  data=load_data(state['coordinates'],'operator_test');metrics,predictions=evaluate_groups(model,data,torch,jvp_diagnostics=True)
  for parent,arrays in predictions.items():np.savez_compressed(dest/f'{parent}.npz',**arrays)
  write(dest/'RESULTS.json',dict(status='COMPLETE',mode=state['mode'],seed=state['seed'],checkpoint=str(args.checkpoint.resolve()),
    cold_checkpoint_and_torch_load_seconds=checkpoint_load,cold_model_build_and_state_load_seconds=cold_load,
    validation_gate=state['validation']['empirical_validation_gate'],metrics=metrics,
    scope='EXPOSED_OPERATOR_HOLDOUT_NOT_LOCKED; saved paid Full references, no new evaluator query',
    source_teacher_geometry_RHS=459712,NN_predictions_not_certificates=True,requires_Full_correction_for_exact_probability_claim=True,
    actual_new_physics_RHS=0,acceleration_claim=False))
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--train',action='store_true');parser.add_argument('--evaluate',action='store_true')
 parser.add_argument('--mode',choices=('value','jvp'),default='value');parser.add_argument('--seed',type=int,choices=SEEDS,default=SEEDS[0]);parser.add_argument('--coordinates',type=int,choices=(64,66),default=66)
 parser.add_argument('--updates',type=int,choices=(500,1000,2000,3000),default=1000);parser.add_argument('--job');parser.add_argument('--checkpoint',type=Path);parser.add_argument('--max-wall',type=float,default=900)
 args=parser.parse_args()
 if sum((args.prepare,args.train,args.evaluate))!=1:parser.error('choose exactly one explicit action')
 if not 0<args.max_wall<=900:parser.error('maximum900 controlled seconds per CPU job')
 if args.prepare:
  manifest=json.loads((PACK/'MANIFEST.json').read_text())
  if manifest['status']!='READY_NO_NN_OPTIMIZATION':raise ValueError('forward/JVP pack not ready')
  write(OUT/'AUDIT/forward_jvp_surrogate/TRAINING_CLI_READY.json',dict(status='PREPARED_NOT_TRAINED',modes=['value','jvp'],seeds=SEEDS,
    parameters66=parameter_count(66),parameters64=parameter_count(64),initial_updates=1000,maximum_updates=3000,device='cpu',new_physics_RHS=0,NN_optimizer_steps=0,inference_wall='NOT_MEASURED_UNTRAINED'))
  print('Forward/JVP training CLI ready; no optimizer or physics launched.');return
 if not args.job:parser.error('explicit unique --job required')
 if args.train:train(args)
 else:
  if args.checkpoint is None:parser.error('--checkpoint required')
  evaluate(args)

if __name__=='__main__':main()
