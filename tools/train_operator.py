#!/usr/bin/env python3
"""Prepare or explicitly train the <=0.8M empirical operator-action ranker.

--prepare performs no NN optimization. --train is a separate future invocation;
validation controls the checkpoint. Operator test parents are opened only by
--evaluate-operator-test with a frozen checkpoint. All labels are paid-J local
GN action errors, never material GT risk or a certified fallback rule.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,argparse,time,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write
from p3r4.operator_surrogate import (FEATURE_NAMES,BASELINES,ERROR_SCALE,ranker_parameter_count,
 build_action_error_ranker,empirical_action_errors)
PACK=OUT/'results/II/operator_pack_v1'
PROTOCOL=OUT/'configs/OPERATOR_MODEL_PROTOCOL.json'
SEEDS=(2026102501,2026102502,2026102503)


def protocol():
 return dict(status='PREPARED_NOT_TRAINED',family='empirical_action_error_refreshranker',
  DINO_replication=False,coordinates_primary=66,coordinates_separate_control=64,
  input_dim=len(FEATURE_NAMES),feature_names=FEATURE_NAMES,baselines=BASELINES,
  target='log1p(absolute local material GN action error/.02)',error_scale=ERROR_SCALE,
  architecture=dict(hidden=[256,128,64],activation='SiLU',normalization='LayerNorm',outputs=4),
  parameters=ranker_parameter_count(),parameter_cap=800000,
  seeds=SEEDS,updates=6000,batch=256,learning_rate=.0003,
  training_sampling='uniform parent then uniform pair; no Heun-parent oversampling advantage',
  scaling='parent-balanced train-only mean/variance, frozen for validation/operator test',
  checkpoint='minimum validation mean of per-parent MSE; include untrained update0',
  train_parents=list(range(16)),validation_parents=list(range(16,24)),operator_test_parents=list(range(24,32)),
  exposure='ALL_EXPOSED_DEVELOPMENT_NOT_LOCKED',test_read_during_training=False,
  deployment_inputs='cached anchor Full/J plus candidate8-RHS Full residual/field and analytic decoder, not J1',
  labels='paid FullJ1 teacher; no GT',predictions='empirical error scores; not certified bounds',
  new_physics_RHS=0,teacher_geometry_acquisition_RHS=459712,
  budget='one explicit future training job per seed, CPU max900 controlled seconds, no GPU by default',
  optimization_authorization='--prepare does not train; main orchestrator must explicitly invoke --train')


def parents(which):
 return tuple(f'R4DEV_{i:05d}' for i in (range(16) if which=='train' else range(16,24) if which=='validation' else range(24,32)))


def load_group(coordinates,which):
 folder=PACK/('probit66' if coordinates==66 else 'field64');data=[]
 for parent in parents(which):
  with np.load(folder/f'{parent}.npz',allow_pickle=False) as z:
   if str(z['split'])!=which or int(z['dimension'])!=coordinates or z['feature_names'].tolist()!=list(FEATURE_NAMES):raise ValueError('pack identity/split/feature order')
   data.append(dict(parent=parent,features=z['features'].copy(),target=z['log_error_target'].astype(np.float32),true_error=z['material_action_error'].copy()))
 return data


def balanced_scaler(training):
 mean=np.mean([d['features'].astype(np.float64).mean(axis=0) for d in training],axis=0)
 second=np.mean([np.mean(d['features'].astype(np.float64)**2,axis=0) for d in training],axis=0)
 scale=np.sqrt(np.maximum(second-mean*mean,1e-12))
 return mean.astype(np.float32),scale.astype(np.float32)


def validate(model,data,mean,scale,torch):
 model.eval();loss=[]
 with torch.no_grad():
  for d in data:
   total=0.;n=0
   for first in range(0,len(d['features']),256):
    x=torch.from_numpy((d['features'][first:first+256]-mean)/scale);y=torch.from_numpy(d['target'][first:first+256])
    residual=model(x)-y;total+=float((residual*residual).sum());n+=residual.numel()
   loss.append(total/n)
 return float(np.mean(loss)),loss


def train(args):
 if not args.job:raise ValueError('explicit unique --job required for NN optimization')
 if json.loads(PROTOCOL.read_text())!=json.loads(json.dumps(protocol())):raise ValueError('model protocol mismatch')
 import torch
 torch.set_num_threads(1);torch.manual_seed(args.seed)
 destination=OUT/'results/II/operator_training'/args.job
 if destination.exists():raise FileExistsError('training result exists')
 b=Book('II',args.job,max_wall=args.max_wall,kind='training',gpu=False)
 status='COMPLETE';error=None;started=time.perf_counter()
 try:
  destination.mkdir(parents=True)
  training=load_group(args.coordinates,'train');validation=load_group(args.coordinates,'validation')
  mean,scale=balanced_scaler(training);model=build_action_error_ranker();optimizer=torch.optim.AdamW(model.parameters(),lr=.0003,weight_decay=.00001)
  rng=np.random.default_rng(args.seed);best=np.inf;best_update=None
  def checkpoint(update):
   nonlocal best,best_update
   loss,byparent=validate(model,validation,mean,scale,torch)
   record=dict(update=update,validation_parent_balanced_MSE=loss,parent_losses=byparent)
   with (destination/'VALIDATION.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
   if loss<best:
    best=loss;best_update=update
    torch.save(dict(model=model.state_dict(),mean=mean,scale=scale,coordinates=args.coordinates,
      feature_names=FEATURE_NAMES,protocol=protocol(),seed=args.seed,selected_update=update,
      validation_MSE=loss,operator_test_opened=False),destination/'BEST.pt')
   model.train()
  checkpoint(0)
  for update in range(1,6001):
   parent_choices=rng.integers(0,len(training),size=256);features=[];target=[]
   for p in parent_choices:
    d=training[int(p)];i=int(rng.integers(len(d['features'])));features.append(d['features'][i]);target.append(d['target'][i])
   x=torch.from_numpy((np.asarray(features)-mean)/scale);y=torch.from_numpy(np.asarray(target))
   optimizer.zero_grad(set_to_none=True);loss=((model(x)-y)**2).mean();loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step()
   if update%200==0:checkpoint(update);print(json.dumps(dict(update=update,loss=float(loss.detach()),best_validation=best)),flush=True)
  write(destination/'RESULTS.json',dict(status='COMPLETE',coordinates=args.coordinates,seed=args.seed,selected_update=best_update,best_validation_MSE=best,
   parameters=sum(p.numel() for p in model.parameters()),new_physics_RHS=0,operator_test_opened=False,teacher_geometry_acquisition_RHS=459712,
   empirical_only=True,posterior_or_GT_risk_claim=False,elapsed_seconds=time.perf_counter()-started))
 except BaseException as exc:
  status='ERROR';error=repr(exc);write(destination/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:b.close(status,error)


def evaluate(args):
 if not args.job or args.checkpoint is None:raise ValueError('explicit evaluation --job and --checkpoint required')
 import torch
 torch.set_num_threads(1);state=torch.load(args.checkpoint,map_location='cpu',weights_only=False)
 if list(state['feature_names'])!=list(FEATURE_NAMES):raise ValueError('checkpoint features mismatch')
 b=Book('II',args.job,max_wall=args.max_wall,analysis_only=True);status='COMPLETE';error=None
 dest=OUT/'results/II/operator_evaluation'/args.job
 try:
  if dest.exists():raise FileExistsError('evaluation result exists')
  dest.mkdir(parents=True);model=build_action_error_ranker();model.load_state_dict(state['model']);model.eval();data=load_group(int(state['coordinates']),'operator_test');rows=[]
  with torch.no_grad():
   for d in data:
    predictions=[]
    for first in range(0,len(d['features']),256):
     x=torch.from_numpy((d['features'][first:first+256]-state['mean'])/state['scale']);predictions.append(model(x).numpy())
    predicted=empirical_action_errors(np.concatenate(predictions));choice=np.argmin(predicted,axis=1)
    selected=d['true_error'][np.arange(len(choice)),choice]
    np.savez_compressed(dest/f"{d['parent']}.npz",predicted_error=predicted,selected_baseline=choice,observed_teacher_error=d['true_error'],selected_observed_error=selected)
    rows.append(dict(parent=d['parent'],pairs=len(choice),selected_observed_error_mean=float(selected.mean()),fixed_baseline_error_mean=d['true_error'].mean(axis=0).tolist()))
  write(dest/'RESULTS.json',dict(status='COMPLETE',scope='EXPOSED_OPERATOR_HOLDOUT_DEVELOPMENT_ONLY',parents=rows,
    source_checkpoint=str(args.checkpoint.resolve()),new_physics_RHS=0,empirical_errors_not_certificates=True,no_GT=True))
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--train',action='store_true');parser.add_argument('--evaluate-operator-test',action='store_true')
 parser.add_argument('--coordinates',type=int,choices=(64,66),default=66);parser.add_argument('--seed',type=int,choices=SEEDS,default=SEEDS[0]);parser.add_argument('--job');parser.add_argument('--checkpoint',type=Path);parser.add_argument('--max-wall',type=float,default=900)
 args=parser.parse_args()
 if sum((args.prepare,args.train,args.evaluate_operator_test))!=1:parser.error('choose exactly one action')
 if not 0<args.max_wall<=900:parser.error('max wall <=900')
 if args.prepare:
  manifest=json.loads((PACK/'MANIFEST.json').read_text())
  if manifest['status']!='READY_NO_NN_OPTIMIZATION':raise ValueError('operator pack not ready')
  if PROTOCOL.exists():raise FileExistsError('frozen model protocol exists')
  write(PROTOCOL,protocol());print(json.dumps(dict(status='PREPARED_NOT_TRAINED',protocol=str(PROTOCOL),parameters=ranker_parameter_count(),features=len(FEATURE_NAMES))));return
 if args.train:train(args)
 else:evaluate(args)

if __name__=='__main__':main()
