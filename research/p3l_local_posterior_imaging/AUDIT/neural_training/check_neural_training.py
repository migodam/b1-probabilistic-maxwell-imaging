"""Three-update analytic dummy check, no campaign arrays or full training."""
from pathlib import Path
import tempfile,sys,json
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
from p3l.neural_training import train_density,load_density_checkpoint,DECLARED_SEEDS
OUT=Path(__file__).resolve().parent
rng=np.random.default_rng(2026102260);checks={}
def check(name,ok,**details):
    checks[name]=dict(passed=bool(ok),**details)
    if not ok:raise AssertionError(name)
n=16;nv=8
scalar=rng.normal(size=(n,6))*np.arange(1,7)+np.arange(6)
validation_scalar=rng.normal(size=(nv,6))+100
labels=rng.normal(size=(n,2))*[2.,3.]+[-5.,8.]
validation_labels=rng.normal(size=(nv,2))*[2.,3.]+[-5.,8.]
image=np.zeros((n,4,24,24));validation_image=np.zeros((nv,4,24,24))
folder=Path(tempfile.mkdtemp(prefix='dummy_',dir=OUT))
callback=[]
result=train_density('mlp_diag',scalar,image,labels,validation_scalar,validation_image,validation_labels,
    output_dir=folder,label_coordinates='SUPPLIED_ANALYTIC_DUMMY_REAL2; NO_PROBIT_OR_PHYSICAL_LABEL',
    seed=DECLARED_SEEDS[0],device='cpu',max_updates=3,batch_size=4,validation_every=1,
    patience_evaluations=10,max_wall_seconds=60,gpu_wall_callback=lambda seconds:callback.append(seconds))
scaling=result['train_scaling']
check('train_only_scaling',np.allclose(scaling['scalar_mean'],scalar.mean(0)) and np.allclose(scaling['label_scale'],labels.std(0)) and np.max(abs(scaling['scalar_mean']-validation_scalar.mean(0)))>90)
check('three_updates_only',result['metadata']['updates_completed']==3 and len(result['logs'])==4,updates=result['metadata']['updates_completed'],complete_validations=len(result['logs']))
model=result['model'];s=torch.tensor(validation_scalar,dtype=torch.float32);im=torch.tensor(validation_image,dtype=torch.float32);v=torch.tensor(validation_labels,dtype=torch.float32)
with torch.no_grad():
    raw=model.log_prob(v,s,im)
    standardized=model.base.log_prob((v-model.label_mean)/model.label_scale,(s-model.scalar_mean)/model.scalar_scale,im)
    difference=standardized-raw
    expected=float(np.log(scaling['label_scale'].astype(float)).sum())
check('affine_label_jacobian_in_raw_nll',np.max(abs(difference.numpy()-expected))<2e-6,expected_log_abs_det_standard_to_raw=expected,maximum_error=float(np.max(abs(difference.numpy()-expected))))
restored=load_density_checkpoint(result['checkpoint_path'],device='cpu')
with torch.no_grad():after=restored['model'].log_prob(v,s,im)
check('best_checkpoint_restores_exact_weights',torch.equal(raw,after))
# Recompute validation using the same batch summation convention as selection.
with torch.no_grad():nll=sum(float((-after[i:i+4]).sum()) for i in range(0,nv,4))/nv
logged=min(row['validation_raw_target_nll'] for row in result['logs'] if row['event']=='COMPLETE_VALIDATION')
check('best_by_raw_validation_nll',abs(nll-logged)<1e-7 and abs(nll-result['metadata']['best_validation_raw_target_nll'])<1e-7,validation_raw_nll=nll,best_update=result['metadata']['best_update'])
check('cpu_not_billed_as_gpu',not callback and result['metadata']['device']=='cpu')
check('coordinates_and_supervision_persist',restored['metadata']['label_coordinates'].startswith('SUPPLIED_ANALYTIC_DUMMY') and restored['metadata']['supervision']=='JOINT_SIMULATION' and restored['metadata']['density_status']=='NORMALIZED_ARCHITECTURE_NOT_POSTERIOR_CERTIFIED')
try:train_density('mlp_diag',scalar,image,labels,validation_scalar,validation_image,validation_labels,output_dir=folder/'blocked',label_coordinates='dummy',supervision='DISTILLATION',max_updates=1)
except ValueError:check('unqualified_distillation_rejected',True)
else:check('unqualified_distillation_rejected',False)
try:train_density('mlp_diag',scalar,image,labels,validation_scalar,validation_image,validation_labels,output_dir=folder/'blocked',label_coordinates='dummy',device='auto',max_updates=1)
except ValueError:check('implicit_device_rejected',True)
else:check('implicit_device_rejected',False)
try:train_density('mlp_diag',scalar,image,labels,validation_scalar,validation_image,validation_labels,output_dir=folder,label_coordinates='dummy',max_updates=1)
except FileExistsError:check('completed_output_not_overwritten',True)
else:check('completed_output_not_overwritten',False)
report=dict(status='PASS_THREE_UPDATE_SYNTHETIC_DUMMY_ONLY',checks=checks,
    dummy_updates=3,campaign_training_runs=0,campaign_data_reads=0,physics_calls=0,
    actual_device='cpu',declared_seeds=list(DECLARED_SEEDS),dummy_output=str(folder.resolve()))
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
