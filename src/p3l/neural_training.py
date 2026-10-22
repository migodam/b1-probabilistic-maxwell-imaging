"""Caller-supplied joint-supervision density training; no campaign data IO.

This module never constructs splits, generates observations, reads held-out or
truth paths, invokes physics, selects a device automatically, or certifies a
posterior. Label coordinates/transforms are explicitly caller-defined. Image
contexts pass through unchanged; scalar and label scaling uses training only.
"""
from pathlib import Path
import json
import time
from copy import deepcopy

import numpy as np
import torch
from torch import nn

from .neural import build_model,parameter_count


DECLARED_SEEDS=(2026102260,2026102261,2026102262)
DEFAULTS=dict(learning_rate=1e-3,batch_size=128,max_updates=3000,
              validation_every=100,patience_evaluations=10,
              gradient_clip_norm=5.,max_wall_seconds=1200.)


def _device(value):
    if value not in ('cpu','mps'):
        raise ValueError("device must explicitly be 'cpu' or 'mps'; no automatic fallback")
    if value=='mps' and not torch.backends.mps.is_available():
        raise RuntimeError('explicit MPS device requested but unavailable')
    return torch.device(value)


def _array(value,name,shape):
    a=np.asarray(value,float)
    if a.ndim!=len(shape) or any(s is not None and a.shape[i]!=s for i,s in enumerate(shape)) or a.shape[0]==0 or not np.all(np.isfinite(a)):
        raise ValueError(f'{name} must have finite nonempty shape {shape}')
    if np.any(np.abs(a)>np.finfo(np.float32).max):
        raise ValueError(f'{name} exceeds the explicit float32 training range')
    return a


def _scaling(train_scalar,train_labels):
    def fit(a):
        mean=np.mean(a,axis=0)
        std=np.std(a,axis=0)
        constant=std==0
        scale=np.where(constant,1.,std)
        # Actual stored float32 affine transform defines the reported density.
        mean,scale=mean.astype(np.float32),scale.astype(np.float32)
        if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(scale)) or np.any(scale<=0):
            raise ValueError('train affine scale not representable in float32')
        return mean,scale,constant
    sm,ss,sc=fit(train_scalar);lm,ls,lc=fit(train_labels)
    return dict(scalar_mean=sm,scalar_scale=ss,label_mean=lm,label_scale=ls,
                scalar_constant=sc,label_constant=lc,
                policy='TRAIN_ONLY_POPULATION_SD; EXACT_ZERO_SD_USES_SCALE_ONE; IMAGES_UNCHANGED',
                log_abs_det_raw_label_to_standard=float(-np.log(ls.astype(float)).sum()))


class ScaledConditionalDensity(nn.Module):
    """Normalized density in caller-supplied raw label coordinates.

    Only scalar context and labels are standardized. No probit, property,
    whitening or other physical transform is performed by this wrapper.
    """
    def __init__(self,base,scaling):
        super().__init__();self.base=base
        for name in ('scalar_mean','scalar_scale','label_mean','label_scale'):
            self.register_buffer(name,torch.as_tensor(np.asarray(scaling[name]),dtype=torch.float32).clone())

    def log_prob(self,values,scalar,image):
        z=(values-self.label_mean)/self.label_scale
        context=(scalar-self.scalar_mean)/self.scalar_scale
        return self.base.log_prob(z,context,image)-torch.log(self.label_scale).sum()

    def sample(self,count,scalar,image):
        context=(scalar-self.scalar_mean)/self.scalar_scale
        return self.base.sample(count,context,image)*self.label_scale+self.label_mean


def _cpu_state(model):
    return {key:value.detach().cpu().clone() for key,value in model.state_dict().items()}


def _scaling_plain(scaling):
    return {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in scaling.items()}


def load_density_checkpoint(path,*,device='cpu'):
    """Restore this module's own checkpoint; no training/data lookup occurs."""
    target_device=_device(device)
    payload=torch.load(Path(path),map_location='cpu',weights_only=True)
    if payload.get('format')!='P3L_SUPPLIED_DATA_DENSITY_V1':
        raise ValueError('unsupported density checkpoint format')
    base=build_model(payload['model_type'],**payload['model_options'])
    scaling=payload['train_scaling']
    model=ScaledConditionalDensity(base,scaling)
    model.load_state_dict(payload['state_dict'],strict=True)
    model.to(target_device);model.eval()
    return dict(model=model,checkpoint=payload,train_scaling=scaling,
                metadata=payload['metadata'])


def train_density(model_type,train_scalar,train_image,train_labels,
                  validation_scalar,validation_image,validation_labels,*,
                  output_dir,label_coordinates,seed=DECLARED_SEEDS[0],device='cpu',
                  supervision='JOINT_SIMULATION',qualified_distillation=False,
                  learning_rate=1e-3,batch_size=128,max_updates=3000,
                  validation_every=100,patience_evaluations=10,
                  gradient_clip_norm=5.,max_wall_seconds=1200.,
                  gpu_wall_callback=None):
    """Train from explicit arrays, select by validation raw-coordinate NLL.

    The three default declared seeds are documented, not three automatic runs.
    An explicitly supplied different integer seed is recorded as an override.
    Only CPU(single thread) or explicit MPS is supported; arrays never transfer
    externally. A nonempty label-coordinate description is mandatory. Default
    supervision is independent joint simulation; caller split/independence
    assertions remain unverified. DISTILLATION requires the explicit qualified
    flag and is labeled as a caller assertion, never internally certified.

    MPS callback receives total synchronized allocated-job wall seconds once,
    including host overhead, not isolated GPU kernel time. Caller bills it.
    Exceptions still invoke this callback when MPS work has begun. Wall is a
    cooperative cap checked between batches/validation batches; a running
    operation/checkpoint can overrun. Partial validation never selects a best.
    """
    target_device=_device(device)
    if not isinstance(label_coordinates,str) or not label_coordinates.strip():
        raise ValueError('explicit nonempty label_coordinates metadata required')
    if supervision not in ('JOINT_SIMULATION','DISTILLATION'):
        raise ValueError('supervision must be JOINT_SIMULATION or DISTILLATION')
    if supervision=='DISTILLATION' and not qualified_distillation:
        raise ValueError('distillation requires explicit qualified_distillation caller assertion')
    if not isinstance(seed,(int,np.integer)) or isinstance(seed,bool) or not 0<=seed<2**63:
        raise ValueError('nonnegative integer seed below2**63 required')
    for key,value in dict(batch_size=batch_size,max_updates=max_updates,validation_every=validation_every,patience_evaluations=patience_evaluations).items():
        if not isinstance(value,int) or isinstance(value,bool) or value<1:raise ValueError(f'positive integer {key} required')
    if any(not np.isfinite(v) or v<=0 for v in (learning_rate,gradient_clip_norm,max_wall_seconds)):
        raise ValueError('positive finite learning rate/gradient clip/wall cap required')
    ts=_array(train_scalar,'train_scalar',(None,None))
    tl=_array(train_labels,'train_labels',(len(ts),None))
    ti=_array(train_image,'train_image',(len(ts),4,24,24))
    vs=_array(validation_scalar,'validation_scalar',(None,ts.shape[1]))
    vl=_array(validation_labels,'validation_labels',(len(vs),tl.shape[1]))
    vi=_array(validation_image,'validation_image',(len(vs),4,24,24))
    if ts.shape[1]<1 or tl.shape[1]<1:raise ValueError('nonempty context/label dimensions required')
    if max_updates>3000 or max_wall_seconds>1200:
        raise ValueError('requested cap exceeds3000updates or1200seconds')
    scaling=_scaling(ts,tl)
    arrays=[np.ascontiguousarray(a,dtype=np.float32) for a in (ts,ti,tl,vs,vi,vl)]
    # CPU tensors back each explicit minibatch; MPS transfer is local and only
    # performed on the requested device. No full data copy to device is hidden.
    tensors=[torch.from_numpy(a) for a in arrays]
    output=Path(output_dir);output.mkdir(parents=True,exist_ok=True)
    if (output/'BEST_DENSITY.pt').exists() or (output/'TRAINING_REPORT.json').exists():
        raise FileExistsError('training output already exists; select a new run directory or explicitly restore checkpoint')
    checkpoint_path=output/'BEST_DENSITY.pt'
    logs=[];best_state=None;best_nll=np.inf;best_update=None;stale=0;updates=0
    start=time.perf_counter();stop_reason='MAX_UPDATES';gpu_wall=0.
    payload=None
    previous_threads=torch.get_num_threads()
    if device=='cpu':torch.set_num_threads(1)
    options=dict(scalar_dim=ts.shape[1],image_shape=(4,24,24),output_dim=tl.shape[1])
    metadata=dict(label_coordinates=label_coordinates,
        raw_target_nll_coordinates='EXACT_CALLER_LABEL_COORDINATES_WITH_AFFINE_STANDARDIZATION_JACOBIAN',
        additional_label_coordinate_jacobians='CALLER_RESPONSIBILITY; NONE_INFERRED_OR_APPLIED',
        supervision=supervision,distillation_qualification='CALLER_ASSERTED_NOT_VERIFIED' if supervision=='DISTILLATION' else 'NOT_DISTILLATION',
        split_provenance='EXPLICIT_CALLER_ARRAYS_NO_SPLIT_CONSTRUCTION_OR_TEST_READ',
        image_treatment='UNCHANGED_CALLER_CONTEXT_NO_GT_TRANSFORM',
        seed=int(seed),declared_seeds=list(DECLARED_SEEDS),seed_override=bool(seed not in DECLARED_SEEDS),
        device=device,dtype='float32',train_count=len(ts),validation_count=len(vs),
        learning_rate=learning_rate,batch_size=batch_size,max_updates=max_updates,
        validation_every=validation_every,patience_evaluations=patience_evaluations,
        gradient_clip_norm=gradient_clip_norm,max_wall_seconds=max_wall_seconds,
        initial_validation_is_candidate=True,checkpoint_selection='LOWEST_COMPLETE_VALIDATION_RAW_TARGET_NLL',
        density_status='NORMALIZED_ARCHITECTURE_NOT_POSTERIOR_CERTIFIED',
        wall_scope='SYNCHRONIZED_ALLOCATED_DEVICE_JOB_WALL_INCLUDING_HOST_OVERHEAD',
        wall_cap_type='COOPERATIVE_BETWEEN_OPERATIONS_NOT_PREEMPTIVE',torch_version=str(torch.__version__))
    def sync():
        if device=='mps':torch.mps.synchronize()
    def elapsed():
        sync();return time.perf_counter()-start
    def expired():return elapsed()>=max_wall_seconds
    def batch(indices,offset):
        return [tensors[offset+i][indices].to(target_device) for i in range(3)]
    try:
        # Preserve caller's CPU RNG state. There is no stochastic model layer
        # or MPS RNG draw during training; initialization occurs on CPU.
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(int(seed));rng=np.random.default_rng(seed)
            model=ScaledConditionalDensity(build_model(model_type,**options),scaling).to(target_device)
            metadata['parameter_count']=parameter_count(model)
            optimizer=torch.optim.Adam(model.parameters(),lr=learning_rate)
            def evaluate(update,training_nll=None,gradient_norm=None):
                nonlocal best_nll,best_state,best_update,stale,payload
                model.eval();total=0.;count=0
                with torch.no_grad():
                    for left in range(0,len(vs),batch_size):
                        if expired():
                            logs.append(dict(update=update,event='VALIDATION_INCOMPLETE_WALL_CAP',evaluated=count,required=len(vs),wall_seconds=elapsed()))
                            return False
                        scalar,image,labels=batch(slice(left,left+batch_size),3)
                        logp=model.log_prob(labels,scalar,image)
                        if not torch.isfinite(logp).all():raise FloatingPointError('nonfinite validation density')
                        total+=float((-logp).sum().cpu());count+=len(labels)
                nll=total/count
                improved=nll<best_nll
                if improved:
                    best_nll=nll;best_state=_cpu_state(model);best_update=update;stale=0
                    payload=dict(format='P3L_SUPPLIED_DATA_DENSITY_V1',model_type=model_type,model_options=options,
                        state_dict=best_state,train_scaling=_scaling_plain(scaling),metadata=deepcopy(metadata),
                        best_update=best_update,best_validation_raw_target_nll=best_nll)
                    temporary=checkpoint_path.with_suffix('.pt.tmp')
                    torch.save(payload,temporary);temporary.replace(checkpoint_path)
                else:stale+=1
                logs.append(dict(update=update,event='COMPLETE_VALIDATION',validation_raw_target_nll=nll,
                    training_batch_raw_target_nll=training_nll,preclip_gradient_norm=gradient_norm,
                    improved=improved,stale_evaluations=stale,wall_seconds=elapsed()))
                return True
            complete=evaluate(0)
            if not complete:stop_reason='WALL_CAP_BEFORE_COMPLETE_INITIAL_VALIDATION'
            else:
                for update in range(1,max_updates+1):
                    if expired():stop_reason='WALL_CAP';break
                    model.train()
                    indices=rng.integers(0,len(ts),size=min(batch_size,len(ts)))
                    scalar,image,labels=batch(indices,0)
                    optimizer.zero_grad(set_to_none=True)
                    loss=-model.log_prob(labels,scalar,image).mean()
                    if not torch.isfinite(loss):raise FloatingPointError('nonfinite training NLL')
                    loss.backward()
                    norm=torch.nn.utils.clip_grad_norm_(model.parameters(),gradient_clip_norm,error_if_nonfinite=True)
                    optimizer.step();updates=update
                    if update%validation_every==0 or update==max_updates:
                        if not evaluate(update,float(loss.detach().cpu()),float(norm.detach().cpu())):
                            stop_reason='WALL_CAP_DURING_VALIDATION';break
                        if stale>=patience_evaluations:stop_reason='VALIDATION_PATIENCE';break
            if best_state is not None:model.load_state_dict(best_state,strict=True)
            model.eval();wall=elapsed()
            metadata.update(updates_completed=updates,stop_reason=stop_reason,
                            best_update=best_update,best_validation_raw_target_nll=None if best_state is None else best_nll,
                            wall_seconds=wall,best_checkpoint_available=best_state is not None)
            # Update metadata on the already selected best checkpoint without
            # altering selected weights or its NLL.
            if best_state is not None:
                payload['metadata']=deepcopy(metadata)
                temporary=checkpoint_path.with_suffix('.pt.tmp')
                torch.save(payload,temporary);temporary.replace(checkpoint_path)
            report=dict(metadata=metadata,train_scaling=_scaling_plain(scaling),logs=logs,
                        checkpoint_path=str(checkpoint_path.resolve()) if best_state is not None else None)
            (output/'TRAINING_REPORT.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            return dict(model=model,state_dict=best_state,train_scaling=scaling,
                        checkpoint_path=report['checkpoint_path'],metadata=metadata,logs=logs)
    except BaseException as error:
        failure=dict(metadata=metadata,updates_completed=updates,error=repr(error),logs=logs,
                     best_checkpoint_available=best_state is not None,
                     status='TRAINING_FAILED_NOT_A_QUALIFIED_DENSITY')
        (output/'TRAINING_FAILURE.json').write_text(json.dumps(failure,indent=2,allow_nan=False)+'\n')
        raise
    finally:
        try:
            sync()
        finally:
            gpu_wall=time.perf_counter()-start if device=='mps' else 0.
            if device=='cpu':torch.set_num_threads(previous_threads)
            if device=='mps' and gpu_wall_callback is not None:gpu_wall_callback(float(gpu_wall))
