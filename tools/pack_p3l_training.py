"""Dedicated supervised-label packing process; DO NOT import in online inference.

This script intentionally reads latent truths only after all public contexts
and generator/context closure receipts have passed preflight. It never imports
Target, Book, a physical backend, an inversion routine, or a trainer. It does
not construct splits, inspect locked/development truth, filter by GT quality,
fit scaling or train a network. Parent launches this separately after online
context generation and bills ordinary packing wall externally if required.

Main files have exactly scalar/image/labels/parent. Strong truth-coordinate
information is written under diagnostics/, outside default training arrays.
"""
from pathlib import Path
import argparse
import json
import os
import re
import time

import numpy as np
from scipy.special import ndtri


ROOT=Path(__file__).resolve().parents[1]
DEFAULT_OUT=ROOT/'research/p3l_local_posterior_imaging'
GROUP_COUNTS={'train':512,'validation':128}
FROZEN_TRAINING_SEEDS=[2026102270,2026102271,2026102272]


def _inside(path,root):
    path=Path(path).resolve()
    if not path.is_relative_to(root):raise ValueError('source/output path escapes supplied campaign root')
    return path


def _json(path,root):
    return json.loads(_inside(path,root).read_text())


def _array(data,name,shape,path):
    if name not in data:raise ValueError(f'missing {name} in {path}')
    value=np.asarray(data[name],float)
    if value.shape!=shape or not np.all(np.isfinite(value)):
        raise ValueError(f'invalid finite {name}{shape} in {path}')
    return value


def _save(path,**arrays):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream,**arrays);stream.flush();os.fsync(stream.fileno())
    temporary.replace(path)


def _write(path,record):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def _complete_job(root,job):
    if not isinstance(job,str) or not re.fullmatch(r'[A-Za-z0-9_.-]+',job) or job in ('.','..'):
        raise ValueError('invalid recorded source job')
    receipt=_json(root/f'D/jobs/{job}/CLOSED.json',root)
    if receipt.get('status')!='COMPLETE' or receipt.get('error') is not None:
        raise RuntimeError(f'source job is not successfully closed: {job}')
    return dict(job=job,closure=str(root/f'D/jobs/{job}/CLOSED.json'),status='COMPLETE')


def pack_training(campaign_root=DEFAULT_OUT):
    """Pack all512train and128validation parents; no other groups accepted.

    This is the ONLY intentional latent-truth access in this tool. All public
    contexts are validated before the first offline_only file is opened. The
    caller must start this as a dedicated process after generators/inference
    exit; successful closure receipts are evidence, not an OS liveness proof.
    No arrays are read on module import. No execution performed by authoring.
    """
    started=time.perf_counter();root=Path(campaign_root).expanduser().resolve()
    if not root.is_dir():raise ValueError('campaign root must already exist')
    destination=_inside(root/'D/training',root)
    outputs=[destination/'TRAIN.npz',destination/'VALIDATION.npz',
             destination/'diagnostics/TRAIN_STRONG.npz',destination/'diagnostics/VALIDATION_STRONG.npz',
             destination/'PACKING_METADATA.json']
    if any(path.exists() for path in outputs):raise FileExistsError('packing output already present; no implicit overwrite or resume')
    protocol=_json(root/'configs/PROTOCOL.json',root)
    seeds=protocol.get('learning',{}).get('seeds')
    if seeds!=FROZEN_TRAINING_SEEDS:
        raise ValueError('frozen protocol learning seeds differ from the declared70/71/72 contract')
    manifest=_json(root/'configs/PARENTS.json',root)
    if not isinstance(manifest,list):raise ValueError('parent manifest must be a list')
    groups={};generator_receipts={};closures={}
    for group,count in GROUP_COUNTS.items():
        rows=[row for row in manifest if row.get('group')==group]
        expected=[f'SBI_{group}_{i:04d}' for i in range(count)]
        actual=[row.get('parent') for row in rows]
        if len(actual)!=count or len(set(actual))!=count or set(actual)!=set(expected):
            raise ValueError(f'exactly the frozen{count}distinct{group}parents are required')
        groups[group]=expected
        receipt_path=root/f'AUDIT/{group}_DATA_GENERATOR_EXIT_REQUIRED.json'
        receipt=_json(receipt_path,root)
        if (receipt.get('group')!=group or receipt.get('generated')!=count
                or receipt.get('GT_inference') is not False
                or receipt.get('no_inference_in_this_process') is not True
                or receipt.get('process')!='dedicated truth generator'):
            raise RuntimeError(f'incomplete or incompatible generator exit-required receipt for{group}')
        closure=_complete_job(root,'dataset_'+group)
        generator_receipts[group]=dict(receipt=str(receipt_path),receipt_fields=receipt,closure=closure)
    frame_path=_inside(root/'D/FIXED_FRAME.npz',root)
    frame_meta=_json(frame_path.with_suffix('.json'),root)
    if (frame_meta.get('training_only') is not True or frame_meta.get('GT_used') is not False
            or frame_meta.get('weak_dimension')!=8):
        raise ValueError('fixed frame lacks required training-only/noGT/8D metadata')
    with np.load(frame_path,allow_pickle=False) as data:
        Vw=_array(data,'Vw',(64,8),frame_path)
        Vs=_array(data,'Vs',(64,56),frame_path)
        Q=_array(data,'Q',(64,64),frame_path)
    frame=np.column_stack((Vw,Vs))
    orthogonal_error=float(np.linalg.norm(frame.T@frame-np.eye(64)))
    determinant_sign,log_abs_det=np.linalg.slogdet(frame)
    if (orthogonal_error>1e-10 or not np.allclose(Q,frame,rtol=0,atol=1e-12)
            or determinant_sign==0 or abs(log_abs_det)>1e-10):
        raise ValueError('fixed frame must be a complete orthogonal Gaussian64 split')
    # Preflight ALL public files and all successful job closures before truth.
    public={}
    for group,parents in groups.items():
        contexts=[]
        for parent in parents:
            if not parent.startswith(f'SBI_{group}_'):raise ValueError('forbidden parent group')
            context_path=_inside(root/f'D/contexts/{parent}.npz',root)
            meta=_json(context_path.with_suffix('.json'),root)
            if meta.get('parent')!=parent or meta.get('GT_used') is not False:
                raise ValueError(f'incompatible noGT public context metadata for{parent}')
            job=meta.get('source_job')
            if not isinstance(job,str) or not job.startswith('contexts_'+group+'_'):
                raise ValueError(f'context source job does not match declared{group}')
            if job not in closures:
                source_receipt=_json(root/f'D/jobs/{job}/CLOSED.json',root)
                if source_receipt.get('status')=='INTERRUPTED_MEMORY_PRESSURE':
                    # Job completion failed, but this parent was fully saved
                    # and independently equation-audited by a completed resume.
                    recovery=_json(context_path.with_suffix('.recovery.json'),root)
                    if recovery.get('parent')!=parent or recovery.get('original_source_job')!=job or recovery.get('status')!='PAID_COMPLETED_CONTEXT_REUSED_AFTER_MEMORY_STOP' or recovery.get('GT_used') is not False:raise RuntimeError('invalid completed-parent recovery receipt')
                    resumed=_complete_job(root,recovery['validated_by_job'])
                    closures[job+'::'+parent]=dict(job=job,original_status='INTERRUPTED_MEMORY_PRESSURE',parent=parent,completed_parent_receipt=str(context_path.with_suffix('.recovery.json')),validation_job=resumed)
                    # Validate every recovered parent, not just the first job.
                else:closures[job]=_complete_job(root,job)
            with np.load(context_path,allow_pickle=False) as data:
                scalar=_array(data,'scalar',(1154,),context_path)
                image=_array(data,'image',(4,24,24),context_path)
                q_context=_array(data,'q',(66,),context_path)
            if np.any((q_context[64:]<0)|(q_context[64:]>1)):
                raise ValueError('online context has illegal bounded properties')
            if not np.array_equal(scalar[1024:1090],q_context):
                raise ValueError('scalar anchor coordinates disagree with separate context q')
            truth_path=_inside(root/f'inputs/offline_only/{parent}.npz',root)
            if not truth_path.is_file():raise FileNotFoundError(truth_path)
            contexts.append(dict(parent=parent,scalar=scalar,image=image,q=q_context,
                                 context_path=context_path,truth_path=truth_path))
        public[group]=contexts
    # Dedicated supervised process boundary: after this line only the exact
    # train/validation files from the frozen manifest are intentionally opened.
    packed={};source_inventory=[]
    for group,contexts in public.items():
        labels=[];strong=[];probit_logj=[]
        for context in contexts:
            truth_path=context['truth_path']
            with np.load(truth_path,allow_pickle=False) as truth:
                # Do not materialize chi/prediction/noise or inspect GT quality.
                q_truth=_array(truth,'q',(66,),truth_path)
            properties=q_truth[64:]
            if np.any((properties<=0)|(properties>=1)):
                raise ValueError(f'{context["parent"]}: exact probit label requires strict0<u<1; no clipping or resampling permitted')
            shift=q_truth[:64]-context['q'][:64]
            eta=ndtri(properties)
            label=np.r_[Vw.T@shift,eta]
            if not np.all(np.isfinite(label)):raise ValueError('nonfinite exact weak/probit label; no repair permitted')
            labels.append(label);strong.append(Vs.T@shift)
            probit_logj.append(float(-.5*(eta@eta)-np.log(2*np.pi)))
            source_inventory.append(dict(parent=context['parent'],group=group,
                context=str(context['context_path']),latent_supervision_source=str(truth_path),
                truth_fields_materialized=['q'],material_truth_read=False))
        packed[group]=dict(scalar=np.stack([c['scalar'] for c in contexts]),
            image=np.stack([c['image'] for c in contexts]),labels=np.stack(labels),
            parent=np.array([c['parent'] for c in contexts]),strong=np.stack(strong),
            log_abs_det_probit_eta_to_u=np.array(probit_logj))
    # Persist supervision and diagnostics separately. Main training arrays have
    # no chi, full qtruth, strong truth coordinates or posterior teacher data.
    destination.mkdir(parents=True,exist_ok=True)
    for group,name in [('train','TRAIN'),('validation','VALIDATION')]:
        arrays=packed[group]
        _save(destination/f'{name}.npz',**{key:arrays[key] for key in ('scalar','image','labels','parent')})
        _save(destination/f'diagnostics/{name}_STRONG.npz',strong=arrays['strong'],
              parent=arrays['parent'],log_abs_det_probit_eta_to_u=arrays['log_abs_det_probit_eta_to_u'])
    metadata=dict(status='COMPLETE_SUPERVISED_LABEL_PACKING_NOT_POSTERIOR_CERTIFICATION',
        process_role='DEDICATED_SUPERVISED_LABEL_PROCESS_AFTER_PUBLIC_CONTEXT_GENERATION',
        caller_must_verify_process_exit=True,receipts_are_not_independent_OS_liveness_proof=True,
        groups=GROUP_COUNTS,forbidden_truth_groups=['development','locked'],
        all_public_contexts_preflighted_before_truth=True,source_inventory=source_inventory,
        generator_receipts=generator_receipts,context_job_closures=list(closures.values()),
        frame=dict(path=str(frame_path),fixed=True,weak_dimension=8,strong_dimension=56,
                   orthogonality_frobenius_error=orthogonal_error,log_abs_determinant=float(log_abs_det),
                   chart_absolute_determinant=1.,source_metadata=frame_meta),
        input_shapes=dict(scalar=[1154],image=[4,24,24],labels=[10]),
        training_file_fields=['scalar','image','labels','parent'],
        diagnostics_not_default_training_input='diagnostics/{TRAIN,VALIDATION}_STRONG.npz',
        label_definition='[Vw.T@(ztruth-qcontext[:64]),ndtri(utruth)]',
        strong_definition='Vs.T@(ztruth-qcontext[:64])',
        label_coordinates='FIXED_GAUSSIAN_WEAK8_ANCHOR_RELATIVE_PLUS_EXACT_PROBIT_PROPERTY2',
        reference_measure='da8 deta2',target='FULL weak8/property2 joint conditional on supplied public context; NOT LOCAL-C',
        strong_coordinates='Marginalized in supervised conditional law; not frozen to zero or supplied to trainer by default',
        original_prior_reference='dz64 du2; original Gaussian64 x uniform2 unchanged',
        complete_frame_measure='dz64=da8 db56 at fixed public context; orthogonal chart absolute determinant1',
        density_probit_relation='log p_(a,u)(a,u|c)=log p_(a,eta)(a,ndtri(u)|c)-sum_j log standard_normal_pdf(eta_j)',
        nll_probit_relation='NLL_(a,u)=NLL_(a,eta)+sum_j log standard_normal_pdf(eta_j)',
        density_affine_relation='Trainer standardizes labels with training-only scales s: log p_(a,eta)=log p_standard-sum log(s)',
        combined_density_relation='log p_(a,u)=log p_standard-sum log(s)-sum log standard_normal_pdf(eta)',
        normalized_density_not_full66_proposal=True,
        support='a in R8; eta in R2 corresponds to strict0<u<1. Exact floating endpoints rejected; no label clipping/resampling.',
        local_selection='No LOCAL-C event, anchor-success/GT-distance filtering, nearest-surface repair, or oracle ranking',
        supervised_extra_information='Latent truth enters labels only in this dedicated offline supervised process; online contexts were generated independently of truth access. This is joint-simulation supervision, not posterior-reference distillation.',
        shared_parent_policy='Exactly one frozen context/label per distinct train/validation parent; no split construction or parent migration',
        frozen_training_seeds=FROZEN_TRAINING_SEEDS,
        trainer_seed_requirement='Explicitly override neural_training defaults with these frozen70/71/72 seeds; do not alter frozen protocol',
        fitted_scaling_here=False,training_runs=0,new_Maxwell_solves=0,new_physical_RHS=0,
        packing_wall_seconds=time.perf_counter()-started,
        physical_cost='Zero added physics; parent accounts ordinary packing wall externally if required')
    _write(destination/'PACKING_METADATA.json',metadata)
    return metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-root',type=Path,default=DEFAULT_OUT)
    args=parser.parse_args()
    record=pack_training(args.campaign_root)
    print(json.dumps(dict(status=record['status'],groups=record['groups'],
                          packing_wall_seconds=record['packing_wall_seconds'],
                          new_physical_RHS=0),allow_nan=False))


if __name__=='__main__':main()
