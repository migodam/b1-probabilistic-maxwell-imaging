#!/usr/bin/env python3
"""Frozen, truth-free locked inference. This script never opens offline labels."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from p3r2.risk_models import load_predictor
from p3r2.residual_risk import load_residual_predictor
from p3r2.density import load_density
from p3r2.fusion import passive_box,prior_material_mean
from p3r2.classical_risk import apply_action

DEFAULT_OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
ONLINE_FIELDS=('observation','q0','chi0','J','prediction','residual')
POINT_METHODS=('exactprior_material_mean','direct','projected_direct','residual','projected_residual','fused','anchor','probit','native')


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)


def save(path,**arrays):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix('.tmp')
    with temp.open('wb') as handle:np.savez_compressed(handle,**arrays)
    temp.replace(path)


def resolve_path(value,out):
    path=Path(value)
    if not path.is_absolute():path=out/path
    path=path.resolve()
    if 'offline_only' in path.parts:raise ValueError('offline labels prohibited in frozen inference inputs')
    return path


def check_matched_receipt(out,parent,frozen_unix):
    folder=out/'L/contexts/details'/parent/'matched'
    point,receipt=folder/'MATCHED.npz',folder/'MATCHED.json'
    for path in (point,receipt):
        if not path.is_file():raise FileNotFoundError('all256 matched online receipts required: '+str(path))
        if path.stat().st_mtime<float(frozen_unix):raise ValueError('matched receipt predates freeze: '+str(path))
    metadata=json.loads(receipt.read_text())
    if metadata.get('GT_used') is not False or metadata.get('extra_fresh_state_cap')!=16 or metadata.get('upper_new_RHS')!=640:
        raise ValueError('matched176 receipt must declare truth-free extra16 states/640RHS: '+parent)


def frozen_contract(out,methods_path):
    m=json.loads(Path(methods_path).read_text())
    if m.get('status')!='FROZEN':raise ValueError('methods marker status must be FROZEN')
    revision=m.get('source_revision')
    if not revision:raise ValueError('source_revision required')
    paths=m.get('frozen_source_paths',['src/p3r2','src/p1_material','tools/*p3r2*.py'])
    if not isinstance(paths,list) or not paths or any(not isinstance(path,str) for path in paths):
        raise ValueError('frozen_source_paths must be a nonempty whitelist')
    required={'src/p3r2','src/p1_material','tools/infer_p3r2_locked.py','tools/evaluate_p3r2_locked.py',
              'tools/run_p3r2_anchors.py','tools/generate_p3r2_data.py','tools/run_p3r2_geometry.py'}
    for path in paths:
        if path in ('src/p3r2','src/p1_material','tools/*p3r2*.py'):continue
        if not re.fullmatch(r'tools/[A-Za-z0-9_]*p3r2[A-Za-z0-9_]*\.py',path) or not (ROOT/path).is_file():
            raise ValueError('invalid frozen source whitelist path: '+path)
    covered=set(paths)
    if 'tools/*p3r2*.py' in covered:covered.update(path for path in required if path.startswith('tools/'))
    if not required.issubset(covered):raise ValueError('frozen source whitelist omits required online/data/evaluation code')
    subprocess.run(['git','diff','--exit-code',revision,'--',*paths],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
    untracked=subprocess.check_output(['git','ls-files','--others','--',*paths],cwd=ROOT,text=True)
    untracked_code=[name for name in untracked.splitlines() if name.endswith('.py')]
    if untracked_code:raise ValueError('untracked frozen source: '+', '.join(untracked_code))
    timestamp=float(m['frozen_unix'])
    parents=m.get('parents',[f'R2_locked_{i:05d}' for i in range(256)])
    if m.get('expected_count',256)!=256 or len(parents)!=256 or len(set(parents))!=256:
        raise ValueError('frozen locked cohort must contain exactly256 unique parents')
    if any(not p.startswith('R2_locked_') or '/' in p for p in parents):raise ValueError('invalid locked parent identifiers')
    if not isinstance(m.get('include_matched_point176',False),bool):raise ValueError('include_matched_point176 must be explicit boolean')
    fusion_value=m.get('fusion_json',m.get('fusion'))
    if fusion_value is None:raise ValueError('fusion_json required')
    artifacts={key:resolve_path(m[key],out) for key in ('direct_checkpoint','residual_checkpoint','density_checkpoint','weak_density_checkpoint')}
    artifacts['fusion_json']=resolve_path(fusion_value,out)
    artifacts['target_config']=(out/'configs/TARGET.json').resolve()
    artifacts['geometry_config']=(ROOT/'research/p1_material_target_20261018/inputs/geometry.json').resolve()
    if m.get('classical_action_json') is not None:
        artifacts['classical_action_json']=resolve_path(m['classical_action_json'],out)
        classical=json.loads(artifacts['classical_action_json'].read_text())['selected']
        if float(classical['sigma'])<0 or not 0<=float(classical['blend'])<=1 or not 0<=float(classical['prior_shrink'])<=1:
            raise ValueError('invalid frozen classical action')
    ablations=m.get('ablation_checkpoints',[])
    if not isinstance(ablations,list):raise ValueError('ablation_checkpoints must be a frozen list')
    seen=set()
    for entry in ablations:
        name=entry.get('name','')
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,63}',name) or name in seen:
            raise ValueError('ablation names must be unique safe identifiers')
        if entry.get('kind') not in ('direct','residual'):raise ValueError('ablation kind must be direct or residual')
        seen.add(name);artifacts['ablation__'+name]=resolve_path(entry['checkpoint'],out)
    fusion=json.loads(artifacts['fusion_json'].read_text())
    selected=fusion['selected']
    artifacts['fusion_checkpoint']=resolve_path(selected['checkpoint'],out)
    if not 0<=float(selected['alpha'])<=1 or selected.get('projection','none') not in ('none','passive_rectangle'):
        raise ValueError('invalid frozen fusion action')
    if not any(selected['name'].startswith(prefix) for prefix in ('direct_','scale2_','residual_')):
        raise ValueError('unknown frozen fusion family')
    records=m.get('frozen_files')
    if not isinstance(records,list):raise ValueError('frozen_files list of path,size_bytes,mtime_ns required')
    registered={str(resolve_path(r['path'],out)):r for r in records}
    for path in set(artifacts.values()):
        stat=path.stat(); record=registered.get(str(path))
        if record is None:raise ValueError('missing frozen artifact registration: '+str(path))
        if stat.st_size!=int(record['size_bytes']) or stat.st_mtime_ns!=int(record['mtime_ns']):
            raise ValueError('frozen artifact bytes/mtime changed: '+str(path))
        if stat.st_mtime>timestamp:raise ValueError('artifact postdates method freeze: '+str(path))
    for parent in parents:
        for path in (out/'inputs/online'/f'{parent}.npz',out/'L/contexts'/f'{parent}.npz'):
            if not path.exists():raise FileNotFoundError(path)
            if path.stat().st_mtime<timestamp:raise ValueError('locked online/context file predates freeze: '+str(path))
        if m.get('include_matched_point176',False):check_matched_receipt(out,parent,timestamp)
    return m,parents,artifacts,selected


def read_online(path):
    # Deliberate per-key reads. No arbitrary dictionary of archive contents.
    with np.load(path,allow_pickle=False) as archive:
        data={key:archive[key] for key in ONLINE_FIELDS}
    expected={'observation':(32,8),'q0':(66,),'chi0':(576,),'J':(512,66),
              'prediction':(32,8),'residual':(32,8)}
    for key,shape in expected.items():
        if data[key].shape!=shape or not np.isfinite(data[key]).all():raise ValueError('invalid online '+key)
    return {key:value[None,...] for key,value in data.items()}


def read_chi(path):
    with np.load(path,allow_pickle=False) as archive:chi=archive['chi'].copy()
    if chi.shape!=(576,) or not np.isfinite(chi).all():raise ValueError('invalid online material point')
    return chi


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_OUT)
    parser.add_argument('--methods',type=Path)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();out=args.root.resolve();methods=args.methods or out/'METHODS_FROZEN_FOR_LOCKED.json'
    torch.set_num_threads(1)
    if torch.get_num_interop_threads()!=1:torch.set_num_interop_threads(1)
    manifest,parents,paths,fusion=frozen_contract(out,methods)
    dest=out/'L/estimates';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE.json').exists():raise FileExistsError('locked inference already complete; no implicit replay')
    cfg=json.loads(paths['target_config'].read_text())
    points=np.asarray(json.loads(paths['geometry_config'].read_text())['points'])
    roi=np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
    if len(roi)!=256 or float(cfg['field_mean'])!=0:raise ValueError('frozen symmetric P1 mean/ROI contract')
    prior_mean=prior_material_mean(cfg,roi)
    classical=json.loads(paths['classical_action_json'].read_text())['selected'] if 'classical_action_json' in paths else None
    start_wall,start_cpu=time.monotonic(),time.process_time();load_wall=time.monotonic()
    models={};load_failures={}
    loaders={'direct_checkpoint':load_predictor,'residual_checkpoint':load_residual_predictor,
             'density_checkpoint':load_density,'weak_density_checkpoint':load_density,
             'fusion_checkpoint':load_residual_predictor if fusion['name'].startswith('residual_') else load_predictor}
    residual_keys={'residual_checkpoint'}
    if fusion['name'].startswith('residual_'):residual_keys.add('fusion_checkpoint')
    cache={}
    for key,loader in loaders.items():
        try:
            if str(paths[key]) not in cache:cache[str(paths[key])]=loader(paths[key],device='cpu')
            models[key]=cache[str(paths[key])]
        except Exception as error:models[key]=None;load_failures[key]=repr(error)
    loading_seconds=time.monotonic()-load_wall
    extra_loading_start=time.monotonic()
    for entry in manifest.get('ablation_checkpoints',[]):
        key='ablation__'+entry['name'];loader=load_residual_predictor if entry['kind']=='residual' else load_predictor
        if entry['kind']=='residual':residual_keys.add(key)
        try:
            if str(paths[key]) not in cache:cache[str(paths[key])]=loader(paths[key],device='cpu')
            models[key]=cache[str(paths[key])]
        except Exception as error:models[key]=None;load_failures[key]=repr(error)
    extra_loading_seconds=time.monotonic()-extra_loading_start
    frame_arrays={}
    for name,key in (('density','density_checkpoint'),('weak_density','weak_density_checkpoint')):
        frame_arrays[name+'_frame']=models[key].frame if models[key] is not None else np.eye(64)
    save(dest/'DENSITY_FRAMES.npz',**frame_arrays)
    write(dest/'INFERENCE_PROTOCOL.json',dict(source_revision=manifest['source_revision'],frozen_unix=manifest['frozen_unix'],
        parents=parents,online_fields=list(ONLINE_FIELDS),truth_read=False,physics_calls=0,density_sample_count=1024,
        density_scope='fixed10 coordinate intervals/proposal samples; no material Bayes mean or full calibration claim',
        model_loading_seconds=loading_seconds,extra_model_loading_seconds=extra_loading_seconds,
        ablation_checkpoints=manifest.get('ablation_checkpoints',[]),load_failures=load_failures,
        include_matched_point176=manifest.get('include_matched_point176',False),
        classical_action=classical,
        frozen_source_paths=manifest.get('frozen_source_paths',['src/p3r2','src/p1_material','tools/*p3r2*.py']),
        fallback='anchor per failed material method; unit prior per failed density method'))
    all_status=[]
    for index,parent in enumerate(parents):
        path=dest/f'{parent}.npz'
        if path.exists():
            if not args.resume:raise FileExistsError(path)
            with np.load(path,allow_pickle=False) as previous:
                if str(previous['source_revision'])!=manifest['source_revision'] or float(previous['frozen_unix'])!=float(manifest['frozen_unix']):
                    raise ValueError('existing estimate from another freeze')
                all_status.append(json.loads(str(previous['status_json'])))
            continue
        before=time.monotonic();status=dict(parent=parent,failures={},baseline_sources={},method_seconds={},neural_wall_seconds=0.,ablation_neural_wall_seconds=0.)
        anchor_source='probit_context'
        try:online=read_online(out/'L/contexts'/f'{parent}.npz');anchor=online['chi0'][0].copy()
        except Exception as error:
            status['failures']['context']=repr(error);online=None
            anchor=None
            for fallback,label in ((out/'L/contexts/cache'/f'{parent}.npz','probit_cache_after_context_failure'),
                                   (out/'L/contexts/details'/parent/'NATIVE.npz','native_fallback')):
                try:anchor=read_chi(fallback);anchor_source=label;break
                except Exception:pass
            if anchor is None:anchor=prior_mean.copy();anchor_source='exactprior_material_mean_fallback'
            status['failures']['anchor']='context failure; retained '+anchor_source
            status['failures']['probit']='context failure; retained '+anchor_source
        status['baseline_sources']['anchor']=anchor_source;status['baseline_sources']['probit']=anchor_source
        arrays={name:anchor.copy() for name in POINT_METHODS}
        arrays['exactprior_material_mean']=prior_mean.copy()
        if classical is not None:
            tick=time.monotonic()
            try:
                arrays['classical_postprocess']=apply_action(anchor,cfg,roi,classical['sigma'],classical['blend'],classical['prior_shrink'])
                if 'probit' in status['failures']:status['failures']['classical_postprocess']='physics baseline failed; action applied to retained '+anchor_source
            except Exception as error:arrays['classical_postprocess']=anchor.copy();status['failures']['classical_postprocess']=repr(error)
            status['method_seconds']['classical_postprocess']=time.monotonic()-tick
        try:
            arrays['native']=read_chi(out/'L/contexts/details'/parent/'NATIVE.npz');status['baseline_sources']['native']='native_point'
        except Exception as error:
            status['failures']['native']=repr(error);status['baseline_sources']['native']='anchor_fallback:'+anchor_source
        def neural(name,key):
            tick=time.monotonic()
            try:
                if online is None or models[key] is None:raise RuntimeError(load_failures.get(key,'online context unavailable'))
                result=models[key].predict_full(online if key in residual_keys else online['observation'])[0]
                if result.shape!=(576,) or not np.isfinite(result).all():raise ValueError('nonfinite or invalid material output')
                return result
            except Exception as error:status['failures'][name]=repr(error);return anchor.copy()
            finally:
                elapsed=time.monotonic()-tick;status['method_seconds'][name]=elapsed;status['neural_wall_seconds']+=elapsed
        arrays['direct']=neural('direct','direct_checkpoint')
        arrays['residual']=neural('residual','residual_checkpoint')
        for name in ('direct','residual'):
            projected='projected_'+name
            if name in status['failures']:
                arrays[projected]=anchor.copy();status['failures'][projected]='source '+name+' failed; anchor fallback'
            else:arrays[projected]=passive_box(arrays[name],cfg,roi)
        if paths['fusion_checkpoint']==paths['direct_checkpoint'] and not fusion['name'].startswith('residual_'):
            source='direct';candidate=arrays['direct']
        elif paths['fusion_checkpoint']==paths['residual_checkpoint'] and fusion['name'].startswith('residual_'):
            source='residual';candidate=arrays['residual']
        else:source='fusion_family';candidate=neural(source,'fusion_checkpoint')
        if source in status['failures']:
            arrays['fused']=anchor.copy();status['failures']['fused']='source '+source+' failed; anchor fallback'
        else:
            if fusion.get('projection','none')=='passive_rectangle':candidate=passive_box(candidate,cfg,roi)
            arrays['fused']=anchor+float(fusion['alpha'])*(candidate-anchor)
        for entry in manifest.get('ablation_checkpoints',[]):
            key='ablation__'+entry['name'];arrays[key]=neural(key,key)
            arrays[key+'_used_anchor_fallback']=np.array(key in status['failures'])
            status['ablation_neural_wall_seconds']+=status['method_seconds'][key]
        for name,key in (('density','density_checkpoint'),('weak_density','weak_density_checkpoint')):
            tick=time.monotonic()
            try:
                if online is None or models[key] is None:raise RuntimeError(load_failures.get(key,'online context unavailable'))
                samples=models[key].sample_labels(online['observation'],1024,seed=2026102300+index+(100000 if name=='weak_density' else 0))[0]
                if samples.shape!=(1024,10) or not np.isfinite(samples).all():raise ValueError('invalid density samples')
                arrays[name+'_used_prior_fallback']=np.array(False)
            except Exception as error:
                status['failures'][name]=repr(error);samples=np.random.default_rng(2026102300+index+(100000 if name=='weak_density' else 0)).normal(size=(1024,10))
                arrays[name+'_used_prior_fallback']=np.array(True)
            arrays[name+'_samples']=samples.astype('float32')
            arrays[name+'_interval90']=np.quantile(samples,[.05,.95],axis=0)
            elapsed=time.monotonic()-tick;status['method_seconds'][name]=elapsed;status['neural_wall_seconds']+=elapsed
        matched=out/'L/contexts/details'/parent/'matched/MATCHED.npz'
        if manifest.get('include_matched_point176',False):
            try:arrays['matched_point176']=read_chi(matched)
            except Exception as error:arrays['matched_point176']=anchor.copy();status['failures']['matched_point176']=repr(error)
        for method in list(POINT_METHODS)+['ablation__'+entry['name'] for entry in manifest.get('ablation_checkpoints',[])]+(['matched_point176'] if manifest.get('include_matched_point176',False) else [])+(['classical_postprocess'] if classical is not None else []):
            arrays[method+'_used_fallback']=np.array(method in status['failures'])
        status['wall_seconds']=time.monotonic()-before
        save(path,**arrays,parent=np.array(parent),roi_indices=roi,status_json=np.array(json.dumps(status,allow_nan=False)),
             source_revision=np.array(manifest['source_revision']),frozen_unix=np.array(manifest['frozen_unix']))
        all_status.append(status)
        print(json.dumps(dict(parent=parent,complete=True,failures=status['failures'],neural_wall_seconds=status['neural_wall_seconds'])),flush=True)
    if len(all_status)!=256:raise RuntimeError('no exclusions: all256 parents required')
    method_names=list(POINT_METHODS)+['ablation__'+entry['name'] for entry in manifest.get('ablation_checkpoints',[])]+(['matched_point176'] if manifest.get('include_matched_point176',False) else [])+(['classical_postprocess'] if classical is not None else [])
    fallback_counts={name:sum(name in row['failures'] for row in all_status) for name in method_names+['context','density','weak_density']}
    baseline_contains_fallbacks=any(fallback_counts[key] for key in ('native','probit','context'))
    write(dest/'COMPLETE.json',dict(status='COMPLETE',parents=parents,count=256,source_revision=manifest['source_revision'],
        frozen_unix=manifest['frozen_unix'],truth_read=False,physics_calls=0,model_loading_seconds=loading_seconds,
        extra_model_loading_seconds=extra_loading_seconds,
        fallback_counts=fallback_counts,primary_baseline_contains_fallbacks=baseline_contains_fallbacks,primary_baseline_nominal_contexts_complete=not baseline_contains_fallbacks,
        primary_baseline_qualification='nominal native96+probit64 comparator contains failures/fallbacks; all256 retained' if baseline_contains_fallbacks else 'all256 nominal native96+probit64 baselines available',
        wall_seconds=time.monotonic()-start_wall,cpu_seconds=time.process_time()-start_cpu,
        neural_wall_seconds=sum(row['neural_wall_seconds'] for row in all_status),
        ablation_neural_wall_seconds=sum(row.get('ablation_neural_wall_seconds',0.) for row in all_status),parent_status=all_status))


if __name__=='__main__':main()
