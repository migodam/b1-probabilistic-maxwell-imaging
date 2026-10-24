#!/usr/bin/env python3
"""Atomic zero-physics packing of completed strong TRAIN/validation contexts.

Never generates contexts or trains a model. Missing artifacts produce WAITING
and no output directory. --check-only inspects availability and parent rosters.
"""
from pathlib import Path
import argparse,json,os,sys,time,tempfile,shutil,resource
import numpy as np
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parents[1]
STAGE=ROOT/'research/p3l_r4_conditional_transport'
ONLINE=('observation','q0','chi0','J','prediction','residual','phi','x','H')
SHAPES={'observation':(32,8),'q0':(66,),'chi0':(576,),'J':(512,66),'prediction':(32,8),'residual':(32,8),'phi':(),'x':(66,),'H':(66,66)}


def atomic_npz(path,arrays):
    path=Path(path);temporary=path.with_suffix('.writing')
    with temporary.open('xb') as f:
        np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
    temporary.replace(path)


def packed_y(observation):
    a=np.asarray(observation)
    return np.concatenate((a.real.transpose(0,2,1),a.imag.transpose(0,2,1)),2).reshape(len(a),512).astype('float32')


def expected_parents(split,count):return [f'R2_{split}_{i:05d}' for i in range(count)]


def inspect_contexts(root,split,parents,target_id):
    missing=[];invalid=[];identities=[]
    for parent in parents:
        path=root/split/f'{parent}.npz';side=path.with_suffix('.json')
        if not path.exists() or not side.exists():missing.append(parent);continue
        try:
            m=json.loads(side.read_text())
            if m.get('parent')!=parent or m.get('GT_read') is not False or (m.get('native_cap'),m.get('probit_cap'))!=(96,64) or m.get('target_id')!=target_id:
                raise ValueError('context sidecar must declare correct parent,target,truth-free96+64')
            identities.append(dict(parent=parent,path=str(path.resolve()),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns,sidecar=str(side.resolve())))
        except Exception as e:invalid.append(dict(parent=parent,error=repr(e)))
    return missing,invalid,identities


def pack_split(joint_path,contexts_root,split,parents,identities,directory,chol_scale):
    # TRAIN/validation labels only. No locked or offline-only path is accepted.
    if 'locked' in str(joint_path).lower() or 'offline_only' in str(joint_path):raise ValueError('training joint path only')
    with np.load(joint_path,allow_pickle=False) as z:
        ids=z['parent'].astype(str)
        if len(ids)!=len(set(ids)) or any(not p.startswith(f'R2_{split}_') for p in ids):raise ValueError('joint fold identity/uniqueness')
        lookup={p:i for i,p in enumerate(ids)}
        try:index=np.array([lookup[p] for p in parents])
        except KeyError as e:raise ValueError('requested context parent absent from matching joint fold') from e
        labels_x=z['x'][index].astype('float32');labels_y=z['y'][index].astype('float32');material=z['material'][index].astype('complex64')
        s=z['s'][index].astype('float32') if 's' in z else None
    n=len(parents)
    if labels_x.shape!=(n,66) or labels_y.shape!=(n,512) or material.shape!=(n,576):raise ValueError('joint labels shape')
    if not all(np.isfinite(a).all() for a in (labels_x,labels_y,material)):raise ValueError('nonfinite training labels')
    arrays={};max_y_roundoff=0.;max_context_probit_error=0.
    identity_map={m['parent']:m for m in identities}
    for i,parent in enumerate(parents):
        path=contexts_root/split/f'{parent}.npz';stat=path.stat();before=identity_map[parent]
        if (stat.st_size,stat.st_mtime_ns)!=(before['bytes'],before['mtime_ns']):raise ValueError('context changed during packing')
        with np.load(path,allow_pickle=False) as z:
            if str(z['parent'].item())!=parent:raise ValueError('internal context parent mismatch')
            row={k:z[k] for k in ONLINE}
        for key,shape in SHAPES.items():
            value=row[key]
            if value.shape!=shape or not np.isfinite(value).all():raise ValueError(f'{parent}:invalid {key}')
            if key not in arrays:arrays[key]=np.empty((n,)+shape,dtype='complex128' if np.iscomplexobj(value) else 'float32')
            arrays[key][i]=value
        if not np.iscomplexobj(row['observation']) or not np.iscomplexobj(row['chi0']):raise ValueError('complex context observation/material')
        yy=packed_y(row['observation'][None])[0]
        yerror=float(np.max(np.abs(yy-labels_y[i])));max_y_roundoff=max(max_y_roundoff,yerror)
        if not np.array_equal(yy,labels_y[i]):raise ValueError('context observation differs from original parent joint noisy y(float32 packing)')
        u=ndtr(row['x'][64:]);error=float(np.max(np.abs(u-row['q0'][64:])))
        max_context_probit_error=max(max_context_probit_error,error)
        if error>2e-6 or np.max(np.abs(row['x'][:64]-row['q0'][:64]))>2e-6:raise ValueError('context x must be Gaussianized q0, not truth x')
        if np.max(np.abs(row['residual']-(row['observation']-row['prediction'])))>1e-12:raise ValueError('context residual provenance mismatch')
    arrays['parent']=np.array(parents)
    arrays['chi_GT']=material
    # q_GT is reconstructed from original float32 Gaussianized label x; do not
    # silently label it original native float64 truth from another source.
    arrays['q_GT']=np.c_[labels_x[:,:64].astype('float64'),ndtr(labels_x[:,64:].astype('float64'))]
    arrays['x_GT']=labels_x
    atomic_npz(directory/f'anchor_{split}.npz',arrays)
    density=dict(x=labels_x,y=labels_y,parent=np.array(parents),base_mean=arrays['x'].astype('float32'),
                 context_x=arrays['x'].astype('float32'),context_q0=arrays['q0'].astype('float32'),
                 context_chi0=arrays['chi0'].astype('complex64'),context_residual=arrays['residual'].astype('complex64'),
                 context_J=arrays['J'].astype('float32'),material=material,
                 label_scope=np.array('TRAIN/validation supervision only; no locked labels'),
                 base_provenance=np.array('truth-free strong native96+probit64 Gaussianized context x'),
                 target_id=np.array('P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018'))
    if s is not None:density['s']=s
    if chol_scale is not None:
        density['base_chol']=np.broadcast_to(np.eye(66,dtype='float32')*chol_scale,(n,66,66)).copy()
    # Independent truth-free companion for future conditioner interfaces; no
    # ambiguous x key (density x is truth; online x0 is context_x).
    atomic_npz(directory/f'online_{split}.npz',{k:v for k,v in arrays.items() if k not in ('chi_GT','q_GT','x_GT')})
    atomic_npz(directory/f'density_{split}.npz',density)
    return dict(split=split,n_parent=n,parent_first=parents[0],parent_last=parents[-1],joint_source=str(joint_path.resolve()),
                context_provenance=identities,maximum_packed_y_difference=max_y_roundoff,maximum_context_CDF_error=max_context_probit_error,
                q_GT_provenance='reconstructed from original R4float32 Gaussianized x label; not untouched original native float64 q',
                anchor_arrays={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in arrays.items()},
                density_arrays={k:dict(shape=list(v.shape),dtype=str(v.dtype)) for k,v in density.items()})


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--train-count',type=int,default=1024);ap.add_argument('--validation-count',type=int,default=256)
    ap.add_argument('--contexts-root',type=Path,default=STAGE/'results/contexts')
    ap.add_argument('--joint-train',type=Path,default=STAGE/'inputs/joint/train32768.npz')
    ap.add_argument('--joint-validation',type=Path,default=STAGE/'inputs/joint/validation.npz')
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--base-chol-scale',type=float)
    ap.add_argument('--check-only',action='store_true');a=ap.parse_args()
    wall,cpu=time.monotonic(),time.process_time()
    if not 1<=a.train_count<=2048 or not 1<=a.validation_count<=512:ap.error('authorized TRAIN<=2048 validation<=512, positive counts')
    if a.base_chol_scale is not None and (not np.isfinite(a.base_chol_scale) or a.base_chol_scale<=0 or not np.isfinite(np.float32(a.base_chol_scale)) or np.float32(a.base_chol_scale)<=0):ap.error('strictly positive finite float32-representable fullrank identity scale required')
    if a.joint_train.resolve()==a.joint_validation.resolve():raise ValueError('joint folds must differ')
    if a.out.exists():raise FileExistsError('fresh output directory required; never overwrite previous packs')
    cfg=json.loads((STAGE/'configs/TARGET.json').read_text());rosters={s:expected_parents(s,n) for s,n in [('train',a.train_count),('validation',a.validation_count)]}
    if set(rosters['train'])&set(rosters['validation']):raise ValueError('fold overlap')
    for split,joint in [('train',a.joint_train),('validation',a.joint_validation)]:
        if 'locked' in str(joint).lower() or 'offline_only' in str(joint):raise ValueError('training joint path only')
        with np.load(joint,allow_pickle=False) as z:
            ids=z['parent'].astype(str)
        if len(ids)!=len(set(ids)) or any(not p.startswith(f'R2_{split}_') for p in ids) or not set(rosters[split]).issubset(set(ids)):
            raise ValueError('joint fold identity or requested parent roster mismatch')
    inspections={s:inspect_contexts(a.contexts_root,s,rosters[s],cfg['target_id']) for s in rosters}
    waiting={s:dict(missing=inspections[s][0],invalid=inspections[s][1],requested=len(rosters[s]),ready=len(inspections[s][2])) for s in rosters}
    ready=all(not v['missing'] and not v['invalid'] for v in waiting.values())
    if a.check_only or not ready:
        print(json.dumps(dict(status='READY' if ready else 'WAITING_FOR_COMPLETE_CONTEXTS',splits=waiting,output_written=False,physics_calls=0,training_updates=0)))
        if not ready and not a.check_only:raise SystemExit(3)
        return
    a.out.parent.mkdir(parents=True,exist_ok=True)
    temporary=Path(tempfile.mkdtemp(prefix='.r4_context_pack_',dir=a.out.parent))
    try:
        rows=[]
        for split,joint in [('train',a.joint_train),('validation',a.joint_validation)]:
            rows.append(pack_split(joint,a.contexts_root,split,rosters[split],inspections[split][2],temporary,a.base_chol_scale))
        metadata=dict(status='PACK_COMPLETE',target_id=cfg['target_id'],physics_calls=0,training_updates=0,
            parent_disjointness='CHECKED by fixed R2_train/R2_validation identity; no augmentation/new parents',
            context_budget=dict(packed_native=96,packed_probit=64,historical_native=48,historical_probit=48,
                                context_shift='new common96+64 contexts; historical512/128 used48+48, do not pool as same estimator'),
            schemas=dict(anchor='supervised packed R2-style: original online fields plus chi_GT,q_GT,x_GT; x remains online Gaussianized x0',
                         density='x is original float32 label truth; y is original raw packed512; base_mean/context_x are truth-free x0',
                         online='truth-free exact field allowlist companion'),
            density_compatibility='existing trainer uses base_mean/base_chol only for CFM --base anchor; spline/Gaussian do not auto-consume context_*',
            base_covariance='omitted: caller must explicitly select positive anchor_scale' if a.base_chol_scale is None else f'fullrank fixed identity factor {a.base_chol_scale}*I66; not posterior covariance',
            splits=rows,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
        (temporary/'PACK_COMPLETE.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
        os.rename(temporary,a.out)
    except BaseException:
        shutil.rmtree(temporary,ignore_errors=True);raise
    print(json.dumps(dict(status='PACK_COMPLETE',out=str(a.out.resolve()),train=a.train_count,validation=a.validation_count,physics_calls=0,training_updates=0)))

if __name__=='__main__':main()
