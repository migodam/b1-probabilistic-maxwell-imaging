"""Independent parent-level risk statistics; no model fitting or physics calls."""
import json,os,time
from pathlib import Path
import numpy as np

ONLINE_FIELDS=('observation','q0','chi0','J','residual')
CI_LEVEL=0.9833333333333333
BOOTSTRAP_SEED=20261025099


def write_new_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def checkpoint_identity(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def read_online_measurement(path,parent=None):
    """Observation-only raw512 loader; no anchor/J/GT array is materialized."""
    path=Path(path)
    if 'offline_only' in path.resolve().parts:raise ValueError('offline_only prohibited online')
    with np.load(path,allow_pickle=False) as z:
        key='observation' if 'observation' in z else 'y'
        value=z[key]
        single=(key=='observation' and value.shape==(32,8)) or (key=='y' and value.shape==(512,))
        if single:value=value[None]
        else:
            if parent is None:raise ValueError('packed measurements require parent identity')
            ids=z['parent'].astype(str);rows=np.flatnonzero(ids==parent)
            if len(rows)!=1:raise ValueError('one matching packed measurement parent required')
            value=value[rows]
    if value.shape==(1,32,8) and np.iscomplexobj(value):
        value=np.concatenate((value.real.transpose(0,2,1),value.imag.transpose(0,2,1)),axis=2).reshape(1,512)
    if value.shape!=(1,512) or np.iscomplexobj(value) or not np.isfinite(value).all():
        raise ValueError('finite source-major measurement512 required')
    return value.astype('float32')


def read_online_context(path,parent=None):
    """Only five online arrays are materialized. Never enumerates/reads GT arrays.

    Single-parent and packed contexts both supported. Packed parent is identity
    metadata; no q_GT/chi_GT/x_GT/material field is accessed.
    """
    with np.load(path,allow_pickle=False) as z:
        obs=z['observation']
        if obs.ndim==2:
            data={k:z[k][None] for k in ONLINE_FIELDS}
        else:
            if parent is None:raise ValueError('packed context requires parent ID')
            ids=z['parent'].astype(str);where=np.flatnonzero(ids==parent)
            if len(where)!=1:raise ValueError('parent must occur once in context')
            i=int(where[0]);data={k:z[k][i:i+1] for k in ONLINE_FIELDS}
    expected={'observation':(1,32,8),'residual':(1,32,8),'q0':(1,66),'chi0':(1,576),'J':(1,512,66)}
    for key,shape in expected.items():
        if data[key].shape!=shape or not np.isfinite(data[key]).all():raise ValueError(f'invalid online {key}')
    if not np.iscomplexobj(data['observation']) or not np.iscomplexobj(data['chi0']):raise ValueError('complex observation and material required')
    return data


def read_anchor_fallback(path,parent=None):
    """Fallback material only; cannot rescue a missing anchor."""
    with np.load(path,allow_pickle=False) as z:
        chi=z['chi0'];q=z['q0']
        if chi.ndim==1:chi,q=chi[None],q[None]
        else:
            idx=np.flatnonzero(z['parent'].astype(str)==parent)
            if len(idx)!=1:raise ValueError('missing unique fallback parent')
            chi,q=chi[idx],q[idx]
    if chi.shape!=(1,576) or q.shape!=(1,66) or not np.isfinite(chi).all() or not np.isfinite(q).all():raise ValueError('invalid fallback anchor')
    return chi,q


def material_losses(prediction,truth,roi,anchor_phase=None,truth_phase=None):
    pred=np.asarray(prediction);truth=np.asarray(truth)
    if pred.shape!=truth.shape or pred.ndim!=2 or pred.shape[1]!=576:raise ValueError('full576 complex fields required')
    error=np.abs(pred-truth);sq=error**2
    result={'full_mse':sq.mean(1),'full_mae':error.mean(1),'roi_mse':sq[:,roi].mean(1),'roi_mae':error[:,roi].mean(1)}
    for prefix,phase in (('anchor',anchor_phase),('truth',truth_phase)):
        if phase is None:continue
        for name,mask in (('host',phase<=-1),('transition',(phase>-1)&(phase<1)),('pure',phase>=1)):
            count=mask.sum(1);den=np.maximum(count,1)
            # Contribution uses full576 denominator and sums to full ROI error.
            result[f'{prefix}_{name}_count']=count
            result[f'{prefix}_{name}_mse']=np.where(count>0,(sq[:,roi]*mask).sum(1)/den,np.nan)
            result[f'{prefix}_{name}_mae']=np.where(count>0,(error[:,roi]*mask).sum(1)/den,np.nan)
            result[f'{prefix}_{name}_mse_contribution']=(sq[:,roi]*mask).sum(1)/576
    return result


def paired_bootstrap(candidate,reference,draws=4000,seed=BOOTSTRAP_SEED,level=CI_LEVEL):
    """Resample independent parents; seed runs must already be mean losses."""
    a=np.asarray(candidate,dtype=float);b=np.asarray(reference,dtype=float)
    if a.shape!=b.shape or a.ndim!=1 or len(a)<2 or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('paired finite parent losses required')
    rng=np.random.default_rng(seed);diff=[];rel=[]
    for start in range(0,draws,100):
        ids=rng.integers(0,len(a),size=(min(100,draws-start),len(a)))
        am=a[ids].mean(1);bm=b[ids].mean(1);diff.extend(bm-am);rel.extend((bm-am)/np.maximum(bm,1e-300))
    alpha=(1-level)/2
    return dict(n_parent=len(a),candidate_mean=float(a.mean()),reference_mean=float(b.mean()),
                improvement_absolute=float(b.mean()-a.mean()),improvement_relative=float((b.mean()-a.mean())/max(b.mean(),1e-300)),
                absolute_CI=np.quantile(diff,[alpha,1-alpha]).tolist(),relative_CI=np.quantile(rel,[alpha,1-alpha]).tolist(),
                CI_level=level,bootstrap_draws=draws,bootstrap_seed=seed,unit='independent parent; seed losses averaged before resampling')


def family_parent_mean(losses,methods,family,action,metric='full_mse'):
    selected=[m for m in methods if m['family']==family]
    if not selected:raise ValueError(f'no family {family}')
    seed=[m['seed'] for m in selected]
    if len(seed)!=len(set(seed)):raise ValueError('duplicate family training seed')
    values=np.stack([losses[m['id']][action][metric] for m in selected])
    return values.mean(0),selected
