"""Metrics/serialization only. No inverse, forward, generator, or truth imports."""
from pathlib import Path
import csv
import json
import math
import numpy as np


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_new(path, value):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')


def stamp(path):
    path=Path(path).resolve()
    s=path.stat()
    return dict(path=str(path),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def unchanged(entry):
    return stamp(entry['path'])==entry


def public_metrics(path):
    with np.load(path,allow_pickle=False) as data:
        keys=('Y','sigma','noise_radius','host','amplitude','material_bounds',
              'smooth_weight','tv_weight','tv_epsilon','prior_width')
        return {key:data[key].copy() for key in keys}


def prediction(path, shape):
    with np.load(path,allow_pickle=False) as data:
        m=data['m'].copy();yp=data['Y_prediction'].copy()
    if m.shape!=(int(np.prod(shape)),) or yp.shape!=(48,12):
        raise ValueError('Invalid saved prediction dimensions')
    if not np.isfinite(m).all() or not np.isfinite(yp).all():
        raise ValueError('Nonfinite saved prediction')
    return m,yp


def otsu(m):
    # Exact frozen R2 score_blind.py implementation.
    if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
    hist,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2
    weight=np.cumsum(hist);moment=np.cumsum(hist*centers);den=weight*(weight[-1]-weight)
    between=np.full(256,-np.inf);ok=den>0
    between[ok]=(moment[-1]*weight[ok]-moment[ok]*weight[-1])**2/den[ok]
    threshold=float(centers[np.argmax(between)])
    return m>threshold,threshold


def iou(p,t):
    n=np.count_nonzero(p|t)
    return float(np.count_nonzero(p&t)/n) if n else 1.


def prior(m,public,shape):
    # Same physical nearest-neighbor prior; shape parameter supports grid32.
    z=m.reshape(shape)
    d=np.r_[np.diff(z,axis=0).ravel(),np.diff(z,axis=1).ravel()]/float(public['prior_width'])
    eps=float(public['tv_epsilon'])
    return float(float(public['smooth_weight'])*np.mean(d*d)+
                 float(public['tv_weight'])*np.mean(np.sqrt(d*d+eps*eps)-eps))


def physical(m,yp,public,shape):
    res=float(np.linalg.norm(yp-public['Y']))
    relative=res/float(np.linalg.norm(public['Y']))
    return dict(relative_residual=relative,whitened_residual=res/float(public['sigma']),
                physical_objective=.5*relative**2+prior(m,public,shape))


def metrics(m,yp,public,truth,support,shape):
    mask,threshold=otsu(m)
    # Frozen amplitude normalization is 12, not max/mean of truth.
    nrmse=float(np.linalg.norm(m-truth)/np.sqrt(truth.size)/12.)
    overlap=iou(mask,support)
    result=physical(m,yp,public,shape)
    lo,hi=map(float,public['material_bounds'])
    violation=float(max(0.,lo-m.min(),m.max()-hi));legal=violation<=1e-12
    mat=nrmse<=.15;image=mat and overlap>=.7
    data=result['whitened_residual']<=float(public['noise_radius'])
    result.update(amplitude_normalized_nrmse=nrmse,RMSE=12*nrmse,otsu_iou=overlap,
        otsu_threshold=threshold,noise_radius99=float(public['noise_radius']),
        material_gate=bool(mat),image_gate=bool(image),noise_gate=bool(data),
        material_and_data_gate=bool(mat and data and legal),joint_gate=bool(image and data and legal),
        material_legal=bool(legal),material_constraint_violation=violation,
        material_min=float(m.min()),material_max=float(m.max()))
    return result


def sum_known(values):
    values=list(values)
    known=[float(x) for x in values if isinstance(x,(int,float)) and math.isfinite(x)]
    return dict(known_subtotal=sum(known),complete=len(known)==len(values),
                missing_count=len(values)-len(known),total=sum(known) if len(known)==len(values) else None)


def write_csv_new(path,rows):
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key:json.dumps(value,allow_nan=False) if isinstance(value,(list,dict)) else value
                             for key,value in row.items()})


def frozen_contract():
    return dict(metric_source='R2 score_blind.py: RMSE/12, prediction-only Otsu256, public 99% CN radius',
        NRMSE_gate=.15,Otsu_IoU_gate=.7,material_box=[.02,20.],box_numeric_tolerance=1e-12,
        source_count=12,receiver_count=48,actual_GPU_tasks=0,actual_GPU_workers=0,
        prior='smooth1e-4/TV1e-4,width19.5,epsilon.01',new_gates=False,hash_checks=False)
