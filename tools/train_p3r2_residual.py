#!/usr/bin/env python3
"""Train anchor residual/geometry/decoder material-risk ablations without physics."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3r2.residual_risk import MODELS,DEFAULT_MODELS,read_anchor,subset,fit_geometry,train_residual
from p3r2.risk_training import write_json


def main(argv=None,default_models=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--train',required=True,type=Path); p.add_argument('--validation',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--target',type=Path,default=ROOT/'research/p3l_local_posterior_imaging/configs/TARGET.json')
    p.add_argument('--geometry',type=Path,default=ROOT/'research/p1_material_target_20261018/inputs/geometry.json')
    p.add_argument('--sizes',nargs='+',type=int,default=[512])
    p.add_argument('--seeds',nargs='+',type=int,default=[11,22,33])
    p.add_argument('--models',nargs='+',choices=MODELS,default=list(DEFAULT_MODELS if default_models is None else default_models))
    p.add_argument('--epochs',type=int,default=100); p.add_argument('--batch-size',type=int,default=128)
    p.add_argument('--patience',type=int,default=15); p.add_argument('--max-seconds',type=float,default=600)
    p.add_argument('--learning-rate',type=float,default=5e-4); p.add_argument('--weight-decay',type=float,default=1e-3)
    p.add_argument('--weak-rank',type=int,default=16); p.add_argument('--strong-rank',type=int,default=12)
    p.add_argument('--j-unwhitened',action='store_true',help='J is raw measurement derivative; whiten it internally')
    p.add_argument('--device',choices=['cpu','mps'],default='cpu'); p.add_argument('--resume',action='store_true')
    a=p.parse_args(argv)
    if not 0<a.epochs<=100 or not 0<a.max_seconds<=600 or min(a.sizes)<1 or a.batch_size<1:
        p.error('positive sizes/batch and epochs<=100, max-seconds<=600 required')
    if min(a.weak_rank,a.strong_rank)<1 or a.weak_rank+a.strong_rank>64:
        p.error('disjoint weak/strong ranks must total at most64')
    tr,va=read_anchor(a.train,True),read_anchor(a.validation,True)
    if set(tr['parent'])&set(va['parent']): raise ValueError('train/validation parent overlap')
    if max(a.sizes)>len(tr['parent']): raise ValueError('not enough independent training parents')
    cfg=json.loads(a.target.read_text()); points=json.loads(a.geometry.read_text())['points']
    a.output.mkdir(parents=True,exist_ok=True)
    write_json(a.output/'RUN_CONTRACT.json',dict(train=str(a.train.resolve()),validation=str(a.validation.resolve()),
        models=a.models,sizes=a.sizes,seeds=a.seeds,physics_calls=0,locked_truth_read=False,
        geometry_frame_fit='train J only, per training size',j_is_whitened=not a.j_unwhitened,
        subset_policy='first N independent train parents; no augmentation counted as new parent',
        criterion='validation actual complex ROI material MSE; all update-zero checkpoints eligible',
        density='NOT_CLAIMED'))
    results=[]
    for n in a.sizes:
        train=subset(tr,n)
        prep_wall,prep_cpu=time.monotonic(),time.process_time()
        geometry=fit_geometry(train,cfg,points,a.weak_rank,a.strong_rank,not a.j_unwhitened)
        geometry_path=a.output/f'n{n}';geometry_path.mkdir(parents=True,exist_ok=True)
        write_json(geometry_path/'GEOMETRY_FIT.json',dict(train_size=n,weak_rank=a.weak_rank,strong_rank=a.strong_rank,
                   cpu_seconds=time.process_time()-prep_cpu,wall_seconds=time.monotonic()-prep_wall,
                   data_used='train J only; prior decoder A from frozen config/geometry',
                   gram_eigenvalues=geometry['eigenvalues'].tolist()))
        if len(geometry['roi'])!=256: raise ValueError('frozen ROI must contain256 cells')
        for kind in a.models:
            for seed in ([a.seeds[0]] if kind in ('anchor','residual_ridge') else a.seeds):
                dest=a.output/f'n{n}'/f'{kind}_s{seed}'
                print(json.dumps(dict(event='start',kind=kind,seed=seed,n=n)),flush=True)
                if a.resume and (dest/'RESULT.json').exists(): result=json.loads((dest/'RESULT.json').read_text())
                else:
                    result=train_residual(kind,seed,train,va,geometry,dest,a.epochs,a.batch_size,a.patience,
                        a.max_seconds,a.learning_rate,a.weight_decay,a.device)
                results.append(dict(path=str(dest.resolve()),**{k:v for k,v in result.items() if k!='history'}))
                write_json(a.output/'RESULTS.json',dict(results=results))
                print(json.dumps(dict(event='completed',path=str(dest),validation=result['validation'],
                                      seconds=result['wall_seconds'])),flush=True)
    best=min(results,key=lambda r:(r['validation']['complex_roi_mse'],r['train_size'],r['kind'],r['seed']))
    write_json(a.output/'SELECTED.json',dict(checkpoint=str(Path(best['path'])/'BEST.pt'),result=best,
               selection='minimum validation complex ROI MSE; locked truth never accessed',density='NOT_CLAIMED'))


if __name__=='__main__': main()
