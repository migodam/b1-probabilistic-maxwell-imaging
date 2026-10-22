#!/usr/bin/env python3
"""Train fixed-coordinate10 proper NLL proposals from independent offline parents."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r2.density import KINDS,load_density_data,coordinate_frame,fit_context,train_density
from p3r2.risk_training import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--train',required=True,type=Path);p.add_argument('--validation',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--sizes',nargs='+',type=int,default=[16384]);p.add_argument('--seeds',nargs='+',type=int,default=[11,22,33])
    p.add_argument('--models',nargs='+',choices=KINDS,default=list(KINDS));p.add_argument('--pca-components',type=int,default=32)
    p.add_argument('--epochs',type=int,default=100);p.add_argument('--batch-size',type=int,default=256)
    p.add_argument('--patience',type=int,default=12);p.add_argument('--max-seconds',type=float,default=300)
    p.add_argument('--learning-rate',type=float,default=2e-4);p.add_argument('--weight-decay',type=float,default=1e-4)
    p.add_argument('--device',choices=['cpu','mps'],default='cpu');p.add_argument('--resume',action='store_true')
    p.add_argument('--train-only-frame',type=Path,help='optional NPZ basis[64,8], frozen from train-J only; default first8 prior coordinates')
    a=p.parse_args()
    if not 0<a.epochs<=100 or not 0<a.max_seconds<=300 or min(a.sizes)<1 or a.batch_size<1 or a.pca_components<1:
        p.error('positive sizes/context/batch and epochs<=100,max-seconds<=300 required')
    obs,q,parent=load_density_data(a.train);vo,vq,vp=load_density_data(a.validation)
    if set(parent)&set(vp): raise ValueError('train/validation parents overlap')
    if max(a.sizes)>len(obs): raise ValueError('insufficient independent train parents')
    basis=None
    if a.train_only_frame:
        with np.load(a.train_only_frame,allow_pickle=False) as z: basis=z['basis']
    frame=coordinate_frame(basis);a.output.mkdir(parents=True,exist_ok=True)
    write_json(a.output/'RUN_CONTRACT.json',dict(train=str(a.train.resolve()),validation=str(a.validation.resolve()),
        models=a.models,sizes=a.sizes,seeds=a.seeds,dimension=10,physics_calls=0,locked_truth_read=False,
        coordinates='first8 fixed prior z coords + probit u2' if basis is None else 'supplied frozen train-only orthonormal z8 + probit u2',
        optional_frame_source=str(a.train_only_frame) if basis is not None else None,
        full66='learned10*prior56 with exact probit Jacobian; no local frame or material mean reconstruction',
        criterion='validation proper NLL dz8 d eta2; update0 eligible',posterior_reference='NOT_ESTABLISHED'))
    results=[]
    for size in a.sizes:
        prep_wall,prep_cpu=time.monotonic(),time.process_time()
        context=fit_context(obs[:size],a.pca_components)
        context_path=a.output/f'n{size}';context_path.mkdir(parents=True,exist_ok=True)
        write_json(context_path/'CONTEXT_FIT.json',dict(train_size=size,components=context['components'],
                   cpu_seconds=time.process_time()-prep_cpu,wall_seconds=time.monotonic()-prep_wall,
                   data_used='train observation only; no validation labels or locked data'))
        for kind in a.models:
            for seed in ([a.seeds[0]] if kind=='unit_prior' else a.seeds):
                dest=a.output/f'n{size}'/f'{kind}_s{seed}'
                print(json.dumps(dict(event='start',kind=kind,size=size,seed=seed)),flush=True)
                if a.resume and (dest/'RESULT.json').exists(): result=json.loads((dest/'RESULT.json').read_text())
                else:
                    result=train_density(kind,seed,obs[:size],q[:size],parent[:size],vo,vq,vp,dest,context,frame,
                        a.epochs,a.batch_size,a.patience,a.max_seconds,a.learning_rate,a.weight_decay,a.device)
                results.append(dict(path=str(dest.resolve()),**{k:v for k,v in result.items() if k!='history'}))
                write_json(a.output/'RESULTS.json',dict(results=results))
                print(json.dumps(dict(event='completed',path=str(dest),nll=result['validation_nll'],seconds=result['wall_seconds'])),flush=True)
    best=min(results,key=lambda r:(r['validation_nll'],r['train_size'],r['kind'],r['seed']))
    write_json(a.output/'SELECTED.json',dict(checkpoint=str(Path(best['path'])/'BEST.pt'),result=best,
        criterion='minimum validation proper NLL over declared runs; full66 density includes independent prior56',
        scientific_status='proposal density only; posterior fidelity NOT_ESTABLISHED'))


if __name__=='__main__':main()
