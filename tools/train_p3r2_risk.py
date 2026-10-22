#!/usr/bin/env python3
"""Train Track A material risk models from already generated, disjoint splits."""
import os
for thread_variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[thread_variable] = '1'
import argparse
import json
import sys
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from p3r2.risk_training import load_data, train_variant, write_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--train', required=True, type=Path)
    p.add_argument('--validation', required=True, type=Path)
    p.add_argument('--target', type=Path, default=ROOT/'research/p3l_local_posterior_imaging/configs/TARGET.json')
    p.add_argument('--geometry', type=Path, default=ROOT/'research/p1_material_target_20261018/inputs/geometry.json')
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--sizes', type=int, nargs='+', default=[4096,16384])
    p.add_argument('--seeds', type=int, nargs='+', default=[11,22,33])
    p.add_argument('--models', nargs='+', choices=['mean','ridge','mlp_small','mlp_medium','cnn','cnn_context'],
                   default=['mean','ridge','mlp_small','mlp_medium','cnn','cnn_context'])
    p.add_argument('--epochs', type=int, default=100)
    p.add_argument('--batch-size', type=int, default=256)
    p.add_argument('--patience', type=int, default=12)
    p.add_argument('--max-seconds', type=float, default=600)
    p.add_argument('--learning-rate', type=float, default=1e-3)
    p.add_argument('--weight-decay', type=float, default=1e-4)
    p.add_argument('--device', choices=['cpu','mps'], default='cpu')
    p.add_argument('--resume', action='store_true', help='skip existing completed variant RESULT.json files')
    args = p.parse_args()
    if args.epochs > 100 or args.max_seconds > 600:
        p.error('frozen safety bound: epochs<=100 and max-seconds<=600 per variant')
    if min(args.sizes) < 1 or not args.seeds or args.epochs < 1 or args.max_seconds <= 0 or args.batch_size < 1:
        p.error('positive sizes, epochs, time limit and batch size required')
    obs, chi, parents = load_data(args.train)
    vo, vc, vp = load_data(args.validation)
    if set(parents) & set(vp):
        raise ValueError('train and validation parents overlap')
    if max(args.sizes)>len(obs):
        raise ValueError('requested training size exceeds available parents')
    cfg = json.loads(args.target.read_text())
    points = np.asarray(json.loads(args.geometry.read_text())['points'])
    if points.shape != (576,2):
        raise ValueError('frozen geometry must have 576 two-dimensional points')
    roi = np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
    if len(roi)!=256:
        raise ValueError(f'frozen ROI requires 256 cells, got {len(roi)}')
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output/'RUN_CONTRACT.json', dict(train=str(args.train.resolve()),validation=str(args.validation.resolve()),
               target=str(args.target.resolve()),geometry=str(args.geometry.resolve()),roi_indices=roi.tolist(),
               sizes=args.sizes,seeds=args.seeds,models=args.models,locked_truth_read=False,physics_calls=0,
               subset_policy='first N independent train parents; identical across model and seed',
               selection='validation only; supervised mean estimator is not a calibrated posterior'))
    results=[]; started=time.monotonic()
    for size in args.sizes:
        for model in args.models:
            # Deterministic baselines do not depend on random seeds.
            for seed in ([args.seeds[0]] if model in ('mean','ridge') else args.seeds):
                dest=args.output/f'n{size}'/f'{model}_s{seed}'
                if args.resume and (dest/'RESULT.json').exists():
                    result=json.loads((dest/'RESULT.json').read_text())
                else:
                    print(json.dumps(dict(event='start',model=model,seed=seed,size=size)),flush=True)
                    result=train_variant(model,seed,obs[:size],chi[:size],parents[:size],vo,vc,vp,cfg,roi,dest,
                           epochs=args.epochs,batch_size=args.batch_size,patience=args.patience,max_seconds=args.max_seconds,
                           learning_rate=args.learning_rate,weight_decay=args.weight_decay,device=args.device)
                results.append(dict(path=str(dest.resolve()),**{k:v for k,v in result.items() if k!='history'}))
                write_json(args.output/'RESULTS.json',dict(results=results,wall_seconds=time.monotonic()-started))
                print(json.dumps(dict(event='completed',path=str(dest),risk=result['validation'],seconds=result['wall_seconds'])),flush=True)
    # Freeze the sole downstream model by validation only, with a deterministic tie break.
    selected=min(results,key=lambda r:(r['validation']['complex_roi_mse'],r['train_size'],r['kind'],r['seed']))
    write_json(args.output/'SELECTED.json',dict(checkpoint=str(Path(selected['path'])/'BEST.pt'),
               criterion='minimum validation complex ROI MSE across declared runs; locked test untouched',
               result=selected,calibrated_posterior='NOT_ESTABLISHED'))


if __name__=='__main__':
    main()
