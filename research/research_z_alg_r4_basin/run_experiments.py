#!/usr/bin/env python3
"""Bounded public-only R4 experiment runner; truth scoring is a separate process."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time
import traceback
import numpy as np

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
HISTORY=ROOT.parent/'research_z_alg_r3_digits'
METHODS=['trf','data_lbfgs','globaldual','basin_bridge','ray_best_only']
V2_METHODS=['uniform_bridge','uniform_best_only','uniform_multistart']
CASES={'holdout7':'obs_c3604cf3555d','holdout9':'obs_42c48e795766','holdout8':'obs_a63df7f20566'}

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def run_one(job):
    from research.research_z_alg_r3_digits.src.solver import load_public
    from research.research_z_alg_r4_basin.src.basin_solver import solve
    label,case,method,public,config,outdir=job
    output=Path(outdir)/(label+'_'+method+'.json')
    begun=time.perf_counter()
    receipt=dict(label=label,case_id=case,method=method,public=str(public),
                 started_utc=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),
                 truth_access='NONE_PUBLIC_INPUT_ONLY')
    dump(output.with_suffix('.started.json'),receipt)
    print(json.dumps(dict(event='START',**receipt)),flush=True)
    try:
        result=solve(load_public(public,HISTORY/'public_data/operators.npz'),method,config)
        result.update(experiment=receipt,scope='EXPOSED_TRUTH_DEVELOPMENT' if not label.startswith('fresh') else 'FRESH_PREDICTION_BEFORE_SCORING')
        dump(output,result)
        trajectory=result['trajectory']
        np.savez_compressed(output.with_suffix('.trajectory.npz'),
            m=np.asarray([r['m'] for r in trajectory],float),
            seconds=np.asarray([r['seconds'] for r in trajectory],float),
            objective=np.asarray([r['objective'] for r in trajectory],float),
            relative_residual=np.asarray([r['relative_residual'] for r in trajectory],float),
            stage=np.asarray([r['stage'] for r in trajectory],int),
            origin=np.asarray([r['origin'] for r in trajectory],str))
        row=dict(label=label,case_id=case,method=method,status=result['status'],
                 wall_seconds=result['wall_seconds'],cpu_seconds=result['cpu_seconds'],
                 relative_residual=None if result['top1'] is None else result['top1']['relative_residual'],
                 objective=None if result['top1'] is None else result['top1']['objective'],
                 candidate_count=len(result['candidates']),trajectory_count=len(trajectory),
                 costcounts=result['costcounts'],accounting=result['accounting'],output=str(output.relative_to(ROOT)))
    except Exception as exc:
        row=dict(label=label,case_id=case,method=method,status='FAILED_EXCEPTION',
                 exception=repr(exc),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-begun,
                 costcounts='UNAVAILABLE_AFTER_EXCEPTION',output=str(output.relative_to(ROOT)))
        dump(output,dict(**row,experiment=receipt))
    print(json.dumps(dict(event='DONE',**row)),flush=True)
    return row

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cases',nargs='+',default=['holdout7'])
    parser.add_argument('--methods',nargs='+',choices=METHODS+V2_METHODS,default=METHODS)
    parser.add_argument('--wall-seconds',type=float,default=60.)
    parser.add_argument('--max-lu',type=int,default=2500)
    parser.add_argument('--max-rhs',type=int)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--run-name',default='development_holdout7_v1')
    parser.add_argument('--factor',action='store_true')
    parser.add_argument('--source-snapshot')
    args=parser.parse_args()
    contract=json.loads((ROOT/'configs/DEVELOPMENT_CONTRACT.json').read_text())
    config=dict(wall_seconds=args.wall_seconds,max_lu=args.max_lu,stages=contract['material_stages'],
                lower=.5,upper=20.,initial_material=contract['initial_material'],
                smooth_weight=contract['regularizer']['smooth'],tv_weight=contract['regularizer']['tv'],
                tv_epsilon=contract['regularizer']['epsilon'])
    if args.max_rhs is not None:config['max_rhs']=args.max_rhs
    outdir=ROOT/'results'/args.run_name
    outdir.mkdir(parents=True,exist_ok=False)
    jobs=[]
    for label in args.cases:
        case=label if args.factor or label.startswith('fresh') else CASES[label]
        public=ROOT/'public_data'/('factor' if args.factor else 'fresh')/(case+'.npz') if args.factor or label.startswith('fresh') else HISTORY/'public_data/holdout'/(case+'.npz')
        for method in args.methods:jobs.append((label,case,method,str(public),config,str(outdir)))
    manifest=dict(schema='r4-run-manifest-v1',run_name=args.run_name,config=config,
                  solver_code='research/research_z_alg_r4_basin/src/basin_solver.py',
                  config_source='research/research_z_alg_r4_basin/configs/DEVELOPMENT_CONTRACT.json',
                  solver_source_snapshot=args.source_snapshot,
                  jobs=[dict(label=j[0],case_id=j[1],method=j[2],public=j[3]) for j in jobs],
                  workers=args.workers,blas_threads=1,truth_access='NONE_PUBLIC_INPUT_ONLY',
                  started_utc=datetime.now(timezone.utc).isoformat(),rows=[])
    dump(outdir/'MANIFEST.json',manifest)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_one,job) for job in jobs]
        for future in as_completed(futures):
            manifest['rows'].append(future.result());dump(outdir/'MANIFEST.json',manifest)
    manifest['completed_utc']=datetime.now(timezone.utc).isoformat()
    dump(outdir/'MANIFEST.json',manifest)

if __name__=='__main__':main()
