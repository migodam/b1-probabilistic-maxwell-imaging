#!/usr/bin/env python3
"""Run one R4 solver using public measurements only."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from research.research_z_alg_r3_digits.src.solver import load_public
from research.research_z_alg_r4_basin.src.basin_solver import solve


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--public',required=True)
    parser.add_argument('--operators')
    parser.add_argument('--method',choices=['data_lbfgs','globaldual','trf','basin_bridge','ray_best_only','uniform_bridge','uniform_best_only','uniform_multistart'],default='basin_bridge')
    parser.add_argument('--config')
    parser.add_argument('--wall-seconds',type=float)
    parser.add_argument('--max-lu',type=int)
    parser.add_argument('--max-rhs',type=int)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    config=json.loads(Path(args.config).read_text()) if args.config else {}
    if args.wall_seconds is not None:config['wall_seconds']=args.wall_seconds
    if args.max_lu is not None:config['max_lu']=args.max_lu
    if args.max_rhs is not None:config['max_rhs']=args.max_rhs
    result=solve(load_public(args.public,args.operators),args.method,config)
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(method=result['method'],status=result['status'],wall_seconds=result['wall_seconds'],
                         objective=None if result['top1'] is None else result['top1']['objective'],
                         candidates=len(result['candidates']),costcounts=result['costcounts'],output=str(output))))

if __name__=='__main__':main()
