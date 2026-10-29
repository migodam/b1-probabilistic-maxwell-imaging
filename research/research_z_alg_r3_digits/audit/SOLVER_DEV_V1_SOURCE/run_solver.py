#!/usr/bin/env python3
"""Run an R3-DIGITS method from public observations/operators only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,json
from pathlib import Path
from src import solve_case


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--public',required=True)
    p.add_argument('--operators')
    p.add_argument('--method',choices=['fullj_native','fullj_reciprocal','globaldual','fullj','reciprocal'],default='fullj_native')
    p.add_argument('--config')
    p.add_argument('--output',required=True)
    a=p.parse_args()
    out=Path(a.output)
    if out.exists():raise FileExistsError('Preserving existing result '+str(out))
    cfg=json.loads(Path(a.config).read_text()) if a.config else {}
    if a.method in cfg and isinstance(cfg[a.method],dict):cfg=cfg[a.method]
    result=solve_case(a.public,a.method,cfg,operators=a.operators)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(method=result['method'],status=result['status'],output=str(out),
                         relative_residual=result['top1']['relative_residual'] if result['top1'] else None,
                         wall_seconds=result['wall_seconds'],cpu_seconds=result['cpu_seconds'],
                         all_unknowns_released=result['all_unknowns_released'],accounting=result['accounting'])))


if __name__=='__main__':main()
