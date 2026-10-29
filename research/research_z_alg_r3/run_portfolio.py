#!/usr/bin/env python3
"""V7 guarded portfolio CLI using immutable public-observation V6 solvers."""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import json
from pathlib import Path
from src.portfolio import solve_portfolio


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public',required=True)
    parser.add_argument('--method',choices=['guarded_globaldual','guarded_reciprocal'],required=True)
    parser.add_argument('--config',help='Optional JSON wrapper configuration')
    parser.add_argument('--output',required=True,help='New JSON path; existing outputs are preserved')
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        raise FileExistsError('Preserving existing result: '+str(out))
    config = json.loads(Path(args.config).read_text()) if args.config else {}
    result = solve_portfolio(args.public,args.method,config)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'],method=result['method'],output=str(out),
                          residual=result['top1']['relative_residual'] if result['top1'] else None,
                          wall_seconds=result['wall_seconds'],accounting=result['accounting'])))


if __name__=='__main__':
    main()
