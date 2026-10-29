#!/usr/bin/env python3
"""Run R3 from a public NPZ containing geometry/operators,y,box only."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import json
from pathlib import Path
from src import solve_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--public', required=True, help='Public NPZ with complex y plus operators or geometry')
    parser.add_argument('--method', choices=['algebraic','fullj','coordinate','globaldual'], default='algebraic')
    parser.add_argument('--config', help='JSON configuration path')
    parser.add_argument('--output', required=True, help='New result JSON path; existing files are preserved')
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        raise FileExistsError('Preserving existing result: '+str(out))
    config = json.loads(Path(args.config).read_text()) if args.config else {}
    result = solve_case(args.public, args.method, config)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], method=result['method'], output=str(out),
                          top1_relative_residual=result['top1']['relative_residual'] if result['top1'] else None,
                          wall_seconds=result['wall_seconds'], accounting=result['accounting'])))


if __name__ == '__main__':
    main()
