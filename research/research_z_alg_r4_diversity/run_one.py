"""Run one archived anonymous observation using the frozen solver API."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse,json,pathlib
import numpy as np
from src.solver import solve

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--observation',required=True)
    ap.add_argument('--method',required=True,choices=('native_trf','reciprocal_trf','data_lbfgs','globaldual'))
    ap.add_argument('--backend',default='cpu',choices=('cpu','cuda'))
    ap.add_argument('--operators',default='public_data/operators.npz')
    ap.add_argument('--contract',default='configs/CONTRACT.json')
    ap.add_argument('--output',required=True)
    a=ap.parse_args()
    with np.load(a.operators,allow_pickle=False) as f:public={k:f[k] for k in f.files}
    with np.load(a.observation,allow_pickle=False) as f:public.update({k:f[k] for k in f.files})
    cfg=json.loads(pathlib.Path(a.contract).read_text())['solver']
    result=solve(public,a.method,cfg,backend=a.backend)
    pathlib.Path(a.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
