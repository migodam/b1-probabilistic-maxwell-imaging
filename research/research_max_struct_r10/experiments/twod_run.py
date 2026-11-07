"""One R10 public task; mature R8 implementation, no GT or NN access."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest',default=str(ROOT/'configs/twod_tasks.json'))
    ap.add_argument('--task-id',required=True)
    ap.add_argument('--config',default=str(ROOT/'configs/twod_solver.json'))
    ap.add_argument('--output-dir')
    ap.add_argument('--r8-repo',default=str(REPO))
    ap.add_argument('--dpcs-rule',default='noise_ball_99',choices=['noise_ball_99'])
    ap.add_argument('--prepare-only',action='store_true')
    ap.add_argument('--cpu',type=int,help='Affinity for this worker only, before importing NumPy')
    ap.add_argument('--vendor-path',help='Existing threadpoolctl vendor, no installation')
    args=ap.parse_args()
    if args.cpu is not None:
        import psutil
        psutil.Process().cpu_affinity([args.cpu])
        if os.name=='nt':psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    if args.vendor_path:sys.path.insert(0,str(Path(args.vendor_path).resolve()))
    cfg=json.loads(Path(args.config).read_text())
    if cfg['status']!='FROZEN_FOR_INVERSION':raise ValueError('Root has not frozen inverse contract')
    sys.path.insert(0,str(ROOT/'src'))
    from twod_adapter import execute
    receipt,out=execute(args)
    print(json.dumps(dict(task_id=args.task_id,status=receipt['status'],output_dir=str(out),
        truth_access='NONE',scoring='NOT_RUN',CPU_seconds=receipt.get('r10_total_cpu_seconds',receipt.get('cpu_seconds')))),flush=True)
    return 1 if receipt['status'].startswith('FAILED') else 0


if __name__=='__main__':raise SystemExit(main())
