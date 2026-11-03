"""Portable single-task replay of published observations; no generator or GT.

This is a post-campaign replay convenience, not the original Windows executor.
Use a new output directory; published data never become blind again.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r4.src.confirmation_solver import solve

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--case',required=True)
    p.add_argument('--method',required=True,choices=['onehop','global','native_extended'])
    p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    with np.load(ROOT/'data_public/confirmation/operators.npz',allow_pickle=False) as f:
        public={k:f[k] for k in f.files}
    with np.load(ROOT/'data_public/confirmation'/(args.case+'.npz'),allow_pickle=False) as f:
        public.update({k:f[k] for k in f.files})
    cfg=json.loads((ROOT/'configs/confirmation_frozen.json').read_text())['methods'][args.method]
    result=solve(public,args.method,cfg)
    np.savez_compressed(args.output/'prediction.npz',**result.pop('trajectory_arrays'))
    result.update(case_id=args.case,classification='PUBLIC_REPLAY_NOT_BLIND',caller_wall_seconds=time.perf_counter()-start)
    (args.output/'RECEIPT.json').write_text(json.dumps(result,indent=2,default=lambda x:x.tolist(),allow_nan=False))
    print(json.dumps({k:result.get(k) for k in ('case_id','method','status','physical_relative_residual','wall_seconds')}))
if __name__=='__main__':main()
