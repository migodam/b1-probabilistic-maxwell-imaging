"""One disclosed development replay on the existing Windows NumPy/SciPy stack.

The purpose is cross-platform repeatability of the critical positive, not a new
sample or algorithm calibration. Public observation only; frozen solver API.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
sys.path.insert(0,str(ROOT.parents[1]/'vendor'))
from research.research_max_struct_r4.src.confirmation_solver import solve

def main():
    p=argparse.ArgumentParser();p.add_argument('--public',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    with np.load(a.public,allow_pickle=False) as f:public={k:f[k] for k in f.files}
    cfg=json.loads(a.config.read_text())['methods']['onehop']
    result=solve(public,'onehop',cfg,backend='cpu')
    arrays=result.pop('trajectory_arrays');np.savez_compressed(a.output/'prediction.npz',**arrays)
    result.update(case_id='obs_b20d5a55677d4eca',classification='DISCLOSED_DEVELOPMENT_CROSS_PLATFORM_VALIDATION',
        caller_wall_seconds=time.perf_counter()-start,reason='Different existing NumPy/SciPy versions in local and remote execution; assess positive repeatability without changing algorithm')
    (a.output/'RECEIPT.json').write_text(json.dumps(result,indent=2,default=lambda x:x.tolist(),allow_nan=False))
    print(json.dumps({k:result.get(k) for k in ('status','case_id','physical_objective','physical_relative_residual','wall_seconds')}),flush=True)
if __name__=='__main__':main()
