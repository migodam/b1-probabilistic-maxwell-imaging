"""Finite local expansion queue, scheduled after the two VAL physics workers."""
from pathlib import Path
import json, os, subprocess, sys, time
S=Path(__file__).resolve().parents[1]
env=os.environ.copy()
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    env[k]='1'
while not (S/'evaluation/WEAKREG45_VAL_SCORES.json').exists():
    time.sleep(30)
start=time.time()
cmd=[sys.executable,str(S/'experiments/run_anchors.py'),'--id-file',str(S/'data/manifests/train_expansion_local.json'),'--workers','2']
with (S/'audit/EXPAND1024_LOCAL.log').open('w') as f:
    p=subprocess.run(cmd,cwd=S,env=env,stdout=f,stderr=subprocess.STDOUT)
(S/'audit/EXPAND1024_LOCAL_DONE.json').write_text(json.dumps(dict(exit_code=p.returncode,start_unix=start,end_unix=time.time(),wall_seconds=time.time()-start,command=cmd)))
