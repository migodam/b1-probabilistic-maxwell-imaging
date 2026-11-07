"""One-shot separately costed weakreg45 VAL comparison after mainVAL finishes."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time
REPO=Path(__file__).resolve().parents[3]
STUDY=Path(__file__).resolve().parents[1]
ids=[r['obs_id']for r in json.loads((STUDY/'data/manifests/online.json').read_text())['records']if r['split']=='val']
if len(ids)!=128:raise ValueError('Expected frozen128VAL observations')
last=None
while True:
    count=sum((STUDY/'anchors'/f'{i}_weakreg1e-5.npz').exists()and(STUDY/'anchors'/f'{i}_weakreg1e-5.json').exists()for i in ids)
    main_done=(STUDY/'evaluation/PILOT256_VAL_SCORES.json').exists()
    status=(count,main_done)
    if status!=last:
        print(json.dumps(dict(phase='WAITING_WEAKREG45_VAL_AND_MAIN_COMPLETION',complete=count,main_done=main_done,unix=time.time())),flush=True);last=status
    if count==128 and main_done:break
    time.sleep(60)
env=os.environ.copy()
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='1'
cli=STUDY/'experiments/run_evaluation.py'
children=[];streams=[]
try:
    for partition in (0,1):
        stream=(STUDY/f'evaluation/WEAKREG45_VAL_EXECUTION_P{partition}.log').open('a');streams.append(stream)
        command=[sys.executable,str(cli),'evaluate','--controls','weakreg45_anchor','--splits','val','--partition-count','2','--partition-index',str(partition)]
        children.append(subprocess.Popen(command,cwd=REPO,env=env,stdout=stream,stderr=subprocess.STDOUT))
    codes=[child.wait()for child in children]
    if any(c!=0 for c in codes):raise RuntimeError(f'Weakreg45 partitions failed:{codes}')
finally:
    for stream in streams:stream.close()
seal=STUDY/'evaluation/WEAKREG45_VAL_SEALED.json';scores=STUDY/'evaluation/WEAKREG45_VAL_SCORES.json'
for command in ([sys.executable,str(cli),'seal','--controls','weakreg45_anchor','--splits','val','--output',str(seal)],
                [sys.executable,str(cli),'score','--seal',str(seal),'--output',str(scores)]):
    subprocess.run(command,cwd=REPO,env=env,check=True)
print(json.dumps(dict(phase='WEAKREG45_VAL_COMPLETE_PARENT_REVIEW_REQUIRED',scores=str(scores),no_test_truth=True)),flush=True)
