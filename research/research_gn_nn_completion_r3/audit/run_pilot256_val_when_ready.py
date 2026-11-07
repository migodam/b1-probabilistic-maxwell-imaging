"""One-shot finite VAL evaluator, waiting for the authorized frozen input cohort."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

REPO=Path(__file__).resolve().parents[3]
STUDY=Path(__file__).resolve().parents[1]
RUNS=[f'pilot256_{a}_seed7'for a in ('v_image','v_physics','j_full','j_null','j_null_current')]
CONTROLS=['continue_trf','weakreg_trf']
ENV=os.environ.copy()
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    ENV[name]='1'
manifest=json.loads((STUDY/'data/manifests/online.json').read_text())
ids=[r['obs_id']for r in manifest['records']if r['split']=='val']
if len(ids)!=128 or len(ids)!=len(set(ids)):
    raise ValueError('Expected fixed128-VAL cohort')
last=None
while True:
    counts={'anchors':sum((STUDY/'anchors'/f'{i}.npz').exists()and(STUDY/'anchors'/f'{i}.json').exists()for i in ids)}
    for run in RUNS:
        directory=STUDY/'predictions'/run
        counts[run]=sum((directory/f'{i}.npz').exists()and(directory/f'{i}.json').exists()for i in ids)
    if counts!=last:
        print(json.dumps(dict(phase='WAITING_FROZEN_VAL_INPUTS',unix=time.time(),counts=counts)),flush=True)
        last=counts
    if all(v==128 for v in counts.values()):
        # A setup receipt is written before predictions, but the entire128-file
        # cohort is already present here. Check provenance/accounting only.
        for run in RUNS:
            directory=STUDY/'predictions'/run
            setup=json.loads((directory/'inference_setup.json').read_text())
            if setup.get('dev_only',False):
                raise ValueError('DEV checkpoint cannot enter pilot256VAL')
            for obs_id in ids:
                receipt=json.loads((directory/f'{obs_id}.json').read_text())
                if not receipt['preparation_cost'].get('complete',False):
                    raise ValueError(f'Incomplete incremental preparation receipt: {run}/{obs_id}')
        break
    time.sleep(60)
cli=STUDY/'experiments/run_evaluation.py'
log=STUDY/'evaluation/PILOT256_VAL_EXECUTION.log'
log.parent.mkdir(parents=True,exist_ok=True)
steps=[
    [sys.executable,str(cli),'seal','--runs',*RUNS,'--controls',*CONTROLS,'--splits','val','--output',str(STUDY/'evaluation/PILOT256_VAL_SEALED.json')],
    [sys.executable,str(cli),'score','--seal',str(STUDY/'evaluation/PILOT256_VAL_SEALED.json'),'--output',str(STUDY/'evaluation/PILOT256_VAL_SCORES.json'),'--oracle'],
]
# Input readiness implies the local anchor cohort is finished. Two finite
# BLAS1 children own disjoint manifest-index partitions; no shared output files.
children=[]
streams=[]
try:
    for partition in (0,1):
        command=[sys.executable,str(cli),'evaluate','--runs',*RUNS,'--controls',*CONTROLS,'--splits','val','--partition-count','2','--partition-index',str(partition)]
        stream=(STUDY/f'evaluation/PILOT256_VAL_EXECUTION_P{partition}.log').open('a')
        streams.append(stream)
        print(json.dumps(dict(phase='RUNNING_TWO_WORKERS',partition=partition,command=command)),flush=True)
        children.append(subprocess.Popen(command,cwd=REPO,env=ENV,stdout=stream,stderr=subprocess.STDOUT))
    codes=[child.wait()for child in children]
    if any(c!=0 for c in codes):
        raise RuntimeError(f'VAL partitions failed: {codes}; preserve logs and all existing outputs')
finally:
    for stream in streams:stream.close()
with log.open('a')as stream:
    for command in steps:
        print(json.dumps(dict(phase='RUNNING',command=command)),flush=True)
        subprocess.run(command,cwd=REPO,env=ENV,stdout=stream,stderr=subprocess.STDOUT,check=True)
sys.path.insert(0,str(REPO))
from research.research_gn_nn_completion_r3.src.score import validation_gate
outputs=validation_gate(STUDY/'evaluation/PILOT256_VAL_SCORES.json',output=STUDY/'evaluation/PILOT256_VAL_GATE.json')
print(json.dumps(dict(phase='VAL_COMPLETE_PARENT_REVIEW_REQUIRED',outputs=[str(p)for p in outputs],log=str(log),no_test_truth=True)),flush=True)
