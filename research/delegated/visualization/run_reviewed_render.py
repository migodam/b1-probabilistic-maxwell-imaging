#!/usr/bin/env python3
"""Measure the renderer process without touching the shared scientific cost ledger."""
import json, os, resource, shutil, subprocess, sys, time
from pathlib import Path
root=Path(__file__).resolve().parents[3]
out=root/'research/delegated/visualization'
if (out/'receipt.json').exists() and not (out/'receipt_before_visual_review.json').exists():
    shutil.copy2(out/'receipt.json',out/'receipt_before_visual_review.json')
env=os.environ.copy()
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    env[k]='1'
exe=root.parents[1]/'.venv_nn/bin/python'
start=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (out/'stdout.log').open('w') as so,(out/'stderr.log').open('w') as se:
    proc=subprocess.run([str(exe),str(out/'make_figures.py')],cwd=root,env=env,stdout=so,stderr=se)
after=resource.getrusage(resource.RUSAGE_CHILDREN)
run={'run_index':1,'wall_seconds':time.perf_counter()-start,
     'cpu_user_seconds':after.ru_utime-before.ru_utime,'cpu_system_seconds':after.ru_stime-before.ru_stime,
     'exit_code':proc.returncode,'physics_calls':0,'training_updates':0}
run['cpu_seconds']=run['cpu_user_seconds']+run['cpu_system_seconds']
history=out/'render_runs.json'
runs=json.loads(history.read_text()) if history.exists() else []
run['run_index']=len(runs)+1;runs.append(run);history.write_text(json.dumps(runs,indent=2)+'\n')
if proc.returncode==0:
    receipt=json.loads((out/'receipt.json').read_text())
    receipt['measured_process']=run
    receipt['visual_review_render_cumulative']={'runs':len(runs),'wall_seconds':sum(r['wall_seconds'] for r in runs),
                                               'cpu_seconds':sum(r['cpu_seconds'] for r in runs),
                                               'physics_calls':0,'cost_ledger_touched':False}
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(run))
if proc.returncode:print((out/'stderr.log').read_text()[-4500:])
sys.exit(proc.returncode)
