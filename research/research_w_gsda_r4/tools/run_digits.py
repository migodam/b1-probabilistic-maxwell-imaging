"""Serial, label-free batch driver; every listed public case gets an attempt.

This launches one process at a time. The outer orchestrator owns shared queue
slots and campaign clock. Each case has identical frozen starts/options/caps.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from gsda.accounting import BudgetBook,dump


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--online-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--options',type=Path,required=True)
    ap.add_argument('--cases',nargs='+',required=True,help='anonymous public case directory names')
    ap.add_argument('--python',default=sys.executable)
    ap.add_argument('--job-prefix',default='fwi_batch')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if len(set(a.cases))!=len(a.cases):raise ValueError('duplicate cases')
    if not a.options.exists():raise ValueError('frozen method options required for batch')
    records=[]
    for case in a.cases:
        online=a.online_root/case;out=a.output/case
        record=dict(case=case,started_unix=time.time(),status='NOT_RUN')
        try:
            book=BudgetBook(root=ROOT,job=f'{a.job_prefix}_{case}');book.check(new_job=True)
            command=[a.python,str(ROOT/'tools/run_fwi.py'),'--config',str(online/'config.json'),
                '--measurement',str(online/'observation.npz'),'--output',str(out),
                '--options',str(a.options),'--job',f'{a.job_prefix}_{case}']
            out.mkdir(parents=True,exist_ok=True)
            with (out/'stdout.log').open('w') as log:
                process=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
            record.update(returncode=process.returncode,status='COMPLETED' if process.returncode==0 else 'FAILED')
            if (out/'ONLINE_SEALED.json').exists():record['online_status']=json.loads((out/'ONLINE_SEALED.json').read_text())['status']
        except Exception as error:
            record['reason']=str(error)
        record['finished_unix']=time.time();records.append(record)
        dump(a.output/'BATCH_STATUS.json',dict(cases=a.cases,attempts=records,
            remaining=[c for c in a.cases if c not in {r['case'] for r in records}],
            same_options=str(a.options),same_budget=True))
    print(json.dumps(dict(output=str(a.output),attempted=len(records),cases=a.cases)),flush=True)

if __name__=='__main__':main()
