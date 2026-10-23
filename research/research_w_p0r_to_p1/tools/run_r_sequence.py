"""Resumable serial R phase; no historical run is restarted."""
import sys,subprocess,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];py=sys.executable
jobs=[]
for seed in [2026102201,2026102202]:
    jobs.append((f'oracle-fit-{seed}-wide',['-m','wr.oracle_fit',str(seed),'wide','TR','14']))
for seed in [2026102201,2026102202]:
    for band in ['narrow','wide']:
        jobs.append((f'online-fit-{seed}-{band}',['-m','wr.online_inference',str(seed),band,'TR','14']))
jobs.append(('online-lbfgsb-2026102201-wide',['-m','wr.online_inference','2026102201','wide','LBFGSB','18']))
for label,args in jobs:
    receipt=root/'results/jobs'/label/'RECEIPT.json'
    if receipt.exists():
        print(json.dumps({'skip_existing_job':label,'status':json.loads(receipt.read_text())['status']}),flush=True);continue
    print(json.dumps({'start_job':label}),flush=True)
    ret=subprocess.run([py,str(root/'tools/run_job.py'),'R',label,py,*args],cwd=root)
    if ret.returncode:
        print(json.dumps({'stopped_after_failure':label,'code':ret.returncode}),flush=True);sys.exit(ret.returncode)
print(json.dumps({'sequence':'COMPLETE'}),flush=True)
