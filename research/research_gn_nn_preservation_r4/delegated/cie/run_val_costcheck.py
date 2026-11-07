"""One finite, public VAL timing pair after M4 preparation; no GT scoring."""
from pathlib import Path
import json,subprocess,sys,time
study=Path(__file__).resolve().parents[2]
request=study/'audit/VAL_EVALUATION_REQUEST.json'
initial=json.loads(request.read_text());planned=initial['runs'][:]
initial.update(runs=[],temporary_parent_pause='Two public VAL CIE timing cases; restore after this finite job')
request.write_text(json.dumps(initial,indent=2)+'\n')
receipt={'status':'WAIT_M4_PREPARATION','truth_access':'NONE','started_unix':time.time()}
dest=study/'audit/CIE_VAL_TIMING_PROCESS.json'
dest.write_text(json.dumps(receipt,indent=2)+'\n')
try:
    for _ in range(90):
        d=json.loads((study/'audit/LOCAL_PREPARATION_PROGRESS.json').read_text())
        if d['ready']==d['expected'] and (d.get('local_process_returncode')==0 or d['status'].startswith('LOCAL_')):break
        time.sleep(2)
    else:raise TimeoutError('M4 timing window not ready; retain main experiment priority')
    records=json.loads((study/'data/manifests/online.json').read_text())['records']
    ids=[r['obs_id']for r in records if r['split']=='val'][:2]
    receipt.update(status='RUNNING_PUBLIC_VAL_PAIR',obs_ids=ids)
    dest.write_text(json.dumps(receipt,indent=2)+'\n')
    result=subprocess.run([sys.executable,'-u',str(study/'experiments/run_cie_reference.py'),
        '--split','val','--obs-id',*ids,'--methods','cie_dual','--workers','2','--execute'],cwd=study)
    receipt.update(status='COMPLETED'if result.returncode==0 else'FAILED',returncode=result.returncode)
finally:
    current=json.loads(request.read_text());current['runs']=list(dict.fromkeys(planned+current.get('runs',[])))
    current.pop('temporary_parent_pause',None);request.write_text(json.dumps(current,indent=2)+'\n')
    receipt['ended_unix']=time.time();dest.write_text(json.dumps(receipt,indent=2)+'\n')
