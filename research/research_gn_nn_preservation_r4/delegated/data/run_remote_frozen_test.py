"""One finite remote TEST launch after this study's CIE and public transfer."""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess,sys,time,os
s=Path(__file__).resolve().parents[2]
out=s/'audit/REMOTE_TEST_PARTITION_PROCESS.json'
end=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()
state=dict(pid=os.getpid(),status='WAIT_OWN_CIE_AND_PUBLIC_PREDICTION_TRANSFER',truth_access='NONE')
out.write_text(json.dumps(state,indent=2)+'\n')
while time.time()<end:
    cie=sorted((s/'audit').glob('REMOTE_TRANSFER_queue_cie_*.json'))
    predictions=sorted((s/'audit').glob('REMOTE_TRANSFER_prediction_upload_*.json'))
    if cie and predictions:
        c=json.loads(cie[-1].read_text());p=json.loads(predictions[-1].read_text())
        if c['status']=='REMOTE_COMMAND_RETURNED' and p['status']=='PUBLIC_REMOTE_PREDICTIONS_UPLOADED':break
        if c['status']=='FAILED_OR_INTERRUPTED' or p['status'].startswith(('FAILED','HOLD')):
            state.update(status='NOT_LAUNCHED_AFTER_RECORDED_FAILURE',cie=c,prediction_transfer=p)
            out.write_text(json.dumps(state,indent=2)+'\n');raise SystemExit(1)
    time.sleep(5)
else:
    state['status']='NOT_RUN_ABSOLUTE_DEADLINE';out.write_text(json.dumps(state,indent=2)+'\n');raise SystemExit(0)
state.update(status='REMOTE_TEST_RUNNING',started_unix=time.time())
out.write_text(json.dumps(state,indent=2)+'\n')
tick=time.perf_counter();r=subprocess.run([sys.executable,'-u',str(s/'src/remote_jobs.py'),'run','--kind','test'],cwd=s)
state.update(status='REMOTE_TEST_COMMAND_RETURNED'if r.returncode==0 else'REMOTE_TEST_FAILED_NO_AUTOMATIC_RETRY',returncode=r.returncode,elapsed_seconds=time.perf_counter()-tick)
out.write_text(json.dumps(state,indent=2)+'\n')
