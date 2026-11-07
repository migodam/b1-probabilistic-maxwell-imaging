"""Finite public VAL depth diagnostic after the existing CPU/MPS lanes close."""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess,sys,time,os
s=Path(__file__).resolve().parents[2]
receipt=s/'audit/RECURRENT_DEPTH_PROCESS.json'
end=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()
state=dict(pid=os.getpid(),status='WAIT_CORE_MPS_AND_VAL_CPU_COMPLETION',truth_access='NONE',tasks=[])
receipt.write_text(json.dumps(state,indent=2)+'\n')
while time.time()<end:
    training=json.loads((s/'audit/CORE_TRAINING_TAIL.json').read_text())['status']
    validation=json.loads((s/'audit/VAL_EVALUATION_CONTROLLER.json').read_text())['status']
    if training=='COMPLETED_CORE_CONFIRMATION_AND_E2' and validation=='CLOSED_ALL_REQUESTED_RUNS_COMPLETED':
        break
    time.sleep(5)
else:
    state['status']='NOT_RUN_ABSOLUTE_DEADLINE';receipt.write_text(json.dumps(state,indent=2)+'\n');raise SystemExit(0)
for mode in ('refreshed_q','feedback1','feedback_hinge'):
    path=s/f'models/recurrent/recurrent_{mode}_seed7_cascade2/best.pt'
    state.update(status='RUNNING_PUBLIC_VAL_DEPTH_DIAGNOSTIC',active=mode)
    receipt.write_text(json.dumps(state,indent=2)+'\n')
    tick=time.perf_counter()
    r=subprocess.run([sys.executable,'-u',str(s/'experiments/recurrent_iterations.py'),'--checkpoint',str(path),'--device','mps'],cwd=s)
    state['tasks'].append(dict(mode=mode,returncode=r.returncode,wall_seconds=time.perf_counter()-tick))
    if r.returncode:
        state['status']='FAILED_NO_AUTOMATIC_RETRY';receipt.write_text(json.dumps(state,indent=2)+'\n');raise SystemExit(r.returncode)
state.update(status='COMPLETED_THREE_VAL_DEPTH_DIAGNOSTICS',active=None)
receipt.write_text(json.dumps(state,indent=2)+'\n')
