"""Wait for the finite own TEST queues, fetch once, then seal without GT scoring."""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess,sys,time
s=Path(__file__).resolve().parents[2]
end=datetime(2026,11,7,10,44,27,tzinfo=timezone.utc).timestamp()
start=time.perf_counter();last=0.
while time.time()<end:
    d=json.loads((s/'audit/REMOTE_TEST_PARTITION_PROCESS.json').read_text())
    if d['status']=='REMOTE_TEST_COMMAND_RETURNED':break
    if d['status'] in ('REMOTE_TEST_FAILED_NO_AUTOMATIC_RETRY','NOT_LAUNCHED_AFTER_RECORDED_FAILURE'):
        raise SystemExit('Remote command failed; inspect retained logs before any retry')
    if time.perf_counter()-last>=30:
        print(json.dumps(dict(status=d['status'],waiting_elapsed_seconds=time.perf_counter()-start,scientific_methods='FROZEN',TEST_truth_scoring='NOT_STARTED')),flush=True)
        last=time.perf_counter()
    time.sleep(5)
else:raise SystemExit('Hard deadline reached; no further automatic task')
for cmd in (
    [sys.executable,'-u','src/remote_jobs.py','fetch-evaluation'],
    [sys.executable,'experiments/score_predictions.py','seal','--split','test','--runs',*json.loads((s/'audit/FINAL_RUNS.json').read_text())['runs'],'--seal','evaluation/TEST_SEALED.json'],
):
    subprocess.run(cmd,cwd=s,check=True)
print(json.dumps(dict(status='TEST_FETCHED_AND_SEALED_NO_GT_SCORING',wall_seconds=time.perf_counter()-start,seal='evaluation/TEST_SEALED.json')),flush=True)
