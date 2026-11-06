"""Root-authorized whole v04 CPU-dense Mac rerun, reuses original public bytes."""
import os,sys,time,subprocess,resource,json
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'delegated/vector/v04_mac_primary';OUT.mkdir(parents=True,exist_ok=False)
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
    if (ROOT/'RESULTS/vector'/('v04_'+arm)).exists():raise RuntimeError('Existing local primary must not be overwritten')
queue=dict(status='ROOT_AUTHORIZED_WHOLE_V04_MAC_CPU_DENSE_PRIMARY_MIGRATION',authorized_task_count=3,
    tasks=[dict(task_id='v04_'+a,public_id='v04_matched12',arm=a) for a in ('sourceCIE','GlobalDual','AplusGlobal')],
    Y_policy='Original Windows-generated public NPZ bytes reused; no new noise',GTscore='HOLD_ALL18')
(ROOT/'configs/vector_v04_mac_primary_queue.json').write_text(json.dumps(queue,indent=2)+'\n')
command=[sys.executable,str(ROOT/'experiments/vector_run.py'),'invert','--queue','configs/vector_v04_mac_primary_queue.json','--max-tasks','3','--backend','cpu-dense','--output','RESULTS/vector/case_v04_mac_primary_process']
with (OUT/'RUN.log').open('w') as log:
    p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
    (OUT/'ACTIVE.json').write_text(json.dumps(dict(pid=p.pid,command=command,started_unix=time.time()),indent=2)+'\n');print(json.dumps(dict(event='START_V04_MAC_PRIMARY',pid=p.pid)),flush=True)
    pid,status,usage=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status)
row=dict(status='COMPLETE' if p.returncode==0 else 'FAILED',exit=p.returncode,platform='Mac',CPU_seconds=usage.ru_utime+usage.ru_stime,wall_seconds=time.perf_counter()-T,
    peak_RSS_bytes=usage.ru_maxrss,additive=True,scope='Whole wait4 v04 scientific child; internal setup/arm timers nested',controller_CPU_seconds=time.process_time()-C)
(OUT/'PROCESS_COST.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
raise SystemExit(p.returncode)
