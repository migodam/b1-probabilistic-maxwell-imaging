"""One root-explicit postseal scorer,with external whole-child Macwait4 cost."""
import time,sys,json,argparse,subprocess,os
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];V=ROOT/'delegated/vector'
p=argparse.ArgumentParser();p.add_argument('--authorization',required=True);p.add_argument('--authorization-note',required=True);a=p.parse_args()
if a.authorization!='ROOT_EXPLICIT_FRESH50_POSTSEAL':raise RuntimeError('FreshGT requires later explicitrootauthorization')
b=json.loads((V/'confirmation_host_union/ALL50_PUBLIC_PREDICTION_AND_PREFIX_SEALS.json').read_text());s=json.loads((V/'confirmation_host_union/SELECTIONS.json').read_text())
if b['global_tasks']!=50 or s['unique_task_ids']!=50:raise RuntimeError('Full50publicdeliverybarrier required')
OUT=V/'fresh_postseal_external';OUT.mkdir(exist_ok=False)
command=[sys.executable,str(ROOT/'experiments/vector_postseal_score.py'),'--cohort','confirmation','--authorization',a.authorization,'--authorization-note',a.authorization_note]
with (OUT/'SCORER.log').open('x') as handle:
 child_start=time.perf_counter();child=subprocess.Popen(command,cwd=ROOT.parents[1],stdout=handle,stderr=subprocess.STDOUT)
 _,status,usage=os.wait4(child.pid,0);child.returncode=os.waitstatus_to_exitcode(status)
record=dict(status='COMPLETE' if child.returncode==0 else 'FAILED_POSTSEAL_METADATA_ONLY',authorization_note=a.authorization_note,command=command,exit=child.returncode,whole_scorer_child_CPU_seconds=usage.ru_utime+usage.ru_stime,whole_scorer_child_wall_seconds=time.perf_counter()-child_start,controller_wall_seconds=time.perf_counter()-ENTRY,controller_own_CPU_seconds=time.process_time(),scope='One externalwait4 scorer includesallmetrics/plots/exports andallchild CPU;replacesnested SCORING_COST;controller ownCPU separate',new_physical_RHS=0,all50publicseals_checked_before_GT=True)
(OUT/'PROCESS_COST.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True)
raise SystemExit(child.returncode)
