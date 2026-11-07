"""Root-authorized own next-arm mkdir guard; never touches current sourceDual."""
import time,json,datetime
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1]
process=ROOT/'delegated/vector/gpu_recovery_attempt1_process';raw=ROOT/'delegated/vector/gpu_recovery_attempt1_inversions'
current=raw/'c17_matched12_sourceDual_CONF';nextarm=raw/'c17_matched12_globalDual_CONF'
if (process/'PROCESS_COST.json').exists() or (current/'SEALED.json').exists():raise RuntimeError('Boundary already passed; refuse late engineering sentinel')
rows=[json.loads(s) for s in (process/'RHS_PROGRESS.jsonl').read_text().splitlines()]
observed=dict(completed_primal_RHS=sum(not x['adjoint'] for x in rows),completed_adjoint_RHS=sum(x['adjoint'] for x in rows),
 completed_primal_batch_labels=len({x['label'] for x in rows if not x['adjoint']}),scope='Only persisted completed RHS rows; not live model counters; in-flight work unknown until finally')
nextarm.mkdir(exist_ok=False)
(nextarm/'CONTROLLED_DRAIN_SENTINEL.json').write_text(json.dumps(dict(status='ROOT_AUTHORIZED_CONTROLLED_ENGINEERING_DRAIN_NOT_ALGORITHM_FAILURE',
 UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),current_task='c17_matched12_sourceDual_CONF',current_directory_untouched=True,
 next_task='c17_matched12_globalDual_CONF',next_task_new_FW=0,next_task_initial_evaluate_prevented_by_existing_exclusive_mkdir=True,
 before_observed_counters=observed,final_counts_source='gpu_recovery_attempt1_process/PROCESS_COST.json after currentarm natural seal and controlled FileExistsError',
 action_CPU_seconds=time.process_time(),action_wall_seconds=time.perf_counter()-ENTRY,GT_access=False),indent=2)+'\n')
print(json.dumps(dict(status='SENTINEL_PLACED_BEFORE_NEXT_TASK_INITIAL_FW',next_path=str(nextarm),before_observed_counters=observed)),flush=True)
