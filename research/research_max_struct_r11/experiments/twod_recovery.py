"""Metadata-only active-attempt/path and immutable-queue resolution; no physics/GT."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/twod/independent48'
def registry():
 p=ROOT/'configs/twod_independent48_recovery_registry.json';return json.loads(p.read_text()) if p.exists() else {}
def active_task(task):
 task=dict(task);info=registry().get('path_overrides',{}).get(task['task_id'])
 if info:task.update(output_dir=info['active_output'],engineering_attempt_id=info['attempt_id'],engineering_original_output=info['original_output'])
 return task
def active_job(job):
 info=registry().get('path_overrides',{}).get(job['task_id'])
 return not info or str(job.get('output_dir','')).replace('\\','/')==info['active_output']
def queue_receipt(path):
 path=Path(str(path).replace('\\','/'))
 recovered=path.parent/'RECOVERED_TERMINAL_PREFIX.json'
 if path.parent.name=='expanded_queue' and recovered.exists():return recovered
 checkpoints=path.parent/'CHECKPOINTS'
 if checkpoints.exists():
  files=sorted(checkpoints.glob('*.json'))
  if files:return files[-1]
 return path
def queue_paths():
 return [queue_receipt(D/name/'RECEIPT.json') for name in ('queue','interleaved_queue','expanded_queue','recovery_queue','final_retry03_queue') if (D/name/'RECEIPT.json').exists() or (D/name/'CHECKPOINTS').exists()]
def affected_objects():return registry().get('affected_objects',[])
def attempt_cost_info(task_id,valid_cpu):
 reg=registry();info=reg.get('path_overrides',{}).get(task_id)
 if not info:return dict(engineering_retry=False,exact_end_to_end_CPU=valid_cpu,end_to_end_CPU_lower_bound=valid_cpu,end_to_end_CPU_upper_bound=valid_cpu,exact_end_to_end_cost_qualified=True)
 original=info['original_attempt'];startup=sum(j['child_process_CPU_seconds'] for j in reg.get('zero_physics_retry02_failures',[]) if j['task_id']==task_id);lower=original.get('CPU_lower_bound_seconds');upper=original.get('CPU_upper_observed_wall_bound_seconds')
 return dict(engineering_retry=True,active_attempt_id=info['attempt_id'],exact_end_to_end_CPU=None,end_to_end_CPU_lower_bound=valid_cpu+startup+(lower if lower is not None else 0) if valid_cpu is not None else None,end_to_end_CPU_upper_bound=valid_cpu+startup+upper if valid_cpu is not None and upper is not None else None,exact_end_to_end_cost_qualified=False,aborted_CPU_lower_bound_recorded=lower,zero_physics_startup_CPU=startup,original_attempt_id=original['attempt_id'],scope='Scientific valid480 prefix known; original engineering-aborted fullCPU unknown. Lower bound omits unmeasured nonnegative CPU, never asserts it was0.')

def deployment_cost_info(prefix,fallback,trigger):
 members=[prefix]+([fallback] if trigger else [])
 lo=[x.get('end_to_end_CPU_lower_bound') for x in members];hi=[x.get('end_to_end_CPU_upper_bound') for x in members];qualified=all(x.get('exact_end_to_end_cost_qualified',True) for x in members)
 return dict(scientific_prefix_deployment_child_CPU_seconds=sum(x['full_child_CPU_seconds'] for x in members) if all(x['full_child_CPU_seconds'] is not None for x in members) else None,whole_deployment_end_to_end_CPU_seconds=sum(x['full_child_CPU_seconds'] for x in members) if qualified and all(x['full_child_CPU_seconds'] is not None for x in members) else None,whole_deployment_CPU_lower_bound=sum(lo) if all(x is not None for x in lo) else None,whole_deployment_CPU_upper_observed_bound=sum(hi) if all(x is not None for x in hi) else None,exact_end_to_end_cost_qualified=qualified,engineering_affected_object=any(x.get('engineering_retry') for x in members))
