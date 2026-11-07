"""Create an explicit recovery manifest/overlay without altering science or old outputs."""
import json,time
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];D=ROOT/'delegated/twod/independent48'
def read(p):return json.loads(p.read_text())
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
master=read(ROOT/'configs/twod_independent48_all_tasks.json');inventory=read(D/'control_abort_20261107/INVENTORY.json');q=read(D/'expanded_queue/RECEIPT.json');done={j['task_id'] for n in ('queue','interleaved_queue','expanded_queue') for j in read(D/n/'RECEIPT.json')['jobs']};aborted={r['task_id']:r for r in inventory['inflight_attempts']}
if len(done)!=132 or len(aborted)!=4 or any(r['completed_valid'] for r in aborted.values()):raise ValueError('Recovery inventory changed; preserve and audit')
registry=dict(status='ENGINEERING_RECOVERY_OVERLAY_SCIENCE_UNCHANGED',version=3,created_utc=datetime.now(timezone.utc).isoformat(),affected_objects=sorted({r['case_id'] for r in aborted.values()}),interrupted_attempts=list(aborted.values()),path_overrides={},extra_attempt_cost_exact=False,primary_all48_retained=True,paired_sensitivity_excludes_affected_objects_only=True,GT_access='NONE')
tasks=[]
for original in master['tasks']:
 if original['task_id'] in done:continue
 t=dict(original)
 if t['task_id'] in aborted:
  original_output=t['output_dir'];new='research/research_max_struct_r11/RESULTS/independent48_engineering_retry02/'+t['task_id']
  if (REPO/new).exists():raise FileExistsError(new)
  registry['path_overrides'][t['task_id']]=dict(original_output=original_output,active_output=new,attempt_id=t['task_id']+'__engineering_retry02',original_attempt=aborted[t['task_id']])
  t.update(output_dir=new,engineering_attempt_id=t['task_id']+'__engineering_retry02',engineering_original_output=original_output)
 elif (REPO/t['output_dir']).exists():raise FileExistsError('Unstarted output unexpectedly exists: '+t['task_id'])
 tasks.append(t)
# Preserve cyclic within-case order, not the original baseline-first ordering.
methods=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT');ordered=[]
for case in sorted({t['case_id'] for t in tasks}):
 k=int(case.rsplit('_',1)[1])%6
 for method in methods[k:]+methods[:k]:
  ordered += [t for t in tasks if t['case_id']==case and t.get('material_readout',{}).get('label',t['arm'])==method]
if len(ordered)!=156:raise ValueError('Expected4 incomplete plus152 unstarted')
for p in (ROOT/'configs/twod_independent48_recovery_registry.json',ROOT/'configs/twod_independent48_recovery_tasks.json',D/'expanded_queue/RECOVERED_TERMINAL_PREFIX.json'):
 if p.exists():raise FileExistsError(p)
q.update(status='INTERRUPTED_CONTROLLER_RETAINED_TERMINAL_PREFIX',running=[],pending=0,not_started_task_ids=[t['task_id'] for t in ordered],interrupted_attempts=list(aborted),original_receipt_preserved='RECEIPT.json',original_tmp_preserved='RECEIPT.json.tmp',external_interrupted_attempt_CPU='UNKNOWN_NOT_ZERO',controller_CPU_scope='Last known controller selfCPU lowerbound; final tail unknown')
dump(D/'expanded_queue/RECOVERED_TERMINAL_PREFIX.json',q);dump(ROOT/'configs/twod_independent48_recovery_registry.json',registry);dump(ROOT/'configs/twod_independent48_recovery_tasks.json',dict(status='ENGINEERING_V3_UNCHANGED156_SCIENTIFIC_TASKS',tasks=ordered));dump(D/'control_abort_20261107/RECOVERY_PLAN.json',dict(retained132=len(done),retry4=list(aborted),unstarted152=152,total_scientific_tasks=288,remaining156=len(ordered),cpus=[0,2,8,10],new_materials=0,GT_access='NONE',registry=registry,CPU_seconds=time.process_time()))
print(json.dumps(dict(status='RECOVERY_READY',retained=132,engineering_retries=4,unstarted=152,remaining=156)))
