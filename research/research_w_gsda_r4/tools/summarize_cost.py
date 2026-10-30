"""Aggregate authoritative R4 reservations and one wrapper receipt per job/PID."""
from pathlib import Path
import json,sys,collections,time
ROOT=Path(__file__).resolve().parents[1]
def main():
 ledger=ROOT/'PHYSICS_LEDGER.jsonl';rows=[json.loads(l) for l in ledger.read_text().splitlines()] if ledger.exists() else []
 receipts=[json.loads(p.read_text()) for p in sorted((ROOT/'receipts').glob('*.json'))]
 unique={}
 for r in receipts:unique[(r.get('job'),r.get('pid'))]=r
 categories=collections.defaultdict(lambda:dict(rhs=0,full_fe=0,records=0))
 dense=dict(LU_attempts=0,LU_completed=0,condition_SVDs=0,dense_cpu_seconds=0.,dense_wall_seconds=0.)
 grouped=collections.defaultdict(lambda:dict(full_fe_reserved_attempts=0,complex_rhs_reserved_attempts=0,events=0,failures=0))
 for r in rows:
  g=grouped[r['job']]
  label=str(r.get('label','')).lower()
  if r['kind']=='reservation':
   cat='GD_Green_actions' if ' gd' in label or ' gdh' in label else ('weighted_adjoint' if 'weighted adjoint' in label else ('receiver_adjoint_GS_or_fullJ' if 'receiver adjoint' in label or 'receiver-adjoint' in label else ('full_factorization' if r.get('factorizations',0) else ('full_primal_or_other_RHS'))))
   categories[cat]['rhs']+=r.get('rhs',0);categories[cat]['full_fe']+=r.get('factorizations',0);categories[cat]['records']+=1
  if r.get('event_kind')=='MODAL_SCHUR_PRECONDITIONER':
   dense['LU_attempts']+=r.get('dense_LU_attempts',r.get('dense_LU_factorizations',0));dense['LU_completed']+=r.get('dense_LU_factorizations',0);dense['condition_SVDs']+=r.get('dense_condition_SVDs',0);dense['dense_cpu_seconds']+=r.get('dense_process_cpu_seconds',0.);dense['dense_wall_seconds']+=r.get('dense_wall_seconds',0.)
  if r['kind']=='reservation':g['full_fe_reserved_attempts']+=r.get('factorizations',0);g['complex_rhs_reserved_attempts']+=r.get('rhs',0)
  elif r['kind']=='failure':g['failures']+=1
  else:g['events']+=1
 for (job,pid),r in unique.items():
  g=grouped[job];g.setdefault('wrapper_receipts',[]).append(r)
 total=dict(full_fe_reserved_attempts=sum(g['full_fe_reserved_attempts'] for g in grouped.values()),complex_rhs_reserved_attempts=sum(g['complex_rhs_reserved_attempts'] for g in grouped.values()),known_completed_wrapper_cpu_seconds=sum(r.get('wrapper_process_cpu_seconds',0) for r in unique.values()),known_completed_wrapper_wall_seconds_serial_jobs=sum(r.get('wall_seconds',0) for r in unique.values()),peak_known_job_rss_bytes=max([r.get('peak_rss_bytes',0) for r in unique.values()] or [0]))
 out=dict(updated_unix=time.time(),scope='new R4 campaign only; historical costs not free deployment',totals=total,by_job=dict(grouped),rhs_categories=dict(categories),reduced_dense_core=dense,receipt_count=len(unique),reservations_counted_once=True,kernel_times_not_added_to_wrappers=True,active_jobs_without_final_receipt=sorted({r['job'] for r in rows}-{r.get('job') for r in unique.values()}),not_measured='interactive coding/reasoning and some read-only file commands; wrapper CPU covers actual scientific job setup/actions/logging/I/O')
 aggregates=[json.loads(p.read_text()) for p in (ROOT/'results').rglob('BATCH_RECEIPT.json') if json.loads(p.read_text()).get('aggregate_only')]
 missing_receipts=out['active_jobs_without_final_receipt']
 out['jobs_without_final_CPU_receipt']=missing_receipts
 out['completed_zero_physics_failures_without_CPU_receipt']=[j for j in missing_receipts if grouped[j]['full_fe_reserved_attempts']==0 and grouped[j]['complex_rhs_reserved_attempts']==0 and grouped[j]['failures']>0]
 out['active_jobs_without_final_receipt']=[j for j in missing_receipts if j not in out['completed_zero_physics_failures_without_CPU_receipt']]
 out['reduced_dense_core']['scope']='Logged modal preconditioner LU/condition-SVD only; not a full profiler of every NumPy/SciPy dense operation. Other small GN/Gram/projected diagnostic solves are included in CPU but lack a unified invocation count.'
 out['authoritative_batch_wrapper_aggregates_not_added_twice']=aggregates
 out['wrapper_cpu_unpartitioned_batch_seconds']=sum(a.get('unpartitioned_process_cpu_seconds',0) for a in aggregates)
 out['totals']['completed_wrapper_cpu_plus_unpartitioned_batch_seconds']=out['totals']['known_completed_wrapper_cpu_seconds']+out['wrapper_cpu_unpartitioned_batch_seconds']
 out['timing_scope_audit']=dict(legacy_receipt_scope='BudgetBook construction or batch-partition start; Python bootstrap/import CPU omitted; original values preserved',strict_whole_campaign_CPU='NOT_FULLY_MEASURED',external_whole_child_records='END_TO_END_CPU.json / EXTERNAL_PROCESS_COST.json where actually measured; do not add entire child to overlapping Book receipts',unmeasured_preimport_failures='six launcher interpreter failures; actual wall logged, CPU not measured',scientific_conclusion='No same-quality or strict end-to-end CPU speedup claim')
 (ROOT/'PHYSICS_COST.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(total))
if __name__=='__main__':main()
