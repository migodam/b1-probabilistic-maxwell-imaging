"""All-27, two-worker, pre-GT budget-fill launcher; originals stay unchanged."""
import json,os,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
def main():
 marker=ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_HOLDOUT.json'
 if not marker.exists() or not (ROOT/'audit/BUDGET_FILL_SOURCE_SNAPSHOT/FREEZE.json').exists():raise SystemExit('Original prediction and supplement source freezes required')
 record=json.loads(marker.read_text());jobs=sorted(record['files'],key=lambda r:(r['observation_id'],r['method']));rows=[];start=time.perf_counter()
 if len(jobs)!=27:raise SystemExit('All27 originals required')
 env=os.environ.copy()
 for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
 def run(job):
  ident=job['observation_id'];method=job['method'];out=ROOT/'results/budget_filled'/method/f'{ident}.json';out.parent.mkdir(parents=True,exist_ok=True)
  if out.exists():raise FileExistsError(out)
  args=[sys.executable,str(ROOT/'budget_fill.py'),'--original',str(ROOT/job['snapshot']),'--public',str(ROOT/'public_data/holdout'/f'{ident}.npz'),'--operators',str(ROOT/'public_data/operators.npz'),'--output',str(out)]
  begun=time.perf_counter();r=subprocess.run(args,cwd=ROOT,env=env,capture_output=True,text=True);out.with_suffix('.log').write_text(r.stdout+'\nSTDERR\n'+r.stderr)
  row=dict(observation_id=ident,method=method,process_exit=r.returncode,process_wall_seconds=time.perf_counter()-begun,output=str(out.relative_to(ROOT)))
  if out.exists():
   result=json.loads(out.read_text());row.update(status=result.get('status'),relative_residual=(result.get('top1') or {}).get('relative_residual'),incremental_wall_seconds=result.get('incremental_wall_seconds'),incremental_cpu_seconds=result.get('incremental_cpu_seconds'),incremental_accounting=result.get('incremental_accounting'))
  else:row['status']='PROCESS_FAILED_NO_OUTPUT'
  print(json.dumps({k:v for k,v in row.items() if k!='incremental_accounting'}),flush=True);return row
 with ThreadPoolExecutor(max_workers=2) as executor:
  for future in as_completed([executor.submit(run,j) for j in jobs]):rows.append(future.result())
 receipt=dict(split='budget_filled',base_split='holdout',workers=2,BLAS_threads_per_worker=1,rows=rows,total_batch_wall_seconds=time.perf_counter()-start,outputs=len(rows),truth_access='NONE_PUBLIC_INPUT_ONLY',primary_originals_unchanged=True,all_processes_succeeded=all(r['process_exit']==0 for r in rows),accounting='Add only incremental_accounting/incremental_wall_seconds/incremental_cpu_seconds to original costs; filled root totals are cumulative')
 (ROOT/'audit/EXPERIMENT_BATCH_BUDGET_FILLED.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(outputs=len(rows),wall_seconds=receipt['total_batch_wall_seconds'],success=receipt['all_processes_succeeded'])),flush=True)
if __name__=='__main__':main()
