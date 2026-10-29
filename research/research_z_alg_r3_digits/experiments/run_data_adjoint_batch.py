"""Nine ordinary data-adjoint control jobs, two workers, public inputs only."""
import json,os,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1]
def main():
 for marker in ['audit/EXPERIMENT_PREDICTIONS_FREEZE_HOLDOUT.json','audit/EXPERIMENT_PREDICTIONS_FREEZE_BUDGET_FILLED.json','audit/DATA_ADJOINT_CONTROL_SOURCE_SNAPSHOT/FREEZE.json']:
  if not (ROOT/marker).exists():raise SystemExit('Original, budget-fill and control-source freezes required')
 cases=sorted((ROOT/'public_data/holdout').glob('*.npz'));rows=[];start=time.perf_counter()
 if len(cases)!=9:raise SystemExit('Nine public observations required')
 env=os.environ.copy()
 for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='1'
 def run(case):
  out=ROOT/'results/data_adjoint/data_adjoint'/case.name.replace('.npz','.json');out.parent.mkdir(parents=True,exist_ok=True)
  if out.exists():raise FileExistsError(out)
  command=[sys.executable,str(ROOT/'data_adjoint_control.py'),'--public',str(case),'--operators',str(ROOT/'public_data/operators.npz'),'--output',str(out)]
  begun=time.perf_counter();r=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True);out.with_suffix('.log').write_text(r.stdout+'\nSTDERR\n'+r.stderr)
  row=dict(observation_id=case.stem,method='data_adjoint',process_exit=r.returncode,process_wall_seconds=time.perf_counter()-begun,output=str(out.relative_to(ROOT)))
  if out.exists():
   result=json.loads(out.read_text());row.update(status=result.get('status'),relative_residual=(result.get('top1') or {}).get('relative_residual'),wall_seconds=result.get('wall_seconds'),cpu_seconds=result.get('cpu_seconds'),accounting=result.get('accounting'),all_unknowns_released=result.get('all_unknowns_released'))
  else:row['status']='PROCESS_FAILED_NO_OUTPUT'
  print(json.dumps({k:v for k,v in row.items() if k!='accounting'}),flush=True);return row
 with ThreadPoolExecutor(max_workers=2) as executor:
  for future in as_completed([executor.submit(run,c) for c in cases]):rows.append(future.result())
 receipt=dict(split='data_adjoint',base_split='holdout',workers=2,BLAS_threads_per_worker=1,rows=rows,total_batch_wall_seconds=time.perf_counter()-start,outputs=len(rows),truth_access='NONE_PUBLIC_INPUT_ONLY',all_processes_succeeded=all(r['process_exit']==0 for r in rows))
 (ROOT/'audit/EXPERIMENT_BATCH_DATA_ADJOINT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(outputs=len(rows),wall_seconds=receipt['total_batch_wall_seconds'],success=receipt['all_processes_succeeded'])),flush=True)
if __name__=='__main__':main()
