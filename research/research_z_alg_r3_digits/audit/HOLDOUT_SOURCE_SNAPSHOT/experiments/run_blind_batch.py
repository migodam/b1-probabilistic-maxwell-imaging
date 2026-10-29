"""Two-worker public-input-only solver batch with retained process failures."""
import argparse,json,os,subprocess,sys,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
def main():
 ap=argparse.ArgumentParser();ap.add_argument('split',choices=['dev','dev_v2','holdout']);ap.add_argument('--config',required=True);ap.add_argument('--parent-go',action='store_true');a=ap.parse_args()
 if a.split=='holdout' and (not a.parent_go or not (ROOT/'audit/HOLDOUT_FREEZE.json').exists()):raise SystemExit('Holdout run requires source/config freeze and parent GO')
 input_split='dev' if a.split.startswith('dev') else a.split
 cases=sorted((ROOT/'public_data'/input_split).glob('*.npz'));expected=2 if a.split.startswith('dev') else 9
 if len(cases)!=expected:raise SystemExit(f'Expected {expected} public cases')
 cfgpath=ROOT/a.config;cfg=json.loads(cfgpath.read_text());env=os.environ.copy()
 for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
 jobs=[(case,method) for case in cases for method in METHODS];started=time.perf_counter();rows=[]
 def run(job):
  case,method=job;out=ROOT/'results'/a.split/method/(case.stem+'.json');out.parent.mkdir(parents=True,exist_ok=True);log=out.with_suffix('.log')
  if out.exists():raise RuntimeError(f'Will not overwrite {out}')
  command=[sys.executable,str(ROOT/'run_solver.py'),'--public',str(case),'--operators',str(ROOT/'public_data/operators.npz'),'--method',method,'--config',str(cfgpath),'--output',str(out)]
  start=time.perf_counter();completed=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True);log.write_text(completed.stdout+'\nSTDERR\n'+completed.stderr)
  row=dict(observation_id=case.stem,method=method,process_exit=completed.returncode,process_wall_seconds=time.perf_counter()-start,output=str(out.relative_to(ROOT)),log=str(log.relative_to(ROOT)))
  if out.exists():
   r=json.loads(out.read_text());row.update(status=r.get('status'),wall_seconds=r.get('wall_seconds'),cpu_seconds=r.get('cpu_seconds'),all_unknowns_released=r.get('all_unknowns_released'),relative_residual=(r.get('top1') or {}).get('relative_residual'),accounting=r.get('accounting'))
  else:row['status']='PROCESS_FAILED_NO_OUTPUT'
  print(json.dumps({k:v for k,v in row.items() if k!='accounting'}),flush=True);return row
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(run,job) for job in jobs]
  for future in as_completed(futures):rows.append(future.result())
 receipt=dict(split=a.split,workers=2,BLAS_threads_per_worker=1,configuration=cfg,configuration_path=str(cfgpath.relative_to(ROOT)),cases=[c.stem for c in cases],methods=METHODS,rows=rows,total_batch_wall_seconds=time.perf_counter()-started,truth_access='NONE_PUBLIC_INPUT_ONLY',all_processes_succeeded=all(r['process_exit']==0 for r in rows))
 p=ROOT/f'audit/EXPERIMENT_BATCH_{a.split.upper()}.json';p.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(receipt=str(p.relative_to(ROOT)),completed=len(rows),wall_seconds=receipt['total_batch_wall_seconds'],success=receipt['all_processes_succeeded'])),flush=True)
if __name__=='__main__':main()
