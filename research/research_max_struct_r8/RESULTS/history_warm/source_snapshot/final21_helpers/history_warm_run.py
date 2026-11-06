"""One active process, oneBLASthread, fixed21 via unchanged R8run_case."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,gc
from types import SimpleNamespace
from research.research_max_struct_r8.experiments.run_case import execute,ROOT,REPO,dump

def main():
    wall,cpu=time.perf_counter(),time.process_time();path=ROOT/'experiments/history_warm_manifest.json';manifest=json.loads(path.read_text());rows=[]
    out=ROOT/'RESULTS/history_warm';out.mkdir(parents=True,exist_ok=True)
    for task in manifest['tasks']:
        dest=REPO/task['output_dir'];receipt=dest/'RECEIPT.json'
        if receipt.exists():
            r=json.loads(receipt.read_text())
            if r['status']=='SEALED_PUBLIC_PREDICTION':rows.append(dict(task_id=task['task_id'],status=r['status'],reused_existing=True,cpu_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],physical_model_cost=r['physical_model_cost']));continue
            raise RuntimeError('Preservefailedattemptdirectorybeforeexplicitrepair: '+str(dest))
        args=SimpleNamespace(manifest=str(path),public_list=None,task_id=task['task_id'],config=None,output_dir=None)
        r,d=execute(args);rows.append(dict(task_id=task['task_id'],status=r['status'],cpu_seconds=r['cpu_seconds'],wall_seconds=r['wall_seconds'],physical_model_cost=r['physical_model_cost']))
        dump(out/'RUN_PROGRESS.json',dict(status='RUNNING',completed=len(rows),planned=len(manifest['tasks']),records=rows,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
        print(json.dumps(rows[-1]),flush=True);gc.collect()
    dump(out/'ALL_PREDICTIONS_SEALED.json',dict(status=f'ALL{len(rows)}_HISTORICAL_PREDICTIONS_SEALED' if all(x['status']=='SEALED_PUBLIC_PREDICTION' for x in rows) else 'COMPLETED_WITH_FAILURES',records=rows,tasks=len(rows),objects=4,cpu_seconds_this_invocation=time.process_time()-cpu,wall_seconds_this_invocation=time.perf_counter()-wall,truth_access='NONE',cold_imports='Shared once withinoneprocess; per-taskrun_case importscachedafterfirst; reusedreceiptcosts notchargedagain; final21 mayresumeafteroriginal20batch'))

if __name__=='__main__':main()
