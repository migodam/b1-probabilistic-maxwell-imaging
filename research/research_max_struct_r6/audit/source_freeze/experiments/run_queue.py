"""Finite four-process CPU queue; scientific failures are results, exceptions stop dispatch."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
import argparse,json,subprocess,sys,time
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--name',default='primary');ap.add_argument('--workers',type=int,default=4)
    args=ap.parse_args();manifest=REPO/args.manifest;data=json.loads(manifest.read_text())
    out=ROOT/'audit'/('queue_'+args.name);out.mkdir(exist_ok=False)
    start=time.perf_counter();result=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),workers=args.workers,
        manifest=args.manifest,truth_access='NONE',tasks=[],failed_attempts=[],batch_wall_seconds=None)
    def write():
        result['elapsed_wall_seconds']=time.perf_counter()-start
        (out/'RECEIPT.json').write_text(json.dumps(result,indent=2))
    def task(row):
        t=time.perf_counter()
        log=out/(row['task_id']+'.log')
        with log.open('w') as stream:
            p=subprocess.run([sys.executable,str(ROOT/'experiments/run_case.py'),'--manifest',str(manifest),'--task-id',row['task_id']],cwd=REPO,stdout=stream,stderr=subprocess.STDOUT)
        receipt_path=REPO/row['output_dir']/'RECEIPT.json'
        rec=json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        return dict(task_id=row['task_id'],returncode=p.returncode,process_wall_seconds=time.perf_counter()-t,
            output_dir=row['output_dir'],status=rec.get('status','NO_RECEIPT'),error=rec.get('error'),
            physical_objective=rec.get('physical_objective'),physical_relative_residual=rec.get('physical_relative_residual'))
    write();tasks=iter(data['tasks']);stopped=False
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={}
        for _ in range(args.workers):
            row=next(tasks,None)
            if row is not None:pending[pool.submit(task,row)]=row
        while pending:
            finished,_=wait(pending,return_when=FIRST_COMPLETED)
            for fut in finished:
                row=pending.pop(fut);entry=fut.result();result['tasks'].append(entry)
                if entry['returncode']:
                    result['failed_attempts'].append(entry);stopped=True
                print(json.dumps(dict(completed=len(result['tasks']),total=len(data['tasks']),**entry)),flush=True)
                if not stopped:
                    row=next(tasks,None)
                    if row is not None:pending[pool.submit(task,row)]=row
                write()
    result.update(status='EXECUTION_EXCEPTION_REQUIRES_REPAIR' if stopped else 'COMPLETED_UNSCORED',
        batch_wall_seconds=time.perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat())
    write();print(json.dumps(dict(status=result['status'],batch_wall_seconds=result['batch_wall_seconds'])),flush=True)

if __name__=='__main__':main()
