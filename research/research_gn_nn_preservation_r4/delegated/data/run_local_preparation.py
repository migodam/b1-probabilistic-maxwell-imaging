"""Finite 2-CPU strict public preparation with periodic canonical pair upload."""
from pathlib import Path
import sys,json,time,subprocess,collections,importlib.util
STUDY=Path(__file__).resolve().parents[2];sys.path.insert(0,str(STUDY.parents[1]))
from research.research_gn_nn_preservation_r4.src.remote_jobs import upload_local_pairs
all_ids=json.loads((STUDY/'delegated/data/local_reserved_ids.json').read_text());uploaded=set();started=time.time();history=[]
def progress(stage,process):
    canonical=collections.Counter();completed=[];nonfull=[];solve=0.
    for obs in all_ids:
        p=STUDY/'anchors'/f'{obs}.json'
        if not p.exists() or not p.with_suffix('.npz').exists():continue
        r=json.loads(p.read_text())
        if r.get('execution_origin',{}).get('machine')!='arm64':continue
        completed.append(obs);canonical.update(r['counts']);solve+=r['solve_wall_seconds']
        if r['final_stage']!=24 or not r['all_unknowns_released']:nonfull.append(obs)
    out=dict(status=stage,local_process_pid=process.pid if process else None,local_process_returncode=process.poll() if process else None,ready=len(completed),expected=len(all_ids),uploaded=len(uploaded),actual_local_counts=dict(canonical),actual_local_solve_seconds=solve,controller_wall_seconds=time.time()-started,nonfull=nonfull,truth_access='NONE',test_truth_status='SEALED_NOT_SCORED')
    p=STUDY/'audit/LOCAL_PREPARATION_PROGRESS.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(out,indent=2)+'\n');tmp.replace(p)
    print(json.dumps(out),flush=True)
def upload():
    try:
        report=upload_local_pairs(all_ids,uploaded);uploaded.update(report['imported']+report['existing']);history.append(dict(unix=time.time(),**report))
        (STUDY/'audit/LOCAL_PREPARATION_UPLOADS.json').write_text(json.dumps(history,indent=2)+'\n')
        if report['held']:print(json.dumps(dict(HELD_CANONICAL_PARTIAL=report['held'])),flush=True)
    except Exception as exc:print(json.dumps(dict(upload_error=repr(exc))),flush=True)
if __name__=='__main__':
    for group in ('val','train'):
        command=[sys.executable,str(STUDY/'experiments/run_anchors.py'),'--id-file',str(STUDY/f'delegated/data/local_{group}_ids.json'),'--workers','2']
        with (STUDY/f'delegated/data/local_{group}_run.log').open('w') as log:
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                upload();progress('RUNNING_'+group.upper(),process);time.sleep(30)
            upload();progress('COMPLETED_'+group.upper(),process)
            if process.returncode:raise SystemExit(process.returncode)
    # Give a transient transfer error finite retries after physics finishes.
    for _ in range(5):
        upload();progress('LOCAL_PHYSICS_COMPLETE',None)
        if len(uploaded)==len(all_ids):break
        time.sleep(30)
    progress('LOCAL_COMPLETE_UPLOADED' if len(uploaded)==len(all_ids) else 'LOCAL_COMPLETE_UPLOAD_INCOMPLETE',None)
