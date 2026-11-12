"""Compose terminal physics attempts into a NEW evidence root; no NPZ reads/solves.

Copies bytes only. Original archives, queues and results are never modified.
A parent invokes this only after the dependency-only recovery has terminated.
"""
import argparse,copy,json,shutil,time
from pathlib import Path

TERMINAL={'COMPLETED','COMPLETED_WITH_ERRORS','COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY'}
ATTEMPT_KEYS={'campaign','parent_attempt','parent_attempt_metadata','recovery_of','recovery_source','parent_campaign','technical_recovery'}


def load(path):return json.loads(path.read_text())

def normalize(value):
    if isinstance(value,dict):return {k:normalize(v) for k,v in value.items() if k not in ATTEMPT_KEYS}
    if isinstance(value,list):return [normalize(v) for v in value]
    return value.replace('\\','/') if isinstance(value,str) else value

def indexed(rows):
    out={}
    for row in rows:
        tid=row['task_id']
        if tid in out:raise ValueError('Duplicate task_id: '+tid)
        if not tid or Path(tid).name!=tid or '/' in tid or '\\' in tid:raise ValueError('Unsafe task_id')
        out[tid]=row
    return out

def terminal(receipt):
    if receipt.get('status') not in TERMINAL or receipt.get('pending',0) or receipt.get('running',[]):raise ValueError('Explicit terminal idle queue receipt required')

def successful(job):return job.get('returncode')==0 and job.get('child_status')=='COMPLETE'

def cpu(job):return (job.get('external') or {}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))

def known_sum(values):
    values=list(values)
    return sum(values) if all(v is not None for v in values) else None

def validate(original,recovery,first,second):
    terminal(first);terminal(second)
    if original['campaign']!='physics_v1' or recovery['campaign']!='physics_recovery_v1':raise ValueError('Unexpected campaign')
    if normalize(original['protocol'])!=normalize(recovery['protocol']):raise ValueError('Numerical protocol changed')
    tasks=indexed(original['tasks']);rt=indexed(recovery['tasks']);a=indexed(first['jobs']);b=indexed(second['jobs'])
    if set(a)!=set(tasks) or len(tasks)!=60:raise ValueError('Original must have all 60 expected attempts')
    failures={tid for tid in tasks if not successful(a[tid])};successes=set(tasks)-failures
    if len(failures)!=42 or len(successes)!=18 or set(rt)!=failures or set(b)!=failures:raise ValueError('Recovery must cover exactly original 42 failures, preserving 18 successes')
    if any(t.get('family') for t in tasks.values()):
        if {tid for tid,t in tasks.items() if t.get('family')=='classic'}!=successes:raise ValueError('Successful original cohort must be the 18 classic Austria tasks')
    for tid in tasks:
        for key in ['case_id','method']:
            if a[tid].get(key)!=tasks[tid][key]:raise ValueError('Original job identity mismatch: '+tid)
    for tid in failures:
        if normalize(tasks[tid])!=normalize(rt[tid]):raise ValueError('Recovery task/input/budget/prior mismatch: '+tid)
        if not successful(b[tid]):raise ValueError('Recovery technical failure remains: '+tid)
        if any(b[tid].get(k)!=tasks[tid][k] for k in ['case_id','method']):raise ValueError('Recovery job identity mismatch: '+tid)
    return tasks,a,b,failures

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--original-manifest',type=Path,required=True);p.add_argument('--recovery-manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--archive-name',default='physics_v1_recovered_full_results_v2.zip');a=p.parse_args()
    root=a.root.resolve();out=a.out.resolve()
    if out.exists() or out==root or out.is_relative_to(root/'RESULTS') or root.is_relative_to(out):raise ValueError('Require a NEW isolated evidence root, never an existing study/RESULTS directory')
    if Path(a.archive_name).name!=a.archive_name or a.archive_name=='physics_v1_full_results.zip' or not a.archive_name.endswith('.zip'):raise ValueError('Require a versioned archive name distinct from the immutable initial archive')
    original=load(a.original_manifest);recovery=load(a.recovery_manifest);q=root/'queues';first=load(q/'physics_v1/RECEIPT.json');second=load(q/'physics_recovery_v1/RECEIPT.json');tasks,initial,retried,failures=validate(original,recovery,first,second)
    # Validate metadata before creating any output; no arrays or scientific values read.
    for tid,t in tasks.items():
        campaign='physics_recovery_v1' if tid in failures else 'physics_v1';folder=root/'RESULTS'/campaign/tid
        r=load(folder/'RECEIPT.json')
        if r.get('status')!='COMPLETE' or normalize(r['task'])!=normalize(t) or normalize(r['protocol'])!=normalize(original['protocol']):raise ValueError('Result receipt task/protocol mismatch: '+tid)
    started,wall=time.process_time(),time.perf_counter();out.mkdir(parents=True)
    def write(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2)+'\n')
    final=[]
    for tid in tasks:
        recovered=tid in failures;campaign='physics_recovery_v1' if recovered else 'physics_v1'
        shutil.copytree(root/'RESULTS'/campaign/tid,out/'RESULTS/physics_v1'/tid)
        job=copy.deepcopy(retried[tid] if recovered else initial[tid]);job['attempt_source_campaign']=campaign;final.append(job)
        if recovered and (root/'RESULTS/physics_v1'/tid).exists():shutil.copytree(root/'RESULTS/physics_v1'/tid,out/'audit/physics_failed_initial'/tid)
    shutil.copytree(q/'physics_recovery_v1',out/'queues/physics_recovery_v1')
    # Keep recovery receipts available for accounting references, not a second primary cohort.
    shutil.copytree(root/'RESULTS/physics_recovery_v1',out/'RESULTS/physics_recovery_v1')
    dst=out/'queues/physics_v1';dst.mkdir(parents=True)
    (dst/'attempts').mkdir();shutil.copy2(q/'physics_v1/RECEIPT.json',dst/'attempts/initial_RECEIPT.json')
    failed=[dict(copy.deepcopy(initial[tid]),attempt_source_campaign='physics_v1_initial',result_receipt_relative=f'audit/physics_failed_initial/{tid}/RECEIPT.json') for tid in tasks if tid in failures]
    merged=dict(status='COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY',campaign='physics_v1',pending=0,running=[],jobs=final,failed_attempts=failed,finished_child_CPU_seconds=known_sum(cpu(j) for j in final),failed_attempt_child_CPU_seconds=known_sum(cpu(j) for j in failed),controller_CPU_seconds=known_sum([first.get('controller_CPU_seconds'),second.get('controller_CPU_seconds')]),wall_seconds=known_sum([first.get('wall_seconds'),second.get('wall_seconds')]),wall_scope='Sum of original and recovery active-controller windows; excludes intervention gap; never summed child wall',recovery_source='../physics_recovery_v1/RECEIPT.json',original_source='attempts/initial_RECEIPT.json',old_successful_runs_repeated=0,accounting='60 primary jobs plus 42 separate failed initial attempts; recovery queue is a nested reference. Controller windows counted once. No scientific quality gate.',initial_archive_immutable=True,archive_filename=a.archive_name)
    write(dst/'RECEIPT.json',merged);write(out/'configs/physics_v1.json',original);write(out/'configs/physics_recovery_v1.json',recovery)
    write(out/'audit/PHYSICS_RECOVERY_COMPOSITION.json',dict(status='COMPOSED_NEW_EVIDENCE_ROOT_NOT_SCORED_OR_ARCHIVED',source_root=str(root),output_root=str(out),planned_archive_name=a.archive_name,primary_jobs=60,initial_successes=18,recovered_jobs=42,failed_initial_attempts=42,truth_access='NONE',new_physics=0,CPU_seconds=time.process_time()-started,wall_seconds=time.perf_counter()-wall,seal_policy='New output root requires fresh seals and parent authorization; copied receipts preserve attempt campaign provenance. Input arrays remain in original source root unless separately restored.'))
    print(json.dumps(dict(status=merged['status'],out=str(out),archive_created=False,primary_jobs=60,failed_attempts=42)))

if __name__=='__main__':main()
