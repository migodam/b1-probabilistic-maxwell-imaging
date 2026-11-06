"""Finite metadata coordinator; waits for original24, then16/5/29seal-score.
No numerical work overlaps the existing two-lane inversion queue.
"""
import time
WALL,CPU=time.perf_counter(),time.process_time()
from pathlib import Path
import os,sys,json,subprocess,argparse,datetime
ROOT=Path(__file__).resolve().parents[2];REPO=ROOT.parents[1];LOGS=ROOT/'delegated/3d_extension';OUT=ROOT/'RESULTS/vector_extension'

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,value):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(p)
def ps_cpu(pid):
    s=subprocess.run(['ps','-p',str(pid),'-o','state=','-o','cputime='],capture_output=True,text=True).stdout.strip()
    if not s:return None,False
    state,s=s.split();parts=s.split(':');return sum(float(v)*60**i for i,v in enumerate(reversed(parts))),'Z' in state

def main():
    p=argparse.ArgumentParser();p.add_argument('--original-parent-pid',required=True,type=int);args=p.parse_args()
    receipt=dict(status='WAITING_ORIGINAL24',started_utc=utc(),original_parent_pid=args.original_parent_pid,
        scientific_change='Root-authorized solegrid16validation and five resolvedhollow controls; original24 unchanged; GT held until29 sealed',jobs=[])
    path=LOGS/'AMENDMENT_COORDINATOR.json';last_cpu=None;last_sample=None
    while True:
        cpu,zombie=ps_cpu(args.original_parent_pid)
        if cpu is None:break
        last_cpu=cpu;last_sample=time.perf_counter();receipt.update(original_parent_observed_CPU_seconds=cpu,wait_wall_seconds=time.perf_counter()-WALL)
        write(path,receipt)
        if zombie:break
        time.sleep(.5)
    queue=read(LOGS/'phase2_queue/RECEIPT.json')
    if queue['status']!='COMPLETE' or len(queue['jobs'])!=24:raise RuntimeError('Original24 queue incomplete; no expansion/score launch')
    holdpath=OUT/'PUBLIC_PREDICTIONS_SEAL.json';hold=read(holdpath)
    if hold.get('status')!='SCORE_HOLD_ROOT_AMENDMENT_29' or hold.get('is_prediction_seal') is not False:raise RuntimeError('Expected explicit pre-GT HOLD, never replace true seal')
    if (OUT/'SCORES.json').exists():raise RuntimeError('Unexpected prior scoring')
    holdpath.rename(LOGS/'SCORE_HOLD_ROOT_AMENDMENT_29.json')
    receipt.update(status='ORIGINAL24_COMPLETE_GT_NOT_SCORED',original24_exit_scope='automatic24 scoring intentionally rejected at existing seal pre-GT guard; no worker interruption',
        original_parent_CPU_observed_lower_bound_seconds=last_cpu,
        original_parent_CPU_unobserved_upper_bound_extra_seconds=(time.perf_counter()-last_sample+.02) if last_sample else None,
        coordinator_cost_note='Original parent was already running before amendment. CPU sampled with ps until exit; retain stated bound instead of claiming an exact postmortem CPU. Validation/generation individual exact CPU also retained; do not add nested costs twice.')
    write(path,receipt)
    for name,script,arguments in [('mie16','mie_refinement.py',[]),('resolved_hollow5','resolved_hollow_control.py',['all']),('finalize29','finalize29.py',[])]:
        command=[sys.executable,str(LOGS/script),*arguments];log=LOGS/(name+'.log');started=time.perf_counter()
        row=dict(stage=name,command=command,started_utc=utc(),log=str(log.relative_to(REPO)));receipt.update(status='RUNNING_'+name);write(path,receipt)
        with log.open('x') as stream:child=subprocess.Popen(command,cwd=REPO,stdout=stream,stderr=subprocess.STDOUT);code=child.wait()
        row.update(exit_code=code,finished_utc=utc(),wall_seconds=time.perf_counter()-started);receipt['jobs'].append(row);write(path,receipt)
        print(json.dumps(row),flush=True)
        if code and name!='mie16':
            receipt.update(status='FAILED_REQUIRES_ROOT_REVIEW',error_stage=name);write(path,receipt);return code
    children=[read(OUT/t['case_id']/t['arm']/'RECEIPT.json') for t in read(ROOT/'configs/vector_extension_all29.json')['tasks']]
    extra_parent=read(LOGS/'resolved_hollow_queue/PARENT_COST.json');mie=read(OUT/'validation_refinement16/RECEIPT.json');final=read(OUT/'FINALIZATION_COST.json')
    known=sum(r['phase2_incremental_CPU_seconds'] for r in children)+extra_parent['parent_CPU_seconds']+mie['CPU_seconds']+final['CPU_seconds']+time.process_time()-CPU+(last_cpu or 0)
    receipt.update(status='COMPLETED',finished_utc=utc(),coordinator_CPU_seconds=time.process_time()-CPU,whole_amendment_wall_seconds=time.perf_counter()-WALL,
        worker_CPU_seconds=sum(r['phase2_incremental_CPU_seconds'] for r in children),known_total_incremental_CPU_lower_bound_seconds=known,
        total_incremental_CPU_upper_bound_seconds=known+(receipt['original_parent_CPU_unobserved_upper_bound_extra_seconds'] or 0),
        total_scope='Originalparent lowerbound +29wholeworkers +newcontrolparent +Mie16 +finalization +metadata coordinator. Originalparent includes originalvalidation/generation; these nestedindividual receipts not added again.',
        arms=29,max_numeric_processes=2,BLAS_threads=1,scoring='ALL29_PUBLIC_SEAL_BEFORE_PRIVATE_SCORE')
    write(path,receipt);write(OUT/'PHASE2_COMPLETE.json',receipt);return 0

if __name__=='__main__':raise SystemExit(main())
