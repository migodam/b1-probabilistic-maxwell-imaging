"""Serial guarded queue; no truth or scores are opened here."""
from pathlib import Path
import argparse,subprocess,sys,time,json,datetime,os,resource
ROOT=Path(__file__).resolve().parents[1]
METHODS=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--noise',type=int,choices=[0,1],default=0);a=ap.parse_args();clock=json.loads((ROOT/'RUN_CLOCK.json').read_text());cut=datetime.datetime.fromisoformat(clock['no_new_heavy_after_utc'].replace('Z','+00:00')).timestamp();log=ROOT/'results'/f'CORE_QUEUE_r{a.noise}.jsonl'
    if not (ROOT/'METHODS_FROZEN.json').exists():raise RuntimeError('Definitive methods must be frozen before blindqueue')
    cases=[f'{base}_r{a.noise}' for base in ['blind_normal','blind_cavity','blind_strong']]
    jobs=([(case,method) for method in METHODS for case in cases] if a.noise==1 else [(case,method) for case in cases for method in METHODS])
    time_stopped=False
    for job_index,(case,method) in enumerate(jobs):
        inp=ROOT/'inputs/online'/case
        out=ROOT/'results/online'/case/method;out.mkdir(parents=True,exist_ok=True)
        if (out/'ONLINE_SEALED.json').exists():continue
        if time.time()>cut:
            time_stopped=True
            for pending_case,pending_method in jobs[job_index:]:
                pending=ROOT/'results/online'/pending_case/pending_method
                if (pending/'ONLINE_SEALED.json').exists():continue
                pending.mkdir(parents=True,exist_ok=True);(pending/'NOT_RUN_TIME.json').write_text(json.dumps(dict(status='NOT_RUN_TIME',reason='Original no-new-heavy cutoff; solver budgets unchanged',time_unix=time.time()),indent=2)+'\n')
            break
        command=[sys.executable,str(ROOT/'tools/guarded_online.py'),'--audit-log',str(out/'FILE_ACCESS.jsonl'),str(ROOT/'tools/run_r3.py'),'run','--measurement',str(inp/'observation.npz'),'--config',str(inp/'config.json'),'--output',str(out),'--method',method]
        started=time.time();row=dict(event='job_start',case=case,method=method,started_unix=started,command=command)
        with log.open('a') as f:f.write(json.dumps(row)+'\n')
        print('START',case,method,flush=True)
        before_usage=resource.getrusage(resource.RUSAGE_CHILDREN);queue_cpu=time.process_time()
        with (out/'OUTPUT.log').open('w') as stdout:result=subprocess.run(command,cwd=ROOT,stdout=stdout,stderr=subprocess.STDOUT,env={**os.environ,'R3_SOURCE_ROOT':str(ROOT/'frozen/src'),'PYTHONDONTWRITEBYTECODE':'1'})
        after_usage=resource.getrusage(resource.RUSAGE_CHILDREN);child_cpu=after_usage.ru_utime+after_usage.ru_stime-before_usage.ru_utime-before_usage.ru_stime;resource_receipt=dict(child_process_cpu_seconds=child_cpu,queue_process_cpu_seconds=time.process_time()-queue_cpu,authority='RUSAGE_CHILDREN before/after single guarded subprocess; full process CPU incl imports',children_max_rss_cumulative=after_usage.ru_maxrss,children_max_rss_scope='cumulative max across exitedchildren; not per-job RSS; use child receipt peak')
        (out/'QUEUE_RESOURCE_RECEIPT.json').write_text(json.dumps(resource_receipt,indent=2)+'\n')
        row=dict(event='job_finish',case=case,method=method,returncode=result.returncode,wall_seconds=time.time()-started,sealed=(out/'ONLINE_SEALED.json').exists(),child_process_cpu_seconds=child_cpu)
        with log.open('a') as f:f.write(json.dumps(row)+'\n')
        print('FINISH',case,method,row,flush=True)
        if result.returncode:print('FAILED retained',out,flush=True)
    (ROOT/'results'/f'CORE_QUEUE_r{a.noise}_FINISHED.json').write_text(json.dumps(dict(time_unix=time.time(),noise=a.noise,stopped_no_new_heavy=time_stopped,methods=METHODS,cases=['blind_normal','blind_cavity','blind_strong']),indent=2))
if __name__=='__main__':main()
