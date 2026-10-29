"""Meter a serial local job once, including subprocess CPU, failures and monitoring."""
from pathlib import Path
import sys,os,time,json,fcntl,resource,subprocess
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from r2.accounting import BudgetBook
def main():
    phase,name,*cmd=sys.argv[1:]
    if not cmd or Path(name).name!=name:raise ValueError('unique safe label required')
    b=BudgetBook();out=ROOT/'results/jobs'/name;out.mkdir(exist_ok=False)
    with (ROOT/'results/job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        left=min(b.cap['process_cpu_seconds']-b.snapshot()['process_cpu_seconds'],b.auth['caps'][phase]['process_cpu_seconds']-b.snapshot(phase)['process_cpu_seconds'])
        if left<1:raise RuntimeError('CPU cap exhausted')
        env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']})
        old=ROOT.parent/'research_w_gq'
        env.update(PYTHONPATH=os.pathsep.join(map(str,[ROOT/'src',ROOT/'reuse/previous_gq/src',ROOT/'reuse/previous_gq/tools',ROOT/'reuse/fixed_vendor',old/'src',old/'vendor',ROOT.parents[2]/'research-w-continuation/.deps'])),MPLCONFIGDIR=str(ROOT/'results/.mpl'),PYTHONDONTWRITEBYTECODE='1')
        def limit():resource.setrlimit(resource.RLIMIT_CPU,(int(left),int(left)))
        t=time.perf_counter();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
        b.event(phase,'job_start',job=name,command=cmd)
        with (out/'OUTPUT.log').open('w') as f:
            p=subprocess.Popen(cmd,env=env,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limit)
            killed=False;peak=0
            while p.poll() is None:
                try:
                    rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(p.pid)],text=True).strip() or 0)*1024;peak=max(peak,rss)
                    if rss>b.auth['rss_bytes']:killed=True;p.kill()
                except (ValueError,subprocess.CalledProcessError):pass
                time.sleep(.25)
        a=resource.getrusage(resource.RUSAGE_CHILDREN)
        result=dict(job=name,phase=phase,command=cmd,exit_code=p.returncode,process_cpu_seconds=time.process_time()-cpu+a.ru_utime+a.ru_stime-before.ru_utime-before.ru_stime,wall_seconds=time.perf_counter()-t,polled_peak_RSS_bytes=peak,memory_killed=killed,authority='wrapper own plus children once; physical reservations include failed requests')
        b.event(phase,'job_completion',**{k:v for k,v in result.items() if k!='phase'})
        (out/'RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
        if p.returncode:b.failure(phase,name,'job failed',receipt=str(out/'RECEIPT.json'))
        print(json.dumps(result));return p.returncode
if __name__=='__main__':raise SystemExit(main())
