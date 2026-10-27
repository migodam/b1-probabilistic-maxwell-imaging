"""Serial BGGA job wrapper. Prep uses no PDE; physics requires an actual grant."""
from pathlib import Path
import fcntl,json,os,resource,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from bgga.accounting import BudgetBook

def main():
    phase,label,*cmd=sys.argv[1:]
    if Path(label).name!=label or not cmd:raise ValueError('safe unique job label and command required')
    book=BudgetBook();cap=None
    if phase!='PREFLIGHT':cap=book.authorization(phase)
    out=ROOT/'results'/'jobs'/label;out.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    # Reuse already installed read-only local dependencies; publication users
    # install requirements.txt. These are libraries, never new source evidence.
    paths=[ROOT/'src',ROOT/'vendor']
    existing=ROOT.parents[2]/'research-w-continuation'/'.deps'
    if existing.is_dir():paths.append(existing)
    env.update(PYTHONPATH=os.pathsep.join(map(str,paths)),BGGA_STAGE=phase,
        PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(ROOT/'results'/'.mpl'))
    with (ROOT/'results'/'job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        t=time.monotonic();self0=time.process_time();r0=resource.getrusage(resource.RUSAGE_CHILDREN)
        book.event(phase,'job_start',job=label,command=cmd)
        def limits():
            if cap:
                remain=max(1,int(cap['cpu_seconds']-book.snapshot(phase)['process_cpu_seconds']))
                resource.setrlimit(resource.RLIMIT_CPU,(remain,remain))
        with (out/'OUTPUT.log').open('w') as log:
            p=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
        r1=resource.getrusage(resource.RUSAGE_CHILDREN)
        rec=dict(job=label,phase=phase,command=cmd,exit_code=p.returncode,
            process_cpu_seconds=time.process_time()-self0+r1.ru_utime+r1.ru_stime-r0.ru_utime-r0.ru_stime,
            wall_seconds=time.monotonic()-t,children_peak_RSS=r1.ru_maxrss,
            RSS_unit='bytes on macOS; high-water, not additive',authority='wrapper own CPU + final child CPU, once',
            new_pde_authorized=bool(cap))
        book.event(phase,'job_completion',**{k:v for k,v in rec.items() if k!='phase'})
        (out/'RECEIPT.json').write_text(json.dumps(rec,indent=2)+'\n')
        if p.returncode:
            with (ROOT/'FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
        print(json.dumps(rec));return p.returncode

if __name__=='__main__':raise SystemExit(main())
