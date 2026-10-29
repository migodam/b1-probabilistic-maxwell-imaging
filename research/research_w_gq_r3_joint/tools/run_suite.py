"""Run exactly one declared task with bounded RSS, one BLAS thread and durable cost receipt.

Physical JointModel methods own .fem.lock, so this wrapper does not recursively
acquire it. The parent schedules generation and inversion serially.
"""
from pathlib import Path
import argparse,os,sys,json,time,subprocess,resource,datetime,fcntl
ROOT=Path(__file__).resolve().parents[1]
PYTHON=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/A25/a25-adaptive-om-feedback/.venv_a25/bin/python")
def main():
    p=argparse.ArgumentParser();p.add_argument('label');p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    if Path(a.label).name!=a.label or not a.command:raise ValueError('safe unique label and command required')
    clock=json.loads((ROOT/'RUN_CLOCK.json').read_text());deadline=datetime.datetime.fromisoformat(clock['deadline_utc'].replace('Z','+00:00')).timestamp()
    dest=ROOT/'results/jobs'/a.label;dest.mkdir(parents=True,exist_ok=False)
    env=os.environ.copy();env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']})
    r2=ROOT.parent/'research_w_gq_r2';deps=ROOT.parents[2]/'research-w-continuation/.deps'
    source=ROOT/'frozen/src' if (ROOT/'frozen/src').exists() else ROOT/'src'
    env['R3_SOURCE_ROOT']=str(source)
    env['PYTHONPATH']=os.pathsep.join(map(str,[source,ROOT/'tools',r2/'src',r2/'reuse/previous_gq/src',r2/'reuse/previous_gq/tools',r2/'reuse/fixed_vendor',deps]))
    env.update(MPLCONFIGDIR=str(ROOT/'results/.mpl'),PYTHONDONTWRITEBYTECODE='1')
    before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.perf_counter();peak=0;killed=None
    command=[str(PYTHON) if x=='PYTHON' else x for x in a.command]
    with (dest/'OUTPUT.log').open('w') as f:
        child=subprocess.Popen(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        while child.poll() is None:
            try:
                rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(child.pid)],text=True).strip() or 0)*1024;peak=max(peak,rss)
                if rss>clock['single_job_rss_bytes']:killed='RSS_LIMIT';child.kill()
                elif time.time()>deadline:killed='CAMPAIGN_DEADLINE';child.kill()
            except (ValueError,subprocess.CalledProcessError):pass
            time.sleep(.25)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt=dict(label=a.label,pid=child.pid,command=command,exit_code=child.returncode,wall_seconds=time.perf_counter()-t,process_cpu_seconds=(after.ru_utime+after.ru_stime)-(before.ru_utime+before.ru_stime),polled_peak_RSS_bytes=peak,killed=killed,finished_unix=time.time(),cost_authority='child process plus its reaped children, once')
    (dest/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    ledger=ROOT/'JOB_LEDGER.jsonl'
    with (ROOT/'.job-ledger.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with ledger.open('a') as f:f.write(json.dumps(receipt)+'\n')
        fcntl.flock(lock,fcntl.LOCK_UN)
    print(json.dumps(receipt));return child.returncode
if __name__=='__main__':raise SystemExit(main())
