"""Serial one-thread scientific job and process-CPU accounting."""
import sys,os,time,json,resource,subprocess,math,shutil
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from wr.accounting import Book,write
stage,label=sys.argv[1:3]; command=sys.argv[3:]; book=Book(stage)
if '/' in command[0] and not os.path.isabs(command[0]):command[0]=str(Path(command[0]).absolute())
dest=root/'results/jobs'/label; dest.mkdir(parents=True,exist_ok=False)
shutil.copytree(root/'src',dest/'SOURCE_SNAPSHOT',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
shutil.copytree(root/'configs',dest/'CONFIG_SNAPSHOT')
dep='/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-continuation/.deps'
hist=root.parent/'research_w_continuation'
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=os.pathsep.join([str(root/'src'),str(hist/'src'),str(hist/'vendor'),dep]),MPLCONFIGDIR=str(root/'results/mpl'))
left=book.config['stages'][stage]['cpu_seconds']-book.totals()[stage]['cpu_seconds']
def limit(): resource.setrlimit(resource.RLIMIT_CPU,(max(1,math.ceil(left)),max(2,math.ceil(left)+1)))
start=time.perf_counter(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
past=[json.loads(x) for x in book.path.read_text().splitlines()] if book.path.exists() else []
wall_left=book.config['stages'][stage]['cpu_seconds']-sum(x.get('wall_seconds',0.) for x in past if x.get('stage')==stage and x.get('event')=='job_cost')
reason=None;peak_sampled=0
with (dest/'OUTPUT.log').open('w') as stream:
    try:
        proc=subprocess.Popen(command,cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=limit)
        while proc.poll() is None:
            if time.perf_counter()-start>wall_left:reason='STAGE_WALL_LIMIT'
            try:
                sample=int(subprocess.check_output(['ps','-o','rss=','-p',str(proc.pid)],text=True).strip())*1024
                peak_sampled=max(peak_sampled,sample)
                if sample>book.config['max_process_RSS_bytes']:reason='RSS_LIMIT'
            except (ValueError,subprocess.CalledProcessError):pass
            if reason:
                proc.terminate()
                try:proc.wait(timeout=10)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()
                break
            time.sleep(.5)
        code=proc.wait()
    except Exception as exc:
        reason=f'LAUNCH_FAILURE: {type(exc).__name__}: {exc}';code=-999
after=resource.getrusage(resource.RUSAGE_CHILDREN)
cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
row=dict(event='job_cost',stage=stage,label=label,command=command,status='COMPLETE' if code==0 else 'FAILED',returncode=code,wall_seconds=time.perf_counter()-start,charged_cpu_seconds=cpu,peak_rss_bytes=after.ru_maxrss,peak_sampled_rss_bytes=peak_sampled,termination_reason=reason)
write(dest/'RECEIPT.json',row); row.pop('stage'); book.event(**row)
if code:
    with (root/'FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
print(json.dumps(row),flush=True); sys.exit(code)
