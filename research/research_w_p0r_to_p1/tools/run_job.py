"""Serial one-thread scientific job and process-CPU accounting."""
import sys,os,time,json,resource,subprocess,math
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from wr.accounting import Book,write
stage,label=sys.argv[1:3]; command=sys.argv[3:]; book=Book(stage)
if '/' in command[0] and not os.path.isabs(command[0]):command[0]=str(Path(command[0]).absolute())
dest=root/'results/jobs'/label; dest.mkdir(parents=True,exist_ok=False)
dep='/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/Research_W/research-w-continuation/.deps'
hist=root.parent/'research_w_continuation'
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=os.pathsep.join([str(root/'src'),str(hist/'src'),str(hist/'vendor'),dep]),MPLCONFIGDIR=str(root/'results/mpl'))
left=book.config['stages'][stage]['cpu_seconds']-book.totals()[stage]['cpu_seconds']
def limit(): resource.setrlimit(resource.RLIMIT_CPU,(max(1,math.ceil(left)),max(2,math.ceil(left)+1)))
start=time.perf_counter(); before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (dest/'OUTPUT.log').open('w') as stream:
    proc=subprocess.run(command,cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,preexec_fn=limit)
after=resource.getrusage(resource.RUSAGE_CHILDREN)
cpu=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
row=dict(event='job_cost',stage=stage,label=label,command=command,status='COMPLETE' if proc.returncode==0 else 'FAILED',returncode=proc.returncode,wall_seconds=time.perf_counter()-start,charged_cpu_seconds=cpu,peak_rss_bytes=after.ru_maxrss)
write(dest/'RECEIPT.json',row); row.pop('stage'); book.event(**row)
if proc.returncode:
    with (root/'FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
print(json.dumps(row),flush=True); sys.exit(proc.returncode)
