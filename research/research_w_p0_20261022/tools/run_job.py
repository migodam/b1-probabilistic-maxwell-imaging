"""Run substantive jobs with BLAS=1 and retain failures and CPU/wall costs."""
import sys,os,time,subprocess,json,resource
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
from research_w.accounting import Book
label=sys.argv[1];cmd=sys.argv[2:];book=Book(root/'results')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=str(root/'src'),MPLCONFIGDIR=str(root/'results/mpl'))
env['W_JOB_DEADLINE_MONOTONIC']=str(time.monotonic()+max(0,7200-book.totals()['cpu_seconds']))
dest=root/'results/jobs'/label;dest.mkdir(parents=True,exist_ok=False)
t=time.perf_counter();r0=resource.getrusage(resource.RUSAGE_CHILDREN)
with (dest/'OUTPUT.log').open('w') as f:
 p=subprocess.run(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
r1=resource.getrusage(resource.RUSAGE_CHILDREN);wall=time.perf_counter()-t;cpu=r1.ru_utime+r1.ru_stime-r0.ru_utime-r0.ru_stime
receipt=dict(event='cost',label=label,command=cmd,status='COMPLETE' if p.returncode==0 else 'FAILED',exit_code=p.returncode,wall_seconds=wall,process_cpu_seconds=cpu,charged_cpu_seconds=max(wall,cpu),peak_rss_bytes=r1.ru_maxrss)
(dest/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));book.event(**receipt)
print(json.dumps(receipt));sys.exit(p.returncode)
