"""Single scientific job: bounded CPU, BLAS=1, persistent stdout and receipts."""
from pathlib import Path
import sys,os,time,resource,subprocess,json,math
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'src'))
from wq.accounting import Book,SearchBook
is_search=os.environ.get('WQ_STAGE')=='p0b'
book=SearchBook(root/'results/p0b') if is_search else Book(root/'results')
label=sys.argv[1]; command=sys.argv[2:]
dest=root/'results/jobs'/label;dest.mkdir(parents=True,exist_ok=False)
left=7200-book.totals()['cpu_seconds']
if left<=1:raise RuntimeError('CPU budget exhausted')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=os.pathsep.join([str(root/'src'),str(root/'vendor'),str(root/'.deps')]),MPLCONFIGDIR=str(root/'results/mpl'),WQ_CPU_START='0')
env['WQ_JOB_ID']=label
def limit():resource.setrlimit(resource.RLIMIT_CPU,(math.ceil(left),math.ceil(left)+1))
t=time.perf_counter();c0=resource.getrusage(resource.RUSAGE_CHILDREN)
with (dest/'OUTPUT.log').open('w') as f:
    p=subprocess.run(command,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,preexec_fn=limit)
c1=resource.getrusage(resource.RUSAGE_CHILDREN);wall=time.perf_counter()-t;cpu=c1.ru_utime+c1.ru_stime-c0.ru_utime-c0.ru_stime
row=dict(event='cost',label=label,command=command,exit_code=p.returncode,status='COMPLETE' if p.returncode==0 else 'FAILED',wall_seconds=wall,process_cpu_seconds=cpu,charged_cpu_seconds=cpu,peak_rss_bytes=c1.ru_maxrss)
(dest/'RECEIPT.json').write_text(json.dumps(row,indent=2));book.event(**row)
if p.returncode:
    with (root/'results/FAILURE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
print(json.dumps(row),flush=True);sys.exit(p.returncode)
