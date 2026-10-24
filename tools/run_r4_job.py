"""Profile a zero-physics training/analysis command under the R4 stage ledger."""
from pathlib import Path
import os,sys,argparse,subprocess,time,signal,resource,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import Book,OUT,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--job',required=True)
 p.add_argument('--max-wall',type=float,default=3600);p.add_argument('--gpu',action='store_true')
 p.add_argument('--training',action='store_true');p.add_argument('--rss-limit-mib',type=float,default=6144)
 p.add_argument('--cpu-parallel',action='store_true',help='released small CPU risk fits only, separate bounded slot')
 p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
 cmd=a.command[1:] if a.command and a.command[0]=='--' else a.command
 if not cmd:p.error('command required')
 if a.cpu_parallel:
  policy=json.loads((OUT/'configs/RESOURCE_POLICY_AMENDMENT_01.json').read_text())
  if not a.training or a.gpu or not policy['released'] or '--device' not in cmd or cmd[cmd.index('--device')+1]!='cpu':raise ValueError('CPU slot policy')
  allowed={'p3r4.train_risk','p3r4.train_risk_bin','p3r4.train_risk_anchorbin','p3r4.train_risk_fusion'}
  if '-m' not in cmd or cmd[cmd.index('-m')+1] not in allowed:raise ValueError('small risk modules only')
  a.rss_limit_mib=min(a.rss_limit_mib,2048)
 b=Book(a.stage,a.job,max_wall=a.max_wall,analysis_only=True,gpu=a.gpu,kind='training_cpu' if a.cpu_parallel else ('training' if a.training else 'analysis'))
 env=os.environ.copy();env.update(PYTHONPATH=str(ROOT/'src'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',
    PYTORCH_MPS_HIGH_WATERMARK_RATIO='0.5',PYTORCH_MPS_LOW_WATERMARK_RATIO='0.4',
    R4_JOBDIR=str(b.jobdir),PYTHONUNBUFFERED='1')
 write(b.jobdir/'COMMAND.json',dict(argv=cmd,cwd=str(ROOT),GPU_allocation='whole MPS job wall' if a.gpu else 'none',physical_calls_allowed=False))
 child=None;error=None;status='FAILED';peak=0;start=time.monotonic()
 try:
  with (b.jobdir/'stdout.log').open('x') as f,(b.jobdir/'stderr.log').open('x') as e:
   child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=f,stderr=e,start_new_session=True)
   write(b.jobdir/'CHILD.json',dict(pid=child.pid,wall_start_unix=time.time()))
   while child.poll() is None:
    if time.monotonic()-start>a.max_wall-15:raise TimeoutError('controlled job reserve minus close margin')
    rss=subprocess.run(['ps','-o','rss=','-p',str(child.pid)],capture_output=True,text=True).stdout.strip()
    if rss:
     peak=max(peak,int(rss)*1024)
     if int(rss)>a.rss_limit_mib*1024:raise MemoryError('bounded RSS policy; child stopped')
    time.sleep(1)
   if child.returncode:raise RuntimeError(f'child exit {child.returncode}; preserve stderr')
   status='COMPLETE'
 except BaseException as e:
  error=repr(e)
  if child is not None and child.poll() is None:
   os.killpg(child.pid,signal.SIGTERM)
   try:child.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
 finally:
  write(b.jobdir/'CHILD_PROFILE.json',dict(returncode=None if child is None else child.returncode,sampled_RSS_peak_bytes=peak,
    child_rusage=resource.getrusage(resource.RUSAGE_CHILDREN),status=status,error=error))
  b.close(status,error)
 print(json.dumps(dict(job=a.job,status=status,error=error,receipt=str(b.jobdir/'CLOSED.json'))))
 if error:raise SystemExit(1)

if __name__=='__main__':main()
