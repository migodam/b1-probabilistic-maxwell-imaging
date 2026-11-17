"""Finite serial executor, subprocess CPU/RSS receipts including startup."""
from common import *
import subprocess,signal
budget=7200.;spent=0.;receipts=[]
for path in [ROOT/'results/debug.json',ROOT/'results/generation_receipt.json']:
 if path.exists():
  d=json.loads(path.read_text());spent+=d.get('receipt',d).get('CPU_seconds',0.)
for case in ['T0','T1','T2','T3','T4','T5']:
 if spent>=budget:break
 log=(ROOT/'results'/(case+'_process.log')).open('x');begin=time.monotonic()
 def limits():resource.setrlimit(resource.RLIMIT_CPU,(int(min(1100,budget-spent)),int(min(1100,budget-spent))+1))
 child=subprocess.Popen([sys.executable,str(ROOT/'src/experiment.py'),case],stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits,start_new_session=True)
 maxrss=0;reason=None
 while True:
  pid,status,usage=os.wait4(child.pid,os.WNOHANG)
  if pid:break
  try:rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(child.pid)],text=True).strip())*1024;maxrss=max(maxrss,rss)
  except Exception:rss=0
  if rss>8*1024**3 or time.monotonic()-begin>1800:
   reason='RSS_OR_WALL_HARD_LIMIT';os.killpg(child.pid,signal.SIGKILL)
  time.sleep(.5)
 cpu=usage.ru_utime+usage.ru_stime;spent+=cpu;rc=os.waitstatus_to_exitcode(status);child.returncode=rc;log.close()
 row=dict(case=case,returncode=rc,CPU_seconds_including_reaped_children=cpu,wall_seconds=time.monotonic()-begin,peak_rss_bytes=max(maxrss,usage.ru_maxrss),cumulative_CPU_seconds=spent,reason=reason)
 receipts.append(row);save(ROOT/'results'/(case+'_process_receipt.json'),row);print(row,flush=True)
 if rc!=0:break # preserve failure and stop for engineering diagnosis; never regenerate targets.
save(ROOT/'results/executor_receipt.json',dict(receipts=receipts,total_CPU_seconds=spent,budget=budget))
