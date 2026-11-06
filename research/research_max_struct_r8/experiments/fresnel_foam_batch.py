"""Corrected-paper-stop Foam5 finite queue, then resume preserved2GHz queue."""
from pathlib import Path
import subprocess,sys,time,json,os,signal,argparse
ap=argparse.ArgumentParser();ap.add_argument('--resume-2ghz-controller-pid',type=int);args=ap.parse_args()
ROOT=Path(__file__).resolve().parents[1];jobs=['cie_dual','a_plus','cie_single','global_dual','a_material'];running=[];receipts=[];t=time.perf_counter()
while jobs or running:
 while jobs and len(running)<2:
  arm=jobs.pop(0);fh=(ROOT/'delegated/fresnel'/f'Foam_n60_{arm}_v3.log').open('w');p=subprocess.Popen([sys.executable,str(ROOT/'experiments/fresnel_foam_invert.py'),'--arm',arm,'--side','60','--attempt','v3'],stdout=fh,stderr=subprocess.STDOUT);running.append((p,fh,arm,time.perf_counter()));print('START',arm,flush=True)
 for row in list(running):
  p,fh,arm,start=row;status=p.poll()
  if status is not None:fh.close();receipts.append(dict(arm=arm,exitcode=status,wall_seconds=time.perf_counter()-start));running.remove(row);print('DONE',arm,status,flush=True)
 if running:time.sleep(1)
(ROOT/'delegated/fresnel/FOAM_BATCH_RECEIPT.json').write_text(json.dumps(dict(jobs=receipts,wall_seconds=time.perf_counter()-t,parallel_limit=2,BLAS_threads_each=1),indent=2))
if args.resume_2ghz_controller_pid is not None:
 try:os.kill(args.resume_2ghz_controller_pid,signal.SIGCONT);print('RESUMED_EXPLICIT_2GHZ_QUEUE',flush=True)
 except ProcessLookupError:print('EXPLICIT_2GHZ_CONTROLLER_NOT_RUNNING',flush=True)
