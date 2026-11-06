"""Task-local finite queue handoff; keep maximum two active physics children.

Uses only explicitly supplied, command-validated controller PIDs. No physics
child is interrupted. Overlap pending 2GHz work with the final sole Foam arm.
"""
from pathlib import Path
import subprocess,sys,time,json,os,signal,argparse
ROOT=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser();ap.add_argument('--foam-controller-pid',type=int,required=True);ap.add_argument('--regression-controller-pid',type=int,required=True);args=ap.parse_args()
def command(pid):
 q=subprocess.run(['ps','-p',str(pid),'-o','command='],capture_output=True,text=True);return q.stdout.strip()
def owns(pid,script):return str(ROOT/'experiments'/script) in command(pid) or command(pid).endswith('research/research_max_struct_r8/experiments/'+script)
first=['cie_dual','a_plus','cie_single','global_dual'];last=ROOT/'RESULTS/fresnel/FoamTwinDielTM_8GHz_n60_a_material_v3';start=time.perf_counter()
while not all((ROOT/'RESULTS/fresnel'/f'FoamTwinDielTM_8GHz_n60_{a}_v3'/'RECEIPT.json').exists() for a in first) or not (last/'FROZEN_CONFIG.json').exists():time.sleep(2)
if not owns(args.foam_controller_pid,'fresnel_foam_batch.py') or not owns(args.regression_controller_pid,'fresnel_batch.py'):raise RuntimeError('Explicit controller command identity no longer matches; no signal sent')
os.kill(args.foam_controller_pid,signal.SIGSTOP)
jobs=[(case,arm,'uniform') for case in ['dielTM_dec8f','twodielTM_8f'] for arm in ['cie_single','cie_dual','global_dual','a_material','a_plus']]
jobs += [(case,arm,'neural') for case in ['dielTM_dec8f','twodielTM_8f'] for arm in ['cie_dual','a_plus']]
def result_folder(case,arm,init):return ROOT/'RESULTS/fresnel'/(f'{case}_{arm}_multipole_full_v2'+('_nn' if init=='neural' else ''))
already=[j for j in jobs if (result_folder(*j)/'RECEIPT.json').exists()];jobs=[j for j in jobs if j not in already];running=[];receipts=[]
print('HANDOFF_PAUSED_CONTROLLERS_ONLY',len(jobs),'remaining',flush=True)
try:
 while jobs or running or not ((last/'RECEIPT.json').exists() or (last/'FAILED.json').exists()):
  foam_active=not ((last/'RECEIPT.json').exists() or (last/'FAILED.json').exists())
  while jobs and len(running)+int(foam_active)<2:
   case,arm,init=jobs.pop(0);fh=(ROOT/'delegated/fresnel'/f'{case}_{arm}_{init}_handoff.log').open('w');p=subprocess.Popen([sys.executable,str(ROOT/'experiments/fresnel_invert.py'),'--case',case,'--arm',arm,'--initialization',init],stdout=fh,stderr=subprocess.STDOUT);running.append((p,fh,case,arm,init,time.perf_counter()));print('START',case,arm,init,flush=True)
  for row in list(running):
   p,fh,case,arm,init,t=row;code=p.poll()
   if code is not None:fh.close();receipts.append(dict(case=case,arm=arm,initialization=init,exitcode=code,wall_seconds=time.perf_counter()-t));running.remove(row);print('DONE',case,arm,init,code,flush=True)
  if jobs or running or foam_active:time.sleep(1)
 (ROOT/'delegated/fresnel/BATCH_RECEIPT.json').write_text(json.dumps(dict(jobs=receipts,already_completed=already,wall_seconds=time.perf_counter()-start,parallel_physics_limit_including_Foam=2,BLAS_threads_each=1,controller_handoff=True),indent=2))
 if not owns(args.regression_controller_pid,'fresnel_batch.py'):raise RuntimeError('Idle original controller identity changed; not signalled')
 # The original regression controller is paused, with no live numeric child.
 # Its completed child receipts are retained; new queue has done every pending job.
 os.kill(args.regression_controller_pid,signal.SIGKILL)
finally:
 if owns(args.foam_controller_pid,'fresnel_foam_batch.py'):os.kill(args.foam_controller_pid,signal.SIGCONT)
(ROOT/'delegated/fresnel/QUEUE_HANDOFF_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETED',jobs=receipts,physics_children_interrupted=0,original_idle_regression_controller_retired=True,wall_seconds=time.perf_counter()-start),indent=2))
