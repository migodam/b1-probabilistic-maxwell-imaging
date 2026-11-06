"""Finite two-process queue; each child fixes BLAS/PyTorch to one thread."""
from pathlib import Path
import subprocess,sys,time,json
ROOT=Path(__file__).resolve().parents[1]
jobs=[(case,arm,'uniform') for case in ['dielTM_dec8f','twodielTM_8f'] for arm in ['cie_single','cie_dual','global_dual','a_material','a_plus'] if not (ROOT/'RESULTS/fresnel'/f'{case}_{arm}_multipole_full_v2'/'RECEIPT.json').exists()]
jobs += [(case,arm,'neural') for case in ['dielTM_dec8f','twodielTM_8f'] for arm in ['cie_dual','a_plus']]
running=[];receipts=[];t=time.perf_counter()
while jobs or running:
 while jobs and len(running)<2:
  case,arm,init=jobs.pop(0);log=ROOT/'delegated/fresnel'/f'{case}_{arm}_{init}.log';fh=log.open('w');p=subprocess.Popen([sys.executable,str(ROOT/'experiments/fresnel_invert.py'),'--case',case,'--arm',arm,'--initialization',init],stdout=fh,stderr=subprocess.STDOUT)
  running.append((p,fh,case,arm,init,time.perf_counter()));print('START',case,arm,init,flush=True)
 for row in list(running):
  p,fh,case,arm,init,start=row;status=p.poll()
  if status is not None:
   fh.close();receipts.append(dict(case=case,arm=arm,initialization=init,exitcode=status,wall_seconds=time.perf_counter()-start));running.remove(row);print('DONE',case,arm,init,status,flush=True)
 if running:time.sleep(1)
(ROOT/'delegated/fresnel/BATCH_RECEIPT.json').write_text(json.dumps(dict(jobs=receipts,wall_seconds=time.perf_counter()-t,parallel_limit=2,BLAS_threads_each=1),indent=2))
