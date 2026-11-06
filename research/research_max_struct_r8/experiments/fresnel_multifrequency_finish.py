"""Run the one frozen6+8 followup only after all previously authorized work."""
from pathlib import Path
import subprocess,time,sys,json
ROOT=Path(__file__).resolve().parents[1];start=time.perf_counter()
while not (ROOT/'delegated/fresnel/FINALIZATION_RECEIPT.json').exists():time.sleep(2)
jobs=[]
for mode in ['control_8','joint_6_8']:
 f=(ROOT/'delegated/fresnel'/f'multifrequency_{mode}.log').open('w');p=subprocess.Popen([sys.executable,str(ROOT/'experiments/fresnel_multifrequency.py'),'--mode',mode],stdout=f,stderr=subprocess.STDOUT);jobs.append((p,f,mode))
records=[]
for p,f,mode in jobs:
 code=p.wait();f.close();records.append(dict(mode=mode,exitcode=code));print(json.dumps(records[-1]),flush=True)
with (ROOT/'delegated/fresnel/MULTIFREQUENCY_REPORT.log').open('w') as f:r=subprocess.run([sys.executable,str(ROOT/'experiments/fresnel_report.py')],stdout=f,stderr=subprocess.STDOUT)
(ROOT/'delegated/fresnel/MULTIFREQUENCY_BATCH_RECEIPT.json').write_text(json.dumps(dict(jobs=records,report_exitcode=r.returncode,wall_seconds_including_predecessor_wait=time.perf_counter()-start,maximum_physics_processes=2),indent=2))
