"""Finish only after both finite physics queues complete; never add a service."""
from pathlib import Path
import subprocess,time,sys,json
ROOT=Path(__file__).resolve().parents[1];t=time.perf_counter()
for name in ['FOAM_BATCH_RECEIPT.json','BATCH_RECEIPT.json']:
 while not (ROOT/'delegated/fresnel'/name).exists():time.sleep(2)
rows=[]
for script in ['fresnel_foam_nie.py','fresnel_illumination_diagnostics.py','fresnel_report.py']:
 start=time.perf_counter();log=ROOT/'delegated/fresnel'/f'{script}.log'
 with log.open('w') as f:r=subprocess.run([sys.executable,str(ROOT/'experiments'/script)],stdout=f,stderr=subprocess.STDOUT)
 rows.append(dict(script=script,exitcode=r.returncode,wall_seconds=time.perf_counter()-start));print(json.dumps(rows[-1]),flush=True)
(ROOT/'delegated/fresnel/FINALIZATION_RECEIPT.json').write_text(json.dumps(dict(stages=rows,wall_seconds_including_queue_wait=time.perf_counter()-t),indent=2))
