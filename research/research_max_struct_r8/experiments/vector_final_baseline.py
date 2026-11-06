from pathlib import Path
import subprocess,sys,json,time
ROOT=Path(__file__).resolve().parents[1];logs=ROOT/'delegated/vector3d'
for name in ['vector_pr240_fourier','vector_report']:
 start=time.perf_counter()
 with (logs/(name+'.log')).open('w') as f:p=subprocess.run([sys.executable,str(ROOT/'experiments'/(name+'.py'))],stdout=f,stderr=subprocess.STDOUT)
 record=dict(job=name,returncode=p.returncode,wall_s=time.perf_counter()-start)
 with (logs/'final_baseline_queue.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
 print(json.dumps(record),flush=True)
