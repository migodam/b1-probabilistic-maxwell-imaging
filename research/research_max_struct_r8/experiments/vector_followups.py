from pathlib import Path
import subprocess,sys,json,time
ROOT=Path(__file__).resolve().parents[1];logs=ROOT/'delegated/vector3d'
# Finite queue, one numeric subprocess at a time. No daemon or service.
for name in ['vector_mie_geometry','vector_original_size','vector_pr240_benchmark','vector_refinement','vector_report']:
    start=time.perf_counter()
    with (logs/(name+'.log')).open('w') as f:
        p=subprocess.run([sys.executable,str(ROOT/'experiments'/(name+'.py'))],stdout=f,stderr=subprocess.STDOUT)
    record=dict(job=name,returncode=p.returncode,wall_s=time.perf_counter()-start)
    with (logs/'followup_queue.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
    print(json.dumps(record),flush=True)
