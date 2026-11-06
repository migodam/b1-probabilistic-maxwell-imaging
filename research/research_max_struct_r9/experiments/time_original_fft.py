"""One parent-authorized <10 CPU-second FFT worker engineering probe."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import sys,time,json,resource
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
from scipy import fft
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel,four_cubes
wall,cpu=time.perf_counter(),time.process_time()
m=four_cubes(60);model=OriginalVectorModel(60,active_indices=np.flatnonzero(m!=0),generation_only=True)
rng=np.random.default_rng(60);v=rng.normal(size=model.M)+1j*rng.normal(size=model.M)
records=[]
for workers in (1,2,4):
    with fft.set_workers(workers):
        start,cstart=time.perf_counter(),time.process_time();model.green(v)
        records.append(dict(workers=workers,wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cstart))
result=dict(scope='ONE_ENGINEERING_PROBE; no solves or inversion, same cached original60 Green action',records=records,
    CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
    setup=model.setup,counts=dict(model.counts),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
path=Path(__file__).resolve().parents[1]/'RESULTS'/'original_vector3d'/'FFT_WORKER_PROBE.json'
path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
