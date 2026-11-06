"""Single authorized full-vs-cropped finite-DDA Green action and timing."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import sys,time,json,resource
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
from scipy import linalg as la
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel,four_cubes
wall,cpu=time.perf_counter(),time.process_time();m=four_cubes(60)
model=OriginalVectorModel(60,active_indices=np.flatnonzero(m!=0),generation_only=True)
rng=np.random.default_rng(61);v=rng.normal(size=model.M)+1j*rng.normal(size=model.M)
times=[]
start,cstart=time.perf_counter(),time.process_time();a=model.green(v)
times.append(dict(embedding='full60',wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cstart))
model.configure_embedding(True)
start,cstart=time.perf_counter(),time.process_time();b=model.green(v)
times.append(dict(embedding='exact_active_boundingbox',wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cstart))
relative=float(la.norm(a-b)/la.norm(a));result=dict(scope='ONE exact zero-exterior embedding check; no solves',
    relative_action_error=relative,PASS=relative<1e-12,times=times,setup=model.setup,counts=dict(model.counts),
    CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    unchanged='original60 h=.05; absolute source/cell coordinates; self; active current unknowns; S and inc')
path=Path(__file__).resolve().parents[1]/'RESULTS'/'original_vector3d'/'CROP_CHECK.json'
path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if not result['PASS']:raise RuntimeError('Exact embedding check failed')
