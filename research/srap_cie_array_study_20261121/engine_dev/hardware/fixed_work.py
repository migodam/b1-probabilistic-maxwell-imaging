"""Use the solver's original thread initialization before importing numerical libraries.
Payload contains exactly the same public 32² model and a fixed, exposed material.
"""
import base64
import io
import json
import os
import platform
import sys
import time

import run_study
from threadpoolctl import threadpool_info
import numpy as np
import scipy
from scipy.fft import get_workers

def run(payload):
    with np.load(io.BytesIO(payload), allow_pickle=False) as f:
        public = dict(f)
    m=public.pop('material')
    model=run_study.la.PaperModel(public,{'side':32},run_study.la.Clock(180.))
    model.G  # outside timing, same exact FFT kernel on both hosts
    arrays=[]
    for i in range(3):
        tc,tw=time.process_time(),time.perf_counter()
        phi,res,_=model.objective(m,public['Y'])
        arrays.append(dict(seconds=time.perf_counter()-tw,cpu_seconds=time.process_time()-tc,phi=phi,residual=res))
    tc,tw=time.process_time(),time.perf_counter()
    for i in range(12):
        out=model.green.apply(model.E)
        back=model.green.adjoint(out)
    ft=dict(seconds=time.perf_counter()-tw,cpu_seconds=time.process_time()-tc,norm=float(np.linalg.norm(back)))
    info={k:os.getenv(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
    print(json.dumps(dict(host=platform.node(),platform=platform.platform(),python=sys.version,executable=sys.executable,
          numpy=np.__version__,scipy=scipy.__version__,threadpools=threadpool_info(),environment=info,
          scipy_fft_workers=get_workers(),objective_3=arrays,green_pairs_12=ft,counts=dict(model.clock.counts)),indent=2))

if __name__=='__main__':
    run(base64.b64decode(PAYLOAD))
