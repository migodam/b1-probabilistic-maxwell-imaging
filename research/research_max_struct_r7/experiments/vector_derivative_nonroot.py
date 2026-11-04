"""Non-root directional derivatives: the dP/dchi term must be active."""
from pathlib import Path
import os,sys,json,time
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r7.src.vector_inverse import VectorInverse
with threadpool_limits(limits=1):
    tick=time.perf_counter();model=VectorInverse();x=np.linspace(.1,1.4,model.N);pred,J,E,_=model.physical(x)
    rng=np.random.default_rng(716);J=1.1*J+.05*(rng.normal(size=J.shape)+1j*rng.normal(size=J.shape));Q=model.inc+model.G@J
    v=rng.normal(size=model.N);v/=la.norm(v);step=2e-6;rows={}
    for kind in ['ordinary_joint','dual_like_joint','A_material_joint']:
        f,g=model.material(x,J,Q,kind);fp=model.material(x+step*v,J,Q,kind)[0];fm=model.material(x-step*v,J,Q,kind)[0];fd=(fp-fm)/(2*step);analytic=float(g@v)
        rows[kind]=dict(loss=f,analytic=analytic,centered_difference=fd,relative_error=float(abs(analytic-fd)/max(abs(fd),abs(analytic),1e-12)))
    receipt=dict(status='COMPLETED',reason='initial root derivative check mostly tested prior; this additional check activates nonzero constitutive residual and derivative of material-dependent A preconditioner',relative_constitutive_residual=float(la.norm(J-np.repeat(x,3)[:,None]*Q)/la.norm(J)),derivatives=rows,counts=model.counts,wall_s=time.perf_counter()-tick)
    (ROOT/'RESULTS/vector_3d/NONROOT_DERIVATIVES.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
    if any(r['relative_error']>2e-5 for r in rows.values()):raise RuntimeError('nonroot derivative failed')
