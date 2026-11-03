"""One disclosed R3 endpoint: does the full-wave feedback supply a direction?"""
from pathlib import Path
import json,sys,time,traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.experiments.run_blind import dump,ALLOW
from research.research_max_struct_r4.src.cooperative_geometry import curvature,inherited_hvp

def run():
    start,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/curvature';out.mkdir(parents=True,exist_ok=True)
    path=out/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve prior result')
    model=None
    receipt=dict(run_id='r3_scaled_face_curvature',case_id='obs_b20d5a55677d4eca',
        status='RUNNING',truth_access='NONE',stage='DEV_DIAGNOSTIC',
        physical_convention='exp(-i omega t); scalar-TM disk-pulse VIE; real contrast',
        material_box=[.02,20.],free_face_tolerance=1e-6,
        source='experiments/curvature_diagnostic.py',hash_checks=False,
        question='Is exact residual feedback curvature negative on the current free face?',
        prior='unchanged inherited smooth/TV; exact differentiated TV Hessian')
    dump(path,receipt)
    try:
        with np.load(REPO/'research/research_max_struct_r2/data_public/obs_b20d5a55677d4eca.npz') as f:
            public={k:f[k] for k in f.files if k in ALLOW}
        m=np.load(REPO/'research/research_max_struct_r3/results/development/scaled_original_restart/prediction.npz')['m']
        model=MaxwellModel(public);Y=public['Y']
        val,grad=model.data_value_gradient(m,Y)
        ids=np.flatnonzero((m>.02+1e-6)&(m<20.-1e-6))
        geom=curvature(model,m,Y,ids)
        with model._physics.operation('CURVATURE_EIGEN', symmetric_eigenproblems=1,eigen_dimension=len(ids)):
            values,vectors=la.eigh(geom['H'],check_finite=False)
        h=np.zeros(model.n);h[ids]=vectors[:,0]
        wired=inherited_hvp(model,geom['state'],Y,h)[ids]
        groups=[]
        for i in range(min(3,len(ids))):
            d=vectors[:,i]
            groups.append(dict(eigenvalue=float(values[i]),GN=float(d@geom['GN']@d),
                               feedback=float(d@geom['feedback']@d),prior=float(d@geom['regularizer']@d),
                               gradient_dot=float(grad[ids]@d)))
        np.savez_compressed(out/'geometry.npz',m=m,gradient=grad,indices=ids,
            eigenvalues=values,eigenvectors=vectors,H=geom['H'],GN=geom['GN'],
            feedback=geom['feedback'],regularizer=geom['regularizer'])
        receipt.update(status='COMPLETED_DIAGNOSTIC',objective=float(val),free_count=len(ids),
            eigenpairs=groups,minimum_eigenvalue=float(values[0]),
            projected_gradient_inf=float(la.norm(m-np.clip(m-grad,.02,20.),np.inf)),
            wiring_absolute_error=float(la.norm(wired-geom['H']@vectors[:,0])),
            wiring_reference='one independent analytic HVP of established R4-BASIN formula; no finite difference',
            physical_residual=float(la.norm(geom['state']['Y']-Y)/la.norm(Y)),
            limitations='fixed free face only; no box-critical-cone/global certificate; no material GT scoring')
    except Exception as e:
        receipt.update(status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(physical_cost=dict(model.cost) if model else {},
            wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
        dump(path,receipt)
    print(json.dumps(receipt))
if __name__=='__main__':
    with threadpool_limits(limits=1):run()
