"""Bounded new adapter checks and independent analytic Mie comparison."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,json
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r13.src.vector_adapter import DenseVectorModel
from research.research_max_struct_r1.src.vector_backend import TensorDDA,sphere_grid
from research.research_max_struct_r1.experiments.qualify_backend import mie
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump


def main():
    start,cpu=time.perf_counter(),time.process_time();model=DenseVectorModel()
    rng=np.random.default_rng(20261113032);m=.5+2*rng.random(model.n)
    Y=model.forward(np.full(model.n,1.1)).copy();st=model.state(m)
    yy,jj,qq=model.legacy.physical(m)
    agreement=dict(prediction=float(np.linalg.norm(st['Y']-yy)/np.linalg.norm(yy)),
       current=float(np.linalg.norm(st['j']-jj)/np.linalg.norm(jj)),
       total_field=float(np.linalg.norm(st['Q']-qq)/np.linalg.norm(qq)))
    f,g=model.data_value_gradient(m,Y);fl,gl=model.legacy.physical(m,Y)
    agreement.update(objective_absolute=abs(f-fl),gradient_relative=float(np.linalg.norm(g-gl)/np.linalg.norm(gl)))
    direction=rng.standard_normal(model.n);direction/=np.linalg.norm(direction)
    step=1e-5;fd=(model.physical_score(m+step*direction,Y)-model.physical_score(m-step*direction,Y))/(2*step)
    analytic=float(g@direction)
    directional=dict(analytic=analytic,central_difference=fd,step=step,
        relative_difference=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-15))
    theta=np.linspace(0,np.pi,181);directions=np.c_[np.sin(theta),np.zeros_like(theta),np.cos(theta)]
    etheta=np.c_[np.cos(theta),np.zeros_like(theta),-np.sin(theta)]
    rows=[]
    for epsilon in (3.,9.):
        _,s2,sigma=mie(.09,2*np.pi,epsilon,np.cos(theta));reference=1j/(2*np.pi)*s2[:,None]*etheta
        for h in (.04,.03):
            points,volume=sphere_grid(.09,h);backend=TensorDDA(points,volume)
            state=backend.solve(np.broadcast_to(epsilon*np.eye(3),(len(points),3,3)),[[0,0,1]],[[1,0,0]])
            far=state.far(directions)['E'].reshape(-1,3)
            rows.append(dict(epsilon=epsilon,radius=.09,h=h,N=len(points),
                cell_phase=float(2*np.pi*h*np.sqrt(epsilon)),
                analytic_sphere_far_relative_error=float(np.linalg.norm(far-reference)/np.linalg.norm(reference)),
                sphere_volume_relative_error=float(abs(len(points)*volume-4*np.pi*.09**3/3)/(4*np.pi*.09**3/3)),
                physical_solve=state.receipt))
    dump(HERE/'QUALIFICATION.json',dict(status='COMPLETE',adapter_legacy_agreement=agreement,
        objective_directional_derivative=directional,initial_physical_residual=st['state_relative_residual'],
        analytic_Mie_rows=rows,counts=dict(model.cost),CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-start,
        interpretation='Bounded independent analytic sphere validation at actual h and contrast; does not certify nonsphere continuum accuracy or material recovery. Crossgrid mismatch is separate.'))
    print(json.dumps(dict(status='COMPLETE',agreement=agreement,directional=directional,Mie=rows)),flush=True)
if __name__=='__main__':main()
