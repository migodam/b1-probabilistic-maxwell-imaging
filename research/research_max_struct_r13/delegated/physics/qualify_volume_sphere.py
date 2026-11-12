"""Independent Mie control with cell-volume material averages, matching targets.

Retains separate historical staircase-sphere control in QUALIFICATION.json.
No change to backend, material model, cohort, acquisitions or noise budget.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,numpy as np
from numpy.polynomial.legendre import leggauss
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r1.src.vector_backend import TensorDDA
from research.research_max_struct_r1.experiments.qualify_backend import mie
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump


def sphere_cell_volumes(points,h,radius,order=48):
    """Integrate exact z intersections over tensor Gauss x/y quadrature."""
    nodes,w=leggauss(order);x=points[:,0,None,None]+h/2*nodes[None,:,None]
    y=points[:,1,None,None]+h/2*nodes[None,None,:]
    z=np.sqrt(np.maximum(radius**2-x*x-y*y,0))
    zlo=points[:,2,None,None]-h/2;zhi=zlo+h
    height=np.maximum(np.minimum(z,zhi)-np.maximum(-z,zlo),0)
    return np.einsum('nij,i,j->n',height,w,w)/(4*h)


def main():
    started,cpu=time.perf_counter(),time.process_time();rows=[]
    theta=np.linspace(0,np.pi,181);directions=np.c_[np.sin(theta),np.zeros_like(theta),np.cos(theta)]
    etheta=np.c_[np.cos(theta),np.zeros_like(theta),-np.sin(theta)]
    for epsilon in (3.,9.):
        _,s2,sigma=mie(.09,2*np.pi,epsilon,np.cos(theta));ref=1j/(2*np.pi)*s2[:,None]*etheta
        for h in (.04,.03):
            side=int(round(.24/h));axis=(np.arange(side)+.5)*h-.12
            points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3)
            fractions=sphere_cell_volumes(points,h,.09);active=fractions>1e-14
            points=points[active];fractions=fractions[active]
            model=TensorDDA(points,h**3);eps=(1+(epsilon-1)*fractions)[:,None,None]*np.eye(3)
            state=model.solve(eps,[[0,0,1]],[[1,0,0]])
            far=state.far(directions)['E'].reshape(-1,3)
            rows.append(dict(epsilon=epsilon,radius=.09,h=h,N=len(points),
                material_average='chi_i=(epsilon-1)*cell_sphere_volume_fraction, 48x48 Gauss plus exact z intersection',
                exact_volume=4*np.pi*.09**3/3,integrated_volume=float(sum(fractions)*h**3),
                cell_phase_max=float(2*np.pi*h*np.sqrt(epsilon)),
                analytic_sphere_far_relative_error=float(np.linalg.norm(far-ref)/np.linalg.norm(ref)),
                physical_solve=state.receipt))
    dump(HERE/'VOLUME_MIE_QUALIFICATION.json',dict(status='COMPLETE',rows=rows,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started,
        qualification='bounded analytic sphere evidence at target cell size and contrast, not a non-sphere continuum certificate; finite cell volume averaging is an approximation'))
    print([(r['epsilon'],r['h'],r['analytic_sphere_far_relative_error']) for r in rows],flush=True)
if __name__=='__main__':main()
