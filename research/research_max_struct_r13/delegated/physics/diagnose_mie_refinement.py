"""One finer occupied-sphere diagnostic at the same fixed physical radius.

Retains original negatives. The sphere is known only in this forward control;
unknown-support inversion remains the entire cube. No inverse code imports it.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,time,numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];sys.path.insert(0,str(REPO))
from research.research_max_struct_r1.src.vector_backend import TensorDDA,sphere_grid
from research.research_max_struct_r1.experiments.qualify_backend import mie
from research.research_max_struct_r13.delegated.physics.generate_inputs import dump
from research.research_max_struct_r13.delegated.physics.qualify_volume_sphere import sphere_cell_volumes

def main():
    start,cpu=time.perf_counter(),time.process_time();rows=[];radius=.09;epsilon=9.
    theta=np.linspace(0,np.pi,181);dirs=np.c_[np.sin(theta),np.zeros_like(theta),np.cos(theta)]
    etheta=np.c_[np.cos(theta),np.zeros_like(theta),-np.sin(theta)]
    _,s2,_=mie(radius,2*np.pi,epsilon,np.cos(theta));ref=1j/(2*np.pi)*s2[:,None]*etheta
    for h in (.02,.015):
        extent=int(np.ceil(radius/h+.5));axis=np.arange(-extent,extent+1)*h
        pts=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3)
        frac=sphere_cell_volumes(pts,h,radius);active=frac>1e-14;pts,frac=pts[active],frac[active]
        model=TensorDDA(pts,h**3);eps=(1+8*frac)[:,None,None]*np.eye(3)
        st=model.solve(eps,[[0,0,1]],[[1,0,0]])
        far=st.far(dirs)['E'].reshape(-1,3)
        row=dict(h=h,radius=radius,epsilon=epsilon,N=len(pts),
            integrated_volume=float(sum(frac)*h**3),analytic_volume=4*np.pi*radius**3/3,
            material_average='same 48x48 Gauss with exact z intersection',
            grid_origin='cell center at0, occupied and fractional boundary cells only',
            cell_phase_max=float(2*np.pi*h*3),
            analytic_sphere_far_relative_error=float(np.linalg.norm(far-ref)/np.linalg.norm(ref)),
            fitted_complex_gain_diagnostic=complex(np.vdot(ref,far)/np.vdot(ref,ref)),
            physical_solve=st.receipt)
        gain=row.pop('fitted_complex_gain_diagnostic');row['fitted_complex_gain_diagnostic']=[gain.real,gain.imag]
        row['shape_after_diagnostic_gain_relative_error']=float(np.linalg.norm(far-gain*ref)/np.linalg.norm(far))
        rows.append(row);print((h,len(pts),row['analytic_sphere_far_relative_error']),flush=True)
    staircase=[]
    for h in (.02,.015):
        pts,volume=sphere_grid(radius,h);model=TensorDDA(pts,volume)
        st=model.solve(np.broadcast_to(epsilon*np.eye(3),(len(pts),3,3)),[[0,0,1]],[[1,0,0]])
        far=st.far(dirs)['E'].reshape(-1,3)
        staircase.append(dict(h=h,radius=radius,epsilon=epsilon,N=len(pts),
            material_rule='center-in-sphere occupied homogeneous epsilon9 voxels',
            integrated_volume=float(len(pts)*volume),analytic_volume=4*np.pi*radius**3/3,
            analytic_sphere_far_relative_error=float(np.linalg.norm(far-ref)/np.linalg.norm(ref)),
            physical_solve=st.receipt))
    dump(HERE/'MIE_REFINEMENT_DIAGNOSIS.json',dict(status='COMPLETE',rows=rows,staircase_control=staircase,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
        diagnostic_gain_used_for_any_prediction=False,
        interpretation='Same Maxwell normalization and selfterm, fixed analytic radius and epsilon; refine h only. Volume averaging and grid-origin effects remain stated. Small solve residual alone is not continuum accuracy.'))
if __name__=='__main__':main()
