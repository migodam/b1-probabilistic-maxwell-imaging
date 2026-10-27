"""Initial G0 fixed-space qualification. Requires an explicit new G0 grant.

This does NOT by itself complete the 12-parent information/posterior stress
panel. It tests the prerequisite discretization before any large port setup.
"""
from pathlib import Path
import gc,json,time
import numpy as np
from scipy.sparse.linalg import splu
from skfem import asm
from bgga.accounting import BudgetBook,_default
from bgga.fixed_mesh import public_fixed_mesh,TargetAssembler
from bgga.operators import relative_backward,FixedSpaceCache
from wslam.fem_backend import FEM,PMLConfig,_contrast_mass,_forcing,line_green
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import pack

ROOT=Path(__file__).resolve().parents[1]

def background_system(p,world):
    eps=p._material(world.coefficient(p.freq,384,include_target=False))
    m=asm(_contrast_mass,p.material_basis,contrast=(eps-1).reshape(p._material_shape))
    L=(p._background_matrix-p.k**2*m[p._free_dofs][:,p._free_dofs]).tocsc()
    inc=line_green(p.material_points,p.tx,p.k)
    b=np.column_stack([asm(_forcing,p.material_basis,
        load=(p.k**2*(eps-1)*inc[:,j]).reshape(p._material_shape))[p._free_dofs] for j in range(len(p.tx))])
    return L,b

def main():
    book=BudgetBook();book.authorization('G0')  # before mesh/assembly/physical work
    out=ROOT/'results'/'g0_qualification';out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((ROOT/'G0_CASES.json').read_text());rows=[]
    tx=np.array([[-1.5,0],[1.5,0]]);rx=np.column_stack((np.linspace(-2,2,16),np.zeros(16)))
    p=FEM(1.,tx,rx,book,h=.17,order=3,mesh=public_fixed_mesh(),phase='G0',residual_tolerance=1e-10)
    ass=TargetAssembler(p.basis,p._free_dofs);R=p._receiver_probe[:,p._free_dofs]
    direct=line_green(rx,tx,p.k).T
    reference=Evaluator(book,'G0',h=.17,n=384,grid_refinement=3)
    for case in cfg['qualification']:
        before=time.process_time()
        w=World(np.array(case['support'],bool),case['epsinf'],case['ell'],case['pose'],
            radius=case['radius'],contrast=case['contrast'],target_loss=case['target_loss'])
        L,b=background_system(p,w)
        a=ass.assemble(w.pose,w.radius,p.k,w.contrast+1j*w.target_loss,
            lambda pts:line_green(pts,tx,p.k))
        full=(L-a.Vtarget).tocsc();rhs=b+a.bX
        book.reserve('G0','qualification fixed-space full LU',factorizations=1,case=case['id'])
        factor=splu(full,permc_spec='COLAMD')
        book.reserve('G0','qualification fixed-space primal',rhs=2,case=case['id'])
        u=factor.solve(rhs);data=direct+(R@u).T;defect=relative_backward(full,u,rhs)
        if defect['backward']>1e-10:raise RuntimeError('G0 backward residual failed')
        ref=reference.evaluate(w)
        white=float(np.linalg.norm(pack(data-ref.data)))
        row=dict(case=case['id'],fixed_vs_interface_fitted_white=white,
            response_relative=float(np.linalg.norm(data-ref.data)/np.linalg.norm(ref.data)),
            algebraic_backward=defect,ndofs=p.ndofs,reference_ndofs=ref.backend.ndofs,
            target_ports=int(len(ass.public_roi_free_indices)),touched_ports=a.touched_rank,
            polygon_area_error=a.metadata['polygon_partition_area_error'],
            measured_kernel_cpu_seconds=time.process_time()-before,
            numerical_qualification=white<=.1,evidence='new discrete full-wave, not continuum bound')
        rows.append(row)
        np.savez_compressed(out/(case['id']+'.npz'),fixed_data=data,reference_data=ref.data,
            fixed_field=u,reference_field=ref.state.scattered_dofs,
            material=ref.state.epsilon,points=ref.state.material_points,
            target_quad_points=a.quadrature.points,target_quad_weights=a.quadrature.weights)
        (out/'RESULTS.json').write_text(json.dumps(rows,default=_default,indent=2)+'\n')
        print(json.dumps(row,default=_default),flush=True)
        del factor,u,ref,full,rhs,a;reference.clear();gc.collect()
        if white>.1:
            # Single bounded repair requires a recorded new design before run.
            (out/'STOP.json').write_text(json.dumps(dict(reason='FIXED_SPACE_NUMERICS_NOT_QUALIFIED',
                case=case['id'],no_cache_speed_claim=True,repair_status='NOT_RUN'),indent=2)+'\n')
            return 0
    # Port setup is a separately declared action, not hidden in the qualification.
    m=len(ass.public_roi_free_indices);needed=16*(len(p._free_dofs)*m+3*m*m)
    (out/'PORT_PLAN.json').write_text(json.dumps(dict(public_ports=m,port_rhs=m,
        owned_array_lower_estimate_bytes=needed,setup_status='NOT_RUN',
        reason='inspect actual qualified solve costs before optional cache setup',
        whole_G0_status='PARTIAL: finite perturbation/Jacobian/12-class stress pending'),indent=2)+'\n')

if __name__=='__main__':main()
