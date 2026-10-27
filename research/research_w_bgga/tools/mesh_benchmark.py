"""No PDE: measure actual public mesh/target assembly size and cost only."""
from pathlib import Path
import json,time
import numpy as np
from skfem import Basis,ElementTriP3
from bgga.fixed_mesh import public_fixed_mesh,TargetAssembler
from wslam.fem_backend import line_green
from bgga.accounting import _default

ROOT=Path(__file__).resolve().parents[1]
def main():
    out=ROOT/'results'/'preflight';out.mkdir(parents=True,exist_ok=True)
    t=time.process_time();mesh=public_fixed_mesh();mesh_cpu=time.process_time()-t
    t=time.process_time();basis=Basis(mesh,ElementTriP3(),intorder=8)
    free=np.setdiff1d(np.arange(basis.N),basis.get_dofs().all());ass=TargetAssembler(basis,free)
    basis_cpu=time.process_time()-t
    tx=np.array([[-1.5,0],[1.5,0]])
    t=time.process_time();a=ass.assemble((.2,2.),.12,2*np.pi,3+.02j,lambda pts:line_green(pts,tx,2*np.pi))
    asm_cpu=time.process_time()-t
    m=len(ass.public_roi_free_indices);n=len(free)
    rec=dict(evidence='ZERO_PDE_MESH_AND_COEFFICIENT_BENCHMARK_ONLY',new_full_fe=0,new_rhs=0,
        ndofs=basis.N,free_dofs=n,triangles=mesh.nelements,public_ports=m,
        example_touched_ports=a.touched_rank,example_V_structural_nnz=a.Vtarget.nnz,
        quadrature_points=a.metadata['quadrature_points'],target_area=a.target_area,
        polygon_area_error=a.metadata['polygon_partition_area_error'],
        mesh_cpu_seconds=mesh_cpu,basis_index_cpu_seconds=basis_cpu,target_assembly_cpu_seconds=asm_cpu,
        full_green_image_bytes=16*n*m,core_bytes=16*m*m,
        cache_owned_array_lower_estimate_bytes=16*(n*m+3*m*m),
        port_setup_rhs=m,direct_two_source_rhs=2,
        historical_only_cpu_per_full_state=9288.02304/2950,
        new_full_solve_cpu='NOT_MEASURED_RESOURCE_GRANT_REQUIRED',
        numerical_forward_accuracy='NOT_TESTED',
        no_speed_claim=True)
    (out/'MESH_BENCHMARK.json').write_text(json.dumps(rec,indent=2,default=_default)+'\n')
    np.savez_compressed(out/'FIXED_MESH.npz',points=mesh.p,triangles=mesh.t,
        public_roi_free_indices=ass.public_roi_free_indices,target_quad_points=a.quadrature.points,
        target_quad_weights=a.quadrature.weights)
    print(json.dumps(rec,default=_default))

if __name__=='__main__':main()
