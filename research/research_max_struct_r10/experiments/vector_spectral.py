"""One background-operator top32 numerical capture; no material/GT/solve."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
from scipy import linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from research.research_max_struct_r10.src.vector_grid import VectorGrid


def main():
    output=ROOT/'RESULTS/vector/operator_spectral';output.mkdir(parents=True,exist_ok=False)
    public=ROOT/'data_public/vector/operator_spectral';public.mkdir(parents=True,exist_ok=False)
    config=json.loads((ROOT/'configs/vector_contract_proposed.json').read_text())
    model=VectorGrid(12,.6,config['acquisition']);M=model.M
    rng=np.random.default_rng(202611070032);v0=rng.standard_normal(M)+1j*rng.standard_normal(M);v0/=la.norm(v0)
    def gram(v):return model.green(model.green(v),True)
    op=LinearOperator((M,M),matvec=gram,dtype=np.complex128)
    status='NUMERICAL_TOP32_CAPTURE';error=None
    try:
        values,V=eigsh(op,k=32,which='LM',v0=v0,ncv=80,maxiter=300,tol=1e-7)
    except ArpackNoConvergence as exc:
        status='PARTIAL_ARPACK_NOT_CONVERGED';values,V=exc.eigenvalues,exc.eigenvectors;error=str(exc)
    if V is None or V.shape[1]==0:
        summary=dict(status=status,error=error,captured_dimension=0,orthogonal_remainder_dimension=M)
    else:
        # Small Rayleigh-Ritz cleanup within the returned space, not another run.
        V,_=la.qr(V,mode='economic');GV=model.green(V)
        values,W=la.eigh(GV.conj().T@GV);order=np.argsort(values)[::-1];values=values[order];V=V@W[:,order];GV=GV@W[:,order]
        residual=model.green(GV,True)-V*values[None,:]
        residual_norm=la.norm(residual,axis=0);relative=residual_norm/np.maximum(abs(values),1e-300)
        orth=float(la.norm(V.conj().T@V-np.eye(V.shape[1])))
        captured=V.shape[1];singular=np.sqrt(np.maximum(values,0))
        numerical_ok=bool(orth<1e-8 and relative.max()<1e-5 and singular.min()>1e-10*singular.max())
        np.savez_compressed(public/'GIN_TOP32_CAPTURE.npz',right_vectors=V,images=GV,singular_values=singular,
            gram_pair_residual_norms=residual_norm)
        summary=dict(status=status,error=error,captured_dimension=captured,orthogonal_remainder_dimension=M-captured,
            domain='complex3N_vector_current',M=M,material_scalar_dimension=model.N,
            singular_values=singular.tolist(),relative_gram_pair_residuals=relative.tolist(),
            orthogonality_error=orth,numerically_supported_returned_range=numerical_ok,
            classification='Computed top32 capture requested by largest Ritz values; all remaining energy is orthogonal remainder/uncaptured',
            global_spectral_order_certificate=False,full_spectrum_computed=False,
            remainder_is_weak_or_kernel=False,threshold_010_allstrong_certified=False,
            scalar_material_direct_projection_permitted=False,
            basis_path=str((public/'GIN_TOP32_CAPTURE.npz').relative_to(ROOT)))
    summary.update(CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,counts=dict(model.counts),
        physical_primal_RHS=0,physical_adjoint_RHS=0,material_GT_access=False,private_data_generation=False,
        platform='Mac',scope='Whole process imports/public geometry/one bounded Lanczos/Ritz/residuals/exports; no physical solve')
    (output/'RECEIPT.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:summary.get(k) for k in ('status','captured_dimension','CPU_seconds','wall_seconds','numerically_supported_returned_range')}))


if __name__=='__main__':main()
