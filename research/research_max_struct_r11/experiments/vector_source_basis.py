"""One physics-motivated public U_src=orth([Einc,G Einc]), no GT/solve."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import time,json,sys,argparse
from pathlib import Path
import numpy as np
from scipy import linalg as la
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_grid import VectorGrid

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    wall,cpu=time.perf_counter(),time.process_time();out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((REPO/'research/research_max_struct_r10/configs/vector_actual_contract.json').read_text())
    model=VectorGrid(12,.6,cfg['acquisition']);seed=np.column_stack((model.inc,model.green(model.inc)))
    U,s,Vh=la.svd(seed,full_matrices=False,check_finite=False);keep=s>=s[0]*1e-10;U=U[:,keep]
    GU=model.green(U)
    np.savez_compressed(out/'BASIS.npz',right_vectors=U,images=GU,singular_values=s[keep])
    receipt=dict(status='PUBLIC_BACKGROUND_SOURCE_BASIS',rank=U.shape[1],M=model.M,N=model.N,
        seed_shape=list(seed.shape),relative_SVD_threshold=1e-10,seed_singular_values=s.tolist(),
        basis_formula='orth([Einc,G Einc]) with actual same60 electric dipole sources',
        seed_singular_values_are_not_Green_spectrum=True,source_arrays_bytes=U.nbytes+GU.nbytes,
        scope='Captures first two background uniform-contrast propagation fields; no completeness for unknown D_m multiplication',
        allstrong_certified=False,GT_access=False,physical_primal_RHS=0,physical_adjoint_RHS=0,
        setup_thin_SVD=1,counts=dict(model.counts),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:receipt[k] for k in ('rank','source_arrays_bytes','CPU_seconds','wall_seconds')}))

if __name__=='__main__':main()
