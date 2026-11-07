"""Only the newly contracted3D z-gradient, one actual-source direction check."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import time,json,sys
from pathlib import Path
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r11.src.vector_denominator import VectorLowRankDenominator
wall,cpu=time.perf_counter(),time.process_time();torch.set_num_threads(1)
out=ROOT/'delegated/vector/z_derivative_check';out.mkdir(parents=True,exist_ok=False)
with np.load(ROOT/'delegated/vector/source_basis/BASIS.npz',allow_pickle=False) as f:U=torch.from_numpy(f['right_vectors']);GU=torch.from_numpy(f['images'])
with np.load(ROOT/'delegated/vector/PUBLIC_V04_FIXED_STATE.npz',allow_pickle=False) as f:J=torch.from_numpy(f['J']);Q=torch.from_numpy(f['Q']);omega=torch.from_numpy(f['weights'])
op=VectorLowRankDenominator(U,GU);beta=3.;rho=1/60;m=torch.full((1728,),.7,dtype=torch.float64);z=m/(1+beta*m);V=Q+beta*J
reference=op.factor_contracted(m,beta);rng=np.random.default_rng(202611070012);direction=torch.from_numpy(rng.standard_normal(1728));direction/=torch.linalg.vector_norm(direction);eps=1e-5
value,g,_,info=op.contracted_value_gradient(J,V,z,beta,omega,rho,reference)
plus=op.contracted_value_gradient(J,V,z+eps*direction,beta,omega,rho,reference,False)[0]
minus=op.contracted_value_gradient(J,V,z-eps*direction,beta,omega,rho,reference,False)[0]
fd=float(((plus-minus)/(2*eps)).item());analytic=float(torch.dot(g,direction).item());error=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-14)
receipt=dict(status='PASS' if error<1e-6 and int(info.item())==0 else 'FAIL',mode='contractedLagged',
    analytic=analytic,centered_difference=fd,epsilon=eps,relative_error=error,rank=120,current_shape=list(J.shape),
    material_coordinates=1728,source_weights='Nonuniform all60 ws',beta=beta,rho=rho,
    GT_access=False,new_physical_RHS=0,counts=dict(op.cost),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
