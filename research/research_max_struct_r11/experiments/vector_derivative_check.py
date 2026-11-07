"""One newly wired realN/complex3N directional derivative check per variant."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import time,json,sys,argparse
from pathlib import Path
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r11.src.vector_denominator import VectorLowRankDenominator
from research.research_max_struct_r10.src.vector_grid import VectorGrid

def main():
    p=argparse.ArgumentParser();p.add_argument('--basis',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    wall,cpu=time.perf_counter(),time.process_time();out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1)
    cfg=json.loads((REPO/'research/research_max_struct_r10/configs/vector_actual_contract.json').read_text())
    model=VectorGrid(12,.6,cfg['acquisition'])
    with np.load(args.basis,allow_pickle=False) as z:U=torch.from_numpy(z['right_vectors']);GU=torch.from_numpy(z['images'])
    old=REPO.parent/'r10_worktree/research/research_max_struct_r10/RESULTS/vector/v04_GlobalDual/stage3_joint.npz'
    with np.load(old,allow_pickle=False) as z:J=torch.from_numpy(z['J']);Q=torch.from_numpy(z['Q'])
    m=torch.full((1728,),.7,dtype=torch.float64);rng=np.random.default_rng(202611070011)
    direction=torch.from_numpy(rng.standard_normal(1728));direction/=torch.linalg.vector_norm(direction)
    omega=torch.from_numpy(1/np.sum(abs(model.inc)**2,axis=0));beta=1.;rho=1/60;eps=1e-4
    op=VectorLowRankDenominator(U,GU);rows=[]
    for mode in ('full','lagged'):
        ref=op.factor(m,beta) if mode=='lagged' else None
        value,g,_,info=op.value_gradient(J,Q,m,beta,omega,rho,ref)
        plus=op.value_gradient(J,Q,m+eps*direction,beta,omega,rho,ref,False)[0]
        minus=op.value_gradient(J,Q,m-eps*direction,beta,omega,rho,ref,False)[0]
        fd=float(((plus-minus)/(2*eps)).item());analytic=float(torch.dot(g,direction).item())
        relative=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-14)
        rows.append(dict(mode=mode,analytic=analytic,centered_difference=fd,epsilon=eps,relative_error=relative,
            passed=relative<1e-6,lu_info=int(info.item()),all60_sources=True,source_weights_nonuniform=True,
            material_shape=list(m.shape),current_shape=list(J.shape),only_new_derivative_checked=True))
    receipt=dict(status='PASS' if all(r['passed'] for r in rows) else 'FAIL',checks=rows,
        basis=str(Path(args.basis).resolve()),state=str(old),GT_access=False,physical_solve_RHS=0,
        counts=dict(op.cost),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n');print(json.dumps(receipt))

if __name__=='__main__':main()
