"""One necessary CPU backend selection: public uniform .25, 1LU/60+60RHS."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import sys,json,resource
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_grid import VectorGrid
from research.research_max_struct_r10.src.vector_cpu_model import CpuDenseVectorModel
from research.research_max_struct_r10.src.vector_physical import PhysicalObjective


def main():
    out=ROOT/'RESULTS/vector/cpu_backend_probe';out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    config=json.loads((ROOT/'configs/vector_contract_proposed.json').read_text())
    model=CpuDenseVectorModel(VectorGrid(12,.6,config['acquisition']))
    chi=np.full(model.N,.25)
    start=time.perf_counter();cput=time.process_time();state=model.state(chi)
    forward=dict(model.last_dense_batch)
    # Public performance-only target, not a private object or inverse dataset.
    # Scaling an already-paid prediction introduces no second forward solve.
    Y=1.1*state['Y'];objective=PhysicalObjective(model,Y,upper=8.)
    start=time.perf_counter();cput=time.process_time();record=objective.evaluate(chi,'public_uniform_CPU_probe',True)
    gradient=dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-start)
    receipt=dict(status='COMPLETED_PUBLIC_PERFORMANCE_ONLY',platform='Mac',
        material='Known public uniform chi=.25 on all1728 cells',target='1.1 times already-paid prediction; no private object/GT, not a scientific reconstruction',
        CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        setup=model.setup,forward=forward,gradient_upgrade_same_LU=gradient,counts=dict(model.counts),
        full_state_shapes=dict(J=list(record['J'].shape),Q=list(record['Q'].shape),material=list(record['m'].shape)),
        Phi=record['Phi'],PG_inf=record['projected_gradient_inf'],
        FFT_constitutive_relative_residual=record['constitutive_relative_residual'],
        maximum_linear_residual=max(r['relative_residual'] for r in model.solve_records),
        additive=True,nested_timers_not_additive=True,new_private_arrays_generated=False,new_inverse_runs=0,
        no_derivative_or_identity_replay=True,backend='Same R1/R9 CM+RR, retained dense LU/all60 source and adjoint RHS; FFT current adapter')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()
