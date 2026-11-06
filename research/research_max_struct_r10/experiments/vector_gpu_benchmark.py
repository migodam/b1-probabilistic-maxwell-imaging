"""One public uniform full forward+adjoint, actual resident CUDA workload.

No GT, inverse solve, identity replay, or backend parameter sweep. Same public
uniform chi=.25 and 1.1*paid-prediction performance target as CPU probe.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import time
T=time.perf_counter()
from pathlib import Path
import sys,json,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process,CooperativeGate
process_config=configure_own_windows_process()
spec=importlib.util.spec_from_file_location('job_helper',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object()
import numpy as np
import torch
import psutil
from research.research_max_struct_r10.src.vector_grid import VectorGrid
from research.research_max_struct_r10.src.vector_cuda_model import ResidentVectorModel
from research.research_max_struct_r10.src.vector_physical import PhysicalObjective


def main():
    out=ROOT/'RESULTS/vector/gpu_backend_probe';out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    config=json.loads((ROOT/'configs/vector_actual_contract.json').read_text())
    gate=CooperativeGate(out/'GPU_BOUNDARIES.jsonl')
    model=ResidentVectorModel(VectorGrid(12,.6,config['acquisition']),gate)
    chi=np.full(model.N,.25);mt=torch.as_tensor(chi,device=model.device,dtype=torch.float64)
    t,c=time.perf_counter(),time.process_time()
    J=model.solve(mt,torch.repeat_interleave(mt,3)[:,None]*model.inc,label='public_uniform_GPU_forward')
    prediction=model.S@J;torch.cuda.synchronize()
    forward=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
    Y=(1.1*prediction).cpu().numpy();objective=PhysicalObjective(model,Y,upper=8.)
    phi=float(torch.mean(torch.sum(abs(prediction-model.array(Y))**2,dim=0)/torch.sum(abs(model.array(Y))**2,dim=0)).item())
    key=objective.key(chi)
    objective.records[key]=dict(id='public_state',m=chi.copy(),prediction=prediction.cpu().numpy(),Phi=phi,gradient=None,
        first_owner='public_uniform_GPU_forward',air_analytic=False,uses=[],initial_paid_event=0)
    objective.pending_current[key]=J
    t,c=time.perf_counter(),time.process_time();record=objective.evaluate(chi,'public_uniform_GPU_gradient',True);torch.cuda.synchronize()
    gradient=dict(CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
    p=psutil.Process();cpu=p.cpu_times()
    receipt=dict(status='COMPLETED_PUBLIC_UNIFORM_FULL_VG_CUDA_PROBE',platform='XINAN',CPU_seconds=cpu.user+cpu.system,
        wall_seconds=time.perf_counter()-T,peak_RSS_bytes=p.memory_info().peak_wset,
        GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated(),GPU_peak_reserved_bytes=torch.cuda.max_memory_reserved(),
        forward=forward,gradient_upgrade_paid_J=gradient,counts=dict(model.counts),solve_records=model.solve_records,
        maximum_original_linear_residual=max(r['relative_residual'] for r in model.solve_records),
        constitutive_relative_residual=record['constitutive_relative_residual'],Phi=record['Phi'],PG_inf=record['projected_gradient_inf'],
        actual_full_state_shapes=dict(J=list(record['J'].shape),Q=list(record['Q'].shape),material=list(record['m'].shape)),
        model_setup=model.setup,process=process_config,job=job,
        gate_helper_CPU_seconds=gate.helper_CPU_seconds,gate_helper_invocations=gate.query_count,
        helper_CPU_scope='Separate from main CPU; count once',
        target='Same public chi=.25 and1.1*paidprediction as CPU performance probe, not a private material dataset',
        new_private_arrays=False,new_inverse_runs=0,new_identity_or_derivative_tests=0,
        algorithm='Frozen complex128 packed16 sourcewise qualified GMRES2e-6; no LU')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','CPU_seconds','wall_seconds','forward','gradient_upgrade_paid_J','GPU_peak_reserved_bytes','maximum_original_linear_residual','constitutive_relative_residual','counts')}),flush=True)


if __name__=='__main__':main()
