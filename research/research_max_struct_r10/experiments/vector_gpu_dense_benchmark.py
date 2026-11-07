"""One public uniform 5184-dim GPU retained-LU/all60+60 RHS assessment.

Uses unchanged R1 dyadic matrix and mature Torch LU/solve primitives catalogued
in A9/device.py. No new physics, factorization kernel, GT, parameter sweep, or
change to the running FFT cohort. Its existing external-PID gate briefly yields
at a safe boundary and keeps all numerical state/budgets intact.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import time
T=time.perf_counter()
from pathlib import Path
import sys,json,importlib.util
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.src.vector_resources import configure_own_windows_process,gpu_process_snapshot
process_config=configure_own_windows_process()
spec=importlib.util.spec_from_file_location('own_job_helper',REPO/'research/research_max_struct_r9/delegated/3d_extension/windows_resolved.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);job=module.job_object()
import psutil
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_grid import VectorGrid
from research.research_max_struct_r10.src.vector_cpu_model import CpuDenseVectorModel


def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    out=ROOT/'RESULTS/vector/gpu_dense_backend_probe';out.mkdir(parents=True,exist_ok=False)
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    status='FAILED';failure=None;counts=dict(LU_attempted=0,LU_completed=0,primal_RHS=0,adjoint_RHS=0)
    query_cpu=0.;cooperator=None;timings={};phi=pg=pr=ar=None
    try:
        active=ROOT/'delegated/vector/windows_v08_gpu_queue/ACTIVE.json'
        if active.exists():
            info=json.loads(active.read_text())
            try:
                other=psutil.Process(info['pid']);cmd=' '.join(other.cmdline())
                if 'vector_run.py' in cmd and 'MAX_STRUCT_R10_VECTOR_20261107' in cmd and '--backend cuda' in cmd:
                    cooperator=other.pid
            except psutil.NoSuchProcess:pass
        before=gpu_process_snapshot();query_cpu+=before['query_helper_CPU_seconds'] or 0.
        unrelated=[r for r in before['processes'] if r['pid'] not in (os.getpid(),cooperator)]
        if unrelated:raise RuntimeError('Unrelated active GPU process, no probe launched: '+str(unrelated))
        # Tiny own context makes the existing cohort gate see this authorized
        # probe, then it yields without clearing history or changing any step.
        token=torch.empty(1,device='cuda',dtype=torch.float64);torch.cuda.synchronize()
        handshake=time.perf_counter();observed_pause=False
        while cooperator is not None:
            log=ROOT/'RESULTS/vector/case_v08_gpu_process/GPU_BOUNDARIES.jsonl'
            if log.exists():
                rows=[json.loads(x) for x in log.read_text().splitlines()[-20:]]
                if any(r.get('own_pid')==cooperator and r.get('action')=='WAIT_EXTERNAL_GPU' and any(p['pid']==os.getpid() for p in r.get('external',[])) for r in rows):
                    observed_pause=True;break
            if not psutil.pid_exists(cooperator):cooperator=None;break
            if time.perf_counter()-handshake>60:raise RuntimeError('Safe existing-gate handoff not observed; no large probe launched')
            time.sleep(.1)
        dump(out/'COOPERATIVE_HANDOFF.json',dict(cooperator_pid=cooperator,probe_pid=os.getpid(),
            existing_gate_yield_observed=observed_pause,wait_seconds=time.perf_counter()-handshake,
            no_numeric_state_reset=True,no_algorithm_budget_changed=True,other_processes_modified=False))
        free,total=torch.cuda.mem_get_info()
        if free<3_000_000_000:raise RuntimeError('Insufficient current GPU headroom for bounded dense probe')
        config=json.loads((ROOT/'configs/vector_actual_contract.json').read_text())
        st,ct=time.perf_counter(),time.process_time();cpu=CpuDenseVectorModel(VectorGrid(12,.6,config['acquisition']))
        G=torch.as_tensor(cpu.G,device='cuda',dtype=torch.complex128)
        E=torch.as_tensor(cpu.cpu.inc,device='cuda',dtype=torch.complex128)
        S=torch.as_tensor(cpu.cpu.S,device='cuda',dtype=torch.complex128)
        torch.cuda.synchronize();timings['same_R1_matrix_and_upload']=dict(CPU_seconds=time.process_time()-ct,wall_seconds=time.perf_counter()-st,G_bytes=cpu.G.nbytes)
        def timed(label,fn):
            start,end=torch.cuda.Event(enable_timing=True),torch.cuda.Event(enable_timing=True)
            wall,cput=time.perf_counter(),time.process_time();start.record();value=fn();end.record();end.synchronize()
            timings[label]=dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,CUDA_event_stream_ms=start.elapsed_time(end))
            return value
        core_start=time.perf_counter();core_cpu=time.process_time()
        K=timed('matrix_assembly',lambda:-.25*G);K.diagonal().add_(1.)
        counts['LU_attempted']=1
        factor=timed('retained_LU',lambda:torch.linalg.lu_factor_ex(K,check_errors=False))
        info=int(factor.info.item())
        if info:raise RuntimeError('Mature GPU LU returned info='+str(info))
        counts['LU_completed']=1;del K
        J=timed('all60_primal_same_LU',lambda:torch.linalg.lu_solve(factor.LU,factor.pivots,.25*E));counts['primal_RHS']=60
        Q=E+G@J;prediction=S@J;Y=1.1*prediction;weights=1/torch.sum(abs(Y)**2,dim=0)
        rhs=S.conj().T@((prediction-Y)*weights)
        P=timed('all60_adjoint_same_LU',lambda:torch.linalg.lu_solve(factor.LU,factor.pivots,rhs,adjoint=True));counts['adjoint_RHS']=60
        gradient=(2/60)*torch.real(torch.sum(P.conj()*Q,dim=1)).reshape(1728,3).sum(1)
        phi=float(torch.mean(torch.sum(abs(prediction-Y)**2,dim=0)*weights).item())
        pg=float(torch.max(abs(.25-torch.clamp(.25-gradient,0,8))).item())
        pr=float((torch.linalg.vector_norm(J-.25*Q)/torch.linalg.vector_norm(.25*E)).item())
        ar=float((torch.linalg.vector_norm(P-G.conj().T@(.25*P)-rhs)/torch.linalg.vector_norm(rhs)).item())
        torch.cuda.synchronize();timings['full_warm_value_gradient_and_residual']=dict(CPU_seconds=time.process_time()-core_cpu,wall_seconds=time.perf_counter()-core_start)
        if max(pr,ar)>2e-6:raise RuntimeError('Original frozen residual qualification failed')
        status='COMPLETED_ONE_PUBLIC_UNIFORM_GPU_DENSE_ASSESSMENT'
    except BaseException as exc:
        import traceback
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
    finally:
        p=psutil.Process();c=p.cpu_times()
        receipt=dict(status=status,platform='XINAN',CPU_seconds=c.user+c.system,wall_seconds=time.perf_counter()-T,
            peak_RSS_bytes=p.memory_info().peak_wset,timings=timings,counts=counts,Phi=phi,PG_inf=pg,
            original_primal_relative_residual=pr,original_adjoint_relative_residual=ar,
            GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_initialized() else None,
            GPU_peak_reserved_bytes=torch.cuda.max_memory_reserved() if torch.cuda.is_initialized() else None,
            query_child_CPU_seconds=query_cpu,query_CPU_not_in_main=True,process=process_config,job=job,
            failure=failure,material='Public uniform chi=.25 on all1728 cells',target='1.1*same paid prediction, identical performance target construction as CPU/FFT probes',
            physical_implementation='Unchanged R1/R9 CM+RR dyadic G and E/S, mature Torch LU/lu_solve default backend',
            linalg_library=str(torch.backends.cuda.preferred_linalg_library()),torch_version=torch.__version__,
            no_GT_read=True,no_new_inverse_run=True,no_parameter_sweep=True,no_new_factorization_kernel=True,
            CUDA_event_scope='Default stream intervals, not active kernel time; event intervals not host CPU charges',
            memory_scope='Dense-probe measured peak; does not alter production FFT1.5GB gate policy or current cohort backend')
        dump(out/'RECEIPT.json',receipt);print(json.dumps(receipt),flush=True)
    return 0 if status.startswith('COMPLETED') else 1


if __name__=='__main__':raise SystemExit(main())
