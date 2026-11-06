"""One paid complete CPU/GPU source59 solve, same original60 DDA arrays."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import sys,time,json,traceback,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
import scipy
from scipy import linalg as la
import psutil,torch
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel,four_cubes
from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d'

def main():
    OUT.mkdir(parents=True,exist_ok=True);psutil.Process().cpu_affinity([12,14]);torch.set_num_threads(1);torch.set_num_interop_threads(1)
    wall,cput=time.perf_counter(),time.process_time();result={};cpu=None;gpu=None
    try:
        m=four_cubes(60);active=np.flatnonzero(m!=0);chi=m[active]
        cpu=OriginalVectorModel(60,active_indices=active,generation_only=True,crop_generation_fft=True)
        rhs=np.repeat(chi,3)[:,None]*cpu.inc[:,59:60]
        start,cc=time.perf_counter(),time.process_time();jc=cpu.solve(chi,rhs,label='benchmark_source59_CPU')[:,0]
        result['CPU_solve']=dict(wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cc,
            solve_records=cpu.solve_records,counts=dict(cpu.counts),setup=cpu.setup)
        gpu=CudaOriginalModel(cpu);torch.cuda.reset_peak_memory_stats()
        start,cc=time.perf_counter(),time.process_time();jg=gpu.solve(chi,rhs,label='benchmark_source59_CUDA')[:,0]
        torch.cuda.synchronize();jg_host=jg.cpu().numpy()
        result['CUDA_solve']=dict(wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cc,
            solve_records=gpu.solve_records,counts=dict(gpu.counts),setup=gpu.setup,
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved())
        # Authoritative cross-backend residual of the completed GPU state, not a kernel-only speed test.
        rg=jg_host-np.repeat(chi,3)*cpu.green(jg_host)-rhs[:,0]
        cross_residual=float(la.norm(rg)/la.norm(rhs[:,0]))
        yc=cpu.S@jc;yg=cpu.S@jg_host
        result['qualification']=dict(original_CPU_matrix_residual_of_GPU_state=cross_residual,
            relative_current_difference=float(la.norm(jc-jg_host)/la.norm(jc)),
            relative_prediction_difference=float(la.norm(yc-yg)/la.norm(yc)),
            wall_speed_ratio=result['CPU_solve']['wall_seconds']/result['CUDA_solve']['wall_seconds'],
            scope='Complete RHS numerical backend comparison: SciPyMGS vs residentCGS2/8stepcheck; not pureGPU-only claim')
        result['PASS']=bool(cross_residual<=2e-6 and gpu.solve_records[0]['relative_residual']<=2e-6)
        if result['PASS']:
            folder=OUT/'sealed_truth'/'generation_checkpoints';folder.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(folder/'material.npz',chi=chi,active=active)
            np.savez_compressed(folder/'rhs059.npz',rhs=rhs[:,0],solution=jg_host,Y_clean_column=yg,
                relative_residual=gpu.solve_records[0]['relative_residual'],iterations=gpu.solve_records[0]['iterations'],
                rtol=2e-6,backend='QualifiedTorchCUDA_source59_benchmark')
        else:raise RuntimeError('Whole-RHS backend qualification failed')
    except BaseException as exc:
        result['failure']=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
    finally:
        result.update(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
            process_affinity=psutil.Process().cpu_affinity(),peak_RSS_bytes=psutil.Process().memory_info().peak_wset,
            platform=platform.platform(),numpy_version=np.__version__,scipy_version=scipy.__version__,torch_version=torch.__version__,
            CPU_counts=dict(cpu.counts) if cpu else {},CUDA_counts=dict(gpu.counts) if gpu else {})
        (OUT/'CUDA_WHOLE_RHS_BENCHMARK.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)

if __name__=='__main__':main()
