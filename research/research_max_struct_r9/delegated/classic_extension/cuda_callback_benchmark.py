"""One sealed public4096-material callback: CPU FFT versus existing CUDA LU.

No GT, optimization, repeated callback/warmup, or CPU16 modifications.
Run only after the parent grants exclusive GPU and host logical P12.
"""
import os
import time
START,CPU_START=time.perf_counter(),time.process_time()
for key in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']='0'
import psutil
psutil.Process().cpu_affinity([12])
from pathlib import Path
import json
import sys
import traceback
REPO=Path(__file__).resolve().parents[4]
ROOT=REPO/'research/research_max_struct_r9'
sys.path.insert(0,str(REPO))
OUT=ROOT/'RESULTS/classic_extension_cuda_callback'
OUT.mkdir(parents=True,exist_ok=True)


def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False))


def main():
    if (OUT/'RECEIPT.json').exists() or (OUT/'FAILURE.json').exists():
        raise FileExistsError('Preserve previous callback benchmark')
    phase='scientific_imports';record={}
    try:
        import numpy as np
        import scipy
        from research.research_max_struct_r8.src.paper_physics import acquisition
        from research.research_max_struct_r9.src.classic_extension import ClassicModel
        Y=np.load(ROOT/'data_public/classic_extension/austria_eps3p5_400MHz/PUBLIC_DATA.npz')['Y']
        m=np.load(ROOT/'data_public/classic_extension/cuda_candidate_material.npz')['m']
        if m.shape!=(4096,) or Y.shape!=(40,20):raise ValueError('Unexpected sealed input shapes')
        record.update(input_material='Mac sealed epsilon3.5 single CIE CHECKPOINT_080; copied public prediction',
            input_observations='unmodified original150² R8 Y',GT_loaded=False,
            physical_grid=[64,64],physical_unknowns=4096,sources=20,receivers=40,
            material_box=[0,20],material_min=float(m.min()),material_max=float(m.max()),
            host_logical_Pcore=psutil.Process().cpu_affinity(),dtype='complex128/float64',
            numpy=np.__version__,scipy=scipy.__version__,python=sys.version,
            science='one fixed-point callback only; CPU16 unchanged; no material recovery assessment')
        phase='CPU_setup';t,c=time.perf_counter(),time.process_time()
        cpu_public=acquisition(64,400e6,dense=False);model=ClassicModel(cpu_public,prior=True)
        record['CPU_setup']=dict(wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c)
        phase='CPU_one_full_value_gradient';t,c=time.perf_counter(),time.process_time()
        value_cpu,gradient_cpu=model.data_value_gradient(m,Y)
        record['CPU_callback']=dict(wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c,
            cost=dict(model.cost),value=float(value_cpu),gradient_norm=float(np.linalg.norm(gradient_cpu)))
        st=model.state(m)
        record['CPU_diagnostics']=dict(physical_relative_data_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),
            state_relative_residual=float(np.linalg.norm(st['j']-m[:,None]*st['Q'])/max(np.linalg.norm(m[:,None]*model.E),1e-300)))
        np.savez_compressed(OUT/'CPU_VALUE_GRADIENT.npz',m=m,value=value_cpu,gradient=gradient_cpu)
        dump(OUT/'PROGRESS.json',dict(phase=phase,**record))
        phase='GPU_cold_setup';gpu_setup_t,gpu_setup_c=time.perf_counter(),time.process_time()
        t,c=time.perf_counter(),time.process_time()
        import torch
        torch.set_num_threads(1);torch.set_num_interop_threads(1)
        torch.backends.cuda.preferred_linalg_library('cusolver')
        torch.cuda.set_device(0);torch.cuda.reset_peak_memory_stats()
        record['Torch_context_library_setup']=dict(wall_seconds=time.perf_counter()-t,
            CPU_seconds=time.process_time()-c,torch=torch.__version__,CUDA=torch.version.cuda,
            selected_linalg='cusolver',device=torch.cuda.get_device_name(0))
        from research.research_z_alg_cuda_pipeline.src.device import SharedOperators,ObjectiveLane
        t,c=time.perf_counter(),time.process_time()
        gpu_public=acquisition(64,400e6,dense=True)
        record['GPU_dense_operator_assembly']=dict(wall_seconds=time.perf_counter()-t,
            CPU_seconds=time.process_time()-c,Green_bytes=gpu_public['G'].nbytes,
            convention='exact same equal-area disk formula including diagonal; not default576')
        t,c=time.perf_counter(),time.process_time()
        shared=SharedOperators(gpu_public)
        lane=ObjectiveLane(shared,{'Y':Y},stage=64,objective='data',
            config=dict(smooth_weight=1e-4,tv_weight=1e-4,tv_epsilon=.01))
        torch.cuda.synchronize()
        record['GPU_resident_setup']=dict(wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c,
            counts=shared.setup_counts,extra_setup='unmodified SharedOperators constructs reference.1LU and unused fullB/DB; paid, not hidden')
        record['GPU_cold_setup_total']=dict(wall_seconds=time.perf_counter()-gpu_setup_t,
            CPU_seconds=time.process_time()-gpu_setup_c,
            note='includes Torch/context, denseG assembly, upload, unused basis and referenceLU; components overlap this total')
        phase='GPU_one_full_value_gradient';dump(OUT/'PROGRESS.json',dict(phase=phase,**record))
        event_start=torch.cuda.Event(enable_timing=True);event_end=torch.cuda.Event(enable_timing=True)
        event_start.record(lane.stream)
        t,c=time.perf_counter(),time.process_time()
        value_gpu,gradient_gpu=lane.evaluate(m)
        callback_wall=time.perf_counter()-t;callback_cpu=time.process_time()-c
        event_end.record(lane.stream);event_end.synchronize()
        record['GPU_callback']=dict(wall_seconds=callback_wall,CPU_seconds=callback_cpu,
            CUDA_event_milliseconds=event_start.elapsed_time(event_end),value=float(value_gpu),
            gradient_norm=float(np.linalg.norm(gradient_gpu)),counts=dict(lane.counts),
            return_semantics='one whole callback incl one host result; no repeated callback or warmup')
        phase='cached_state_diagnostics';t,c=time.perf_counter(),time.process_time()
        endpoint=lane.close(m) # Cached same state: no extra LU or physical RHS.
        record['GPU_cached_state_diagnostics']=dict(wall_seconds=time.perf_counter()-t,
            CPU_seconds=time.process_time()-c,physical_relative_data_residual=endpoint['relative_residual'],
            state_relative_residual=endpoint['state_relative_residual'],LU_info=endpoint['lu_info'],
            counts=endpoint['counts'],extra_LU=0,extra_physical_RHS=0)
        record['comparison']=dict(value_absolute_difference=abs(float(value_gpu-value_cpu)),
            value_relative_difference=abs(float(value_gpu-value_cpu))/max(abs(float(value_cpu)),1e-300),
            gradient_relative_L2=float(np.linalg.norm(gradient_gpu-gradient_cpu)/max(np.linalg.norm(gradient_cpu),1e-300)),
            gradient_absolute_Linf=float(np.max(abs(gradient_gpu-gradient_cpu))),
            gradient_relative_Linf=float(np.max(abs(gradient_gpu-gradient_cpu))/max(np.max(abs(gradient_cpu)),1e-300)),
            CPU_callback_wall_over_GPU_callback_wall=record['CPU_callback']['wall_seconds']/callback_wall,
            interpretation='one public state; does not establish whole-optimizer throughput or material recovery')
        record['GPU_peak_allocated_bytes']=torch.cuda.max_memory_allocated()
        record['GPU_peak_reserved_bytes']=torch.cuda.max_memory_reserved()
        np.savez_compressed(OUT/'VALUE_GRADIENT.npz',m=m,value_CPU=value_cpu,value_GPU=value_gpu,
                            gradient_CPU=gradient_cpu,gradient_GPU=gradient_gpu)
        record.update(status='COMPLETED',valid_CPU_callbacks=1,valid_GPU_callbacks=1,
            full_total_wall_seconds=time.perf_counter()-START,
            full_total_internal_process_CPU_seconds=time.process_time()-CPU_START,
            total_setup_LU=shared.setup_counts['LU'],GPU_callback_LU=lane.counts['LU'],
            GPU_total_source_RHS=shared.setup_counts['source_RHS']+lane.counts['source_RHS'],
            GPU_total_adjoint_RHS=lane.counts['adjoint_RHS'],
            costs_scope='all recorded actual setup/callback/diagnostic work; startup/output tail not externally harvested')
        dump(OUT/'RECEIPT.json',record)
        print(json.dumps(dict(status=record['status'],comparison=record['comparison'],
            CPU_callback=record['CPU_callback'],GPU_callback=record['GPU_callback'],
            GPU_cold_setup_total=record['GPU_cold_setup_total'],GPU_diagnostics=record['GPU_cached_state_diagnostics'],
            GPU_peak_allocated_bytes=record['GPU_peak_allocated_bytes'])),flush=True)
    except BaseException as exc:
        dump(OUT/'FAILURE.json',dict(status='FAILED',phase=phase,error=repr(exc),
            wall_seconds=time.perf_counter()-START,internal_process_CPU_seconds=time.process_time()-CPU_START,
            partial_record=record,GT_loaded=False,traceback=traceback.format_exc()))
        raise


if __name__=='__main__':main()
