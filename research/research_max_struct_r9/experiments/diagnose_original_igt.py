"""Predeclared source59 only: original3lambda same target30/60 IGT diagnosis."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import argparse,json,sys,time,traceback,platform
from pathlib import Path
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
import scipy
from scipy import linalg as la
import psutil
from research.research_max_struct_r9.src.vector_original_igt import CellIntegratedOriginalModel
from research.research_max_struct_r9.src.vector_original_extension import four_cubes
BASE=Path(__file__).resolve().parents[1];PARENT=BASE/'RESULTS'/'original_vector3d';OUT=PARENT/'igt_source59'

def dump(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def prepare():
    for side in (30,60):
        dest=OUT/str(side)
        if (dest/'PREPARED.json').exists():raise RuntimeError('Refusing repeated CPU assembly')
        dest.mkdir(parents=True,exist_ok=True);wall,cput=time.perf_counter(),time.process_time();model=None;failure=None
        try:
            chi=four_cubes(side);active=np.flatnonzero(chi!=0) if side==60 else None
            model=CellIntegratedOriginalModel(side,source_indices=(59,),active_indices=active,
                generation_only=side==60,crop_generation_fft=side==60)
            material=chi[model.active]
            for key,value in dict(kernel=model.kernel,inc=model.inc,S=model.S,chi=material,active=model.active,
                green_rows=model.green_rows,local_offdiagonal_Ablock=model.local_offdiagonal_block()).items():np.save(dest/(key+'.npy'),value)
            dump(dest/'PREPARED.json',dict(side=side,N=model.N,M=model.M,P=1,full_N=model.full_N,h=model.h,k=model.k,
                s=[model.s.real,model.s.imag],source_index_zero_based=59,convolution_shape=model.convolution_shape,pad=model.pad,
                convolution_N=model.convolution_N,block_rhs=4,rtol=2e-6,restart=60,cycles=20,generation_only=side==60,
                source_position=model.sources[59].tolist(),source_polarization=model.polarizations[59].tolist(),setup=model.setup))
            print(json.dumps(dict(status='CPU_PREPARED',side=side,setup=model.setup)),flush=True)
        except BaseException as exc:
            failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
        finally:
            dump(dest/'PREPARE_COST.json',dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
                counts=dict(model.counts) if model else{},physics_RHS=0,failure=failure,
                peak_RSS_bytes=psutil.Process().memory_info().peak_wset))
        del model

def solve():
    # Parent separately authorizes exclusive GPU release before this mode starts.
    import torch
    from research.research_max_struct_r9.src.vector_original_cuda import CudaOriginalModel
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    for side in (30,60):
        dest=OUT/str(side);meta=json.loads((dest/'PREPARED.json').read_text())
        if (dest/'RESULT.json').exists():raise RuntimeError('Only one predeclared RHS per side')
        wall,cput=time.perf_counter(),time.process_time();model=None;failure=None
        try:
            cpu=SimpleNamespace(**{k:meta[k] for k in ['side','N','M','P','full_N','block_rhs','rtol','restart','cycles','convolution_N','generation_only']})
            cpu.s=complex(*meta['s']);cpu.convolution_shape=tuple(meta['convolution_shape']);cpu.pad=tuple(meta['pad']);cpu.setup=meta['setup']
            for k in ('kernel','inc','S','active','green_rows'):setattr(cpu,k,np.load(dest/(k+'.npy')))
            chi=np.load(dest/'chi.npy');model=CudaOriginalModel(cpu,fft_layout='packed',block_rhs=16);model.P=1
            rhs=np.repeat(chi,3)[:,None]*cpu.inc
            def progress(r):
                with (dest/'progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(dict(side=side,**r)),flush=True)
            J=model.solve(chi,rhs,label='IGT_original_source59_side'+str(side),progress=progress)
            y=(model.S@J)[:,0].cpu().numpy();j=J[:,0].cpu().numpy()
            # Independent original-CPU-matrix residual for this completed state.
            from scipy import fft
            full=np.zeros((3*cpu.convolution_N,1),complex);full[cpu.green_rows]=j[:,None]
            cube=full.reshape(cpu.convolution_shape+(3,1));hat=fft.fftn(cube,s=cpu.pad,axes=(0,1,2))
            prod=np.einsum('...ij,...jb->...ib',cpu.kernel,hat,optimize=True)
            gj=fft.ifftn(prod,axes=(0,1,2))[tuple(slice(0,d) for d in cpu.convolution_shape)].reshape(3*cpu.convolution_N,1)[cpu.green_rows,0]
            residual=float(la.norm(j-np.repeat(chi,3)*gj-rhs[:,0])/la.norm(rhs[:,0]));ycpu=cpu.S@j
            result=dict(side=side,source_index_zero_based=59,original_CPU_matrix_relative_residual=residual,
                CPU_GPU_prediction_relative_difference=float(la.norm(y-ycpu)/max(la.norm(ycpu),1e-300)),
                PASS=bool(residual<=2e-6),receiver_vector_rows=180,no_noise=True,no_inverse=True,
                physical_model='fixed8/3/16source-cellIGT;pulse/collocation; incidentcenterfield',solve_record=model.solve_records[0])
            np.savez_compressed(dest/'prediction.npz',prediction=y,current=j);dump(dest/'RESULT.json',result)
            if not result['PASS']:raise RuntimeError('Original CPU IGT residual failed')
        except BaseException as exc:
            failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
        finally:
            if model:torch.cuda.synchronize()
            dump(dest/'SOLVE_COST.json',dict(CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
                counts=dict(model.counts) if model else{},solve_records=model.solve_records if model else[],
                CPU_original_matrix_residual_verification_Green_RHS=1 if failure is None else None,
                peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if model else 0,
                peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if model else 0,
                peak_RSS_bytes=psutil.Process().memory_info().peak_wset,failure=failure))
        del model,cpu,J
        torch.cuda.empty_cache()
    a=np.load(OUT/'30/prediction.npz')['prediction'];b=np.load(OUT/'60/prediction.npz')['prediction']
    cm30=np.load(PARENT/'sealed_truth/coarse_truth_response.npz')['prediction'][:,59]
    cm60=np.load(PARENT/'sealed_truth/truth.npz')['Y_clean'][:,59]
    result=dict(scope='PREDECLARED_SOURCE59_ONLY; two clean RHS; same original3lambda aligned-cube target; no imaging/GT pixel scoring',
        IGT_relative_coarse_fine=float(la.norm(a-b)/la.norm(b)),CM_relative_coarse_fine=float(la.norm(cm30-cm60)/la.norm(cm60)),
        IGT_absolute_gap_norm=float(la.norm(a-b)),CM_absolute_gap_norm=float(la.norm(cm30-cm60)),
        IGT30_prediction_norm=float(la.norm(a)),IGT60_prediction_norm=float(la.norm(b)),CM30_prediction_norm=float(la.norm(cm30)),CM60_prediction_norm=float(la.norm(cm60)),
        IGT_vs_CM30_relative=float(la.norm(a-cm30)/la.norm(cm30)),IGT_vs_CM60_relative=float(la.norm(b-cm60)/la.norm(cm60)),
        fixed_quadrature=dict(near=8,far_external=3,self=16,near_radius=2),no_gain_rescaling=True,no_new_60_source_package=True,no_inversion=True)
    dump(OUT/'PAIRED_DIAGNOSTIC.json',result);print(json.dumps(result),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','solve']);args=ap.parse_args()
    psutil.Process().cpu_affinity([12,14]);OUT.mkdir(parents=True,exist_ok=True)
    dump(OUT/'DECLARED_CONTRACT.json',dict(scope='source59_only original3lambda same alignedfourcubes30/60 IGT propagation diagnostic',
        source_index_zero_based=59,sides=[30,60],quadrature=dict(near=8,far_external=3,self=16,near_radius=2),
        inverse_unknown_cells_unchanged=27000,only_fine_generation_known_support_crop=True,noise=None,
        GMRES=dict(restart=60,max_cycles=20,original_relative_tolerance=2e-6),no_inverse=True,no_new_60_source_data=True,
        prior_CM_frozen=True,incident='center_collocation',external_S='source_cell_integral_no_extra_h3',
        GPU_authorization='Parent explicitly released classic GPU benchmark before solve; no automatic later 60source/5arm execution'))
    wall,cput=time.perf_counter(),time.process_time();failure=None
    try:prepare() if args.mode=='prepare' else solve()
    except BaseException as exc:failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
    finally:
        ct=psutil.Process().cpu_times()
        dump(OUT/(args.mode+'_PROCESS_COST.json'),dict(CPU_seconds=ct.user+ct.system,wall_seconds=time.time()-psutil.Process().create_time(),
            main_body_CPU_seconds=time.process_time()-cput,main_body_wall_seconds=time.perf_counter()-wall,
            process_affinity=psutil.Process().cpu_affinity(),numpy_version=np.__version__,scipy_version=scipy.__version__,platform=platform.platform(),failure=failure))

if __name__=='__main__':main()
