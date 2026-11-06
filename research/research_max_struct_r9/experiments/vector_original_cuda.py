"""Finite foreground CUDA jobs; one shared resident operator for all five arms."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import argparse,json,sys,time,traceback,platform,gc
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[3]))
import numpy as np
import scipy
from scipy import linalg as la
import psutil,torch
from research.research_max_struct_r9.src.vector_original_extension import OriginalVectorModel,four_cubes,LocalArcReadout,material_readout
from research.research_max_struct_r9.src.vector_original_cuda import (CudaOriginalModel,CudaReceiverSubspace,CudaImplicitFourier,
    CudaLocalArcReadout,cuda_material_readout,cuda_auxiliary_solve)
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'RESULTS'/'original_vector3d';CFG=BASE/'configs'/'original_vector3d.json'

class ControlledEngineeringStop(SystemExit):pass

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)

def delta_counts(model,before):return {k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}

def wiring(model,receiver,Y):
    """One actual-public-operator GPU-specific integrated check; no GT/solves."""
    wall,cput=time.perf_counter(),time.process_time();before=dict(model.counts)
    B=CudaImplicitFourier(model,receiver,6);rng=np.random.default_rng(2026110620)
    def rand(shape):return model.array(rng.normal(size=shape)+1j*rng.normal(size=shape))
    g=rand((B.K,2));x=rand((model.M,2));bg=B.apply(g)
    dotb=torch.vdot(bg.flatten(),x.flatten());dotba=torch.vdot(g.flatten(),B.adjoint(x).flatten())
    berr=float((abs(dotb-dotba)/torch.clamp(abs(dotb),min=1)).item())
    dotg=torch.vdot(model.green(bg).flatten(),x.flatten());dotga=torch.vdot(bg.flatten(),model.green(x,True).flatten())
    gerr=float((abs(dotg-dotga)/torch.clamp(abs(dotg),min=1)).item())
    J0=receiver.Jd[:,:2];E=model.inc[:,:2];target=Y[:,:2]
    wd=1/torch.sum(abs(target)**2,dim=0);ws=1/torch.sum(abs(E)**2,dim=0)
    m=torch.full((model.N,),.1,device='cuda',dtype=torch.float64);chi=torch.repeat_interleave(m,3)[:,None];p=6/(1+6*chi)
    Q=E+model.green(J0);D=model.S@J0-target;v=p.conj()*p*(J0-chi*Q)*ws
    grad=B.adjoint(model.S.conj().T@(D*wd)+v-model.green(chi*v,True))
    direction=rand(grad.shape);direction=grad+.1*torch.linalg.vector_norm(grad)*direction/torch.linalg.vector_norm(direction)
    direction=direction/torch.linalg.vector_norm(direction)
    def value(gamma):
        j=J0+B.apply(gamma);q=E+model.green(j);r=p*(j-chi*q);d=model.S@j-target
        return float((torch.sum(abs(d)**2*wd)+torch.sum(abs(r)**2*ws)).item())
    step=1e-3;analytic=float((2*torch.vdot(grad.flatten(),direction.flatten()).real).item())
    fd=(value(step*direction)-value(-step*direction))/(2*step)
    jhost=J0.cpu().numpy();qhost=Q.cpu().numpy();whost=ws.cpu().numpy();mhost=m.cpu().numpy()
    closed_gpu=cuda_material_readout(J0,Q,ws,6).cpu().numpy();closed_cpu=material_readout(jhost,qhost,whost,6)
    arc_cpu,info_cpu=LocalArcReadout(model.cpu)(jhost,qhost,mhost,whost,6)
    arc_gpu,info_gpu=CudaLocalArcReadout(model)(J0,Q,m,ws,6);arc_gpu=arc_gpu.cpu().numpy()
    relgrad=abs(fd-analytic)/max(abs(fd),abs(analytic),1e-30)
    arc_metric=abs(info_gpu['fixed_current_loss_after']-info_cpu['fixed_current_loss_after'])/max(abs(info_cpu['fixed_current_loss_after']),1e-30)
    result=dict(scope='ONE GPU-specific public-state wiring check; two actual source columns; no GT and no full solves',
        implicit_adjoint_relative=berr,Green_adjoint_relative=gerr,current_gradient_analytic=analytic,
        current_gradient_finite_difference=fd,current_gradient_relative_error=relgrad,
        closed_material_max_absolute_difference=float(np.max(abs(closed_gpu-closed_cpu))),
        arc_material_max_absolute_difference=float(np.max(abs(arc_gpu-arc_cpu))),arc_metric_relative_difference=arc_metric,
        arc_CPU=info_cpu,arc_GPU=info_gpu,counts=delta_counts(model,before),physics_solves=0,
        CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall)
    result['PASS']=bool(berr<1e-10 and gerr<1e-10 and relgrad<1e-6 and result['closed_material_max_absolute_difference']<1e-9 and arc_metric<1e-8)
    dump(OUT/'CUDA_WIRING.json',result)
    if not result['PASS']:raise RuntimeError('GPU-specific wiring failed')
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['generate','qualify-discretization','invert-five']);args=ap.parse_args()
    psutil.Process().cpu_affinity([12,14]);torch.set_num_threads(1);torch.set_num_interop_threads(1)
    config=json.loads(CFG.read_text());OUT.mkdir(parents=True,exist_ok=True)
    dump(OUT/'ACTIVE_CUDA_JOB.json',dict(pid=os.getpid(),mode=args.mode,affinity=psutil.Process().cpu_affinity(),
        process_started_epoch=psutil.Process().create_time(),control='STOP_CURRENT_JOB at an auxiliary iteration boundary; saves paid material and exact counters'))
    wall,cput=time.perf_counter(),time.process_time();model=None;status='FAILED';failure=None;arm_failures=[]
    try:
        if args.mode=='generate':
            if (OUT/'public'/'observations.npz').exists():raise RuntimeError('Public data already frozen')
            if not json.loads((OUT/'CUDA_WHOLE_RHS_BENCHMARK.json').read_text())['PASS']:raise RuntimeError('Backend not qualified')
            m=four_cubes(60);active=np.flatnonzero(m!=0)
            cpu=OriginalVectorModel(60,active_indices=active,generation_only=True,crop_generation_fft=True)
            model=CudaOriginalModel(cpu,fft_layout=config.get('CUDA_fft_layout','strided'),block_rhs=config.get('CUDA_block_rhs',4));dump(OUT/'GENERATION_SETUP_CUDA.json',model.setup)
            def progress(r):
                with (OUT/'generation_CUDA_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(r),flush=True)
            Yc=model.physical(m[active],label='fine60_CUDA_resume',progress=progress,checkpoint_dir=OUT/'sealed_truth'/'generation_checkpoints')
            sigma=.1*la.norm(Yc)/np.sqrt(Yc.size);rng=np.random.default_rng(config['noise_seed'])
            noise=sigma/np.sqrt(2)*(rng.normal(size=Yc.shape)+1j*rng.normal(size=Yc.shape));Y=Yc+noise
            np.savez_compressed(OUT/'sealed_truth'/'truth.npz',fine_material=m,inverse_material=four_cubes(30),Y_clean=Yc,noise=noise)
            public=OUT/'public';public.mkdir(exist_ok=True)
            np.savez_compressed(public/'observations.npz',Y=Y,sigma=sigma,noise_radius=sigma*np.sqrt(Y.size),
                inverse_side=30,width_lambda=3.,material_box=np.array([0.,4.]))
            dump(public/'manifest.json',dict(scope=config['scope'],config='../../../configs/original_vector3d.json',observation_shape=list(Y.shape),
                sigma=float(sigma),noise_rule=config['noise_implementation'],receiver_rows='60 receiver points x3 Cartesian electric components',
                sources='60 electric dipoles; literal three offset radius3 circles',physics='independent CM+RR DDA,60^3 generation/30^3 inverse; full3lambdaROI',
                backend='Mac SciPy qualified checkpoints + Windows resident CUDA unfinished sources; each resumed state requalified on original matrix',truth_sealed_separately=True))
        elif args.mode=='qualify-discretization':
            if not (OUT/'public'/'observations.npz').exists():raise RuntimeError('Public data must be frozen first')
            if (OUT/'DISCRETIZATION_QUALIFICATION.json').exists():raise RuntimeError('No repeated grid qualification')
            truth=np.load(OUT/'sealed_truth'/'truth.npz');m=truth['inverse_material'];Yc=truth['Y_clean']
            model=CudaOriginalModel(OriginalVectorModel(30),fft_layout=config.get('CUDA_fft_layout','strided'),block_rhs=config.get('CUDA_block_rhs',4))
            def progress(r):
                with (OUT/'discretization_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(r),flush=True)
            coarse=model.physical(m,label='offline_coarse30_GT_qualification',progress=progress)
            np.savez_compressed(OUT/'sealed_truth'/'coarse_truth_response.npz',prediction=coarse)
            dump(OUT/'DISCRETIZATION_QUALIFICATION.json',dict(scope='OFFLINE_REFERENCE_ONLY; no inverse input/tuning/grid sweep',
                fine_side=60,inverse_side=30,width_lambda=3.,relative_global_L2=float(la.norm(coarse-Yc)/la.norm(Yc)),
                relative_source_RMS=float(np.sqrt(np.mean(np.sum(abs(coarse-Yc)**2,axis=0)/np.sum(abs(Yc)**2,axis=0)))),
                caveat='Independent CM+RR DDA refinement discrepancy; not continuum or author-integration error certificate'))
        else:
            # This branch reads only public observations/config; no truth loading.
            if (OUT/'STOP_CURRENT_JOB').exists():raise ControlledEngineeringStop('Explicit stop marker exists before inverse start')
            Yhost=np.load(OUT/'public'/'observations.npz')['Y']
            model=CudaOriginalModel(OriginalVectorModel(30),fft_layout=config.get('CUDA_fft_layout','strided'),block_rhs=config.get('CUDA_block_rhs',4));Y=model.array(Yhost)
            receiver=CudaReceiverSubspace(model,Y,config['dpcs_relative_cutoff'])
            shared=dict(CPU_seconds=time.process_time(),wall_seconds=time.time()-psutil.Process().create_time(),
                counts=dict(model.counts),setup=model.setup,receiver=receiver.info,arms=5,amortization='one fifth to each fixed arm',
                scope='One shared resident G/E/S/rank/Jd/Qd/Dd precompute, including Python imports/startup CPU')
            dump(OUT/'SHARED_CUDA_PRECOMPUTE_COST.json',shared);wiring(model,receiver,Y)
            for arm in config['arms']:
                dest=OUT/'inversions'/arm
                if (dest/'SEALED.json').exists():raise RuntimeError('Refusing duplicate arm '+arm)
                dest.mkdir(parents=True,exist_ok=True);before=dict(model.counts);nrecords=len(model.solve_records)
                awall,acpu=time.perf_counter(),time.process_time();endpoints=[];armstatus='FAILED';error=None
                dump(dest/'FROZEN_CONFIG.json',dict(config,arm=arm,backend='sharedResidentCUDAcomplex128',truth_access=False,
                    shared_precompute='../../SHARED_CUDA_PRECOMPUTE_COST.json',Aplus_variant='R9vectorAplus_with_CIE_seed_proposal'))
                try:
                    def callback(r,m):
                        r.update(cumulative_CPU_seconds=time.process_time()-acpu+shared['CPU_seconds']/5,
                            cumulative_wall_seconds=time.perf_counter()-awall+shared['wall_seconds']/5,
                            process_CPU_seconds=time.process_time(),runtime_counts=dict(model.counts))
                        with (dest/'history.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                        print(json.dumps(dict(arm=arm,**r)),flush=True)
                        if (OUT/'STOP_CURRENT_JOB').exists():
                            np.savez_compressed(dest/'controlled_stop_material.npz',material=m.cpu().numpy())
                            raise ControlledEngineeringStop('Explicit engineering stop at fully accounted auxiliary boundary')
                    def endpoint(name,m,stage):
                        np.savez_compressed(dest/(name+'.npz'),material=m.cpu().numpy())
                        def progress(r):
                            with (dest/'physical_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                            if r['source']%10==0:print(json.dumps(dict(arm=arm,**r)),flush=True)
                        value=model.physical(m,Y,label=arm+'_'+name,progress=progress)
                        np.savez_compressed(dest/(name+'_prediction.npz'),prediction=value.pop('prediction'))
                        endpoints.append(dict(name=name,stage=stage,**value,cumulative_CPU_seconds=time.process_time()-acpu+shared['CPU_seconds']/5,
                            cumulative_wall_seconds=time.perf_counter()-awall+shared['wall_seconds']/5))
                        dump(dest/'ENDPOINTS.json',endpoints)
                    result=cuda_auxiliary_solve(model,receiver,Y,arm,config,callback,endpoint)
                    monotonic=[]
                    if arm in ('sourceCIE','LSIEsource'):
                        for s in range(2):
                            values=np.array([h['auxiliary_sum'] for h in result['history'] if h['stage']==s]);diff=np.diff(values)
                            monotonic.append(dict(stage=s,max_absolute_increase=float(max(0,diff.max())) if len(diff) else 0.,
                                maximum_signed_change=float(diff.max()) if len(diff) else 0.,
                                max_relative_increase=float(max(0,np.max(diff/np.maximum(abs(values[:-1]),1e-300)))) if len(diff) else 0.))
                    dump(dest/'SUMMARY.json',dict(stages=result['stages'],receiver=receiver.info,final=endpoints[-1],
                        endpoint_selection='Final stage primary; all endpoints retained; no GT selection',prior=0,rho=1,
                        within_stage_monotonicity_diagnostic=monotonic,monotonicity_scope='Required for exact sourceCIE/LSIE only; not Dual/Aplus'))
                    dump(dest/'SEALED.json',dict(status='PREDICTIONS_SEALED_BEFORE_OFFLINE_SCORING',arm=arm,stages=result['stages']));armstatus='COMPLETE'
                except BaseException as exc:
                    error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                    arm_failures.append(dict(arm=arm,**error))
                    # Algorithm failure is retained and remaining independently frozen arms continue.
                    if isinstance(exc,(KeyboardInterrupt,SystemExit)):raise
                finally:
                    torch.cuda.synchronize()
                    dump(OUT/('invert_'+arm+'_COST.json'),dict(status=armstatus,backend='matchedResidentCUDAcomplex128',
                        CPU_seconds=time.process_time()-acpu+shared['CPU_seconds']/5,wall_seconds=time.perf_counter()-awall+shared['wall_seconds']/5,
                        raw_arm_CPU_seconds=time.process_time()-acpu,raw_arm_wall_seconds=time.perf_counter()-awall,
                        shared_precompute_CPU_share=shared['CPU_seconds']/5,shared_precompute_wall_share=shared['wall_seconds']/5,
                        shared_precompute_counted_once_globally=True,counts=delta_counts(model,before),solve_records=model.solve_records[nrecords:],
                        max_physical_linear_residual=max([r.get('relative_residual',0.) for r in model.solve_records[nrecords:]] or[0.]),
                        peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),failure=error))
                gc.collect()
        status='COMPLETE_WITH_ARM_FAILURES' if arm_failures else 'COMPLETE'
    except BaseException as exc:
        if isinstance(exc,ControlledEngineeringStop):status='CONTROLLED_ENGINEERING_STOP'
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
    finally:
        if torch.cuda.is_available():torch.cuda.synchronize()
        dump(OUT/(args.mode+'_CUDA_COST.json'),dict(status=status,mode=args.mode,
            CPU_seconds=time.process_time(),wall_seconds=time.time()-psutil.Process().create_time(),
            main_body_CPU_seconds=time.process_time()-cput,main_body_wall_seconds=time.perf_counter()-wall,
            timer_scope='Whole Python process including imports; arm shares count common precompute once',
            platform=platform.platform(),numpy_version=np.__version__,scipy_version=scipy.__version__,torch_version=torch.__version__,
            process_affinity=psutil.Process().cpu_affinity(),peak_RSS_bytes=psutil.Process().memory_info().peak_wset,
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated() if model else 0,
            peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved() if model else 0,
            counts=dict(model.counts) if model else{},setup=model.setup if model else None,
            solve_records=model.solve_records if model else[],failure=failure,arm_failures=arm_failures))

if __name__=='__main__':main()
