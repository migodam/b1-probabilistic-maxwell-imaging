"""Original-scale independent-DDA runner. One numeric process; BLAS threads=1.

Truth is accessed only by generate/offline-score. Invert reads public observations
and the frozen config, never clean fields or material truth.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,json,resource,sys,time,traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy import linalg as la
from research.research_max_struct_r9.src.vector_original_extension import (
    OriginalVectorModel,ReceiverSubspace,ImplicitFourier,LocalArcReadout,material_readout,
    four_cubes,auxiliary_solve)
BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'RESULTS'/'original_vector3d'
CFG=BASE/'configs'/'original_vector3d.json'

def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)

def model_options(config):
    return {k:config[k] for k in ('block_rhs','gmres_rtol','gmres_restart','gmres_cycles')}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['generate','wire','invert','qualify-discretization','offline-score'])
    ap.add_argument('--arm',choices=['sourceCIE','sourceDual','GlobalDual','AplusGlobal','LSIEsource'])
    args=ap.parse_args();config=json.loads(CFG.read_text());OUT.mkdir(parents=True,exist_ok=True)
    wall,cpu=time.perf_counter(),time.process_time();model=None;status='FAILED';failure=None
    receipt=OUT/(args.mode+('_'+args.arm if args.arm else '')+'_COST.json')
    try:
        if args.mode=='generate':
            if (OUT/'public'/'observations.npz').exists():raise RuntimeError('Refusing to overwrite paid public observations')
            truth=four_cubes(config['generation_side']);active=np.flatnonzero(truth!=0)
            model=OriginalVectorModel(config['generation_side'],active_indices=active,generation_only=True,
                crop_generation_fft=config['generation_crop_fft'],**model_options(config))
            dump(OUT/'GENERATION_SETUP.json',model.setup)
            def progress(r):
                with (OUT/'generation_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(r),flush=True)
            Yc=model.physical(truth[active],label='fine60_generation',progress=progress,
                checkpoint_dir=OUT/'sealed_truth'/'generation_checkpoints')
            sigma=config['noise_relative_expected_norm']*la.norm(Yc)/np.sqrt(Yc.size)
            rng=np.random.default_rng(config['noise_seed']);noise=sigma/np.sqrt(2)*(rng.normal(size=Yc.shape)+1j*rng.normal(size=Yc.shape))
            Y=Yc+noise
            sealed=OUT/'sealed_truth';sealed.mkdir(exist_ok=True)
            np.savez_compressed(sealed/'truth.npz',fine_material=truth,inverse_material=four_cubes(config['inverse_side']),Y_clean=Yc,noise=noise)
            public=OUT/'public';public.mkdir(exist_ok=True)
            np.savez_compressed(public/'observations.npz',Y=Y,sigma=sigma,noise_radius=sigma*np.sqrt(Y.size),
                inverse_side=config['inverse_side'],width_lambda=3.,material_box=np.array([0.,4.]))
            dump(public/'manifest.json',dict(scope=config['scope'],config='../../../configs/original_vector3d.json',
                observation_shape=list(Y.shape),sigma=float(sigma),noise_rule=config['noise_implementation'],
                receiver_rows='60 receiver points x3 Cartesian electric components',sources='60 electric dipoles; literal three offset radius3 circles',
                physics='independent CM+RR DDA,60^3 generation/30^3 inverse; full3lambdaROI',truth_sealed_separately=True))
        elif args.mode=='wire':
            # One paid wiring check on the new implementation; not an imaging run.
            model=OriginalVectorModel(6,**model_options(config));rng=np.random.default_rng(119)
            Y=rng.normal(size=(180,60))+1j*rng.normal(size=(180,60))
            r=ReceiverSubspace(model,Y,.1);B=ImplicitFourier(model,r,2)
            gamma=rng.normal(size=(B.K,2))+1j*rng.normal(size=(B.K,2));d=rng.normal(size=gamma.shape)+1j*rng.normal(size=gamma.shape)
            x=rng.normal(size=(model.M,2))+1j*rng.normal(size=(model.M,2))
            bg=B.apply(gamma);bd=B.apply(d);gx=model.green(x)
            dotB=abs(np.vdot(bg,x)-np.vdot(gamma,B.adjoint(x)))/max(abs(np.vdot(bg,x)),1)
            dotG=abs(np.vdot(model.green(bg),x)-np.vdot(bg,model.green(x,True)))/max(abs(np.vdot(model.green(bg),x)),1)
            m=np.linspace(.1,2,model.N);chi=np.repeat(m,3)[:,None];p=6/(1+6*chi)
            target=Y[:,:2];E=model.inc[:,:2];wd=1/np.sum(abs(target)**2,axis=0);ws=1/np.sum(abs(E)**2,axis=0)
            def value(g):
                j=B.apply(g);q=E+model.green(j);rr=p*(j-chi*q);dd=model.S@j-target
                return float(np.sum(abs(dd)**2*wd)+np.sum(abs(rr)**2*ws))
            q=E+model.green(bg);v=p.conj()*p*(bg-chi*q)*ws
            grad=B.adjoint(model.S.conj().T@((model.S@bg-target)*wd)+v-model.green(chi*v,True))
            analytic=2*np.vdot(grad,d).real;step=1e-5
            fd=(value(gamma+step*d)-value(gamma-step*d))/(2*step)
            # New vectorized coordinate arc is also charged in this single check.
            jf=rng.normal(size=(model.M,60))+1j*rng.normal(size=(model.M,60));qf=model.inc+model.green(jf)
            ww=1/np.sum(abs(model.inc)**2,axis=0);old=np.linspace(0,4,model.N)
            mm,ai=LocalArcReadout(model)(jf,qf,old,ww,6)
            result=dict(implicit_adjoint_relative=float(dotB),Green_adjoint_relative=float(dotG),
                current_gradient_analytic=float(analytic),current_gradient_finite_difference=float(fd),
                current_gradient_relative_error=float(abs(fd-analytic)/max(abs(fd),abs(analytic),1e-30)),arc=ai,
                box_min=float(mm.min()),box_max=float(mm.max()),physics_solves=0)
            result['PASS']=bool(dotB<1e-10 and dotG<1e-10 and result['current_gradient_relative_error']<1e-6 and ai['fixed_current_loss_after']<=ai['fixed_current_loss_before']*(1+1e-12))
            dump(OUT/'WIRING.json',result)
            if not result['PASS']:raise RuntimeError('New implicit wiring failed')
        elif args.mode=='qualify-discretization':
            if not (OUT/'public'/'observations.npz').exists():raise RuntimeError('Public data must be frozen before offline truth qualification')
            if (OUT/'DISCRETIZATION_QUALIFICATION.json').exists():raise RuntimeError('No repeated qualification or grid sweep')
            truth=np.load(OUT/'sealed_truth'/'truth.npz');m=truth['inverse_material'];Yclean=truth['Y_clean']
            model=OriginalVectorModel(config['inverse_side'],**model_options(config))
            def progress(r):
                with (OUT/'discretization_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(r),flush=True)
            coarse=model.physical(m,label='offline_coarse30_GT_qualification',progress=progress)
            np.savez_compressed(OUT/'sealed_truth'/'coarse_truth_response.npz',prediction=coarse)
            dump(OUT/'DISCRETIZATION_QUALIFICATION.json',dict(scope='OFFLINE_REFERENCE_ONLY; no inverse inputs, tuning, or grid sweep',
                fine_side=60,inverse_side=30,width_lambda=3.,relative_global_L2=float(la.norm(coarse-Yclean)/la.norm(Yclean)),
                relative_source_RMS=float(np.sqrt(np.mean(np.sum(abs(coarse-Yclean)**2,axis=0)/np.sum(abs(Yclean)**2,axis=0)))),
                caveat='Independent CM+RR DDA refinement discrepancy; not a continuum or author-integration error certificate'))
        elif args.mode=='invert':
            if not args.arm:raise ValueError('--arm required')
            dest=OUT/'inversions'/args.arm
            if (dest/'SEALED.json').exists():raise RuntimeError('Refusing to overwrite sealed predictions')
            dest.mkdir(parents=True,exist_ok=True)
            if not json.loads((OUT/'WIRING.json').read_text())['PASS']:raise RuntimeError('Wiring not qualified')
            Y=np.load(OUT/'public'/'observations.npz')['Y'];model=OriginalVectorModel(config['inverse_side'],**model_options(config))
            dump(dest/'FROZEN_CONFIG.json',dict(config,arm=args.arm,truth_access=False))
            endpoints=[]
            def callback(r):
                r.update(cumulative_CPU_seconds=time.process_time()-cpu,cumulative_wall_seconds=time.perf_counter()-wall)
                with (dest/'history.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                print(json.dumps(r),flush=True)
            def endpoint(name,m,stage):
                np.savez_compressed(dest/(name+'.npz'),material=m)
                def physical_progress(r):
                    with (dest/'physical_progress.jsonl').open('a') as f:f.write(json.dumps(r)+'\n')
                    print(json.dumps(r),flush=True)
                value=model.physical(m,Y,label=name,progress=physical_progress)
                np.savez_compressed(dest/(name+'_prediction.npz'),prediction=value.pop('prediction'))
                endpoints.append(dict(name=name,stage=stage,**value,
                    cumulative_CPU_seconds=time.process_time()-cpu,cumulative_wall_seconds=time.perf_counter()-wall))
                dump(dest/'ENDPOINTS.json',endpoints)
            result=auxiliary_solve(model,Y,args.arm,config,callback,endpoint)
            dump(dest/'SUMMARY.json',dict(stages=result['stages'],receiver=result['receiver'],endpoint_selection='All stage endpoints retained; final stage is primary, no GT selection',
                final=endpoints[-1],prior=0,rho=1))
            dump(dest/'SEALED.json',dict(status='PREDICTIONS_SEALED_BEFORE_OFFLINE_SCORING',arm=args.arm,stages=result['stages']))
        else:
            truth=np.load(OUT/'sealed_truth'/'truth.npz')['inverse_material'];scores=[]
            for dest in sorted((OUT/'inversions').iterdir()):
                if not (dest/'SEALED.json').exists():continue
                for ep in json.loads((dest/'ENDPOINTS.json').read_text()):
                    m=np.load(dest/(ep['name']+'.npz'))['material']
                    scores.append(dict(arm=dest.name,endpoint=ep['name'],**ep,
                        paper_Err=float(np.sqrt(np.mean(abs(m-truth)**2/abs(1+truth)**2))),
                        contrast_NRMSE=float(la.norm(m-truth)/la.norm(truth)),peak_contrast=float(m.max()),
                        lower_bound_fraction=float(np.mean(m==0)),upper_bound_fraction=float(np.mean(m==4))))
            dump(OUT/'OFFLINE_SCORES.json',dict(scope='OFFLINE_AFTER_PREDICTION_SEAL; no inverse feedback',scores=scores))
        status='COMPLETE'
    except BaseException as exc:
        failure=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc());raise
    finally:
        dump(receipt,dict(status=status,mode=args.mode,arm=args.arm,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            counts=dict(model.counts) if model else {},setup=model.setup if model else None,
            max_physical_linear_residual=model.max_residual if model else None,
            solve_records=model.solve_records if model else [],failure=failure))

if __name__=='__main__':main()
