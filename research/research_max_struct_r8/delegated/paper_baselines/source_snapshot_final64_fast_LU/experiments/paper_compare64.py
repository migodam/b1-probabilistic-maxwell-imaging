"""Original150/64 eps4 acquisition, five prescribed full-pixel arms.

Paper coefficient stopping with500 ceiling; matched80 physicalcorrector.
Writes paid stage/current progress for reliable long-run observability.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json,argparse
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write
from research.research_max_struct_r8.src.paper_physics import acquisition,MatrixFreeModel
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_max_struct_r8.src.paper_nie import run_paper_tsom
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout

def run(only=None):
    case='austria_eps4_400MHz';Y=np.load(OUT/case/'PUBLIC_DATA.npz')['Y'];public=acquisition(64,400e6,dense=False)
    for arm in ('cie_single','cie_dual','a_plus','a_plus_sourcewise','paper_nie'):
        if only and arm!=only:continue
        dest=OUT/'original64'/case/arm
        if (dest/'PREDICTION_SEAL.json').exists():continue
        dest.mkdir(parents=True,exist_ok=True);model=MatrixFreeModel(public,prior=arm!='paper_nie');wall,cpu=time.perf_counter(),time.process_time()
        def progress(extra):
            write(dest/'PROGRESS.json',dict(phase='auxiliary',wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,physical_cost=dict(model.cost),**extra))
        stages=[dict(mf=mf,beta_material=b,beta_current=b if arm=='cie_single' else 8,maxiter=500,coefficient_tol=tol**2,paper_delta=tol) for mf,b,tol in [(5,8,.01),(8,2,.001),(10,.5,.001)]]
        m0=np.zeros(model.n)
        try:
            recovered=dest/'AUXILIARY_RECEIPT_RECOVERED.json'
            if recovered.exists():
                m=np.load(dest/'AUXILIARY_PREDICTION.npz')['m'];aux=json.loads(recovered.read_text());result=dict(history=[dict(event='RESUMED_PAID_AUXILIARY',source=aux['original_attempt'],history_gap=aux['auxiliary_operator_counter_gap'])])
            elif arm=='paper_nie':
                result=run_paper_tsom(model,Y,m0,stages,maxiter=500);cand=[s['m'] for s in result['stage_states']];m=min(cand,key=lambda z:model.physical_score(z,Y));aux=result['receipt']
            else:
                cfg=dict(arm='a_plus' if arm=='a_plus_sourcewise' else arm,stages=stages,dpcs_rank=15,material_box=(0,20),miniter=2,sourcewise_cg=arm in ('cie_single','cie_dual','a_plus_sourcewise'),regularizer=lambda z:model.reg.evaluate(z)[:2],fft_geometry=model.paper_geometry,physical_sigma=np.linalg.norm(Y),aplus_regularizer=model.reg,aplus_sweeps=2,aplus_arc_candidates=True,cache_material_value_gradient=True,iteration_callback=lambda st:progress(dict(stage=st['stage'],iteration=st['iteration'])) if st['iteration']%25==0 else None)
                if arm.startswith('a_plus'):cfg['material_readout']=mobius_material_readout
                result=run_cie_fft(model,Y,m0,cfg);m=result['m'];aux=result['execution_receipt']
            auxm=m.copy();history=[];best=[m.copy(),model.physical_score(m,Y)]
            np.savez_compressed(dest/'AUXILIARY_PREDICTION.npz',m=auxm,Y_prediction=model.state(auxm)['Y']);write(dest/'AUXILIARY_SEAL.json',dict(status='SEALED_NO_GT',GT_loaded=False));write(dest/'AUXILIARY_RECEIPT.json',aux);write(dest/'AUXILIARY_HISTORY.json',result['history']);progress(dict(phase_detail='auxiliary_complete'))
            def fg(z):
                try:v,g=model.data_value_gradient(z,Y)
                except RuntimeError as exc:
                    if not str(exc).startswith('Matrix-free physical solve failed'):raise
                    history.append(dict(evaluation=len(history)+1,status='REJECTED_UNSOLVED_LINE_CANDIDATE',error=str(exc),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost)))
                    return np.inf,np.zeros_like(z)
                if v<best[1]:best[:]=[z.copy(),v]
                history.append(dict(evaluation=len(history)+1,Phi=v,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,complex_RHS=model.cost['actual_complex_RHS']))
                if len(history)%10==0:write(dest/'PROGRESS.json',dict(phase='physical_corrector',evaluation=len(history),wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost)))
                return 1e6*v,1e6*g
            opt=minimize(fg,m,jac=True,method='L-BFGS-B',bounds=[(0,20)]*model.n,options=dict(maxiter=80,maxfun=320,ftol=1e-12,gtol=.01,maxls=30));final=best[0];st=model.state(final)
            np.savez_compressed(dest/'PREDICTION.npz',m=final,m_auxiliary=auxm,Y_prediction=st['Y']);write(dest/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_GT',time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),GT_loaded=False))
            reused_cpu=float(aux.get('cpu_seconds',0)) if recovered.exists() else 0.;reused_wall=float(aux.get('wall_seconds',0)) if recovered.exists() else 0.
            physical_cost=dict(model.cost)
            if recovered.exists():
                for k,v in aux['physical_cost'].items():physical_cost[k]=max(physical_cost.get(k,0),v) if k.startswith('maximum_') else physical_cost.get(k,0)+v
            rejections=sum(x.get('status')=='REJECTED_UNSOLVED_LINE_CANDIDATE' for x in history)
            write(dest/'HISTORY.json',dict(auxiliary=result['history'],corrector=history));write(dest/'RECEIPT.json',dict(status='COMPLETED',arm=arm,wall_seconds=time.perf_counter()-wall+reused_wall,cpu_seconds=time.process_time()-cpu+reused_cpu,reused_auxiliary_cpu_seconds=reused_cpu,reused_auxiliary_wall_seconds=reused_wall,corrector_history='freshL-BFGSfrompaidauxiliary' if recovered.exists() else 'freshL-BFGSfromselectedauxiliary',physical_cost=physical_cost,auxiliary=aux,relative_data_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),corrector_iterations=int(opt.nit),corrector_evaluations=int(opt.nfev),stop=str(opt.message),optimizer_numerical_rejections=rejections,convergence_qualified=bool(opt.success and rejections==0),forward_grid=[150,150],inverse_grid=[64,64],all_material_pixels_free=True,GT_loaded=False))
            print('ORIGINAL64',case,arm,time.perf_counter()-wall,flush=True)
        except Exception as e:
            write(dest/'FAILURE.json',dict(status='FAILED',error=repr(e),wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost)));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--arm');a=p.parse_args();run(a.arm)
