"""One passive derivative check and exactly two parent-authorized controls."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,sys,time,traceback
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];R2=REPO/'research/research_max_struct_r2'
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r3.src.passive_continuation import PassiveModel
from research.research_max_struct_r2.experiments.run_blind import dump,compact,metric

CASE='obs_b20d5a55677d4eca'

def check_passive_gradient():
    path=ROOT/'audit/PASSIVE_GRADIENT_WIRING.json'
    if path.exists():raise FileExistsError('Preserve one-direction check')
    start=time.perf_counter();cpu=time.process_time();models=[]
    rec=dict(status='RUNNING',classification='ONE_NEW_ADAPTER_DIRECTIONAL_WIRING_CHECK',gamma=1.,
        h=1e-5,directions=1,inverse_runs=0,truth_access='disclosed constant4x4 root only',
        derivative_review='For fixedgamma dc/dm=1; inherited real adjoint gradient and realm prior are unchanged',
        cache_review='gamma fixed per instance; changinggamma aftercache is unsupported',hash_checks=False)
    try:
        with np.load(REPO/'research/research_max_struct_r1/data_public/dynamics_operators_4x4.npz',allow_pickle=False)as f:
            ops={k:f[k]for k in ('G','E0','S','grid_shape')}
        root=MaxwellModel(ops);models.append(root);Y=root.state(np.full(root.n,.55))['Y'].copy()
        model=PassiveModel(ops,1.);models.append(model)
        m=.2+.03*np.cos(np.arange(model.n));direction=np.cos(np.arange(model.n)+.37);direction/=np.linalg.norm(direction)
        value,g=model.data_value_gradient(m,Y);h=1e-5
        plus=model.physical_score(m+h*direction,Y);minus=model.physical_score(m-h*direction,Y)
        analytic=float(g@direction);finite=float((plus-minus)/(2*h))
        rec.update(status='COMPLETED',analytic=analytic,finite=finite,absolute_error=abs(analytic-finite),value=float(value),
                   relative_error=abs(analytic-finite)/max(abs(analytic),abs(finite),1e-300))
    except Exception as exc:rec.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        costs=Counter()
        for model in models:costs.update(model.cost)
        rec.update(physical_model_cost=dict(costs),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
                   timer_scope='public load,1rootforward,3passiveforwards,1analyticadjoint; ownreceiptwrite excluded')
        dump(path,rec)
    print(json.dumps(rec),flush=True)

def scaled_original(model,Y,m0):
    start=time.perf_counter();cpu=time.process_time();before=Counter(model.cost);scale=1e6
    history=[];best=[m0.copy(),np.inf]
    def objective(m):
        value,g=model.data_value_gradient(m,Y)
        history.append(dict(original_objective=float(value),scaled_objective=float(scale*value),
                            gradient_inf=float(np.linalg.norm(g,np.inf)),wall_seconds=time.perf_counter()-start))
        if value<best[1]:best[:]=[np.asarray(m).copy(),float(value)]
        return scale*value,scale*g
    result=minimize(objective,np.asarray(m0),jac=True,method='L-BFGS-B',bounds=[(.02,20.)]*model.n,
        options=dict(maxiter=1750,maxfun=7000,ftol=1e-12,gtol=1e-8*scale,maxls=30))
    return dict(m=best[0],history=history,cost=dict(Counter(model.cost)-before),stop=str(result.message),
                optimizer_nit=int(result.nit),optimizer_nfev=int(result.nfev),scale=scale,
                gradient_tolerance_scaled=.01,ftol=1e-12,
                wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)

def run(arm):
    directory=ROOT/'results/development'/arm;directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json';prediction=directory/'prediction.npz'
    if path.exists()or prediction.exists():raise FileExistsError('Preserve authorized single disposition')
    start=time.perf_counter();cpu=time.process_time();model=None
    rec=dict(status='RUNNING',run_id=arm,case_id=CASE,classification='POST_HOC_DEVELOPMENT_PUBLISHED_R2',
        started_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',code='experiments/known_controls.py',
        initialization='identical R2 corrected native material',prior='unchanged historical prior',material_box=[.02,20.],
        online_selection='original normalized physicalPhi including sameprior only, includesinitialseed',
        threads=1,device='CPU',hash_checks=False,inverse_runs=1)
    dump(path,rec)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:public={k:f[k]for k in f.files}
        with np.load(R2/'predictions'/CASE/'native.npz',allow_pickle=False)as f:m0=f['m']
        model=MaxwellModel(public);Y=public['Y'];initial=metric(model,Y,m0)
        if arm=='globaldual_restart':
            rec['schedule']='known GlobalDual1000 +same original physical L-BFGS-Bcorrector1000'
            aux=corrector(model,Y,m0,dict(method='globaldual',max_nfev=1000,box=(.02,20.)))
            seed=min([initial,metric(model,Y,aux['m'])],key=lambda c:c['objective'])
            corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
            candidates=[initial,seed,metric(model,Y,corr['m'])]
            dump(directory/'AUXILIARY.json',compact(aux));dump(directory/'CORRECTOR.json',compact(corr))
            rec.update(auxiliary_stop=aux['stop'],corrector_stop=corr['stop'],auxiliary_cost=aux['cost'],corrector_cost=corr['cost'])
        elif arm=='scaled_original_restart':
            rec['schedule']='physicalPhi andgradient×1e6 L-BFGS-B1750;gtol×1e6,ftolunchanged;noextracorrector'
            aux=scaled_original(model,Y,m0);candidates=[initial,metric(model,Y,aux['m'])]
            dump(directory/'OPTIMIZER.json',compact(aux))
            rec.update(optimizer_stop=aux['stop'],optimizer_nit=aux['optimizer_nit'],optimizer_nfev=aux['optimizer_nfev'],
                       objective_scale=1e6,gtol_scaled=.01,gtol_original_equivalent=1e-8,auxiliary_cost=aux['cost'])
        else:raise ValueError(arm)
        best=min(candidates,key=lambda c:c['objective']);state=model.state(best['m']);value,g=model.data_value_gradient(best['m'],Y)
        pg=best['m']-np.clip(best['m']-g,.02,20.)
        np.savez_compressed(prediction,m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str(prediction.relative_to(ROOT)),physical_objective=float(value),
            physical_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),projected_gradient_l2=float(np.linalg.norm(pg)),
            physical_linear_residual=state['state_relative_residual'],candidate_count=len(candidates),
            limitations='same disclosedseed; unequal auxiliary budgets exposed; solver stopping differs; no confirmation or basin certificate')
    except Exception as exc:rec.update(status='FAILED_DEVELOPMENT_DISPOSITION',error=repr(exc),traceback=traceback.format_exc())
    finally:
        rec.update(physical_model_cost=dict(model.cost)if model else {},physical_operation_wall=dict(model.timings)if model else {},
                   physical_operation_cpu=dict(model._physics.cpu_timings)if model else {},wall_seconds=time.perf_counter()-start,
                   cpu_seconds=time.process_time()-cpu,FFT=0,GPU=0,
                   timer_scope='load,optimization,physicalselection,finalgradient,serialization; finalreceiptwrite excluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['run_id','status','physical_residual','projected_gradient_inf','wall_seconds','cpu_seconds','physical_model_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        check_passive_gradient()
        for arm in ('globaldual_restart','scaled_original_restart'):run(arm)
