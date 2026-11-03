"""One fixed-eta lambda0/lambda1 comparison; public observations only."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import json,sys,time,traceback,argparse
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.experiments.run_blind import dump,metric
from research.research_max_struct_r4.experiments.cooperative_paths import refine
from research.research_max_struct_r4.src.port_profile import PortProfile,fixed_eta_reference


def utc():return datetime.now(timezone.utc).isoformat()


def wiring():
    started,cpu=time.perf_counter(),time.process_time()
    path=ROOT/'results/port_profile_setup/GRADIENT_WIRING.json'
    if path.exists():raise FileExistsError('Preserve one authorized gradient direction; no repetition')
    models=[]
    rec=dict(status='RUNNING',started_utc=utc(),purpose='ONE_NEW_ENVELOPE_GRADIENT_WIRING_DIRECTION',
             side=4,lambda_value=1,directions=1,step=1e-5,inverse_runs=0,
             truth_access='Public disclosed constant0.55 root only; no material/support truth file',
             eta_policy='same fixed host0.1 reference for both central perturbations',hash_checks=False)
    dump(path,rec)
    try:
        with np.load(REPO/'research/research_max_struct_r1/data_public/dynamics_operators_4x4.npz',allow_pickle=False)as f:
            public={key:f[key].copy()for key in ('G','E0','S','grid_shape')}
        root=MaxwellModel(public);models.append(root)
        Y=root.state(np.full(root.n,.55))['Y'].copy()
        model=MaxwellModel(public);models.append(model)
        eta,mref,Cref=fixed_eta_reference(model,.1)
        objective=PortProfile(model,Y,eta,1,reference=(mref,Cref))
        m=.2+.03*np.cos(np.arange(model.n));h=np.cos(np.arange(model.n)+.37);h/=la.norm(h)
        value,g,_,_=objective.evaluate(m)
        plus=objective.evaluate(m+1e-5*h,gradient=False)[0]
        minus=objective.evaluate(m-1e-5*h,gradient=False)[0]
        analytic=float(g@h);finite=float((plus-minus)/2e-5)
        rec.update(status='COMPLETED_ONE_DIRECTION',eta=eta,analytic=analytic,finite=finite,
                   absolute_error=abs(analytic-finite),
                   relative_error=abs(analytic-finite)/max(abs(analytic),abs(finite),np.finfo(float).tiny),
                   auxiliary_objective=value)
    except Exception as exc:
        rec.update(status='FAILED_PAID_WIRING',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        costs=Counter()
        for model in models:costs.update(model.cost)
        rec.update(physical_cost=dict(costs),wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
                   timer_scope='public operator load,constant observation,reference,one gradient/two central values;imports/finalreceiptwrite excluded')
        dump(path,rec)
    print(json.dumps(rec),flush=True)


def run_arm(public,reference,lam,cfg):
    started,cpu=time.perf_counter(),time.process_time()
    directory=ROOT/f'results/port_profile_{lam}';path=directory/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve fixed-lambda result; no rerun')
    model=None;history=[];materials=[];best=[None,np.inf,None,None]
    eta,mref,Cref=reference
    rec=dict(status='RUNNING',run_id=f'port_profile_{lam}',started_utc=utc(),
             classification='FIXED_FINITE_ETA_DEVELOPMENT_CONTROL',truth_access='NONE',configuration=cfg,
             lambda_value=lam,eta=eta,shared_reference_receipt='../port_profile_setup/REFERENCE.json',
             label=cfg['lambda0_label']if lam==0 else cfg['lambda1_label'],
             source='experiments/port_profile.py',candidate_selection='original physical Phi, includes public host',
             corrector_calls=0,hash_checks=False)
    dump(path,rec);dump(directory/'FROZEN_CONFIG.json',cfg)
    try:
        model=MaxwellModel(public);Y=public['Y']
        profile=PortProfile(model,Y,eta,lam,reference=(mref,Cref))
        m0=np.full(model.n,cfg['material_initial'])
        scale=cfg['objective_scale']
        def objective(m):
            value,g,detail,state=profile.evaluate(m)
            row=dict(evaluation=len(history),**detail,auxiliary_gradient_inf=float(la.norm(g,np.inf)))
            history.append(row);materials.append(np.asarray(m).copy())
            if detail['physical_objective']<best[1]:
                best[:]=[np.asarray(m).copy(),detail['physical_objective'],state['Y'].copy(),row.copy()]
            if len(history)%50==0:
                rec.update(auxiliary_evaluations=len(history),last_evaluation=row,
                           selected_original_Phi=float(best[1]),physical_cost=dict(model.cost))
                dump(path,rec)
                print(json.dumps(dict(arm=lam,evaluations=len(history),auxiliary=value,
                    physical=detail['physical_objective'],selected_Phi=best[1])),flush=True)
            return scale*value,scale*g
        optimizer=minimize(objective,m0,method='L-BFGS-B',jac=True,bounds=[tuple(cfg['material_box'])]*model.n,
                           options=dict(maxiter=cfg['auxiliary_maxiter'],maxfun=cfg['auxiliary_maxfun'],
                                        gtol=cfg['gtol_scaled'],ftol=cfg['ftol'],maxls=cfg['maxls']))
        aux_cost=dict(model.cost)
        np.savez_compressed(directory/'AUXILIARY_EVALUATIONS.npz',m=np.asarray(materials),
                            auxiliary_objective=np.asarray([row['auxiliary_objective']for row in history]),
                            physical_objective=np.asarray([row['physical_objective']for row in history]))
        dump(directory/'AUXILIARY_HISTORY.json',history)
        np.savez_compressed(directory/'AUXILIARY_OPTIMIZER_ENDPOINT.npz',m=optimizer.x,gradient=optimizer.jac/scale)
        optimizer_info=dict(stop=str(optimizer.message),nit=int(optimizer.nit),nfev=int(optimizer.nfev),
                            auxiliary_objective=float(optimizer.fun/scale),
                            original_equivalent_gtol=cfg['gtol_scaled']/scale)
        dump(directory/'AUXILIARY_OPTIMIZER.json',optimizer_info)
        seed=best[0]
        np.savez_compressed(directory/'SELECTED_BEFORE_CORRECTOR.npz',m=seed,Y_prediction=best[2])
        rec.update(status='CORRECTOR_RUNNING',auxiliary_optimizer=optimizer_info,
                   auxiliary_cost=aux_cost,auxiliary_evaluations=len(history),
                   selected_before_corrector=best[3],corrector_calls=1)
        dump(path,rec)
        final,corr=refine(model,Y,seed,maximum=cfg['corrector_maxiter'])
        dump(directory/'CORRECTOR.json',corr)
        final_score=metric(model,Y,final)
        value,g=model.data_value_gradient(final,Y)
        projected=final-np.clip(final-g,*cfg['material_box'])
        state=model.state(final)
        np.savez_compressed(directory/'prediction.npz',m=final,Y_prediction=state['Y'],gradient=g,projected_gradient=projected)
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=utc(),prediction='prediction.npz',
                   physical_objective=float(value),physical_relative_residual=final_score['relative_residual'],
                   projected_gradient_inf=float(la.norm(projected,np.inf)),corrector_stop=corr['stop'],
                   corrector_iterations=corr['iterations'],corrector_evaluations=corr['evaluations'],
                   corrector_cost=corr['cost'],material_distance_from_selected=float(la.norm(final-seed)))
    except Exception as exc:
        rec.update(status='FAILED_PAID_CONTROL',error=repr(exc),traceback=traceback.format_exc())
        if history:
            dump(directory/'FAILED_PARTIAL_HISTORY.json',history)
            np.savez_compressed(directory/'FAILED_PARTIAL_MATERIALS.npz',m=np.asarray(materials))
        raise
    finally:
        rec.update(physical_cost=dict(model.cost)if model else {},
                   operation_wall_seconds=dict(model.timings)if model else {},
                   operation_cpu_seconds=dict(model._physics.cpu_timings)if model else {},
                   wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,FFT=0,GPU=0,
                   cost_semantics='One physical LU per new material; retained receiver48/current12/adjoint12 RHS per profile gradient; small48-source-system work separately counted',
                   timer_scope='arm model setup,aux optimization/candidate archive,onecorrector,finalgradient,serialization;imports/sharedreference/finalreceiptwrite excluded')
        dump(path,rec)
    print(json.dumps({key:rec[key]for key in ('run_id','status','eta','physical_objective','physical_relative_residual',
                      'projected_gradient_inf','corrector_stop','wall_seconds','physical_cost')}),flush=True)
    return rec


def run():
    cfg=json.loads((ROOT/'configs/port_profile.json').read_text())
    if cfg['lambdas']!=[0,1]:raise ValueError('Only the two frozen arms are authorized')
    check=json.loads((ROOT/'results/port_profile_setup/GRADIENT_WIRING.json').read_text())
    if check['status']!='COMPLETED_ONE_DIRECTION':raise RuntimeError('Completed one-direction wiring required')
    started,cpu=time.perf_counter(),time.process_time()
    path=ROOT/'results/port_profile_setup/REFERENCE.json'
    if path.exists():raise FileExistsError('Preserve shared reference; no eta tuning')
    model=None
    rec=dict(status='RUNNING',started_utc=utc(),truth_access='NONE',configuration=cfg,
             reference_material=.1,eta_policy='trace(Cref Cref*)/nr, fixed across both arms and all material evaluations')
    dump(path,rec)
    try:
        with np.load(ROOT.parent/'research_max_struct_r2/data_public'/f"{cfg['case_id']}.npz",allow_pickle=False)as f:
            public={key:f[key].copy()for key in ('G','E0','S','Y','grid_shape')}
        model=MaxwellModel(public)
        reference=fixed_eta_reference(model,cfg['material_initial'])
        eta,mref,Cref=reference
        np.savez_compressed(path.parent/'REFERENCE.npz',m=mref,C=Cref,eta=eta)
        rec.update(status='FIXED_REFERENCE_SEALED',eta=eta,n=model.n,nr=model.nr,ns=model.ns,
                   arrays='REFERENCE.npz',source_fields=list(public))
    except Exception as exc:
        rec.update(status='FAILED_PAID_REFERENCE',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        rec.update(physical_cost=dict(model.cost)if model else {},wall_seconds=time.perf_counter()-started,
                   cpu_seconds=time.process_time()-cpu,timer_scope='public load,host forward,receiver transpose solve,eta Gram,reference arrays;imports/finalreceiptwrite excluded')
        dump(path,rec)
    Cref.flags.writeable=False
    for lam in cfg['lambdas']:
        run_arm(public,reference,lam,cfg)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('operation',choices=['wiring','run']);args=parser.parse_args()
    with threadpool_limits(limits=1):
        wiring()if args.operation=='wiring'else run()
