"""Exactly two mechanism controls on disclosed R2 curved bar; no new blind GT."""
from pathlib import Path
from datetime import datetime,timezone
import json,sys,time,traceback
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];R2=REPO/'research/research_max_struct_r2'
sys.path.insert(0,str(REPO));sys.path.insert(0,str(R2/'experiments'))
from run_blind import native_endpoints,metric,compact,dump
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import run_cie
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
from research.research_max_struct_r3.src.cie_strong_release import run_strong_release

CASE='obs_b20d5a55677d4eca'
STAGES=[dict(mf=3,beta=1.,maxiter=500),dict(mf=6,beta_current=1.,beta_material=.2,maxiter=500),
        dict(mf=10,beta_current=1.,beta_material=.2,maxiter=500)]

def run(arm):
    directory=ROOT/'results/development'/arm;directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json';prediction=ROOT/'predictions'/(arm+'.npz')
    if path.exists()or prediction.exists():raise FileExistsError('Preserve single authorized disposition')
    start=time.perf_counter();cpu=time.process_time();model=None
    rec=dict(status='RUNNING',classification='OFFLINE_DIAGNOSTIC_POSTHOC_R2_DEVELOPMENT',case=CASE,arm=arm,
        started_utc=datetime.now(timezone.utc).isoformat(),code='experiments/development_controls.py',
        truth_policy='No GT material or support used by online solver; clean arm reads only already-published Y_clean',
        stages=STAGES,source_metrics='wd/ws frozen from original noisy publicY/E, even for clean target',
        coefficient_tol=0.,miniter=50,material_inner_maxiter=60,dpcs_rtol=1e-3,
        initialization='uniform public host .1',prior='same historical smooth1e-4/TV1e-4,width19.5,epsilon.01',
        material_box=[.02,20.],corrector=dict(method='original physical L-BFGS-B',maxiter=1000,maxfun=4000),
        online_selection='original physical objective+same prior, always includes uniforminitial; never GT score',
        threads=1,device='CPU',hash_checks=False,inverse_runs=1)
    dump(path,rec)
    try:
        with np.load(R2/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:
            public={k:f[k]for k in ('G','E0','S','Y','grid_shape','sigma','noise_radius','host')}
        Yobserved=public['Y'];Y=Yobserved
        if arm=='clean_native':
            with np.load(R2/'sealed_truth'/(CASE+'.npz'),allow_pickle=False)as f:Y=f['Y_clean']
            rec['data_change']='published stored clean measurement replaces noisy target; no regeneration'
        elif arm=='noisy_strong_release':
            rec['data_change']='NONE; same noisy publicY'
            rec['representation_change']='only finalmf10 current=[VD,Q],strongcoeff initial0; firsttwo stagesunchanged'
        else:raise ValueError(arm)
        model=MaxwellModel(public);m0=np.full(model.n,float(public['host']));initial=metric(model,Y,m0)
        cfg=dict(stages=STAGES,dpcs_rtol=1e-3,coefficient_tol=0.,miniter=50,material_inner_maxiter=60,
            regularizer=model.regularizer,material_box=(.02,20.),physical_sigma=float(np.linalg.norm(Y)),sourcewise_cg=True,
            data_weights=1/np.sum(abs(Yobserved)**2,axis=0),state_weights=1/np.sum(abs(model.E)**2,axis=0))
        with native_endpoints(directory)as calls:
            aux=(run_strong_release if arm=='noisy_strong_release'else run_cie)(model,Y,m0,cfg)
        scored=[initial]+[metric(model,Y,c['m'])for c in aux['candidates']];seed=min(scored,key=lambda c:c['objective'])
        dump(directory/'AUXILIARY.json',compact(aux))
        corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
        best=min([initial,seed,metric(model,Y,corr['m'])],key=lambda c:c['objective'])
        state=model.state(best['m']);value,g=model.data_value_gradient(best['m'],Y);pg=best['m']-np.clip(best['m']-g,.02,20.)
        prediction.parent.mkdir(exist_ok=True)
        np.savez_compressed(prediction,m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        dump(directory/'CORRECTOR.json',compact(corr))
        ends=[h for h in aux['history']if h.get('event')=='stage_end']
        rec.update(status='SEALED_DEVELOPMENT_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str(prediction.relative_to(ROOT)),physical_objective=float(value),
            physical_relative_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            residual_to_original_noisyY=float(np.linalg.norm(state['Y']-Yobserved)/np.linalg.norm(Yobserved)),
            physical_whitened_residual_original_sigma=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),projected_gradient_l2=float(np.linalg.norm(pg)),
            original_physical_linear_residual=state['state_relative_residual'],corrector_stop=corr['stop'],
            alternating_iterations=sum('iteration'in h for h in aux['history']),readout_calls=calls[0],stage_ends=ends,
            auxiliary_cost=aux['cost'],corrector_cost=corr['cost'],
            free_material_real_dof=model.n,physical_state_complex_entries=model.n*model.ns,
            fixed_DPCS_rank=aux['receiver_rank'],auxiliary_selected_physical_objective=seed['objective'],
            scientific_scope='mechanism development afterGT disclosure, not fresh confirmation or standalone basin proof')
    except Exception as exc:rec.update(status='FAILED_DEVELOPMENT_DISPOSITION',error=repr(exc),traceback=traceback.format_exc())
    finally:
        rec.update(physical_model_cost=dict(model.cost)if model else {},physical_operation_wall=dict(model.timings)if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings)if model else {},wall_seconds=time.perf_counter()-start,
            cpu_seconds=time.process_time()-cpu,FFT=0,GPU=0,
            timer_scope='public/cleandata load,aux1500steps,physicalreranking,corrector,finalgradient,serialization; finalreceiptwrite excluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['arm','status','physical_relative_residual','projected_gradient_inf','wall_seconds','cpu_seconds','physical_model_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):
        for arm in ('noisy_strong_release','clean_native'):run(arm)
