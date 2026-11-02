"""Single POST_HOC_DEVELOPMENT native mf12 control; last inverse run.

Case identity selected after GT release by parent; this is not new blind data.
No GT is opened here. Historical/native source and all old artifacts untouched.
"""
from pathlib import Path
from datetime import datetime,timezone
import json,sys,time,traceback
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO));sys.path.insert(0,str(ROOT/'experiments'))
from run_blind import native_endpoints,metric,compact,dump,ALLOW
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import run_cie
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector

CASE='obs_b20d5a55677d4eca'
STAGES=[dict(mf=3,beta=1.,maxiter=500),
        dict(mf=6,beta_current=1.,beta_material=.2,maxiter=500),
        dict(mf=12,beta_current=1.,beta_material=.2,maxiter=500)]

def main():
    directory=ROOT/'results/capacity_control'/CASE;directory.mkdir(parents=True,exist_ok=True)
    path=directory/'RECEIPT.json';prediction=ROOT/'predictions_capacity'/(CASE+'_native_mf12.npz')
    if path.exists()or prediction.exists():raise FileExistsError('Single control already disposed; no rerun')
    start=time.perf_counter();cpu=time.process_time();model=None
    rec=dict(status='RUNNING',classification='POST_HOC_DEVELOPMENT_AFTER_GT_RELEASE',case_id=CASE,
        started_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE_BY_THIS_RUNNER',
        selection='Parent selected previously failed case after GT release; not a new blind success test',
        frozen_change='Native Fourier stages3,6,12 rather than3,6,10; no beta/iteration/prior/init/corrector change',
        code='experiments/capacity_control.py',hash_checks=False,stages=STAGES,coefficient_tol=0.,miniter=50,
        dpcs_rtol=1e-3,material_inner_maxiter=60,initialization='same uniform public host .1',
        corrector=dict(method='original physical L-BFGS-B',maxiter=1000,maxfun=4000,wallcap=None),
        objective='same original normalized physical data + historical smooth/TV prior',
        threads=1,device='CPU',inverse_runs_authorized=1,last_inverse_run=True,
        scope='Tests available full APCS capacity under fixed measured DPCS; more matrix work is paid and not iteration-cost matched')
    dump(path,rec)
    try:
        with np.load(ROOT/'data_public'/(CASE+'.npz'),allow_pickle=False)as f:
            if not set(f.files)<=ALLOW:raise ValueError('Unexpected public key')
            public={k:f[k]for k in f.files}
        model=MaxwellModel(public);Y=public['Y'];m0=np.full(model.n,float(public['host']))
        initial=metric(model,Y,m0)
        cfg=dict(stages=STAGES,dpcs_rtol=1e-3,coefficient_tol=0.,miniter=50,
            material_inner_maxiter=60,regularizer=model.regularizer,material_box=(.02,20.),
            physical_sigma=float(np.linalg.norm(Y)),sourcewise_cg=True)
        with native_endpoints(directory)as calls:aux=run_cie(model,Y,m0,cfg)
        scored=[initial]+[metric(model,Y,c['m'])for c in aux['candidates']]
        seed=min(scored,key=lambda c:c['objective'])
        dump(directory/'AUXILIARY.json',compact(aux))
        corr=corrector(model,Y,seed['m'],dict(method='lbfgs',max_nfev=1000,box=(.02,20.)))
        candidates=[initial,seed,metric(model,Y,corr['m'])];best=min(candidates,key=lambda c:c['objective'])
        state=model.state(best['m']);value,g=model.data_value_gradient(best['m'],Y)
        pg=best['m']-np.clip(best['m']-g,.02,20.)
        prediction.parent.mkdir(exist_ok=True)
        np.savez_compressed(prediction,m=best['m'],Y_prediction=state['Y'],j=state['j'],Q=state['Q'],gradient=g,projected_gradient=pg)
        dump(directory/'CORRECTOR.json',compact(corr))
        ends=[h for h in aux['history']if h.get('event')=='stage_end']
        capacities=[]
        for h in ends:
            capacities.append(dict(stage=h['stage'],mf=h['mf'],Fourier_complex_coefficients=h['fourier_columns']*model.ns,
                Fourier_real_coefficients=2*h['fourier_columns']*model.ns,
                independent_APCS_spatial_complex_rank=h['apcs_rank'],
                independent_free_current_real_dof=2*h['apcs_rank']*model.ns,
                fixed_DPCS_spatial_complex_rank=aux['receiver_rank'],
                fixed_DPCS_data_determined_real_coordinates=2*aux['receiver_rank']*model.ns,
                total_DPCS_APCS_linear_real_dof=2*(aux['receiver_rank']+h['apcs_rank'])*model.ns,
                material_real_dof=model.n))
        rec.update(status='SEALED_POST_HOC_PREDICTION',finished_utc=datetime.now(timezone.utc).isoformat(),
            prediction=str(prediction.relative_to(ROOT)),physical_objective=float(value),
            physical_relative_residual=best['relative_residual'],physical_whitened_residual=float(np.linalg.norm(state['Y']-Y)/float(public['sigma'])),
            noise_radius99=float(public['noise_radius']),physical_linear_residual=state['state_relative_residual'],
            projected_gradient_mapping_l2=float(np.linalg.norm(pg)),projected_gradient_mapping_inf=float(np.linalg.norm(pg,np.inf)),
            active_lower_exact=int(np.count_nonzero(best['m']==.02)),active_upper_exact=int(np.count_nonzero(best['m']==20.)),
            corrector_stop=corr['stop'],all_alternating_iterations=sum('iteration'in h for h in aux['history']),
            native_forwarding_readout_calls=calls[0],stage_capacity=capacities,
            auxiliary_cost=aux['cost'],corrector_cost=corr['cost'],auxiliary_status=aux['status'],
            corrector_initial_original_objective=seed['objective'],
            limitations='Full APCS still fixes DPCS coefficients from noisy receiver inversion; no GT recovery claim before separate scoring; more matrix work than mf10')
    except Exception as exc:rec.update(status='FAILED_POST_HOC_DISPOSITION',error=repr(exc),traceback=traceback.format_exc())
    finally:
        rec.update(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
            physical_model_cost=dict(model.cost)if model else {},
            physical_operation_wall=dict(model.timings)if model else {},
            physical_operation_cpu=dict(model._physics.cpu_timings)if model else {},
            FFT=0,hashes_computed=0,
            timer_scope='public load+native preparation/1500steps+physical candidates+corrector+final gradient+serialization; receipt final write excluded')
        dump(path,rec)
    print(json.dumps({k:rec.get(k)for k in ['status','physical_relative_residual','physical_whitened_residual','projected_gradient_mapping_inf','wall_seconds','cpu_seconds','physical_model_cost','error']}),flush=True)

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
