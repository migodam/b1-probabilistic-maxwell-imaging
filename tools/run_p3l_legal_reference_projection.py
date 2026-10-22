"""Pure decoder projection of fixed, UNQUALIFIED DEV Gaussian means.

No truth, observations, Target, physical backend or locked data are read.
The result is a prior-legal local action for the supplied wrong/unqualified
reference mean, not a Full-posterior Bayes action or data-qualified image.
"""
from pathlib import Path
import json
import shutil
import sys
import time
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,BudgetStop,write,save
from p1_material.prior import MaterialPrior
from p3l.estimators import project_legal_bayes_action

PARENTS=[f'DEV_{i:03d}' for i in range(32)]
DEST=OUT/'B/legal_reference_projection'
JOB='legal_reference_projection_DEV32'


def array(path,key,shape):
    with zipfile.ZipFile(path) as archive:
        if archive.getinfo(key+'.npy').file_size>128*1024*1024:
            raise MemoryError('array exceeds128MiB cap')
    with np.load(path,allow_pickle=False) as archive:value=archive[key].copy()
    if value.shape!=shape or not np.isfinite(value).all():raise ValueError('input layout/finite contract')
    return value


def diagnostics(prior,q,mean):
    material=prior.decode(q)[0]
    delta=material-mean
    gradient=2*np.real(prior.derivative(q).conj().T@delta)/576
    lower=np.r_[np.full(64,-np.inf),[0.,0.]]
    upper=np.r_[np.full(64,np.inf),[1.,1.]]
    projected=q-np.clip(q-gradient,lower,upper)
    return dict(material=material,gradient=gradient,projected_gradient=projected,
        gradient_inf=float(np.max(np.abs(gradient))),
        projected_gradient_inf=float(np.max(np.abs(projected))),
        supplied_mean_fit_MSE=float(np.mean(np.abs(delta)**2)),
        supplied_mean_fit_RMS=float(np.sqrt(np.mean(np.abs(delta)**2))))


def main():
    if DEST.exists():raise FileExistsError('projection evidence exists; no implicit overwrite/replay')
    protocol=dict(status='FROZEN_BEFORE_PROJECTION',parents=PARENTS,job=JOB,stage='B',
        analysis_only=True,max_wall_seconds=1200,
        input_mean='B/gaussian_images/{parent}.npz:gaussian_material_mean',
        fixed_start='A/anchors/{parent}.npz:q',starts=1,maximum_starts=1,
        max_nfev=256,tolerance=1e-9,posterior_latent_samples=None,
        objective='full-grid||original_decoder(q)-supplied_mean||²/576; no prior penalty',
        selection='one predeclared anchor start; retain start if optimizer fit worsens',
        original_prior='Gaussian64 and native uniform2 box; original MaterialPrior unchanged',
        input_qualification='WRONG_OR_UNQUALIFIED_FIXED_PROPERTY_GN_GAUSSIAN_REFERENCE_MEAN_NOT_FULL_POSTERIOR',
        output_qualification='PRIOR_LEGAL_LOCAL_REFERENCE_PROJECTION_NOT_GLOBAL_ARGMIN_NOT_DATA_QUALIFIED',
        GT_read=False,observations_read=False,locked_read=False,Full_RHS=0,physics_states=0,LU=0,GPU_seconds=0,
        scoring='NOT_RUN; parent may separately specify offline scoring and paid likelihood validation',
        locked_methods_changed=False)
    DEST.mkdir(parents=True)
    write(DEST/'FROZEN_PROTOCOL.json',protocol)
    shutil.copyfile(Path(__file__),DEST/'FROZEN_RUNNER.py')
    shutil.copyfile(ROOT/'src/p3l/estimators.py',DEST/'FROZEN_ESTIMATORS.py')
    book=Book('B',JOB,max_wall=1200,analysis_only=True)
    status='COMPLETE';error=None;rows=[];started=time.perf_counter()
    try:
        cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        geometry=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())
        prior=MaterialPrior(cfg,geometry['points'])
        write(DEST/'TARGET_CONFIG_SNAPSHOT.json',cfg)
        for parent in PARENTS:
            folder=DEST/parent;folder.mkdir();start=time.perf_counter()
            q0=mean=start_material=None
            try:
                q0=array(OUT/f'A/anchors/{parent}.npz','q',(66,))
                if not prior.legal(q0):raise ValueError('fixed anchor is outside original prior support')
                start_material=prior.decode(q0)[0]
                mean=array(OUT/f'B/gaussian_images/{parent}.npz','gaussian_material_mean',(576,))
                save(folder/'INPUTS.npz',q_start=q0,material_start=start_material,supplied_mean=mean)
                result=project_legal_bayes_action(mean,prior,q0,posterior_latent_samples=None,
                    maximum_starts=1,max_nfev=256,tolerance=1e-9)
                selected=np.asarray(result['q']).copy()
                if not prior.legal(selected):raise ValueError('projection returned illegal coordinates')
                check=diagnostics(prior,selected,mean)
                initial=diagnostics(prior,q0,mean)
                fallback=check['supplied_mean_fit_MSE']>initial['supplied_mean_fit_MSE']
                if fallback:selected=q0.copy();check=initial
                trials=result['trials']
                raw=dict(q_start=q0,material_start=start_material,supplied_mean=mean,
                    q=selected,material=check['material'],gradient=check['gradient'],
                    projected_gradient=check['projected_gradient'])
                for i,trial in enumerate(trials):
                    for key,value in trial.items():
                        if isinstance(value,np.ndarray):raw[f'trial{i}_{key}']=value
                save(folder/'RESULT.npz',**raw)
                row=dict(parent=parent,status='PRIOR_LEGAL_LOCAL_REFERENCE_PROJECTION',
                    prior_legal=True,input_reference_qualified=False,posterior_qualified=False,
                    global_argmin_certified=False,Full_likelihood_qualified='NOT_EVALUATED',
                    original_start_preserved_if_worse=bool(fallback or not result['selected_optimized_result']),
                    selected_optimized_result=bool(result['selected_optimized_result'] and not fallback),
                    optimization_success=bool(result['optimization_success'] and not fallback),
                    optimization_status=result['optimization_status'],
                    initial_supplied_mean_fit_MSE=initial['supplied_mean_fit_MSE'],
                    supplied_mean_fit_MSE=check['supplied_mean_fit_MSE'],
                    supplied_mean_fit_RMS=check['supplied_mean_fit_RMS'],
                    gradient_inf=check['gradient_inf'],projected_gradient_inf=check['projected_gradient_inf'],
                    trials=[{key:value for key,value in trial.items() if not isinstance(value,np.ndarray)} for trial in trials],
                    wall_seconds=time.perf_counter()-start,new_physical_RHS=0)
            except BudgetStop:raise
            except Exception as exc:
                # Never invent a legal material if the fixed input itself is
                # unavailable. Empty raw vectors explicitly preserve that case.
                save(folder/'RESULT.npz',q_start=np.array([]) if q0 is None else q0,
                    material_start=np.array([],complex) if start_material is None else start_material,
                    supplied_mean=np.array([],complex) if mean is None else mean,
                    q=np.array([]) if q0 is None else q0,
                    material=np.array([],complex) if start_material is None else start_material)
                row=dict(parent=parent,status='FAILED_PARENT_RETAINED_START_PRESERVED_WHEN_AVAILABLE',
                    error=repr(exc),prior_legal=False if q0 is None else prior.legal(q0),
                    posterior_qualified=False,Full_likelihood_qualified='NOT_EVALUATED',
                    wall_seconds=time.perf_counter()-start,new_physical_RHS=0)
            write(folder/'RESULT.json',row);rows.append(row)
            write(DEST/'PROGRESS.json',dict(rows=rows,expected_parents=PARENTS,
                all_failures_retained=True,new_physical_RHS=0,GT_read=False))
        summary=dict(status='COMPLETE_DEV32_REFERENCE_PROJECTION_NOT_FULL_POSTERIOR_BAYES_ACTION',
            parents=PARENTS,rows=rows,parent_count=32,
            projection_outputs=sum(r['status']=='PRIOR_LEGAL_LOCAL_REFERENCE_PROJECTION' for r in rows),
            optimizer_successes=sum(r.get('optimization_success',False) for r in rows),
            failures=sum(r['status']!='PRIOR_LEGAL_LOCAL_REFERENCE_PROJECTION' for r in rows),
            selected_start_fallbacks=sum(r.get('original_start_preserved_if_worse',False) for r in rows),
            wall_seconds=time.perf_counter()-started,GT_read=False,observations_read=False,locked_read=False,
            new_physical_RHS=0,new_physics_states=0,new_LU=0,GPU_seconds=0,analysis_only=True,
            input_reference='UNQUALIFIED_FIXED_PROPERTY_GN_GAUSSIAN_REFERENCE',
            output='PRIOR_LEGAL_LOCAL_PROJECTION_NOT_GLOBAL_ARGMIN',
            risk_evaluation='NOT_RUN',Full_likelihood_validation='NOT_RUN',locked_methods_changed=False)
        write(DEST/'SUMMARY.json',summary)
        lines=['# DEV32 supplied-reference legal projection: executed evidence',
            '', 'This analysis projects a supplied wrong/unqualified fixed-property GN Gaussian material mean onto the original decoder support. It does not supply a Full-posterior Bayes action, global optimum, GT-risk result or data/noise-qualified reconstruction.',
            '',f"Exactly32 predeclared DEV parents attempted; {summary['projection_outputs']} legal projection outputs, {summary['optimizer_successes']} optimizer successes, {summary['failures']} failures. One fixed anchor start each; max_nfev256, tolerance1e-9. Failed/iteration-capped outcomes remain visible.",
            '', '| Parent | Status | Initial fit MSE | Selected fit MSE | Projected gradient inf | Optimizer success |',
            '| --- | --- | --- | --- | --- | --- |']
        for row in rows:
            lines.append('| '+row['parent']+' | '+row['status']+' | '+str(row.get('initial_supplied_mean_fit_MSE'))+' | '+str(row.get('supplied_mean_fit_MSE'))+' | '+str(row.get('projected_gradient_inf'))+' | '+str(row.get('optimization_success',False))+' |')
        lines+=['','All reported fit losses are distance to the supplied reference mean, never GT risk. Raw q/material/start/mean and trial arrays are saved for every parent; failed input evidence uses explicit empty vectors rather than fabricated materials.',
            '', 'No truth, observations, locked arrays, physical solver or Target was read/executed. New RHS/state/LU/GPU are zero. Actual pure-analysis controlled wall is recorded in [the Book closure](../jobs/'+JOB+'/CLOSED.json). The parent owns any future isolated risk scoring and paid Full validation. Frozen locked methods were not changed.',
            '', '[Frozen protocol](FROZEN_PROTOCOL.json) · [Complete summary](SUMMARY.json)']
        (DEST/'EVIDENCE.md').write_text('\n'.join(lines)+'\n')
        print(json.dumps({key:summary[key] for key in ['status','parent_count','projection_outputs','optimizer_successes','failures','wall_seconds','new_physical_RHS']}))
    except BaseException as exc:
        status='ERROR';error=repr(exc)
        write(DEST/'FAILURE.json',dict(error=error,completed_parents=len(rows),all_failures_retained=True,new_physical_RHS=0))
        raise
    finally:book.close(status,error)


if __name__=='__main__':main()
