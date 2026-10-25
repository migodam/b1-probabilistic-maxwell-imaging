"""Finite-prior E4 approximation comparison; import performs no experiment.

No PDE calls, no continuous quadrature certificate. Offline projections use
exact posterior moments and must not be interpreted as deployable algorithms.
"""
from pathlib import Path
import argparse
import csv
import json
import time

import numpy as np

from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.finite_bayes import (PRIOR_ID, joint_coordinates, posterior_metrics,
    importance_sample, tempered_smc, update_logweights)


def normalize_log(logweights):
    logweights=np.asarray(logweights,float)
    if np.any(np.isnan(logweights)) or np.any(np.isposinf(logweights)):
        raise ValueError('Invalid log mass')
    maximum=logweights.max()
    if not np.isfinite(maximum):raise ValueError('No finite legal mass')
    weights=np.exp(logweights-maximum)
    return weights/weights.sum()


def normalized(weights):
    weights=np.asarray(weights,float)
    if weights.ndim!=1 or not np.isfinite(weights).all() or np.any(weights<0) or weights.sum()<=0:
        raise ValueError('Expected nonnegative finite mass vector')
    return weights/weights.sum()


def product_marginals(reference,coordinates):
    logs=np.zeros(len(reference))
    with np.errstate(divide='ignore'):
        for j in range(coordinates.shape[1]):
            _,inverse=np.unique(coordinates[:,j],return_inverse=True)
            mass=np.bincount(inverse,weights=reference)
            logs+=np.log(mass[inverse])
    return normalize_log(logs)


def moment_kernel(points,weights):
    """Gaussian log kernel with exact weighted moments; strict SPD, no repair."""
    weights=normalized(weights)
    mean=weights@points
    centered=points-mean
    covariance=(centered*weights[:,None]).T@centered
    covariance=(covariance+covariance.T)/2
    eigenvalues=np.linalg.eigvalsh(covariance)
    if not np.isfinite(eigenvalues).all() or eigenvalues[0]<=0:
        raise ValueError('NOT_DEFINED: moment covariance is not strictly SPD; no jitter/pseudoinverse')
    try:
        chol=np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError('NOT_DEFINED: moment covariance Cholesky failed; no jitter/pseudoinverse') from exc
    scaled=np.linalg.solve(chol,centered.T)
    log_kernel=-.5*np.sum(scaled*scaled,axis=0)-np.log(np.diag(chol)).sum()
    return log_kernel,dict(mean=mean,covariance=covariance,eigenvalues=eigenvalues,
                           covariance_repair=False,parameters='EXACT_REFERENCE_OFFLINE_PROJECTION_ONLY')


def gaussian_projection(reference,grid,coordinates):
    logs,parameters=moment_kernel(coordinates,reference)
    with np.errstate(divide='ignore'):
        weights=normalize_log(logs+np.log(grid['prior_weight']))
    return weights,dict(label='FINITE_PRIOR_WEIGHTED_10D_MOMENT_GAUSSIAN_KERNEL_OFFLINE_PROJECTION_ONLY',
        parameters=parameters,normalization='Gaussian kernel times finite prior mass; legal atoms only',
        moment_match_after_discrete_truncation_not_guaranteed=True,deployable=False)


def occupancy_mixture(reference,grid):
    result=np.zeros(len(reference));parameters=[]
    ids=np.asarray(grid['support_id'])
    for sid in range(len(grid['legal_supports'])):
        mask=ids==sid
        mass=float(reference[mask].sum())
        if mass==0:
            parameters.append(dict(support_id=sid,mass=0.,status='ZERO_REFERENCE_MASS'))
            continue
        support=np.asarray(grid['legal_supports'][sid])
        pose=np.asarray(grid['pose'])[mask]
        # The empty branch has no identifiable wall material coordinates;
        # unused material latents are integrated, never jittered into evidence.
        points=np.column_stack((grid['epsinf'][mask],grid['ell'][mask],pose)) if support.any() else pose
        logs,moments=moment_kernel(points,reference[mask]/mass)
        base=np.asarray(grid['prior_weight'])[mask]
        with np.errstate(divide='ignore'):
            result[mask]=mass*normalize_log(logs+np.log(base/base.sum()))
        parameters.append(dict(support_id=sid,mass=mass,
            coordinates=['epsinf','ell','px','py'] if support.any() else ['px','py'],**moments))
    return normalized(result),dict(label='OCCUPANCY_MIXTURE_CONDITIONAL_MOMENT_KERNEL_OFFLINE_PROJECTION_ONLY',
        parameters=parameters,deployable=False,empty_wall_material_integrated=True,
        normalization='Exact support mass times finite conditional-prior-weighted Gaussian kernel')


def resolve_likelihood(reference,prior,history):
    """Disambiguate incremental versus cumulative arrays by exact identity.

    Never add cumulative likelihoods repeatedly. The selected interpretation
    must reproduce the stored reference; no posterior/prior density-ratio
    substitution or likelihood tempering is used to rescue an inconsistent file.
    """
    choices={'cumulative_final_entry':history[-1],
             'incremental_entries_summed_once':np.sum(history,axis=0)}
    matches=[];errors={}
    for name,ll in choices.items():
        candidate=np.exp(update_logweights(prior,[('complete_final_data',ll)]))
        error=float(np.max(np.abs(candidate-reference)))
        errors[name]=error
        if error<=1e-10:matches.append(name)
    if not matches:raise ValueError('Stored likelihood history does not reproduce final exact posterior: '+str(errors))
    name=matches[0]
    return choices[name],dict(chosen=name,also_matching=matches[1:],normalization_errors=errors,
        tolerance=1e-10,tolerance_scope='File algebra consistency only, not scientific accuracy gate')


def write_json(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def persist(out,report,weight_rows):
    write_json(out/'RESULTS.json',report)
    columns=['parent_id','repeat_id','method','particles','seed','status','coordinate_cdf_max',
        'target_cdf_max','coordinate_cdf_within_002','event_error','important_mode_mass_error','action_regret','task_tolerance',
        'cdf_within_002','event_within_002','mode_within_002','task_within_tolerance',
        'all_declared_tolerances_met','wall_seconds','process_cpu_seconds','error']
    temporary=out/'FULL_TABLE.csv.tmp'
    with temporary.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(report['rows'])
    temporary.replace(out/'FULL_TABLE.csv')
    temporary=out/'WEIGHTS.tmp.npz'
    if weight_rows:
        np.savez_compressed(temporary,weights=np.stack(weight_rows),
            parent_id=np.array([r['parent_id'] for r in report['rows']]),
            repeat_id=np.array([r['repeat_id'] for r in report['rows']]),
            method=np.array([r['method'] for r in report['rows']]),
            status=np.array([r['status'] for r in report['rows']]),
            prior_id=np.array(PRIOR_ID),undefined_weights_are_nan=np.array(True))
        temporary.replace(out/'WEIGHTS.npz')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();root=args.root.resolve()
    bankroot=root/'results/e2e3/finite_bank';online=bankroot/'online'
    with np.load(online/'FINITE_GRID.npz',allow_pickle=False) as raw:grid={k:raw[k].copy() for k in raw.files}
    with np.load(online/'FINITE_BANK.npz',allow_pickle=False) as raw:
        completed,total=int(raw['completed']),int(raw['total'])
        bank_prior=str(raw['prior_id'])
    if completed!=total or bank_prior!=PRIOR_ID:raise RuntimeError('Complete matching finite bank required')
    with np.load(online/'OBSERVATIONS.npz',allow_pickle=False) as raw:
        observed_parent=raw['parent_id'].copy();observed_repeat=raw['repeat_id'].copy()
    with np.load(root/'results/e2e3/inference/POSTERIORS.npz',allow_pickle=False) as raw:
        if not bool(raw['complete']) or not np.all(raw['completed_mask']):
            raise RuntimeError('Incomplete posterior checkpoint cannot support E4')
        exact=raw['weights'].copy();likelihoods=raw['log_likelihoods'].copy()
        parent_ids=raw['parent_id'].copy();repeat_ids=raw['repeat_id'].copy()
    if exact.ndim!=4 or exact.shape!=likelihoods.shape or exact.shape[-1]!=len(grid['prior_weight']):
        raise ValueError('Expected (parent,repeat,time,finite_atom) posterior and likelihood arrays')
    if not np.array_equal(parent_ids,observed_parent) or not np.array_equal(repeat_ids,observed_repeat):
        raise ValueError('Inference/observation parent or repeat identities differ')
    coordinates=joint_coordinates(grid)
    out=root/'results/e4/comparison';out.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(out/'REFERENCE_FINAL.npz',weights=exact[:,:,-1],parent_id=parent_ids,
                        repeat_id=repeat_ids,prior_weight=grid['prior_weight'],prior_id=np.array(PRIOR_ID))
    book=BudgetBook(root/'results/budget_ledger.json')
    setup_path=online/'JOB_COST_FINAL.json'
    setup=json.loads(setup_path.read_text())
    if setup.get('complete') is not True:raise RuntimeError('Final measured bank setup cost required')
    report=dict(scope='EXPLICIT_FINITE_PRIOR_APPROXIMATION_COMPARISON',status='RUNNING',rows=[],
        prior_id=PRIOR_ID,prior_definition='3x3 material nodes and 5x5 target-pose nodes DEFINE the finite prior; not an accuracy estimate for a continuous prior',
        continuous_quadrature_certificate=False,global_posterior_certificate=False,
        shared_bank_setup_attribution=setup,shared_bank_source_stage='E2E3',
        e4_physical_rhs=0,e4_physical_factorizations=0,
        setup_cost_warning='Shared bank is not free: construction charged once in E2E3; independent deployment must account for that setup',
        projection_warning='Exact-reference moment and marginal parameters are OFFLINE_PROJECTION_ONLY, not deployable discovery/inference',
        important_mode_rule='Existing posterior_metrics: reference support mass >= .05',
        declared_tolerances=dict(target_cdf=.02,event=.02,important_mode_mass=.02,task_fraction=.05,task_floor=.005**2),
        scientific_decision='PARENT_REVIEW_REQUIRED',budget_start=book.totals('E4'))
    weights_out=[];stopped=False
    methods=[('MAP_POINT',None,None),('PRODUCT_10_COORDINATE_MARGINALS',None,None),
        ('OCCUPANCY_MIXTURE_OFFLINE_PROJECTION',None,None),('MOMENT_GAUSSIAN_OFFLINE_PROJECTION',None,None)]
    methods += [(method,n,seed) for method in ('TEMPERED_SMC','PRIOR_IMPORTANCE') for n in (256,1024) for seed in (0,1,2)]
    for ip,pid in enumerate(parent_ids):
        for ir,rid in enumerate(repeat_ids):
            reference=normalized(exact[ip,ir,-1])
            try:
                ll,ll_audit=resolve_likelihood(reference,grid['prior_weight'],likelihoods[ip,ir])
            except Exception as exc:
                ll=None;ll_audit=dict(error=str(exc))
            for method,n,seed in methods:
                row=dict(parent_id=int(pid),repeat_id=int(rid),method=method,particles=n,seed=seed,
                    status='NOT_RUN',likelihood_audit=ll_audit)
                weight=np.full(len(reference),np.nan)
                start=time.perf_counter();cpu=time.process_time()
                try:
                    if stopped:raise BudgetExceeded('E4 stage stopped; remaining comparisons NOT_RUN')
                    if book.snapshot('E4')['E4']['remaining']['process_cpu_seconds']<=0:
                        raise BudgetExceeded('E4 shared CPU budget exhausted')
                    metadata={}
                    if method=='MAP_POINT':
                        weight=np.zeros_like(reference);weight[int(np.argmax(reference))]=1.
                    elif method=='PRODUCT_10_COORDINATE_MARGINALS':
                        weight=product_marginals(reference,coordinates)
                        metadata=dict(label='PRODUCT_COORDINATE_MASS_RESTRICTED_TO_LEGAL_FINITE_ATOMS',
                            independent_voxel_full_posterior=False,deployable=False,parameters='OFFLINE_PROJECTION_ONLY')
                    elif method=='OCCUPANCY_MIXTURE_OFFLINE_PROJECTION':
                        weight,metadata=occupancy_mixture(reference,grid)
                    elif method=='MOMENT_GAUSSIAN_OFFLINE_PROJECTION':
                        weight,metadata=gaussian_projection(reference,grid,coordinates)
                    else:
                        if ll is None:raise ValueError('Cannot sample: inconsistent or unresolved likelihood history')
                        sampler=tempered_smc if method=='TEMPERED_SMC' else importance_sample
                        sampled=sampler(grid,ll,n_particles=n,seed=seed)
                        weight=sampled['weights']
                        metadata={k:v for k,v in sampled.items() if k!='weights'}
                    metrics=posterior_metrics(weight,reference,grid)
                    target_cdf=float(np.max(metrics['coordinate_cdf_sup'][-2:]))
                    event=metrics['event_px_ge_0_error'];mode=metrics['mode_mass_max_error_ref_ge_005']
                    tolerance=.05*metrics['reference_bayes_risk']+.005**2
                    checks=[target_cdf<=.02,event<=.02,mode<=.02,metrics['action_regret']<=tolerance]
                    row.update(status='EVALUATED',metadata=metadata,metrics=metrics,
                        coordinate_cdf_max=float(np.max(metrics['coordinate_cdf_sup'])),target_cdf_max=target_cdf,
                        coordinate_cdf_within_002=bool(np.max(metrics['coordinate_cdf_sup'])<=.02),
                        event_error=event,important_mode_mass_error=mode,action_regret=metrics['action_regret'],
                        task_tolerance=tolerance,cdf_within_002=checks[0],event_within_002=checks[1],
                        mode_within_002=checks[2],task_within_tolerance=checks[3],all_declared_tolerances_met=all(checks))
                except BudgetExceeded as exc:
                    stopped=True;row.update(status='NOT_RUN',error=str(exc))
                except Exception as exc:
                    row.update(status='NOT_DEFINED' if 'NOT_DEFINED' in str(exc) else 'FAILED',
                        error_type=type(exc).__name__,error=str(exc))
                    weight=np.full(len(reference),np.nan)
                row.update(wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu)
                report['rows'].append(row);weights_out.append(weight);persist(out,report,weights_out)
    report['status']='STOPPED_E4_BUDGET' if stopped else 'FINITE_COMPARISONS_COMPLETE_PARENT_REVIEW_REQUIRED'
    report['budget_end']=book.totals('E4')
    report['all_row_counts']={status:sum(r['status']==status for r in report['rows']) for status in sorted({r['status'] for r in report['rows']})}
    report['method_totals_including_failures']={method:dict(
        planned=sum(r['method']==method for r in report['rows']),
        evaluated=sum(r['method']==method and r['status']=='EVALUATED' for r in report['rows']),
        within_all_tolerances=sum(r['method']==method and r.get('all_declared_tolerances_met',False) for r in report['rows']),
        failed_undefined_or_unrun=sum(r['method']==method and r['status']!='EVALUATED' for r in report['rows'])) for method,_,_ in methods}
    persist(out,report,weights_out)
    write_json(out/'gate.json',dict(decision='INCONCLUSIVE',execution_status=report['status'],
        scope='Declared finite-prior approximation errors only; no continuous or global certificate',
        scientific_review='Parent owns interpretation; all failed and undefined comparisons retained'))


if __name__=='__main__':
    main()
