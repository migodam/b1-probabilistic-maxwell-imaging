"""Exact finite-prior wall map and explicitly approximate freeze/reopen controls.

No PDE and no truth reads in inference. Offline scoring starts only after every
online posterior/action has been durably written. Not continuous quadrature.
"""
from pathlib import Path
import json,csv,time
import numpy as np
from scipy.special import logsumexp
from scipy.stats import chi2
from wslam.accounting import BudgetBook,json_default
from wslam.finite_bayes import (summarize,update_logweights,conditional_target_only,
    posterior_metrics,joint_coordinates,PRIOR_ID)
from wslam.factors import FactorRegistry,refresh_logweights
from wslam.information import local_gaussian_audit
from wslam.physics import acquisition

ROOT=Path(__file__).resolve().parents[1]
METHODS=['HARD_FROZEN_WALL','ITERATIVE_LAYER_POINT','BATCH_JOINT_MAP',
 'INCREMENTAL_JOINT_SMOOTHER','PERIODIC_GLOBAL_REFRESH','SOFT_FREEZE',
 'SPECTRAL_ONLY','INNOVATION_ONLY','TASK_WEIGHTED','QUALIFIED_REOPEN']

def jsave(p,obj):
    p.write_text(json.dumps(obj,indent=2,default=json_default,allow_nan=False)+'\n')

def normalized(logw):return np.exp(logw-logsumexp(logw))

def delta(n,i):
    w=np.zeros(n);w[i]=1.;return w

def likelihood(y,means):
    difference=(means-y)/.01
    return -.5*np.sum(difference.real**2+difference.imag**2,axis=-1)

def gaussian_reopen(old,coords,means,y,reference_tau):
    """Moment regression surrogate; D6-8 valid only for this named Gaussian.

    Uses cached candidate responses available to all arms. It is not the local
    nonlinear Maxwell Jacobian and is not a true-posterior certificate.
    """
    scale=np.r_[np.full(4,5.),np.full(4,.25),2.,1.]
    x=coords/scale;xm=old@x;dx=x-xm
    cov=(dx*old[:,None]).T@dx
    pred=np.column_stack((means.real,means.imag))/.01
    ym=old@pred;dy=pred-ym
    eigen=np.linalg.eigvalsh((cov+cov.T)/2)
    if eigen[0]<=max(eigen[-1],1.)*1e-13:
        return dict(defined=False,reason='singular/concentrated moment covariance; conservative reopen',
            covariance_eigenvalues=eigen,spectral_reopen=True,innovation_reopen=True,task_reopen=True)
    prior=np.linalg.solve(cov,np.eye(10))
    j=np.linalg.solve(cov,(dx*old[:,None]).T@dy).T
    residual=np.r_[y.real,y.imag]/.01-ym
    h=prior+j.T@j;g=j.T@residual
    task=np.zeros((2,10));task[:,8:]=np.diag(scale[8:])
    audit=local_gaussian_audit(h,g,task,np.arange(8))
    innovation_cov=np.eye(len(ym))+(dy*old[:,None]).T@dy
    innovation=float(residual@np.linalg.solve(innovation_cov,residual))
    # Compare marginal wall information, eliminating the same shared pose.
    old_wall=np.linalg.solve(cov[:8,:8],np.eye(8))
    new_wall=h[:8,:8]-h[:8,8:]@np.linalg.solve(h[8:,8:],h[8:,:8])
    spectral=float(np.trace(cov[:8,:8]@(new_wall-old_wall)))
    heuristic_tau=float(reference_tau)
    task_bound=float(audit['task_regret']+np.trace(audit['omitted_task_covariance']))
    return dict(defined=True,scope='MOMENT_REGRESSION_GAUSSIAN_ONLY',
        true_posterior_certificate=False,covariance_eigenvalues=eigen,
        spectral_trace=spectral,spectral_reopen=spectral>.1*8,
        innovation=innovation,innovation_threshold=float(chi2.ppf(.99,len(ym))),
        innovation_reopen=innovation>chi2.ppf(.99,len(ym)),
        heuristic_tau=heuristic_tau,task_regret=audit['task_regret'],
        omitted_task_covariance=audit['omitted_task_covariance'],task_bound=task_bound,
        task_reopen=task_bound>heuristic_tau,
        dual_interval=audit['dual_interval'],dual_residuals=audit['dual_residuals'],
        lower_bound=audit['lower_bound_source'])

def main():
    bankroot=ROOT/'results/e2e3/finite_bank';online=bankroot/'online'
    independent_check=json.loads((ROOT/'results/e2e3/bank_validation/RESULTS.json').read_text())
    if independent_check['decision']!='PASS':raise RuntimeError('Finite-bank fixed independent-grid audit failed')
    metadata=json.loads((bankroot/'COMPLETE.json').read_text())
    if metadata.get('status')!='COMPLETE':raise RuntimeError('Missing terminal physical-bank receipt')
    with np.load(online/'FINITE_GRID.npz',allow_pickle=False) as z:grid={k:z[k] for k in z.files}
    grid.update(label=PRIOR_ID,discrete_prior=True,continuous_uniform_integration_certificate=False)
    with np.load(online/'FINITE_BANK.npz',allow_pickle=False) as z:
        if int(z['completed'])!=2950 or int(z['total'])!=2950 or str(z['prior_id'])!=PRIOR_ID:
            raise RuntimeError('Incomplete or wrong-prior physical reference')
        means=z['data'].copy();backend_residual=float(z['residuals'].max())
        tx,rx=acquisition(outer=True)
        if not np.array_equal(z['tx'],tx) or not np.array_equal(z['rx'],rx) or float(z['frequency'])!=1.:
            raise RuntimeError('Bank acquisition or frequency mismatch')
    with np.load(online/'OBSERVATIONS.npz',allow_pickle=False) as z:
        observed=z['data'].copy();identifiers=z['sequential_factor_id'].copy()
        if complex(z['known_gamma'])!=1 or float(z['sigma_real_imag'])!=.01:
            raise RuntimeError('Finite-prior calibration/noise mismatch')
    protocol=json.loads((ROOT/'FINITE_PILOT_SUPPLEMENT.json').read_text())
    out=ROOT/'results/e2e3/inference';out.mkdir(parents=True,exist_ok=False)
    prior=grid['prior_weight'];coords=joint_coordinates(grid);N=len(prior)
    if (means.shape!=(2950,4,16) or observed.shape!=(12,2,4,16) or identifiers.shape!=(12,2,3)
            or grid['support'].shape!=(2950,4) or grid['pose'].shape!=(2950,2)
            or len(grid['legal_supports'])!=14 or N!=2950 or np.any(prior<=0)
            or not np.isclose(prior.sum(),1) or not np.isfinite(means).all()
            or not np.isfinite(observed).all() or len(np.unique(identifiers))!=72):
        raise RuntimeError('Incomplete or invalid finite physical panel/prior/factors')
    shape=observed.shape[:2]
    references=np.zeros((*shape,3,N));lls=np.zeros_like(references)
    completed_mask=np.zeros((*shape,3),bool)
    allweights=np.zeros((*shape,3,len(METHODS),N));logs=[];summaries=[]
    seq=[([0],[8]),([1],[8]),([0,1],[i for i in range(16) if i!=8])]
    worst_refresh=0.;likelihood_work=0
    for parent in range(shape[0]):
        for repeat in range(shape[1]):
            registry=FactorRegistry();registry.register('prior','prior',payload=PRIOR_ID)
            refs=[];old_methods={};unique={};incremental_log=np.log(prior)
            for t,(sources,receivers) in enumerate(seq):
                newmeans=means[:,sources][:,:,receivers].reshape(N,-1)
                y=observed[parent,repeat][sources][:,receivers].ravel()
                fid=str(identifiers[parent,repeat,t])
                start=time.process_time();ll=likelihood(y,newmeans);likelihood_cpu=time.process_time()-start
                likelihood_work+=N*len(y)
                lls[parent,repeat,t]=ll
                registry.register(fid,'likelihood',observation_id=fid,payload=ll)
                unique[fid]=ll
                start=time.process_time();logref=update_logweights(prior,unique);ref=np.exp(logref)
                batch_cpu=time.process_time()-start
                references[parent,repeat,t]=ref;refs.append(ref)
                # Computational refresh replaces this factor; it is no new data.
                registry.register(fid,'likelihood',version=2,observation_id=fid,payload=ll)
                refreshed=refresh_logweights(logref,ll,ll)
                worst_refresh=max(worst_refresh,float(np.max(np.abs(np.exp(refreshed)-ref))))
                if t==0:
                    early=ref.copy();wall_mass=np.bincount(grid['wall_group_id'],weights=early)
                    frozen_wall=int(np.argmax(wall_mass))
                for mi,method in enumerate(METHODS):
                    started=time.process_time();action='KEEP';criterion=None
                    old=old_methods.get(method,early)
                    if method=='INCREMENTAL_JOINT_SMOOTHER':
                        # Exact update of the retained joint state, not repeated evidence.
                        incremental_log=incremental_log+ll
                        incremental_log=incremental_log-logsumexp(incremental_log)
                        w=np.exp(incremental_log)
                        action='GLOBAL_RELINEARIZE' if t else 'KEEP'
                    elif method=='BATCH_JOINT_MAP':
                        w=delta(N,int(np.argmax(logref)));action='GLOBAL_RELINEARIZE'
                    elif method=='HARD_FROZEN_WALL':
                        selected=grid['wall_group_id']==frozen_wall
                        masked=np.where(selected,logref,-np.inf);w=normalized(masked)
                        action='HARD_FREEZE_WALL'
                    elif method=='ITERATIVE_LAYER_POINT':
                        # Ordinary revisitable block MAP on this finite class:
                        # two wall-then-pose sweeps, retaining all old factors.
                        point=int(np.argmax(prior)) if t==0 else int(np.argmax(old))
                        for sweep in range(2):
                            same_pose=np.all(grid['pose']==grid['pose'][point],axis=1)
                            point=int(np.argmax(np.where(same_pose,logref,-np.inf)))
                            same_wall=grid['wall_group_id']==grid['wall_group_id'][point]
                            point=int(np.argmax(np.where(same_wall,logref,-np.inf)))
                        w=delta(N,point);action='ORDINARY_BLOCK_MAP_WALL_THEN_POSE'
                    elif t==0:
                        w=ref.copy();action='SOFT_FREEZE'
                    else:
                        stale=conditional_target_only(old,ll,grid['wall_group_id'])['weights']
                        if method=='SOFT_FREEZE':w=stale;action='SOFT_FREEZE'
                        elif method=='PERIODIC_GLOBAL_REFRESH':
                            reopen=(t%2==1);w=ref.copy() if reopen else stale
                            action='GLOBAL_RELINEARIZE' if reopen else 'SOFT_FREEZE'
                        else:
                            pose=grid['pose'];mean_pose=ref@pose
                            reference_risk=float(ref@np.sum((pose-mean_pose)**2,axis=1))
                            criterion=gaussian_reopen(old,coords,newmeans,y,.05*reference_risk+.005**2)
                            if method=='TASK_WEIGHTED':
                                criterion['reference_threshold_requires_full_update']=True
                                criterion['reference_cpu_seconds_attribution']=batch_cpu+likelihood_cpu
                            key={'SPECTRAL_ONLY':'spectral_reopen','INNOVATION_ONLY':'innovation_reopen',
                                 'TASK_WEIGHTED':'task_reopen','QUALIFIED_REOPEN':'task_reopen'}[method]
                            reopen=bool(criterion[key])
                            if method=='QUALIFIED_REOPEN':
                                # Full finite reference audit is explicitly paid computation.
                                check=posterior_metrics(stale,ref,grid)
                                tau=.05*check['reference_bayes_risk']+.005**2
                                reopen|=(max(check['coordinate_cdf_sup'][8:])>.02 or
                                    check['mode_mass_max_error_ref_ge_005']>.02 or
                                    check['event_px_ge_0_error']>.02 or check['action_regret']>tau)
                                criterion['full_reference_audit']=check
                                criterion['audit_reference_cpu_seconds']=batch_cpu+likelihood_cpu
                            w=ref.copy() if reopen else stale
                            action='REOPEN' if reopen else 'SOFT_FREEZE'
                    cpu=time.process_time()-started
                    allweights[parent,repeat,t,mi]=w;old_methods[method]=w
                    metrics=posterior_metrics(w,ref,grid)
                    tau=.05*metrics['reference_bayes_risk']+.005**2
                    probability_ok=bool(max(metrics['coordinate_cdf_sup'][8:])<=.02 and
                         metrics['event_px_ge_0_error']<=.02 and metrics['mode_mass_max_error_ref_ge_005']<=.02)
                    logs.append(dict(parent=parent,repeat=repeat,time=t+1,method=method,action=action,
                        reopened_regions=[0,1,2,3] if action in ('REOPEN','GLOBAL_RELINEARIZE') else [],
                        posterior_metrics=metrics,probability_tolerance_met=probability_ok,
                        task_tolerance=tau,task_tolerance_met=bool(metrics['action_regret']<=tau),
                        cpu_seconds=cpu,shared_likelihood_cpu_seconds=likelihood_cpu,
                        standalone_cpu_seconds=cpu+likelihood_cpu+(batch_cpu if method in
                            ('BATCH_JOINT_MAP','HARD_FROZEN_WALL','ITERATIVE_LAYER_POINT',
                             'PERIODIC_GLOBAL_REFRESH','TASK_WEIGHTED','QUALIFIED_REOPEN') else 0.),
                        shared_batch_reference_cpu_seconds=batch_cpu,
                        actual_new_maxwell_rhs=0,setup_rhs_attribution=11800,
                        criterion=criterion,factor_ids=list(unique),old_factor_refresh_is_new_data=False))
                summary=summarize(ref,grid)
                completed_mask[parent,repeat,t]=True
                summaries.append(dict(parent=parent,repeat=repeat,time=t+1,reference=summary,
                    frozen_wall_group=frozen_wall,factor_count=len(unique),registry_factor_count=len(registry)))
            np.savez_compressed(out/'POSTERIORS.npz',weights=references,log_likelihoods=lls,
                method_weights=allweights,methods=np.asarray(METHODS),completed_parent=parent,
                completed_mask=completed_mask,complete=bool(completed_mask.all()),
                parent_id=np.arange(shape[0]),repeat_id=np.arange(shape[1]))
            jsave(out/'METHOD_RESULTS.json',logs);jsave(out/'BELIEFS.json',summaries)
            print(json.dumps(dict(parent=parent,repeat=repeat,status='ONLINE_SAVED')),flush=True)
    # All decisions are now saved. Truth is only an offline scoring input here.
    offline=np.load(bankroot/'offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz',allow_pickle=False)
    truth=offline['truth_index'];score=[]
    for record in logs:
        p,r,t,mi=record['parent'],record['repeat'],record['time']-1,METHODS.index(record['method'])
        w=allweights[p,r,t,mi];s=summarize(w,grid);gt=int(truth[p])
        error=float(np.linalg.norm(s['pose_mean']-grid['pose'][gt]))
        posecdf=s['coordinate_distributions'][8:]
        covered=all(d['quantiles'][0]<=v<=d['quantiles'][-1] for d,v in zip(posecdf,grid['pose'][gt]))
        score.append(dict(parent=p,repeat=r,time=t+1,method=record['method'],pose_error=error,
            truth_mode_mass=float(s['mode_probability'][grid['support_id'][gt]]),
            truth_wall_group_mass=float(w[grid['wall_group_id']==grid['wall_group_id'][gt]].sum()),
            occupancy_brier=float(np.mean((s['occupancy']-grid['support'][gt])**2)),
            marginal95_both_pose_covered=covered,offline_truth_only=True))
    jsave(out/'OFFLINE_SCORES.json',score)
    with (out/'OFFLINE_SCORES.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(score[0]));writer.writeheader();writer.writerows(score)
    reliable=bool(backend_residual<1e-9 and len(prior)==2950 and worst_refresh<1e-12)
    gate=dict(decision='PASS' if reliable else 'FAIL',prior_id=PRIOR_ID,reliable_actual_reference=reliable,
        numerical_identity_scope='exact normalization of declared finite physical prior, not continuous posterior convergence',
        continuous_probability_map='NOT_ESTABLISHED',max_factor_refresh_difference=worst_refresh,
        complex_rhs_setup=11800,fe_factorizations_setup=2950,likelihood_complex_channel_evaluations=likelihood_work,
        ordinary_smoother_uses_same_bank=True,scheduling_acceleration='NOT_ESTABLISHED',
        E3_exact_nonlinear_nonrecoverable_negative_control='NOT_ESTABLISHED',
        conditions='known target intrinsic and calibrated instrument, 2x2 public unknown occupancy',
        stage_cost_snapshot=BudgetBook().totals('E2E3'))
    jsave(out/'gate.json',gate);print(json.dumps(gate,default=json_default))

if __name__=='__main__':main()
