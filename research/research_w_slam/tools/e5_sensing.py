#!/usr/bin/env python3
"""Conditional finite-prior active sensing. Importing performs no experiment."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.finite_bayes import PRIOR_ID,summarize,update_logweights
from wslam.forward import Evaluator
from wslam.physics import acquisition
from wslam.world import World

ROOT=Path(__file__).resolve().parents[1]
POLICIES=('Random','Uniform','EIG','SingularValue','PrincipalAngleSum','GeneralizedRetention','TaskEVSI')
SCALES=np.array([4.5,.25,2.,1.])
SIGMA=.01


def write_json(path,data):
    path=Path(path);temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def save_npz(path,**data):
    path=Path(path);temporary=path.with_suffix('.tmp.npz')
    np.savez_compressed(temporary,**data);temporary.replace(path)


def loglikelihood(prediction,observation):
    difference=(np.asarray(prediction)-np.asarray(observation))/SIGMA
    return -.5*np.sum(difference.real**2+difference.imag**2,axis=tuple(range(1,difference.ndim)))


def weights_from_log(prior,ll):
    return np.exp(update_logweights(prior,[ll]))


def pose_risk(w,pose):
    mean=w@pose
    return float(w@np.sum((pose-mean)**2,axis=1))


def _range(matrix):
    u,s,_=np.linalg.svd(matrix,full_matrices=False)
    threshold=max(matrix.shape)*np.finfo(float).eps*(s[0] if s.size else 0.)
    rank=int(np.count_nonzero(s>threshold))
    return u[:,:rank],s,rank,float(threshold)


def local_action_metrics(jacobian):
    """No inverse, pseudo-inverse, statistical jitter, or truth input."""
    j=np.asarray(jacobian).reshape(4,16,-1)[:,:,:4]
    result=[]
    for action in range(4):
        physical=np.vstack((j[action].real,j[action].imag))/SIGMA
        a,b=physical[:,:2],physical[:,2:4]
        scaled=physical*SCALES
        _,singular,rank,threshold=_range(scaled)
        ua,wall_s,wall_rank,wall_threshold=_range(a)
        ub,pose_s,pose_rank,pose_threshold=_range(b)
        projection=ub@(ub.T@a)
        residual=a-projection
        information=residual.T@residual
        # Retention lives on observable wall response only. No inverse needed.
        cross=ub.T@ua
        retention=np.eye(wall_rank)-cross.T@cross
        retention=(retention+retention.T)/2
        retention_values=np.linalg.eigvalsh(retention)
        cosines=np.linalg.svd(cross,compute_uv=False)
        cosines=np.r_[cosines,np.zeros(wall_rank-len(cosines))]
        angles=np.arccos(np.clip(cosines,0.,1.))
        result.append(dict(action=action,physical_jacobian=physical,physical_A=a,physical_B=b,
            parameter_scales=SCALES,scaled_joint_singular_values=singular,joint_rank=rank,
            rank_threshold=threshold,max_singular=float(singular[0]) if len(singular) else 0.,
            min_nonzero_singular=float(singular[rank-1]) if rank else 0.,
            # A rank-deficient four-parameter local model must not be rewarded
            # for dropping its invisible dimensions (e.g. empty MAP wall).
            singular_value_score=float(singular[-1]**2) if rank==4 else 0.,
            wall_rank=wall_rank,pose_rank=pose_rank,wall_rank_threshold=wall_threshold,
            pose_rank_threshold=pose_threshold,physical_wall_singular_values=wall_s,
            physical_pose_singular_values=pose_s,principal_angles_radians=angles,
            principal_angle_sum_score=float(np.sum(np.sin(angles)**2)),
            retention_eigenvalues=retention_values,retention_score=float(retention_values.sum()),
            physical_absolute_wall_information=information,
            rank_status='NUMERICAL_TOLERANCE_DEPENDENT',zero_wall_response='EXCLUDED_UNDEFINED_RETENTION'))
    return result


def predictive_scores(posterior,bank,pose,parent,repeat):
    """Hypothetical outcomes depend only on current posterior and public bank."""
    rng=np.random.default_rng(np.random.SeedSequence([2026102506,parent,repeat,1]))
    indices=rng.choice(len(posterior),size=32,p=posterior)
    # Identical latent draws and standardized noise across candidate actions.
    noise=rng.normal(size=(32,16))+1j*rng.normal(size=(32,16))
    outcomes=bank[indices]+SIGMA*noise[:,None,:]
    eig=np.empty((4,32));risks=np.empty((4,32))
    for action in range(4):
        for draw in range(32):
            ll=loglikelihood(bank[:,action,:],outcomes[draw,action])
            updated=weights_from_log(posterior,ll)
            positive=updated>0
            eig[action,draw]=np.sum(updated[positive]*(np.log(updated[positive])-np.log(posterior[positive])))
            risks[action,draw]=pose_risk(updated,pose)
    return dict(eig=eig.mean(axis=1),task_evsi=pose_risk(posterior,pose)-risks.mean(axis=1),
        draw_eig=eig,draw_posterior_risk=risks,simulated_indices=indices,
        standardized_noise=noise,outcomes=outcomes,draws=32,
        scope='32 predictive Monte Carlo outcomes; no global utility certificate')


def choose_without_future(posterior,bank,grid,local,parent,repeat):
    """No truth index, truth pose, or realized future measurement is accepted."""
    timings={};wall_timings={};choices={};scores={}
    started=time.process_time();wall_started=time.perf_counter()
    rng=np.random.default_rng(np.random.SeedSequence([2026102506,parent,repeat,0]))
    choices['Random']=int(rng.integers(4));timings['Random']=time.process_time()-started
    wall_timings['Random']=time.perf_counter()-wall_started
    started=time.process_time();wall_started=time.perf_counter()
    choices['Uniform']=(parent*2+repeat)%4;timings['Uniform']=time.process_time()-started
    wall_timings['Uniform']=time.perf_counter()-wall_started
    started=time.process_time();wall_started=time.perf_counter();mc=predictive_scores(posterior,bank,grid['pose'],parent,repeat)
    mc_cpu=time.process_time()-started;mc_wall=time.perf_counter()-wall_started
    for policy,key in [('EIG','eig'),('TaskEVSI','task_evsi')]:
        choices[policy]=int(np.argmax(mc[key]));scores[policy]=mc[key]
        wall_timings[policy]=mc_wall
        timings[policy]=mc_cpu  # Standalone attribution for shared exact updates.
    for policy,key in [('SingularValue','singular_value_score'),('PrincipalAngleSum','principal_angle_sum_score'),('GeneralizedRetention','retention_score')]:
        started=time.process_time();wall_started=time.perf_counter();values=np.array([row[key] for row in local])
        choices[policy]=int(np.argmax(values));scores[policy]=values
        timings[policy]=time.process_time()-started
        wall_timings[policy]=time.perf_counter()-wall_started
    return dict(choices=choices,scores=scores,selection_cpu_seconds=timings,selection_wall_seconds=wall_timings,predictive=mc,
        tie_rule='first action index under exact ties',future_measurement_access=False,
        truth_access=False,shared_mc_cpu_seconds=mc_cpu)


def future_measurements_after_selection(bank,truth_index,parent,repeat):
    """Offline generator, called only after durable action selection."""
    rng=np.random.default_rng(np.random.SeedSequence([2026102507,parent,repeat]))
    noise=SIGMA*(rng.normal(size=(4,16))+1j*rng.normal(size=(4,16)))
    return bank[int(truth_index)]+noise,noise


def run(args):
    gate=json.loads(Path(args.reference_gate).read_text())
    if gate.get('prior_id')!=PRIOR_ID or gate.get('reliable_actual_reference') is not True or gate.get('decision')!='PASS':
        raise RuntimeError('E2 reliable actual finite reference gate required before E5')
    supplement=json.loads((ROOT/'FINITE_PILOT_SUPPLEMENT.json').read_text())
    if supplement.get('prior_id')!=PRIOR_ID:
        raise RuntimeError('finite pilot prior mismatch')
    source=Path(args.bank)
    with np.load(source/'online/FINITE_GRID.npz',allow_pickle=False) as archive:
        grid={name:archive[name] for name in archive.files}
    grid['label']=PRIOR_ID
    with np.load(source/'online/FINITE_BANK.npz',allow_pickle=False) as archive:
        if int(archive['completed'])!=2950 or str(archive['prior_id'])!=PRIOR_ID:
            raise RuntimeError('complete finite prior bank required')
        bank=np.array(archive['data']);tx=np.array(archive['tx']);rx=np.array(archive['rx'])
    with np.load(source/'online/OBSERVATIONS.npz',allow_pickle=False) as archive:
        observations=np.array(archive['data'])
        known_gamma=complex(archive['known_gamma'])
        t1_ids=np.array(archive['sequential_factor_id'])[:,:,0]
    if known_gamma!=1 or observations.shape!=(12,2,4,16):
        raise RuntimeError('known gamma1 and 12 parent × 2 repeat acquisition required')
    if bank.shape!=(2950,4,16) or not np.isfinite(bank).all():
        raise RuntimeError('invalid finite bank')
    declared_tx,declared_rx=acquisition(outer=True)
    if not np.array_equal(tx,declared_tx) or not np.array_equal(rx,declared_rx):
        raise RuntimeError('bank acquisition differs from paid local controls')
    bank_cost_path=source/'online/JOB_COST_FINAL.json'
    bank_cost_receipt=json.loads(bank_cost_path.read_text())
    if bank_cost_receipt.get('complete') is not True or bank_cost_receipt.get('prior_id')!=PRIOR_ID:
        raise RuntimeError('terminal JOB_COST_FINAL.json required; in-flight bank CPU receipts are incomplete')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    offline=out/'offline';offline.mkdir()
    book=BudgetBook(args.ledger)
    manifest=dict(status='RUNNING',prior_id=PRIOR_ID,reference_gate=gate,finite_pilot_supplement=supplement,
        conditional_finite_prior_only=True,continuous_probability_certificate=False,scientific_claim='NOT_RUN',
        parents=12,noise_repeats=2,statistical_unit='parent; the 24 noise cases are not 24 independent parents',
        known_target=dict(radius=.12,contrast=3.,loss=.02),known_gamma=[1.,0.],
        parameter_scales=dict(wall=[4.5,.25],pose=[2.,1.]),noise_sigma_real_imag=SIGMA,
        planned_E5_rhs=480,planned_E5_factorizations=24,initial_costs=book.totals('E5'),
        shared_bank_setup=dict(actual_charged_stage='E2E3',rhs=11800,factorizations=2950,raw_bank_receipt=bank_cost_receipt,scope='full finite bank shared once; not cheap online work'),
        completed_cases=0)
    write_json(out/'MANIFEST.json',manifest)
    online_results=[];gt_results=[]
    try:
        for parent in range(12):
            for repeat in range(2):
                case=out/f'parent_{parent:02d}_repeat_{repeat}';case.mkdir()
                t1_id=str(t1_ids[parent,repeat])
                old_ll=loglikelihood(bank[:,0,8:9],observations[parent,repeat,0,8:9])
                early=np.exp(update_logweights(grid['prior_weight'],[(t1_id,old_ll)]))
                map_index=int(np.argmax(early))
                map_world=World(grid['support'][map_index].reshape(2,2),float(grid['epsinf'][map_index]),float(grid['ell'][map_index]),grid['pose'][map_index],radius=.12,contrast=3.,target_loss=.02)
                solver=Evaluator(book,'E5',h=.17,n=256,order=3,grid_refinement=2,tx=tx,rx=rx)
                before=book.totals('E5');j_cpu=time.process_time();j_wall=time.perf_counter()
                try:
                    ev=solver.evaluate(map_world,f=1.,jacobian=True,kind='shared')
                    j_time=time.process_time()-j_cpu;j_elapsed=time.perf_counter()-j_wall
                    jac=np.array(ev.jacobian);after=book.totals('E5')
                    save_npz(case/'PAID_LOCAL_J.npz',J=jac,names=np.asarray(ev.names),map_index=np.array(map_index),
                        data=ev.data,physical_parameter_scales=SCALES,residual=np.array(ev.state.residual),
                        receipt_json=np.array(json.dumps(dict(before=before,after=after,metadata=ev.metadata),default=json_default)))
                    del ev
                finally:
                    solver.clear()
                audit_start=time.process_time();audit_wall_start=time.perf_counter()
                local=local_action_metrics(jac);audit_cpu=time.process_time()-audit_start;audit_wall=time.perf_counter()-audit_wall_start
                selection=choose_without_future(early,bank,grid,local,parent,repeat)
                write_json(case/'ACTION_SELECTION.json',dict(parent=parent,repeat=repeat,map_index=map_index,
                    t1_factor_id=t1_id,**selection,local_control_metrics=local,
                    shared_local_j_cpu_seconds=j_time,shared_local_j_wall_seconds=j_elapsed,local_audit_cpu_seconds=audit_cpu))
                # Only now may offline truth indices or actual future noise be read.
                with np.load(source/'offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz',allow_pickle=False) as truth:
                    truth_index=int(truth['truth_index'][parent])
                future,noise=future_measurements_after_selection(bank,truth_index,parent,repeat)
                future_ids=np.array([f'E5-parent-{parent:02d}-repeat-{repeat}-action-{a}-new-full16' for a in range(4)])
                save_npz(offline/f'FUTURE_parent_{parent:02d}_repeat_{repeat}.npz',measurements=future,noise=noise,
                    factor_id=future_ids,truth_index=np.array(truth_index),seed=np.array(2026102507))
                raw_weights={};case_online=[];case_gt=[]
                for policy in POLICIES:
                    action=selection['choices'][policy]
                    future_id=str(future_ids[action])
                    new_ll=loglikelihood(bank[:,action,:],future[action])
                    began=time.process_time();inference_wall_start=time.perf_counter()
                    posterior=np.exp(update_logweights(grid['prior_weight'],[(t1_id,old_ll),(future_id,new_ll)]))
                    summary=summarize(posterior,grid)
                    inference_cpu=time.process_time()-began;inference_wall=time.perf_counter()-inference_wall_start
                    raw_weights[policy]=posterior
                    is_local=policy in ('SingularValue','PrincipalAngleSum','GeneralizedRetention')
                    selection_cpu=selection['selection_cpu_seconds'][policy]+(j_time+audit_cpu if is_local else 0.)
                    row=dict(parent=parent,repeat=repeat,policy=policy,action=action,
                        factor_ids=[t1_id,future_id],summary=summary,target_bayes_risk=pose_risk(posterior,grid['pose']),
                        initial_target_bayes_risk=pose_risk(early,grid['pose']),
                        selection_cpu_seconds_standalone=selection_cpu,inference_cpu_seconds=inference_cpu,inference_wall_seconds=inference_wall,
                        selection_wall_seconds_standalone=selection['selection_wall_seconds'][policy]+(j_elapsed+audit_wall if is_local else 0.),
                        standalone_selection_rhs=20 if is_local else 0,standalone_selection_factorizations=1 if is_local else 0,
                        shared_paid_local_J_actual_rhs=after['rhs']-before['rhs'],
                        shared_paid_local_J_actual_factorizations=after['factorizations']-before['factorizations'],
                        bank_setup_selection_required=policy not in ('Random','Uniform'),
                        standalone_common_posterior_bank_rhs=11800,standalone_common_posterior_bank_factorizations=2950,
                        cost_warning='Bank preparation charged once in E2E3 and attributed explicitly; policy rows are standalone alternatives, not additive actual charges.')
                    case_online.append(row)
                    truth_pose=grid['pose'][truth_index]
                    case_gt.append(dict(parent=parent,repeat=repeat,policy=policy,action=action,
                        truth_index=truth_index,truth_pose=truth_pose,
                        target_squared_error=float(np.sum((summary['pose_mean']-truth_pose)**2))))
                save_npz(case/'POSTERIORS.npz',early=early,**raw_weights)
                write_json(case/'ONLINE_RESULTS.json',case_online)
                write_json(offline/f'GT_parent_{parent:02d}_repeat_{repeat}.json',case_gt)
                online_results.extend(case_online);gt_results.extend(case_gt)
                manifest.update(completed_cases=len(online_results)//len(POLICIES),actual_E5_costs=book.totals('E5'))
                write_json(out/'MANIFEST.json',manifest)
        grouped=[]
        for parent in range(12):
            for policy in POLICIES:
                rows=[r for r in gt_results if r['parent']==parent and r['policy']==policy]
                grouped.append(dict(parent=parent,policy=policy,noise_repeats=2,
                    mean_target_squared_error=float(np.mean([r['target_squared_error'] for r in rows]))))
        write_json(offline/'PARENT_GROUPED_GT_SCORES.json',grouped)
        write_json(out/'ONLINE_RESULTS.json',online_results)
        manifest.update(status='COMPLETE',actual_E5_costs=book.totals('E5'),scientific_claim='NOT_RUN',
            prior_reference_scope='Reliable actual reference on explicit finite prior only',
            policy_comparison='Raw conditional experiment; no performance or cost advantage claim')
        write_json(out/'MANIFEST.json',manifest)
        write_json(out/'READY.json',dict(status='RAW_RESULTS_READY',prior_id=PRIOR_ID,manifest=str(out/'MANIFEST.json'),scientific_claim='NOT_RUN'))
    except BaseException as exc:
        manifest.update(status='PARTIAL_FAILURE',error=repr(exc),actual_E5_costs=book.totals('E5'))
        write_json(out/'MANIFEST.json',manifest)
        raise


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank',default=str(ROOT/'results/e2e3/finite_bank'))
    parser.add_argument('--reference-gate',default=str(ROOT/'results/e2e3/inference/gate.json'))
    parser.add_argument('--output',default=str(ROOT/'results/e5/finite_sensing'))
    parser.add_argument('--ledger',default=None)
    run(parser.parse_args(argv))


if __name__=='__main__':
    main()
