#!/usr/bin/env python3
"""Frozen ordinary optimization diagnostic; importing executes no experiment."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.optimize import LOWER,UPPER,X_SCALE,GAIN_PRIOR_SIGMA
from wslam.physics import acquisition,pack
from wslam.world import World
from a0_intrinsic import write_json,save_npz,decode_prediction,objective,material_pose_metrics,holdout_prediction

ROOT=Path(__file__).resolve().parents[1]
WIDTH=UPPER-LOWER
PARENTS=(0,3,6,9)


class CallCeiling(RuntimeError):pass
class GradientAuditFailure(RuntimeError):pass


def select_main(records,index,y):
    candidates=[]
    for row in records:
        if int(row['observation_id'])!=index or row.get('selection') is None:continue
        fit=row['selection'];theta=np.asarray(fit['theta'],float);support=np.asarray(fit['support'],bool)
        if theta.shape!=(6,) or np.any(theta<LOWER) or np.any(theta>UPPER):raise ValueError('saved start outside frozen public box')
        pred=decode_prediction(fit['pred'])
        candidates.append(dict(theta=theta,support=support,pred=pred,objective=objective(theta,pred,y,support),
            main_method=row['method'],main_fit_success=fit.get('success'),main_fit_status=fit.get('fit_status')))
    candidates.sort(key=lambda row:(row['objective'],row['main_method']))
    return (candidates[0] if candidates else None),candidates


def polish(y,start,book,directory,audit):
    from wslam.physics import weighted_material_gradient
    solver=Evaluator(book,'A0',h=.17,n=256,order=3,grid_refinement=2)
    support=start['support'];combined_cache={};attempts=[];unique_count=0
    selected=dict(start,origin='SAVED_MAIN_NOOP');paid_before=book.totals('A0')
    began=time.perf_counter();cpu=time.process_time();audit_record=None
    def evaluate(theta):
        nonlocal selected,unique_count,audit_record
        theta=np.asarray(theta,float);key=tuple(theta)
        attempt=dict(attempt=len(attempts)+1,theta=theta.copy(),status='ATTEMPTED',cached=key in combined_cache)
        attempts.append(attempt);write_json(directory/'ATTEMPTS.json',attempts)
        if key in combined_cache:
            attempt.update(status='CACHED_COMBINED_VALUE_GRADIENT');write_json(directory/'ATTEMPTS.json',attempts)
            return combined_cache[key]
        if unique_count>=20:
            attempt.update(status='UNPAID_CALL_CEILING_STOP');write_json(directory/'ATTEMPTS.json',attempts)
            raise CallCeiling('20 unique joint-theta evaluations exhausted; no further PDE permitted')
        unique_count+=1;before=book.totals('A0')
        pending=dict(theta=theta,unique_evaluation=unique_count,cost_before=before)
        write_json(directory/'PENDING_PHYSICS.json',pending)
        world=World(support,theta[0],theta[1],theta[2:4],radius=.12,contrast=3.,target_loss=.02)
        ev=solver.evaluate(world,jacobian=False)
        f=np.asarray(ev.data).reshape(-1);phase=np.exp(1j*theta[5]);gamma=(1+theta[4])*phase
        raw=gamma*f-y.reshape(-1)
        wall_pose,metadata=weighted_material_gradient(ev,world,raw.reshape(ev.data.shape),gamma)
        wall_pose=np.asarray(wall_pose,float)
        if wall_pose.shape!=(4,) or not np.isfinite(wall_pose).all():raise RuntimeError('weighted gradient requires finite four-vector')
        gain=float(np.vdot(phase*f,raw).real/.01**2+theta[4]/GAIN_PRIOR_SIGMA**2)
        phase_grad=float(np.vdot(1j*gamma*f,raw).real/.01**2+theta[5]/GAIN_PRIOR_SIGMA**2)
        gradient=np.r_[wall_pose,gain,phase_grad]
        prediction=(gamma*f).reshape(2,16);value=objective(theta,prediction,y,support)
        # Persist the native weighted gradient before optional independent J audit.
        save_npz(directory/f'EVAL_{unique_count:02d}.npz',theta=theta,F=f,prediction=prediction,raw_residual=raw,
            gradient=gradient,gradient_normalized=gradient*WIDTH,normalized_coordinates=(theta-LOWER)/WIDTH,objective=np.array(value),support=support,
            weighted_metadata_json=np.array(json.dumps(metadata,default=json_default)),
            costs_before_json=np.array(json.dumps(before)),costs_after_weighted_json=np.array(json.dumps(book.totals('A0'))))
        if audit and unique_count==1:
            upgraded=solver.evaluate(world,jacobian=True,kind='shared')
            direct=np.real((gamma*np.asarray(upgraded.jacobian)[:,:4]).conj().T@raw)/.01**2
            scale=X_SCALE[:4];difference=(wall_pose-direct)*scale
            denominator=max(1.,float(np.max(abs(direct*scale))),float(np.max(abs(wall_pose*scale))))
            error=float(np.max(abs(difference))/denominator)
            audit_record=dict(weighted=wall_pose,direct=direct,physical_scales=scale,relative_scaled_max_error=error,
                threshold=1e-9,passed=bool(error<=1e-9),full_J_metadata=upgraded.metadata)
            save_npz(directory/'NATIVE_FULL_J_AUDIT.npz',theta=theta,J=upgraded.jacobian,F=f,raw_residual=raw,
                direct_gradient=direct,weighted_gradient=wall_pose,scaled_error=np.array(error))
            write_json(directory/'GRADIENT_AUDIT.json',audit_record)
            del upgraded
            if not audit_record['passed']:
                attempt.update(status='WEIGHTED_GRADIENT_AUDIT_FAILED');write_json(directory/'ATTEMPTS.json',attempts)
                raise GradientAuditFailure(f'weighted material gradient audit failed: scaled error {error}')
        combined_cache[key]=(value,gradient)
        attempt.update(status='COMPLETE',objective=value,gradient=gradient,cost_after=book.totals('A0'))
        if value<selected['objective']:
            selected=dict(theta=theta.copy(),support=support.copy(),pred=prediction.copy(),objective=value,
                main_method=start['main_method'],origin='PAID_ADJOINT_POLISH_TRIAL')
        write_json(directory/'CURRENT_ONLINE_SELECTION.json',selected)
        write_json(directory/'ATTEMPTS.json',attempts);(directory/'PENDING_PHYSICS.json').unlink()
        del ev
        return value,gradient
    failure=None;result=None;status='FAILED';fatal=None
    def normalized_evaluate(u):
        value,physical_gradient=evaluate(LOWER+np.asarray(u)*WIDTH)
        return value,physical_gradient*WIDTH
    try:
        result=minimize(normalized_evaluate,(start['theta']-LOWER)/WIDTH,jac=True,method='L-BFGS-B',bounds=[(0.,1.)]*6,
            options=dict(maxfun=20,maxiter=20,maxls=20))
        status='CONVERGED' if result.success else 'MAX_CALLS_OR_NONCONVERGENCE'
    except CallCeiling as exc:status='CALL_CEILING_STOP_NONCONVERGED';failure=repr(exc)
    except GradientAuditFailure as exc:status='WEIGHTED_FORMULA_AUDIT_FAILED';failure=repr(exc);fatal=exc
    except Exception as exc:status='FAILED';failure=repr(exc);fatal=exc
    finally:solver.clear()
    report=dict(status=status,error=failure,selection=selected,saved_main_noop=start,unique_joint_evaluations=unique_count,
        attempts=len(attempts),gradient_audit=audit_record,cost_before=paid_before,cost_after=book.totals('A0'),
        wall_seconds=time.perf_counter()-began,process_cpu_seconds=time.process_time()-cpu,
        solver=None if result is None else dict(success=bool(result.success),status=int(result.status),message=str(result.message),
            nfev=int(result.nfev),njev=int(result.njev),nit=int(result.nit),theta=LOWER+result.x*WIDTH,normalized_coordinates=result.x,gradient_normalized=result.jac,gradient_physical=result.jac/WIDTH))
    write_json(directory/'ONLINE_RESULT.json',report)
    if fatal is not None:raise fatal
    return report


def run(args):
    root=Path(args.root).resolve();out=Path(args.output or root/'results/a0/adjoint_polish').resolve()
    supplement=json.loads((root/'A0_ADJOINT_POLISH_SUPPLEMENT.json').read_text())
    if supplement['parents']!=list(PARENTS) or supplement['optimizer']['max_unique_joint_theta']!=20:raise RuntimeError('frozen polish definition changed')
    out.mkdir(parents=True,exist_ok=False);book=BudgetBook(args.ledger or root/'results/budget_ledger.json')
    remaining=book.snapshot('A0')['A0']['remaining']
    report=dict(status='INITIALIZED',supplement=supplement,online=[],offline_read=False,main_A0_unchanged=True,
        budget_start=book.totals('A0'),starting_headroom=remaining,scientific_gate='NO_NEW_GO_GATE',
        scope='ORDINARY_OPTIMIZATION_ERROR_DIAGNOSTIC_NOT_MATCHED_SPEED_WIN',
        setup_attribution='Starts depend on all three main arms and their shared paid bank; retain those costs in any comparison.')
    write_json(out/'RESULTS.json',report)
    if remaining['rhs']<700 or remaining['factorizations']<180 or remaining['process_cpu_seconds']<=0:
        report.update(status='NOT_RUN',reason='INSUFFICIENT_SHARED_A0_HEADROOM_NO_BORROWING');write_json(out/'RESULTS.json',report);return
    main=json.loads((root/'results/a0/optimization/RESULTS.json').read_text())
    if main.get('status') in (None,'RUNNING'):raise RuntimeError('main online fits must be terminal')
    with np.load(root/'inputs/online/A0_OBSERVATIONS.npz',allow_pickle=False) as archive:
        base=np.array(archive['base']);parents=np.array(archive['parent_id']);repeats=np.array(archive['noise_repeat'])
    conditions=[(p,r,int(np.flatnonzero((parents==p)&(repeats==r))[0])) for p in PARENTS for r in (0,1)]
    try:
        report['status']='ONLINE_RUNNING';write_json(out/'RESULTS.json',report)
        for case_index,(parent,repeat,index) in enumerate(conditions):
            directory=out/f'parent_{parent:02d}_repeat_{repeat}';directory.mkdir()
            initial,candidates=select_main(main['records'],index,base[index])
            write_json(directory/'INITIALIZATION.json',dict(chosen=initial,all_main_candidates=candidates,truth_access=False))
            if initial is None:raise RuntimeError('fixed panel missing saved main candidate; no alternative parent selection')
            online=polish(base[index],initial,book,directory,audit=case_index==0)
            report['online'].append(dict(parent_id=parent,noise_repeat=repeat,observation_id=index,**online));write_json(out/'RESULTS.json',report)
    except BaseException as exc:
        report.update(status='ONLINE_STOPPED_NO_OFFLINE_ACCESS',error=repr(exc),budget_end=book.totals('A0'))
        completed={(r['parent_id'],r['noise_repeat']) for r in report['online']}
        report['remaining_conditions']=[dict(parent_id=p,noise_repeat=r,status='NOT_RUN_AFTER_ONLINE_FAILURE') for p,r,_ in conditions if (p,r) not in completed]
        write_json(out/'RESULTS.json',report);raise
    report.update(status='ALL_ONLINE_TERMINAL_DURABLE',budget_after_online=book.totals('A0'))
    write_json(out/'RESULTS.json',report);write_json(out/'ONLINE_COMPLETE.json',report)
    # Previously exposed development holdout; no claim of a fresh blind test.
    rule=json.loads((root/'results/a0/preparation/DISCREPANCY_RULE.json').read_text());upper=float(rule['upper_noise_energy'])
    if rule.get('dimension')!=64 or not np.isfinite(upper) or upper<=0:raise ValueError('unchanged frozen noise envelope required')
    labels=json.loads((root/'inputs/offline/CONTINUOUS_TRUTHS.json').read_text());labels={int(v['parent_id']):v for v in labels}
    with np.load(root/'inputs/offline/A0_HOLDOUT.npz',allow_pickle=False) as archive:
        holdout=np.array(archive['data']);tx=np.array(archive['tx']);rx=np.array(archive['rx'])
    alltx,allrx=acquisition(outer=True)
    if not np.array_equal(tx,alltx[2:]) or not np.array_equal(rx,allrx):raise ValueError('independent outer acquisition mismatch')
    report.update(offline_read=True,holdout_scope='EXPOSED_DEVELOPMENT_REEVALUATION_NOT_NEW_BLIND_TEST');write_json(out/'RESULTS.json',report)
    scores=[]
    for row in report['online']:
        selected=row['selection'];theta6=np.asarray(selected['theta']);theta=np.r_[theta6,3.,.12]
        support=np.asarray(selected['support']);index=row['observation_id'];directory=out/f"parent_{row['parent_id']:02d}_repeat_{row['noise_repeat']}"
        score=dict(parent_id=row['parent_id'],noise_repeat=row['noise_repeat'],observation_id=index,fit_status=row['status'],status='NOT_RUN')
        try:
            prediction,metadata=holdout_prediction(theta,support,book,tx,rx)
            white=pack(prediction-holdout[index]);training=pack(np.asarray(selected['pred'])-base[index])
            energy=float(white@white);trainenergy=float(training@training)
            save_npz(directory/'PRIMARY_HOLDOUT.npz',theta=theta6,support=support,prediction=prediction,white=white,metadata_json=np.array(json.dumps(metadata,default=json_default)))
            score.update(material_pose_metrics(theta,support,labels[row['parent_id']]),status='SCORED',training_noise_energy=trainenergy,holdout_noise_energy=energy,
                training_noise_envelope_compatible=trainenergy<=upper,holdout_noise_envelope_compatible=energy<=upper,
                both_noise_envelopes_compatible=trainenergy<=upper and energy<=upper)
            if row['parent_id']==0:
                try:
                    finer,meta=holdout_prediction(theta,support,book,tx,rx,finer=True)
                    difference=pack(finer-prediction);fwhite=pack(finer-holdout[index])
                    save_npz(directory/'FINER_HOLDOUT.npz',prediction=finer,difference_white=difference,white=fwhite,metadata_json=np.array(json.dumps(meta,default=json_default)))
                    score.update(finer_status='SCORED',finer_prediction_difference_white=float(np.linalg.norm(difference)),finer_holdout_noise_energy=float(fwhite@fwhite))
                except Exception as exc:score.update(finer_status='FAILED',finer_error=repr(exc))
        except Exception as exc:score.update(status='FAILED',error=repr(exc))
        scores.append(score);write_json(out/'OFFLINE_SCORES.json',scores)
    report.update(status='OPTIMIZER_DIAGNOSTIC_COMPLETE_PARENT_REVIEW_REQUIRED',budget_end=book.totals('A0'),matched_speed_claim=False,posterior_claim=False)
    write_json(out/'RESULTS.json',report)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',default=str(ROOT));parser.add_argument('--output');parser.add_argument('--ledger')
    run(parser.parse_args(argv))


if __name__=='__main__':main()
