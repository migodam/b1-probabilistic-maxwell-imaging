#!/usr/bin/env python3
"""Frozen target-intrinsic nuisance stress; no experiment runs on import."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
from wslam.accounting import BudgetBook,BudgetExceeded,json_default
from wslam.forward import Evaluator
from wslam.geometry import validate_geometry,CELL
from wslam.optimize import support_mass,GAIN_PRIOR_SIGMA
from wslam.physics import acquisition,pack
from wslam.world import World

ROOT=Path(__file__).resolve().parents[1]
PARENTS=(0,3,6,9)
LOWER=np.array([1.5,0.,-1.,1.6,-.01,-.01,2.,.08])
UPPER=np.array([6.,.25,1.,2.6,.01,.01,5.,.18])
SCALE=np.array([4.5,.25,2.,1.,.02,.02,3.,.1])


def write_json(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False,default=json_default)+'\n');temporary.replace(path)


def save_npz(path,**data):
    temporary=path.with_suffix('.tmp.npz');np.savez_compressed(temporary,**data);temporary.replace(path)


def decode_prediction(raw):
    a=np.asarray(raw)
    if a.shape not in ((32,2),(2,16,2)):raise ValueError('invalid saved complex prediction')
    return (a[...,0]+1j*a[...,1]).reshape(2,16)


def objective(theta,prediction,y,support):
    white=pack(np.asarray(prediction)-y)
    data_gain=float(.5*white@white+.5*np.sum((theta[4:6]/GAIN_PRIOR_SIGMA)**2))
    return data_gain-float(np.log(support_mass(support)))+float(np.log(2.))


def world_for(theta,support):
    return World(support,float(theta[0]),float(theta[1]),theta[2:4],contrast=float(theta[6]),radius=float(theta[7]),target_loss=.02)


def choose_main(records,observation_id,y):
    options=[]
    for record in records:
        if int(record['observation_id'])!=observation_id or record.get('selection') is None:continue
        selected=record['selection'];theta=np.r_[np.asarray(selected['theta'],float),3.,.12]
        if theta.shape!=(8,) or np.any(theta<LOWER) or np.any(theta>UPPER):raise ValueError('saved main candidate violates fixed public box')
        prediction=decode_prediction(selected['pred']);support=np.asarray(selected['support'],bool)
        value=objective(theta,prediction,y,support)
        options.append(dict(theta=theta,support=support,pred=prediction,objective=value,
            main_method=record['method'],main_fit_success=selected.get('success'),main_fit_status=selected.get('fit_status'),
            source_saved_hierarchical_objective=selected.get('hierarchical_map_objective')))
    if not options:return None,[]
    options.sort(key=lambda r:(r['objective'],r['main_method']))
    return options[0],options


def fit_one(y,initial,book,directory):
    support=initial['support'];solver=Evaluator(book,'A0',h=.17,n=256,order=3,grid_refinement=2)
    cache={};trials=[];candidate=dict(initial,selection_origin='SAVED_MAIN_NOOP')
    paid_count=0;before=book.totals('A0');began=time.perf_counter();cpu=time.process_time()
    def persist():
        arrays={}
        for i,row in enumerate(trials):
            for key,value in row.items():arrays[f'trial_{i:02d}_{key}']=value
        save_npz(directory/'ALL_TRIALS.npz',**arrays)
    def get(theta):
        nonlocal paid_count
        key=tuple(float(theta[i]) for i in (0,1,2,3,6,7))
        if key not in cache:
            if paid_count>=6:raise RuntimeError('frozen six paid-point ceiling reached')
            paid_count+=1
            write_json(directory/'PENDING_PHYSICS.json',dict(theta=theta,paid_point=paid_count,before=book.totals('A0')))
            ev=solver.evaluate(world_for(theta,support),jacobian=True,kind='shared')
            item=dict(F=np.asarray(ev.data).reshape(-1).copy(),J=np.asarray(ev.jacobian).copy(),metadata=ev.metadata)
            cache[key]=item
            save_npz(directory/f'PAID_{paid_count:02d}.npz',theta=theta,F=item['F'],J=item['J'],support=support,
                names=np.asarray(ev.names),metadata_json=np.array(json.dumps(ev.metadata,default=json_default)),
                costs_json=np.array(json.dumps(book.totals('A0'))))
            del ev
            (directory/'PENDING_PHYSICS.json').unlink()
        return cache[key]
    def fun(theta):
        nonlocal candidate
        attempt=dict(kind=np.array('RESIDUAL_ATTEMPT'),theta=theta.copy());trials.append(attempt);persist()
        item=get(theta);gamma=(1+theta[4])*np.exp(1j*theta[5])
        prediction=(gamma*item['F']).reshape(2,16);white=pack(prediction-y)
        residual=np.r_[white,theta[4:6]/GAIN_PRIOR_SIGMA]
        value=objective(theta,prediction,y,support)
        attempt.update(kind=np.array('RESIDUAL_COMPLETE'),pred=prediction,F=item['F'],J=item['J'],residual=residual,objective=np.array(value));persist()
        if value<candidate['objective']:
            candidate=dict(theta=theta.copy(),support=support.copy(),pred=prediction.copy(),objective=value,
                selection_origin='PAID_EXTENSION_RESIDUAL_TRIAL',main_method=initial['main_method'])
        write_json(directory/'CURRENT_ONLINE_SELECTION.json',candidate)
        return residual
    def jac(theta):
        item=get(theta);phase=np.exp(1j*theta[5]);gamma=(1+theta[4])*phase
        physical=item['J']
        if physical.shape!=(32,7):raise ValueError('expected shared physical J eps,ell,px,py,contrast,loss,radius')
        cj=np.column_stack((gamma*physical[:,:4],phase*item['F'],1j*gamma*item['F'],gamma*physical[:,4],gamma*physical[:,6]))
        real=np.vstack((cj.real,cj.imag))/.01
        prior=np.zeros((2,8));prior[:,4:6]=np.eye(2)/GAIN_PRIOR_SIGMA
        result=np.vstack((real,prior))
        trials.append(dict(kind=np.array('JACOBIAN'),theta=theta.copy(),residual_J=result));persist()
        return result
    fit_status='FAILED';error=None;solution_info=None
    try:
        solution=least_squares(fun,initial['theta'],jac=jac,method='trf',bounds=(LOWER,UPPER),x_scale=SCALE,max_nfev=6)
        fit_status='CONVERGED' if solution.success else 'MAX_NFEV_OR_NONCONVERGENCE'
        solution_info=dict(success=bool(solution.success),status=int(solution.status),message=str(solution.message),nfev=int(solution.nfev),njev=int(solution.njev or 0),theta=solution.x)
    except Exception as exc:
        error=dict(type=type(exc).__name__,message=str(exc))
        if isinstance(exc,BudgetExceeded):fit_status='STOPPED_BUDGET'
    finally:
        persist();solver.clear()
    result=dict(fit_status=fit_status,error=error,solver_result=solution_info,selection=candidate,
        saved_main_noop=initial,paid_points_attempted=paid_count,all_nonconvergence_preserved=True,
        costs_before=before,costs_after=book.totals('A0'),wall_seconds=time.perf_counter()-began,cpu_seconds=time.process_time()-cpu)
    write_json(directory/'ONLINE_RESULT.json',result)
    return result


def material_pose_metrics(theta,support,label):
    predicted=validate_geometry(support,theta[0],theta[1],theta[2:4],theta[7])
    actual=validate_geometry(label['support'],label['epsinf'],label['ell'],label['target_pose'],label['target_intrinsic']['radius'])
    z=predicted['support'];truth=actual['support'];union=np.count_nonzero(z|truth)
    delta=theta[2:4]-np.asarray(label['target_pose'])
    return dict(geometry_feasible=predicted['legal'],geometry_reasons=predicted['reasons'],
        pose_error_lambda0=float(np.linalg.norm(delta)),pose_squared_error_lambda0_squared=float(delta@delta),
        material_eps_contrast_L2=float(np.sqrt(CELL**2*np.sum((predicted['epsinf']-actual['epsinf'])**2))),
        material_ell_L2=float(np.sqrt(CELL**2*np.sum((predicted['ell']-actual['ell'])**2))),
        support_iou=float(np.count_nonzero(z&truth)/union) if union else 1.,
        support_symmetric_difference_area=float(CELL**2*np.count_nonzero(z^truth)),
        target_contrast_error=float(theta[6]-label['target_intrinsic']['contrast']),
        target_radius_error=float(theta[7]-label['target_intrinsic']['radius']))


def holdout_prediction(theta,support,book,tx,rx,finer=False):
    evaluator=Evaluator(book,'A0',h=.17,n=384 if finer else 256,grid_refinement=3 if finer else 2,order=3,tx=tx,rx=rx)
    try:
        ev=evaluator.evaluate(world_for(theta,support),jacobian=False)
        prediction=(1+theta[4])*np.exp(1j*theta[5])*np.asarray(ev.data)
        metadata=ev.metadata
        return prediction,metadata
    finally:evaluator.clear()


def run(args):
    root=Path(args.root).resolve();out=Path(args.output or root/'results/a0/intrinsic_stress').resolve()
    supplement=json.loads((root/'A0_INTRINSIC_SUPPLEMENT.json').read_text())
    if supplement['parents']!=list(PARENTS) or supplement['conditions']!=8:raise RuntimeError('frozen intrinsic panel mismatch')
    out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(args.ledger or root/'results/budget_ledger.json');headroom=book.snapshot('A0')['A0']['remaining']
    report=dict(status='INITIALIZED',scope=supplement['scope'],supplement=supplement,online=[],offline_read=False,
        start_headroom=headroom,budget_start=book.totals('A0'),main_A0_modified=False,scientific_claim='NOT_RUN')
    if headroom['rhs']<900 or headroom['factorizations']<60 or headroom['process_cpu_seconds']<=0:
        report.update(status='NOT_RUN',reason='INSUFFICIENT_SHARED_A0_HEADROOM_NO_BORROWING')
        write_json(out/'RESULTS.json',report);return
    fit_report=json.loads((root/'results/a0/optimization/RESULTS.json').read_text())
    if fit_report.get('status') in (None,'RUNNING'):raise RuntimeError('saved main A0 must be terminal')
    rule=json.loads((root/'results/a0/preparation/DISCREPANCY_RULE.json').read_text());upper=float(rule['upper_noise_energy'])
    if rule.get('dimension')!=64 or not np.isfinite(upper) or upper<=0:raise ValueError('unchanged frozen 64-real noise envelope required')
    with np.load(root/'inputs/online/A0_OBSERVATIONS.npz',allow_pickle=False) as archive:
        base=np.array(archive['base']);parents=np.array(archive['parent_id']);repeats=np.array(archive['noise_repeat'])
    conditions=[(p,r,int(np.flatnonzero((parents==p)&(repeats==r))[0])) for p in PARENTS for r in (0,1)]
    report.update(status='ONLINE_RUNNING',noise_rule=rule);write_json(out/'RESULTS.json',report)
    for parent,repeat,index in conditions:
        directory=out/f'parent_{parent:02d}_repeat_{repeat}';directory.mkdir()
        initial,options=choose_main(fit_report['records'],index,base[index])
        write_json(directory/'ONLINE_INITIALIZATION.json',dict(options=options,chosen=initial,truth_access=False))
        if initial is None:result=dict(fit_status='NOT_RUN',selection=None,error='NO_SAVED_MAIN_CANDIDATE')
        else:result=fit_one(base[index],initial,book,directory)
        report['online'].append(dict(parent_id=parent,noise_repeat=repeat,observation_id=index,**result))
        write_json(out/'RESULTS.json',report)
    report.update(status='ALL_ONLINE_TERMINAL_DURABLE',budget_after_online=book.totals('A0'))
    write_json(out/'RESULTS.json',report)
    write_json(out/'ONLINE_COMPLETE.json',report)
    # Offline files are first opened only after all eight online rows are durable.
    labels=json.loads((root/'inputs/offline/CONTINUOUS_TRUTHS.json').read_text())
    labels={int(label['parent_id']):label for label in labels}
    with np.load(root/'inputs/offline/A0_HOLDOUT.npz',allow_pickle=False) as archive:
        holdout=np.array(archive['data']);tx=np.array(archive['tx']);rx=np.array(archive['rx'])
    expected_tx,expected_rx=acquisition(outer=True)
    if not np.array_equal(tx,expected_tx[2:]) or not np.array_equal(rx,expected_rx):raise ValueError('frozen independent outer acquisition required')
    score=[];report.update(offline_read=True,status='OFFLINE_SCORING');write_json(out/'RESULTS.json',report)
    for row in report['online']:
        selected=row['selection'];directory=out/f"parent_{row['parent_id']:02d}_repeat_{row['noise_repeat']}"
        result=dict(parent_id=row['parent_id'],noise_repeat=row['noise_repeat'],observation_id=row['observation_id'],fit_status=row['fit_status'],status='NOT_RUN')
        if selected is not None:
            theta=np.asarray(selected['theta']);support=np.asarray(selected['support']);index=row['observation_id']
            try:
                prediction,metadata=holdout_prediction(theta,support,book,tx,rx)
                white=pack(prediction-holdout[index]);train=pack(np.asarray(selected['pred'])-base[index])
                energy=float(white@white);training_energy=float(train@train)
                save_npz(directory/'PRIMARY_HOLDOUT.npz',theta=theta,support=support,prediction=prediction,observation=holdout[index],white=white,metadata_json=np.array(json.dumps(metadata,default=json_default)))
                result.update(material_pose_metrics(theta,support,labels[row['parent_id']]),status='SCORED',
                    training_noise_energy=training_energy,holdout_noise_energy=energy,
                    training_noise_envelope_compatible=training_energy<=upper,holdout_noise_envelope_compatible=energy<=upper,
                    both_noise_envelopes_compatible=training_energy<=upper and energy<=upper)
                if row['parent_id']==0:
                    fine,fine_meta=holdout_prediction(theta,support,book,tx,rx,finer=True)
                    difference=pack(fine-prediction);finewhite=pack(fine-holdout[index])
                    save_npz(directory/'FINER_HOLDOUT.npz',prediction=fine,primary=prediction,difference_white=difference,white=finewhite,metadata_json=np.array(json.dumps(fine_meta,default=json_default)))
                    result.update(finer_status='SCORED',finer_prediction_difference_white=float(np.linalg.norm(difference)),finer_holdout_noise_energy=float(finewhite@finewhite))
            except Exception as exc:result.update(status='FAILED',error_type=type(exc).__name__,error=str(exc))
        score.append(result);write_json(out/'OFFLINE_SCORES.json',score)
    report.update(status='STRESS_EXECUTION_COMPLETE_PARENT_REVIEW_REQUIRED',budget_end=book.totals('A0'),offline_score_path=str(out/'OFFLINE_SCORES.json'),posterior_claim=False,unknown_intrinsic_prior_SBC=False)
    write_json(out/'RESULTS.json',report)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',default=str(ROOT));parser.add_argument('--output');parser.add_argument('--ledger')
    run(parser.parse_args(argv))


if __name__=='__main__':main()
