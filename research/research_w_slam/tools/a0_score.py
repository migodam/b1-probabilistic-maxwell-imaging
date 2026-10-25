"""Separate post-fit A0 scorer. No physical calls occur on import.

Offline labels and holdout are opened only after the saved online optimizer
reaches a terminal status. They are never passed to fitting or selection.
"""
from pathlib import Path
import argparse
import csv
import gc
import json
import time

import numpy as np

from wslam.accounting import BudgetBook, BudgetExceeded, json_default
from wslam.forward import Evaluator
from wslam.geometry import CELL, SHAPE, validate_geometry
from wslam.physics import acquisition, pack
from wslam.world import World

FINE_PARENT_IDS = (0,3,6,9)


def write_json(path, value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,default=json_default,allow_nan=False)+'\n')
    temporary.replace(path)


def complex_json(value):
    """Decode accounting.json_default's [real, imaginary] encoding."""
    raw=np.asarray(value)
    if raw.shape not in ((32,2),(2,16,2)):
        raise ValueError('Unexpected serialized complex base prediction')
    return (raw[...,0]+1j*raw[...,1]).reshape(2,16)


def edges(mask):
    """Public-grid material/air boundary segments, including ROI boundary."""
    mask=np.asarray(mask,bool)
    horizontal=np.zeros((mask.shape[0]+1,mask.shape[1]),bool)
    vertical=np.zeros((mask.shape[0],mask.shape[1]+1),bool)
    horizontal[:-1]^=mask;horizontal[1:]^=mask
    vertical[:,:-1]^=mask;vertical[:,1:]^=mask
    return np.r_[horizontal.ravel(),vertical.ravel()]


def offline_metrics(selection, label):
    theta=np.asarray(selection['theta'],float)
    support=np.asarray(selection['support'])
    if theta.shape!=(6,) or not np.isfinite(theta).all():
        raise ValueError('Invalid selected six-parameter fit')
    candidate=validate_geometry(support,theta[0],theta[1],theta[2:4],.12)
    intrinsic=label['target_intrinsic']
    actual=validate_geometry(label['support'],label['epsinf'],label['ell'],
                             label['target_pose'],intrinsic['radius'])
    if not actual['legal']:
        raise ValueError('Stored evaluation label violates frozen geometry')
    cm,tm=candidate['support'],actual['support']
    intersection=int(np.count_nonzero(cm&tm));union=int(np.count_nonzero(cm|tm))
    ce,te=edges(cm),edges(tm)
    common=int(np.count_nonzero(ce&te));nc,nt=int(ce.sum()),int(te.sum())
    precision=common/nc if nc else (1. if not nt else 0.)
    recall=common/nt if nt else (1. if not nc else 0.)
    pose_difference=theta[2:4]-np.asarray(label['target_pose'],float)
    # Contrast is zero in air. All errors use the same public 24x6 grid and
    # area weights, so 2x2 and 3x3 representations are compared consistently.
    eps_difference=(candidate['epsinf']-1)-(actual['epsinf']-1)
    ell_difference=candidate['ell']-actual['ell']
    area=CELL**2
    return dict(geometry_feasible=bool(candidate['legal']),geometry_reasons=candidate['reasons'],
        material_eps_contrast_L2=float(np.sqrt(area*np.sum(eps_difference**2))),
        material_ell_L2=float(np.sqrt(area*np.sum(ell_difference**2))),
        material_eps_contrast_RMSE=float(np.sqrt(np.mean(eps_difference**2))),
        material_ell_RMSE=float(np.sqrt(np.mean(ell_difference**2))),
        public_material_grid=list(SHAPE),material_error_area_unit='lambda0_squared',
        pose_error_lambda0=float(np.linalg.norm(pose_difference)),
        pose_squared_error_lambda0_squared=float(pose_difference@pose_difference),
        support_iou=intersection/union if union else 1.,
        support_symmetric_difference_area=float(area*np.count_nonzero(cm^tm)),
        edge_precision=precision,edge_recall=recall,
        edge_f1=2*precision*recall/(precision+recall) if precision+recall else 0.,
        edge_symmetric_difference_length=float(CELL*np.count_nonzero(ce^te)),
        candidate_components=candidate['components'],label_components=actual['components'],
        candidate_holes=candidate['holes'],label_holes=actual['holes'],
        topology_exact_match=(candidate['components']==actual['components'] and candidate['holes']==actual['holes']),
        topology_metric_scope='component_and_hole_counts_only',
        target_intrinsic_fixed_in_fitting=dict(radius=.12,contrast=3.,loss=.02),
        intrinsic_label_matches_fitting=all(np.isclose(intrinsic[k],v) for k,v in [('radius',.12),('contrast',3.),('loss',.02)]))


def paid_prediction(world,theta,tx,rx,book,*,finer=False):
    """Serial qualified-family forward only; clear state before next solve."""
    solver=Evaluator(book,'A0',tx=tx,rx=rx,h=.17,n=384 if finer else 256,
                     grid_refinement=3 if finer else 2)
    try:
        ev=solver.evaluate(world,jacobian=False)
        data=np.asarray(ev.data).copy()
        metadata=dict(ev.metadata)
        return (1+theta[4])*np.exp(1j*theta[5])*data,metadata
    finally:
        solver.clear()
        gc.collect()


TABLE_COLUMNS=('observation_id','parent_id','noise_repeat','method','fit_status',
    'selected_fit_success','score_status','finer_status','training_noise_energy',
    'holdout_noise_energy','training_noise_envelope_compatible','holdout_noise_envelope_compatible',
    'both_noise_envelopes_compatible','geometry_feasible','pose_error_lambda0',
    'pose_squared_error_lambda0_squared','material_eps_contrast_L2','material_ell_L2',
    'support_iou','edge_f1','edge_symmetric_difference_length','candidate_components',
    'candidate_holes','topology_exact_match','finer_prediction_difference_white',
    'finer_holdout_noise_energy','actual_nfev','actual_njev','error_type','error')


def persist(out,report):
    write_json(out/'RESULTS.json',report)
    temporary=out/'FULL_TABLE.csv.tmp'
    with temporary.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=TABLE_COLUMNS,extrasaction='ignore')
        writer.writeheader();writer.writerows(report['rows'])
    temporary.replace(out/'FULL_TABLE.csv')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args();root=args.root.resolve()
    # A live RUNNING document must not be scored: prevent label access before
    # the online selection process is terminal and durable.
    fit_report=json.loads((root/'results/a0/optimization/RESULTS.json').read_text())
    if fit_report.get('status') in (None,'RUNNING'):
        raise RuntimeError('Online optimization has not reached a saved terminal state')
    if json.loads((root/'results/e0/gate.json').read_text()).get('decision')!='PASS':
        raise RuntimeError('E0 PASS required before holdout forward evaluation')
    frozen=json.loads((root/'EXPERIMENT_CONTRACT.json').read_text())
    rule=json.loads((root/'results/a0/preparation/DISCREPANCY_RULE.json').read_text())
    upper=float(rule['upper_noise_energy'])
    if rule.get('dimension')!=64 or not np.isfinite(upper) or upper<=0:
        raise ValueError('Expected frozen 64-real-component conservative noise envelope')
    with np.load(root/'inputs/online/A0_OBSERVATIONS.npz',allow_pickle=False) as raw:
        base=raw['base'].copy();parent_ids=raw['parent_id'].copy();repeats=raw['noise_repeat'].copy()
    # From this point labels are allowed solely in this post-fit process.
    labels=json.loads((root/'inputs/offline/CONTINUOUS_TRUTHS.json').read_text())
    label_by_parent={int(row['parent_id']):row for row in labels}
    with np.load(root/'inputs/offline/A0_HOLDOUT.npz',allow_pickle=False) as raw:
        holdout=raw['data'].copy();tx=raw['tx'].copy();rx=raw['rx'].copy()
    _,expected_rx=acquisition(outer=True)
    expected_tx=acquisition(outer=True)[0][2:]
    if base.shape!=(24,2,16) or holdout.shape!=base.shape:
        raise ValueError('Unexpected paired base/holdout dimensions')
    if not np.array_equal(tx,expected_tx) or not np.array_equal(rx,expected_rx):
        raise ValueError('Holdout source/receiver positions differ from frozen outer acquisition')
    out=root/'results/a0/scoring';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook(root/'results/budget_ledger.json')
    began=time.perf_counter()
    rows=[]
    for record in fit_report['records']:
        selected=record.get('selection')
        rows.append(dict(observation_id=int(record['observation_id']),parent_id=int(record['parent_id']),
            noise_repeat=int(record['noise_repeat']),method=record['method'],
            online_record_status=record['status'],fit_status=selected.get('fit_status') if selected else record['status'],
            selected_fit_success=selected.get('success') if selected else None,
            actual_nfev=record.get('actual_nfev'),actual_njev=record.get('actual_njev'),
            all_fit_statuses=[f.get('fit_status') for f in record.get('fits',[])],
            score_status='NOT_RUN' if selected else 'NO_SELECTED_CANDIDATE',
            finer_status='NOT_RUN' if selected and int(record['parent_id']) in FINE_PARENT_IDS else 'NOT_SELECTED_BY_FIXED_PARENT_RULE',
            online_error=record.get('error'),error=record.get('error'),error_type=record.get('error_type'),
            online_failure_preserved=not bool(selected)))
    report=dict(scope='SEPARATE_OFFLINE_SCORING_NOT_FITTING',status='RUNNING',rows=rows,
        online_terminal_status=fit_report['status'],noise_discrepancy_rule=rule,
        finer_parent_rule=dict(parent_ids=list(FINE_PARENT_IDS),noise_repeats='all_saved_repeats',
            methods='all_arms',selection='fixed_parent_ID_no_truth_or_error_selection',
            primary=dict(h=.17,n=256,grid=2),independent=dict(h=.17,n=384,grid=3),
            channels='independent_outer_two_Tx_only',threshold='NO_NEW_GATE_DESCRIPTIVE_CHECK_ONLY'),
        noise_gate_scope='Frozen one-sided conservative noise envelope, not fitted dof or calibrated posterior',
        truth_access='POST_FIT_SCORER_ONLY_NO_SELECTION_FEEDBACK',posterior_claim=False,
        performance_claim=False,scientific_decision='PARENT_REVIEW_REQUIRED',
        budget_start=book.totals('A0'))
    persist(out,report)
    stopped=False
    # First obtain the complete primary holdout table before optional finer
    # audits consume more of the same immutable campaign budget.
    for record,row in zip(fit_report['records'],rows):
        selected=record.get('selection')
        if selected is None:continue
        if stopped:
            row['score_stop_reason']='SHARED_A0_BUDGET_EXHAUSTED';continue
        i=row['observation_id'];before=book.totals('A0')
        try:
            if not 0<=i<len(base) or int(parent_ids[i])!=row['parent_id'] or int(repeats[i])!=row['noise_repeat']:
                raise ValueError('Selected record/online observation identity mismatch')
            theta=np.asarray(selected['theta'],float)
            world=World(selected['support'],theta[0],theta[1],theta[2:4],radius=.12,contrast=3.,target_loss=.02)
            metrics=offline_metrics(selected,label_by_parent[row['parent_id']])
            train_pred=complex_json(selected['pred'])
            train_white=pack(train_pred-base[i]);train_energy=float(train_white@train_white)
            prediction,metadata=paid_prediction(world,theta,tx,rx,book)
            holdout_white=pack(prediction-holdout[i]);energy=float(holdout_white@holdout_white)
            stem=f"observation_{i:02d}_{row['method']}"
            path=out/(stem+'_primary.npz')
            np.savez_compressed(path,theta=theta,support=np.asarray(selected['support']),
                train_prediction=train_pred,train_observation=base[i],train_white=train_white,
                holdout_prediction=prediction,holdout_observation=holdout[i],holdout_white=holdout_white,
                tx=tx,rx=rx,metadata_json=np.array(json.dumps(metadata,default=json_default)))
            row.update(metrics,score_status='SCORED',prediction_path=str(path),
                training_noise_energy=train_energy,holdout_noise_energy=energy,
                training_noise_envelope_compatible=train_energy<=upper,
                holdout_noise_envelope_compatible=energy<=upper,
                both_noise_envelopes_compatible=train_energy<=upper and energy<=upper,
                saved_training_energy_difference=train_energy-float(np.asarray(selected['rawtrainwhite'])@np.asarray(selected['rawtrainwhite'])),
                primary_metadata=metadata,selected_theta=theta,selected_support=selected['support'])
        except BudgetExceeded as exc:
            row.update(score_status='STOPPED_BUDGET',error=str(exc));stopped=True
            report['status']='STOPPED_A0_BUDGET'
        except Exception as exc:
            row.update(score_status='FAILED',error_type=type(exc).__name__,error=str(exc))
        finally:
            row['primary_score_stage_cost']={k:book.totals('A0')[k]-before[k] for k in before}
            persist(out,report)
    for row in rows:
        if row['finer_status']!='NOT_RUN':continue
        if stopped:
            row['finer_stop_reason']='SHARED_A0_BUDGET_EXHAUSTED';continue
        if row['score_status']!='SCORED':
            row['finer_stop_reason']='PRIMARY_SCORE_UNAVAILABLE';continue
        before=book.totals('A0')
        try:
            theta=np.asarray(row['selected_theta'],float)
            world=World(row['selected_support'],theta[0],theta[1],theta[2:4],radius=.12,contrast=3.,target_loss=.02)
            fine,metadata=paid_prediction(world,theta,tx,rx,book,finer=True)
            with np.load(row['prediction_path'],allow_pickle=False) as raw:primary=raw['holdout_prediction'].copy()
            difference=pack(fine-primary)
            white=pack(fine-holdout[row['observation_id']])
            path=out/f"observation_{row['observation_id']:02d}_{row['method']}_independent.npz"
            np.savez_compressed(path,theta=theta,support=np.asarray(row['selected_support']),
                primary_prediction=primary,independent_prediction=fine,prediction_difference_white=difference,
                holdout_white=white,tx=tx,rx=rx,metadata_json=np.array(json.dumps(metadata,default=json_default)))
            row.update(finer_status='SCORED',finer_prediction_difference_white=float(np.linalg.norm(difference)),
                finer_holdout_noise_energy=float(white@white),finer_prediction_path=str(path),finer_metadata=metadata)
        except BudgetExceeded as exc:
            row.update(finer_status='STOPPED_BUDGET',error=str(exc));stopped=True
            report['status']='STOPPED_A0_BUDGET'
        except Exception as exc:
            row.update(finer_status='FAILED',error_type=type(exc).__name__,error=str(exc))
        finally:
            row['finer_score_stage_cost']={k:book.totals('A0')[k]-before[k] for k in before}
            persist(out,report)
    if not stopped:report['status']='SCORING_COMPLETE_PARENT_REVIEW_REQUIRED'
    report.update(budget_end=book.totals('A0'),wall_seconds=time.perf_counter()-began,
        row_counts={status:sum(r['score_status']==status for r in rows) for status in sorted({r['score_status'] for r in rows})},
        nonconverged_selected_count=sum(r['selected_fit_success'] is False for r in rows))
    persist(out,report)
    write_json(out/'gate.json',dict(decision='INCONCLUSIVE',execution_status=report['status'],
        reason='Paired train/holdout envelope, offline errors and fixed numerical audits supplied for parent judgment; no automatic scientific acceptance',
        posterior_claim=False,all_failed_and_unrun_rows_preserved=True))


if __name__=='__main__':
    main()
