#!/usr/bin/env python3
"""Post-seal selection geometry with cached coarse GT reference; ZERO new forward.

Requires the new264-task approved campaign seal and completed ALL_VIEWS.json.
Reuses old_endpoint_geometry cached m_gt/Y_gt_coarse only. No direct truth NPZ
loading, candidate re-selection, threshold changes, online predictions or SVD.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
START_CPU=time.process_time();START_WALL=time.perf_counter()
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1];WORKTREE=ROOT.parents[1]
sys.path[:0]=[str(Path(__file__).resolve().parent),str(WORKTREE)]
import analyze_campaign as campaign
base=campaign.base


SCOPE='OFFLINE_EXISTENCE_DIAGNOSTIC'
LIMIT=('Two distinct archived legal/eligible predictions inside the same finite discrete noise ball '
       'only witness finite discretized noise ambiguity at that radius. They do not establish '
       'continuous-Maxwell or GT nonidentifiability, or an online candidate selector.')


def number(value):return float(value) if base.finite(value) else None


def divide(n,d):
    if not base.finite(n) or not base.finite(d):return None,'MISSING_OR_NONFINITE_VALUE'
    if d<=0:return None,'NONPOSITIVE_DENOMINATOR'
    return n/d,None


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--seal',type=Path,required=True,help='Same new264-task parent-approved campaign seal used for main analysis')
    parser.add_argument('--out',type=Path)
    args=parser.parse_args();root=args.root.resolve();seal_path=args.seal.resolve()
    # All authorization and prediction inventory checks precede any GT-derived load.
    seal,cfg,tasks,account,jobs=campaign.authorize(root,seal_path)
    analysis=root/'analysis'/campaign.CAMPAIGN
    if not (analysis/'ALL_VIEWS.json').is_file():raise FileNotFoundError('Completed sealed main ALL_VIEWS.json required')
    main_receipt=base.read(analysis/'RECEIPT.json')
    if main_receipt.get('campaign')!=campaign.CAMPAIGN or main_receipt.get('expected_tasks')!=264 or Path(main_receipt.get('seal','')).resolve()!=seal_path:raise ValueError('Main analysis must refer to this approved264-task seal')
    rows=base.read(analysis/'ALL_VIEWS.json')
    lookup={}
    for r in rows:
        if r['view'] in base.VIEWS:
            k=(r['case_id'],r['solver'],r['view'])
            if k in lookup:raise ValueError('Duplicate selected main-analysis row: '+str(k))
            lookup[k]=r
    out=(args.out or root/'analysis/selection_geometry_v1').resolve()
    if out!=root/'analysis/selection_geometry_v1' or out.exists():raise ValueError('New canonical analysis/selection_geometry_v1 output required')
    out.mkdir(parents=True)
    receipt=dict(status='RUNNING',scope='Posthoc exposed24 frozen selections; cached offline reference only',
        authorization=str(seal_path),main_analysis=str(analysis),expected_tasks=264,expected_cases=24,
        truth_source='Existing authorized audit/old_endpoint_geometry_v1/cases/*/COARSE_GT_REFERENCE.npz; no direct sealed_truth access',
        new_forward=0,new_LU=0,new_complex_RHS=0,new_SVD=0,new_optimization=0,new_candidate_selection=False,
        interpretation_limit=LIMIT,BLAS_threads=1,errors=[])
    refs=[];comparisons=[];views=[];pairs=[];summaries=[];issues=[];model_costs=[]
    try:
        import numpy as np
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        def norm(x):return float(np.linalg.norm(x))
        def squared(x):return float(np.vdot(x,x).real)
        for cid in campaign.CASES:
            case_tasks=[t for t in tasks if t['case_id']==cid]
            reference=dict(case_id=cid,status='UNAVAILABLE',D=None,R=None,Phi=None,delta=None,sigma=None,
                reference_scope='Coarse inverse-grid material state only; not replacement of finer-grid acquisition GT',
                crossgrid_reference_is_noise_certificate=False)
            public=None;model=None;public_error=None;regularizer_evaluations=0
            cached=root/'audit/old_endpoint_geometry_v1/cases'/cid/'COARSE_GT_REFERENCE.npz'
            public_path=root/case_tasks[0]['input_file']
            try:
                with np.load(public_path,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
                model=MaxwellModel(public)  # constructor/import only; no state/forward/gradient call
                Y=public['Y'];sigma=float(public['sigma']);delta=float(public['delta']);Y2=squared(Y)
                reference.update(public_input=str(public_path),cached_reference=str(cached),delta=delta,sigma=sigma,
                    inverse_side=int(public.get('inverse_side',tuple(public['grid_shape'])[0])),
                    forward_side=int(public.get('forward_side',tuple(public['grid_shape'])[0])),
                    observation_contract=case_tasks[0].get('observation_contract'),public_Y_shape=list(Y.shape))
            except Exception as exc:
                public_error=repr(exc);issues.append(dict(case_id=cid,issue='PUBLIC_INPUT_OR_REGULARIZER_UNAVAILABLE',error=public_error))
            if public_error is None:
                try:
                    # This is the only GT-derived array read. No new GT forward.
                    with np.load(cached,allow_pickle=False) as z:mgt=z['m_gt'].copy();Ygt=z['Y_gt_coarse'].copy()
                    if mgt.shape!=(model.n,) or Ygt.shape!=Y.shape:raise ValueError('Cached coarse reference dimensions differ from public observation')
                    if not np.isfinite(mgt).all() or not np.isfinite(Ygt).all():raise ValueError('Nonfinite cached coarse reference')
                    residual=norm(Ygt-Y);D,Dreason=divide(.5*residual**2,Y2)
                    regularizer_evaluations+=1
                    prior=float(model.reg.evaluate(mgt)[0])
                    residual_over_delta,delta_reason=divide(residual,delta)
                    reference.update(status='AVAILABLE',D=D,D_null_reason=Dreason,R=prior,Phi=D+prior if D is not None else None,
                        Phi_null_reason=Dreason,reference_raw_residual_L2=residual,reference_residual_over_delta=residual_over_delta,
                        delta_ratio_null_reason=delta_reason,reference_whitened_residual_L2=divide(residual,sigma)[0],
                        coarse_reference_inside_original_discrete_noise_ball=bool(residual<=delta) if delta>0 else None,
                        ball_membership_null_reason=None if delta>0 else 'NONPOSITIVE_DELTA',
                        qualification='Numerical membership of this cached inverse-grid reference only; crossgrid mismatch is not noise-certified',
                        regularizer_source='Original MaxwellModel.reg.evaluate default smooth=tv=1e-4 epsilon=.01')
                except Exception as exc:
                    reference['unavailable_reason']=repr(exc);issues.append(dict(case_id=cid,issue='CACHED_COARSE_GT_REFERENCE_UNAVAILABLE',error=repr(exc)))
            else:reference['unavailable_reason']=public_error
            refs.append(reference)
            for task in case_tasks:
                method=task['method'];ident=task['task_id'];folder=root/'RESULTS'/campaign.CAMPAIGN/ident
                secondary=lookup.get((cid,method,'secondary'))
                for view in ('primary','raw','corrected','secondary'):
                    selected=lookup.get((cid,method,view));record=dict(case_id=cid,solver=method,view=view,task_id=ident,
                        cohort_denominator=24,total_task_denominator=264,selected_row_present=selected is not None,
                        GT_reference_status=reference['status'],GT_D=reference['D'],GT_R=reference['R'],GT_Phi=reference['Phi'],
                        reference_scope=reference['reference_scope'],noise_certificate_claim=False)
                    if selected:
                        record.update(run_status=selected.get('run_status'),available=selected.get('available'),selected_candidate=selected.get('selection_candidate'),
                            source_file=selected.get('source_file'),source_view=selected.get('source_view'),secondary_fallback=selected.get('secondary_fallback'),
                            selected_legal=selected.get('legal'),selected_eligible=selected.get('eligible'),selected_noise_feasible=selected.get('noise_feasible'))
                    for metric in ('D','R','Phi'):
                        selected_value=number(selected.get(metric)) if selected else None;refvalue=number(reference.get(metric))
                        record['selected_'+metric]=selected_value
                        record['delta_'+metric+'_selected_minus_GT_reference']=selected_value-refvalue if selected_value is not None and refvalue is not None else None
                        record['delta_'+metric+'_null_reason']=None if selected_value is not None and refvalue is not None else 'MISSING_SELECTED_ROW_METRIC_OR_CACHED_REFERENCE'
                    comparisons.append(record)
                unique={};conflicts=set();task_views=[]
                for view in base.VIEWS:
                    source_view=view;fallback=None;source=folder/(view+'.npz');metadata=folder/(view+'.json')
                    # Archived secondary deployment fallback repeats the same primary ID,
                    # never creates another pair or changes an online selection.
                    if view=='secondary' and not source.is_file() and secondary and secondary.get('secondary_fallback'):
                        source_view='primary';source=folder/'primary.npz';metadata=folder/'primary.json';fallback=secondary['secondary_fallback']
                    record=dict(case_id=cid,solver=method,task_id=ident,view=view,source_view=source_view,source_file=str(source),
                        secondary_fallback=fallback,source_exists=source.is_file(),metadata_exists=metadata.is_file(),
                        scope='Only five fixed sealed views; not all paid candidates',candidate=None,pair_eligible=False,exclusion_reason=None,
                        material_R_from_saved_m=None,actual_residual_L2=None,actual_residual_over_delta=None,actual_inside_original_ball=None)
                    try:
                        if not source.is_file() or not metadata.is_file():raise FileNotFoundError('Saved view/metadata absent')
                        meta=base.read(metadata);candidate=meta.get('candidate')
                        record.update(candidate=candidate,archived_legal=meta.get('legal'),archived_eligible=meta.get('eligible'),
                            archived_noise_feasible=meta.get('noise_feasible'),archived_CPU_seconds=meta.get('CPU_seconds'))
                        with np.load(source,allow_pickle=False) as z:
                            m=z['m'].copy() if 'm' in z else None
                            prediction=z['Y_prediction'].copy() if 'Y_prediction' in z else None
                        record['missing_array_fields']=[k for k,v in [('m',m),('Y_prediction',prediction)] if v is None]
                        if m is None or prediction is None:raise ValueError('Required actual m/Y_prediction missing; no forward reconstruction allowed')
                        if public_error is not None:raise ValueError('Public Y/delta unavailable: '+public_error)
                        m=np.asarray(m).reshape(-1)
                        if m.shape!=(model.n,) or prediction.shape!=Y.shape or np.iscomplexobj(m) or not np.isfinite(m).all() or not np.isfinite(prediction).all():raise ValueError('Invalid saved material/data array shape or value')
                        residual=norm(prediction-Y);ratio,why=divide(residual,delta)
                        actual_legal=bool(m.min()>=float(public['lower'])-1e-12 and m.max()<=float(public['upper'])+1e-12)
                        actual_inside=bool(residual<=delta) if delta>0 else None
                        regularizer_evaluations+=1
                        record.update(material_R_from_saved_m=float(model.reg.evaluate(m)[0]),actual_residual_L2=residual,
                            actual_residual_over_delta=ratio,delta_ratio_null_reason=why,actual_inside_original_ball=actual_inside,
                            actual_material_legal=actual_legal)
                        # Same original legality/eligible120CPU rule and same original delta.
                        # These restrict the existence statement, not the campaign/algorithm.
                        allowed=base.eligible(meta) and actual_legal and actual_inside is True
                        record['pair_eligible']=allowed
                        record['exclusion_reason']=None if allowed else 'NOT_ORIGINAL_LEGAL_ELIGIBLE_AND_ACTUALLY_INSIDE_SAVED_NOISE_BALL'
                        if allowed:
                            if candidate in unique:
                                old=unique[candidate]
                                old['view_labels'].append(view)
                                if not np.array_equal(m,old['m']) or not np.array_equal(prediction,old['prediction']):
                                    conflicts.add(candidate);issues.append(dict(task_id=ident,candidate=candidate,issue='DUPLICATE_CANDIDATE_ID_WITH_DIFFERENT_SAVED_ARRAYS'))
                            else:unique[candidate]=dict(candidate=candidate,m=m,prediction=prediction,view_labels=[view],source_file=str(source))
                    except Exception as exc:
                        record.update(pair_eligible=False,exclusion_reason=repr(exc))
                        issues.append(dict(task_id=ident,view=view,issue='SAVED_VIEW_GEOMETRY_UNAVAILABLE',error=repr(exc)))
                    task_views.append(record)
                for record in task_views:
                    if isinstance(record['candidate'],int) and record['candidate'] in conflicts:record.update(pair_eligible=False,exclusion_reason='DUPLICATE_CANDIDATE_ARRAY_CONFLICT')
                views+=task_views
                distinct=[entry for candidate,entry in sorted(unique.items()) if candidate not in conflicts]
                task_pairs=[]
                for a,b in combinations(distinct,2):
                    data_distance=norm(a['prediction']-b['prediction']);distance_ratio,why=divide(data_distance,float(public['delta']))
                    material_rms=float(np.sqrt(np.mean((a['m']-b['m'])**2)))
                    record=dict(case_id=cid,solver=method,task_id=ident,candidate_a=a['candidate'],candidate_b=b['candidate'],
                        views_a=a['view_labels'],views_b=b['view_labels'],material_RMS_distance=material_rms,data_L2_distance=data_distance,
                        data_distance_over_delta=distance_ratio,data_ratio_null_reason=why,delta=float(public['delta']),
                        both_in_same_original_discrete_ball=True,scope=SCOPE,pair_uses_GT=False,online_selector=False,interpretation_limit=LIMIT)
                    task_pairs.append(record);pairs.append(record)
                maximum=max(task_pairs,key=lambda r:(r['material_RMS_distance'],-r['candidate_a'],-r['candidate_b'])) if task_pairs else None
                data_maximum=max(task_pairs,key=lambda r:(r['data_L2_distance'],-r['candidate_a'],-r['candidate_b'])) if task_pairs else None
                summaries.append(dict(case_id=cid,solver=method,task_id=ident,run_status=account[ident]['status'],cohort_denominator=24,total_task_denominator=264,
                    fixed_views_inspected=5,pair_eligible_view_count=sum(r['pair_eligible'] for r in task_views),distinct_eligible_inball_candidate_count=len(distinct),
                    distinct_inball_candidate_ids=[r['candidate'] for r in distinct],enumerated_candidate_pairs=len(task_pairs),
                    missing_or_excluded_views=[dict(view=r['view'],reason=r['exclusion_reason']) for r in task_views if not r['pair_eligible']],
                    secondary_fallback=secondary.get('secondary_fallback') if secondary else None,
                    max_material_distance_pair=maximum,max_data_distance_pair=data_maximum,
                    maximum_null_reason=None if maximum else 'FEWER_THAN_TWO_DISTINCT_SAVED_LEGAL_ELIGIBLE_INBALL_CANDIDATES',
                    maximization_rule='Maximum material RMS distance; additionally maximum data distance. Deterministic candidate-ID tie break; neither uses GT or selects an inverse output.',scope=SCOPE,interpretation_limit=LIMIT))
            if model is not None:model_costs.append(dict(case_id=cid,cost=dict(model.cost),regularizer_evaluations=regularizer_evaluations,scope='Model constructor plus reg.evaluate only; no physical state or gradient call'))
        for name,data in [('GT_REFERENCE24',refs),('SELECTED_VS_GT_REFERENCE',comparisons),('SEALED_VIEW_GEOMETRY',views),('FEASIBLE_SAVED_VIEW_PAIRS',pairs),('TASK_MAX_DISTANCE_PAIRS264',summaries),('MODEL_CONSTRUCTION_COST',model_costs),('ISSUES',issues)]:base.table(out,name,data)
        lu=sum(r['cost'].get('LU',0) for r in model_costs);rhs=sum(r['cost'].get('actual_complex_RHS',0) for r in model_costs)
        receipt.update(status='COMPLETE_WITH_RETAINED_ISSUES' if issues else 'COMPLETE',reference_rows=len(refs),selected_comparison_rows=len(comparisons),
            task_summary_rows=len(summaries),sealed_view_rows=len(views),distinct_candidate_pair_rows=len(pairs),issues=len(issues),
            actual_new_LU=lu,actual_new_complex_RHS=rhs,regularizer_evaluations=sum(r['regularizer_evaluations'] for r in model_costs),model_construction_cost=model_costs,
            eligibility_scope='Only original five sealed views and original legal/eligible120CPU rule plus actual saved Y_prediction residual <= original public delta. Missing fields remain null with reason. No added success threshold or campaign gate.',
            reference_scope='Cached inverse-grid coarse GT D/R/Phi and raw ball membership; no crossgrid noise certificate or new GT forward',
            candidate_pair_scope=LIMIT,complete_task_denominator=264,source_predictions_modified=False)
        report=['# Frozen-view selection geometry','',
            'Post-seal, already exposed same-cohort offline description. All264 task rows and24 cached coarse reference rows are retained. No candidate is changed or reselected.','',
            f"Selected primary/raw/corrected/secondary comparison rows: {len(comparisons)}. Fixed-view inspection rows: {len(views)}. Distinct eligible in-ball candidate pairs: {len(pairs)}.",
            '', 'GT_REFERENCE24 contains inverse-grid reference D=.5||Ygt_coarse−Y||²/||Y||², original prior R and Phi=D+R, and residual/delta membership. The reference is reused from old_endpoint_geometry; it is not a recomputed GT forward and does not convert crossgrid mismatch into a noise certificate.',
            '', 'SELECTED_VS_GT_REFERENCE uses the main analysis archived selections unchanged, including explicit secondary→primary fallback. Signed deltas are selected minus cached inverse-grid reference; missing values have reasons.',
            '', 'FEASIBLE_SAVED_VIEW_PAIRS enumerates every distinct candidate-ID pair in the five fixed sealed views only. It uses actual m/Y_prediction, original legality/eligibility and the original public ball radius. Duplicated view aliases/fallbacks do not create extra candidates. No GT is used to choose the maximum-distance pair.',
            '', LIMIT,'',
            'Maximum material-distance and maximum data-distance pairs are OFFLINE_EXISTENCE_DIAGNOSTIC records, never online selectors. Absence of two stored in-ball views is only absence of this particular archived witness, not a uniqueness certificate.',
            '', f'Actual new physical LU={lu}; actual new physical complexRHS={rhs}; no state/data_value_gradient/forward/SVD/optimizer calls. CPU and wall are recorded in RECEIPT.json. All scientific interpretation remains with the parent.','']
        (out/'REPORT.md').write_text('\n'.join(report),encoding='utf8')
    except BaseException as exc:
        receipt.update(status='FAILED',error=repr(exc));raise
    finally:
        receipt.update(CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,whole_process_CPU_seconds_at_receipt=time.process_time(),created_utc=datetime.now(timezone.utc).isoformat())
        base.write_json(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k] for k in ('status','reference_rows','selected_comparison_rows','task_summary_rows','distinct_candidate_pair_rows','actual_new_LU','actual_new_complex_RHS','CPU_seconds','wall_seconds')},indent=2))


if __name__=='__main__':main()
