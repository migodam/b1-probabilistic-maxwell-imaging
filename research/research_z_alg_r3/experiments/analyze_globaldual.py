"""Post-freeze raw dual endpoints versus full-J correction; no solver reruns.

Raw endpoint candidates already independently validated are reused exactly.
An endpoint merged with a nearby stored candidate is independently solved
instead, with its additional NumPy full RHS charged in a separate receipt.
Truth is used only for the explicit offline RMS scoring after this freeze.
"""
import argparse
import datetime
import json
from pathlib import Path
import time
from collections import Counter
import numpy as np
from validate_predictions import frozen_files, validate

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--freeze-manifest',type=Path,required=True)
    p.add_argument('--validation',type=Path,required=True)
    p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--extra-validation-output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.extra_validation_output.exists():raise FileExistsError('Preserving analysis')
    valid=json.loads(a.validation.read_text())
    scored=json.loads(a.scores.read_text())
    if not scored.get('truth_read_after_prediction_freeze'):raise RuntimeError('Frozen offline scores required')
    scoremap={r['file']:r for r in scored['rows']}
    validmap={r['file']:{v['index']:v for v in r['rows']} for r in valid['rows']}
    started=time.perf_counter();rows=[];extra=[]
    for path,result in frozen_files(a.freeze_manifest):
        if result['method']!='globaldual':continue
        filename=str(path.relative_to(ROOT));score=scoremap[filename];case=score['case_id']
        events=[e for e in result['algebraic_events'] if e.get('method')=='GLOBALDUAL_VARPRO']
        counts=result['accounting']['counts'];fitcounts=Counter()
        for fit in result['fits']:fitcounts.update(fit.get('counts',{}))
        terminal_incomplete=any(f.get('status')=='BUDGET_STOP' for f in result['fits'])
        endpoint_count=sum('candidate' in e for e in events)
        endpoint_full_rhs=counts.get('actual_complex_RHS',0)-fitcounts['actual_complex_RHS']
        endpoint_lu=counts.get('LU',0)-fitcounts['LU']
        exact=(not terminal_incomplete and endpoint_full_rhs==3*endpoint_count and endpoint_lu==endpoint_count)
        with np.load(ROOT/'private_truth'/'holdout'/(case+'.npz'),allow_pickle=False) as f:truth=f['m']
        endpoint_rows=[]
        for event in events:
            if 'endpoint' not in event:continue
            material=np.asarray(event['endpoint']);ci=event.get('candidate')
            endpoint=dict(start=event['start'],dual_objective=event.get('dual_objective'),
                          material_RMS=float(np.sqrt(np.mean((material-truth)**2))),
                          candidate_index=ci)
            if ci is not None and np.array_equal(material,np.asarray(result['candidates'][ci]['m'])):
                check=validmap[filename][ci]
                endpoint.update(relative_residual=check['relative_residual'],validation_reused=True,
                                validation_pass=check['validation_pass'])
            else:
                # The provisional residual is replaced after the independent
                # solve; disagreement against a merged candidate is not a gate.
                stub=dict(candidates=[dict(m=material.tolist(),relative_residual=0.)],top1_index=0,method='globaldual_raw_endpoint',case_id=case)
                checked=validate(path,stub,'holdout','all')
                check=checked['rows'][0]
                passed=bool(check['box_admissible'] and check.get('resolvent_relative_gap',0)>1e-12 and check.get('state_equation_relative',1)<1e-10 and check.get('constitutive_relative',1)<1e-10)
                check['validation_pass']=passed;check['stored_residual_abs_disagreement']=None
                checked.update(globaldual_start=event['start'],all_validated_pass=passed)
                extra.append(checked)
                endpoint.update(relative_residual=check.get('relative_residual'),validation_reused=False,validation_pass=passed)
            fits=[f for f in result['fits'] if f['origin']=='globaldual_fullj_correction_'+str(event['start'])]
            if fits:
                fit=fits[0];corrected=fit.get('candidate')
                endpoint['correction']=dict(status=fit['status'],wall_seconds=fit['wall_seconds'],counts=fit.get('counts'),candidate_index=corrected,
                                            archived_representative_relative_residual=result['candidates'][corrected]['relative_residual'] if corrected is not None else None,
                                            archived_representative_material_RMS=score['all_candidate_RMS'][corrected] if corrected is not None else None,
                                            archive_boundary='A near-duplicate fit endpoint may merge with an earlier candidate; use the exact best-paid-trial trace for actual global improvements')
                improvements=[t for t in result.get('best_so_far',[]) if t['origin']==fit['origin'] and 'm' in t]
                if improvements:
                    paid=min(improvements,key=lambda t:t['relative_residual'])
                    endpoint['correction']['best_paid_global_improvement']=dict(relative_residual=paid['relative_residual'],seconds=paid['seconds'],material_RMS=float(np.sqrt(np.mean((np.asarray(paid['m'])-truth)**2))))
            endpoint_rows.append(endpoint)
        best_raw=min(endpoint_rows,key=lambda e:e['relative_residual']) if endpoint_rows else None
        rows.append(dict(file=filename,case_id=case,started_dual_fits=len(events),completed_raw_endpoints=len(endpoint_rows),
                         phase_count_attribution_exact=exact,
                         raw_endpoint_full_cost=dict(full_complex_rhs=endpoint_full_rhs,full_factorizations=endpoint_lu),
                         fullj_correction_cost=dict(full_complex_rhs=fitcounts['actual_complex_RHS'],full_factorizations=fitcounts['LU'],material_J_builds=counts.get('material_J_builds',0),fit_count=len(result['fits'])),
                         variable_projection_inner_cost=dict(QR_factorizations=counts.get('algebraic_QR',0),complex_LS_RHS=counts.get('algebraic_complex_rhs',0),matrix_shape=[12,9]),
                         correction_fit_wall_seconds=sum(f['wall_seconds'] for f in result['fits']),
                         proposal_setup_endpoints_and_overhead_wall_seconds=result['wall_seconds']-sum(f['wall_seconds'] for f in result['fits']),
                         pure_endpoint_original_residual_top1=best_raw,
                         final_top1_RMS=score['top1_RMS'],final_top1_relative_residual=score['top1_relative_residual'],
                         raw_endpoint_rows=endpoint_rows))
    extra_payload=dict(validated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),freeze_manifest=str(a.freeze_manifest),rows=extra,truth_read=False,
                       purpose='Additional NumPy validation only for exact raw endpoints that differ from their archived near-duplicate; no solver or candidate changes',
                       full_complex_rhs=sum(r['cost']['full_complex_rhs'] for r in extra),full_factorizations=sum(r['cost']['full_factorizations'] for r in extra),wall_seconds=sum(r['cost']['wall_seconds'] for r in extra))
    a.extra_validation_output.write_text(json.dumps(extra_payload,indent=2)+'\n')
    payload=dict(scored_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),freeze_manifest=str(a.freeze_manifest),scoring_only=True,new_physics_rhs=0,
                 extra_physics_receipt=str(a.extra_validation_output),scoring_wall_seconds=max(time.perf_counter()-started-extra_payload['wall_seconds'],0),rows=rows,
                 truth_read_after_prediction_freeze=True,GT_values_disclosed=False,
                 metric_boundary='Pure dual endpoint selection uses independent original-data residual; final ranking includes paid FullJ corrections; inner dual objective is proposal-only')
    a.output.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),phase_counts_exact=all(r['phase_count_attribution_exact'] for r in rows),extra_validation_full_rhs=extra_payload['full_complex_rhs'])))


if __name__=='__main__':main()
