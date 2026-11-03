"""Independent exposed-GT continuation scoring; no forward/gradient/algorithm import."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json
import numpy as np
from offline_score_common import (read,write_new,stamp,unchanged,public_metrics,
    prediction,physical,metrics,write_csv_new,frozen_contract)

ROOT=Path(__file__).resolve().parents[1]
CASE='obs_42cb233e4f224bb3'
METHODS=('onehop','native_extended')


def utc():return datetime.now(timezone.utc).isoformat()


def seal(args):
    public_path=ROOT/'data_public/confirmation'/(CASE+'.npz')
    public=public_metrics(public_path);entries=[]
    for method in METHODS:
        folder=ROOT/'results/failed_endpoint_continuation'/method
        receipt_path=folder/'RECEIPT.json';prediction_path=folder/'prediction.npz'
        rec=read(receipt_path)
        if rec['status']!='SEALED_POST_CONFIRMATION_DIAGNOSTIC':raise RuntimeError('Unclosed continuation: '+method)
        m,yp=prediction(prediction_path,(24,24))
        entries.append(dict(method=method,receipt=stamp(receipt_path),prediction=stamp(prediction_path),
            status=rec['status'],original_Phi_pretruth=physical(m,yp,public,(24,24))['physical_objective']))
    result=dict(status='TWO_CONTINUATION_DIAGNOSTICS_SEALED',sealed_utc=utc(),entries=entries,
        case_id=CASE,classification='POST_CONFIRMATION_BUDGET_DIAGNOSTIC',
        public=stamp(public_path),original18=stamp(ROOT/'results/CONFIRMATION_SCORES.json'),
        old13=stamp(ROOT/'results/DEVELOPMENT_SCORES.json'),
        sources=[stamp(Path(__file__)),stamp(Path(__file__).with_name('offline_score_common.py'))],
        truth_opened_this_pass=False,truth_already_exposed=True,append_to_confirmation=False,
        new_physical_LU=0,new_complex_RHS=0,hash_checks=False)
    write_new(args.seal,result)
    print(json.dumps(dict(status=result['status'],endpoints=2,truth_opened=False)))


def score(args):
    if not args.authorize_exposed_truth:raise PermissionError('Parent must authorize these2 sealed continuation endpoints')
    outputs=[args.output.with_suffix('.'+ext) for ext in ('json','csv','md')]
    if any(path.exists() for path in outputs):raise FileExistsError('Preserve prior continuation score; use a new version')
    closed=read(args.seal)
    if closed['status']!='TWO_CONTINUATION_DIAGNOSTICS_SEALED' or len(closed['entries'])!=2:raise RuntimeError('Complete2-endpoint seal required')
    evidence=[closed['public'],closed['original18'],closed['old13']]+closed['sources']
    evidence += [item for entry in closed['entries'] for item in (entry['receipt'],entry['prediction'])]
    if not all(unchanged(item) for item in evidence):raise RuntimeError('Sealed evidence changed')
    opened=utc()
    with np.load(ROOT/'sealed_truth/confirmation'/(CASE+'.npz'),allow_pickle=False) as data:
        truth=data['m'].copy();support=data['mask'].ravel().astype(bool)
    public=public_metrics(closed['public']['path'])
    original=read(ROOT/'results/CONFIRMATION_SCORES.json');rows=[]
    for entry in closed['entries']:
        rec=read(entry['receipt']['path']);method=entry['method']
        m,yp=prediction(entry['prediction']['path'],(24,24))
        base=next(row for row in original['rows'] if row['observation_id']==CASE and row['method']==method)
        cfg=rec['configuration'];pg=float(rec['projected_gradient_inf']);target=float(cfg['original_gtol'])
        row=dict(method=method,observation_id=CASE,classification='POST_CONFIRMATION_BUDGET_DIAGNOSTIC',
            **metrics(m,yp,public,truth,support,(24,24)))
        row.update(original_confirmation_nrmse=base['amplitude_normalized_nrmse'],
            original_confirmation_iou=base['otsu_iou'],original_confirmation_relative_residual=base['relative_residual'],
            original_confirmation_Phi=base['physical_objective'],
            delta_nrmse_from_original=row['amplitude_normalized_nrmse']-base['amplitude_normalized_nrmse'],
            delta_iou_from_original=row['otsu_iou']-base['otsu_iou'],
            initial_original_Phi=rec['initial']['physical_Phi'],reported_final_original_Phi=rec['physical_objective'],
            original_Phi_decrease=rec['initial']['physical_Phi']-row['physical_objective'],
            physical_objective_receipt_delta=row['physical_objective']-rec['physical_objective'],
            relative_residual_receipt_delta=row['relative_residual']-rec['physical_relative_residual'],
            initial_projected_gradient_inf=rec['initial']['projected_gradient_inf'],
            projected_gradient_inf=pg,original_projected_gradient_target=target,
            original_gradient_precision_met=pg<=target,
            gradient_evidence='Receipt and already stored gradient vector; no gradient recomputation by scorer',
            extra_corrector_allowance=cfg['extra_corrector_maxiter'],
            actual_extra_corrector_iterations=rec['corrector']['iterations'],
            actual_extra_corrector_evaluations=rec['corrector']['evaluations'],
            stop=rec['corrector']['stop'],physical_cost=rec['physical_cost'],
            wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],failed_attempts=rec['failed_attempts'],
            receipt=entry['receipt']['path'],prediction=entry['prediction']['path'],
            qualification='Original bound-projected first-order tolerance attained; no local-minimum/uniqueness certificate' if pg<=target else 'Original bound-projected first-order tolerance not attained; function-reduction stop is not KKT qualification')
        rows.append(row)
    result=dict(status='DEVELOPMENT_CONTINUATION_SCORED',classification='POST_CONFIRMATION_BUDGET_DIAGNOSTIC_NOT_BLIND',
        seal=str(args.seal),first_exposed_truth_opened_this_pass_utc=opened,rows=rows,contract=frozen_contract(),
        original18_preserved=unchanged(closed['original18']),old13_preserved=unchanged(closed['old13']),
        original_confirmation_denominator=6,original_confirmation_tasks=18,original_joint_counts_unchanged=True,
        append_to_confirmation=False,no_parameter_tuning=True,
        budget='At most4000 extra original corrector iterations per assigned endpoint, may stop earlier on original rules; no further extension',
        extra_original_run_physical_cost={key:sum(row['physical_cost'].get(key,0) for row in rows) for key in ('LU','source_RHS','adjoint_RHS','actual_complex_RHS')},
        scoring_new_physical_LU=0,scoring_new_complex_RHS=0,hash_checks=False,
        scientific_scope='Single exposed case and frozen additional allowance, descriptive stop/first-order precision/material metrics; no general convergence or local-minimum claim',
        subsequent_scientific_experiments='STOPPED_BY_PARENT_REQUEST',scientific_judgment='Reserved for parent Codex')
    write_new(outputs[0],result);write_csv_new(outputs[1],rows)
    lines=['# Post-confirmation continuation scoring V1','',
        '**POST_CONFIRMATION_BUDGET_DIAGNOSTIC / exposed case42cb, not blind.** Both endpoints were sealed before this scoring pass reopened the already exposed GT. Original13/18 score tables and all original confirmation predictions/4-of-6 denominators remain unchanged.','',
        'Each assigned endpoint had an allowance of **at most4000** extra original corrector iterations; original FTOL/PGTOL rules could stop it earlier. No added parameter search or further allowance.','',
        '| Endpoint | RMSE/12 | IoU | Relative residual | White residual | Original Phi | Joint |',
        '|---|---:|---:|---:|---:|---:|---|']
    for row in rows:lines.append(f"|{row['method']}|{row['amplitude_normalized_nrmse']:.9f}|{row['otsu_iou']:.9f}|{row['relative_residual']:.12f}|{row['whitened_residual']:.6f}|{row['physical_objective']:.12g}|{'PASS' if row['joint_gate'] else 'FAIL'}|")
    lines+=['','| Endpoint | m min / max | Box violation | Extra iterations / allowance | Stop | Projected gradient inf | Original target met |',
        '|---|---:|---:|---:|---|---:|---|']
    for row in rows:
        stop='FTOL relative reduction' if 'RELATIVE REDUCTION' in row['stop'] else 'PGTOL projected gradient'
        lines.append(f"|{row['method']}|{row['material_min']:.9f} / {row['material_max']:.9f}|{row['material_constraint_violation']:.3g}|{row['actual_extra_corrector_iterations']} / {row['extra_corrector_allowance']}|{stop}|{row['projected_gradient_inf']:.12g}|{row['original_gradient_precision_met']} (1e-8)|")
    lines+=['','| Endpoint | Original Phi before → after | NRMSE before → after | Extra LU | Extra physical RHS | Wall / CPU s |',
        '|---|---:|---:|---:|---:|---:|']
    for row in rows:
        lines.append(f"|{row['method']}|{row['initial_original_Phi']:.12g} → {row['physical_objective']:.12g}|{row['original_confirmation_nrmse']:.9f} → {row['amplitude_normalized_nrmse']:.9f}|{row['physical_cost']['LU']}|{row['physical_cost']['actual_complex_RHS']}|{row['wall_seconds']:.6f} / {row['cpu_seconds']:.6f}|")
    lines+=['',
        'onehop stops on relative function reduction before attaining the original1e-8 projected-gradient precision; that stop does not qualify KKT accuracy. native_extended attains that numerical bound-projected first-order tolerance. This is a precision statement at a single endpoint, not a local-minimum, uniqueness, or general algorithm qualification. All material/data/image gates remain independently reported above.','',
        'Projected gradients and stop reasons come from the existing numerical receipts; the scorer performs no new forward or derivative. Extra-run physical totals use each receipt top-level physical_cost once, including its endpoint diagnostics; nested corrector costs and previous6-LU stationarity diagnostics are not added again. The scoring itself adds0 LU/0 RHS.','',
        'No new result replaces any of the original18 confirmation endpoints or secondary choices. This exposed budget diagnostic cannot be counted as independent confirmation. The assigned scientific experiment sequence ends here; final scientific interpretation belongs to parent Codex.']
    with outputs[2].open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],rows=[{key:row[key]for key in ('method','amplitude_normalized_nrmse','otsu_iou','relative_residual','whitened_residual','physical_objective','material_constraint_violation','projected_gradient_inf','original_gradient_precision_met','actual_extra_corrector_iterations','stop','joint_gate')}for row in rows],extra_physical_cost=result['extra_original_run_physical_cost'],original18_preserved=result['original18_preserved'],scoring_new_physical_LU=0)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['seal','score'])
    parser.add_argument('--seal',type=Path,default=ROOT/'audit/DEVELOPMENT_CONTINUATION_SEALED_V1.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/DEVELOPMENT_CONTINUATION_SCORES_V1')
    parser.add_argument('--authorize-exposed-truth',action='store_true')
    args=parser.parse_args()
    if args.operation=='seal':seal(args)
    else:score(args)
