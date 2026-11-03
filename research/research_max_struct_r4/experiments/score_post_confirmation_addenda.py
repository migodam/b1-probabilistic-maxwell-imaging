"""Separate exposed-development addenda; no physics and no old-table mutation."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json
import numpy as np
from offline_score_common import (read,write_new,stamp,unchanged,public_metrics,
    prediction,physical,metrics,write_csv_new,frozen_contract)

ROOT=Path(__file__).resolve().parents[1]
R2=ROOT.parent/'research_max_struct_r2'
CASE='obs_b20d5a55677d4eca'
FRESH='obs_42cb233e4f224bb3'
SPECS=(
    ('partition_average_R2',CASE,'results/partition_average/'+CASE,'POST_CONFIRMATION_DEVELOPMENT',
     R2/'data_public'/(CASE+'.npz'),R2/'sealed_truth'/(CASE+'.npz')),
    ('partition_average_42cb',FRESH,'results/partition_average/'+FRESH,'POST_CONFIRMATION_DEVELOPMENT',
     ROOT/'data_public/confirmation'/(FRESH+'.npz'),ROOT/'sealed_truth/confirmation'/(FRESH+'.npz')),
    ('onehop_Windows_R2',CASE,'results/cross_platform_remote/development_windows','DISCLOSED_DEVELOPMENT_CROSS_PLATFORM_VALIDATION',
     R2/'data_public'/(CASE+'.npz'),R2/'sealed_truth'/(CASE+'.npz')))


def utc():return datetime.now(timezone.utc).isoformat()


def seal(args):
    entries=[]
    for method,ident,folder,classification,public_path,truth_path in SPECS:
        receipt_path=ROOT/folder/'RECEIPT.json';prediction_path=ROOT/folder/'prediction.npz'
        rec=read(receipt_path)
        if rec['status'] not in ('SEALED_DEVELOPMENT_PREDICTION','SEALED_PREDICTION'):
            raise RuntimeError('Unclosed development output: '+method)
        m,yp=prediction(prediction_path,(24,24));public=public_metrics(public_path)
        entries.append(dict(method=method,observation_id=ident,classification=classification,
            receipt=stamp(receipt_path),prediction=stamp(prediction_path),public=stamp(public_path),
            status=rec['status'],original_Phi_pretruth=physical(m,yp,public,(24,24))['physical_objective']))
    result=dict(status='THREE_EXPOSED_DEVELOPMENT_ADDENDA_SEALED',sealed_utc=utc(),entries=entries,
        old13_table=stamp(ROOT/'results/DEVELOPMENT_SCORES.json'),
        original_confirmation_table=stamp(ROOT/'results/CONFIRMATION_SCORES.json'),
        sources=[stamp(Path(__file__)),stamp(Path(__file__).with_name('offline_score_common.py'))],
        truth_opened_by_this_pass=False,truth_already_exposed=True,
        original_confirmation_denominator=6,append_to_confirmation=False,
        no_candidate_reselection=True,new_physical_LU=0,new_complex_RHS=0,hash_checks=False)
    write_new(args.seal,result)
    print(json.dumps(dict(status=result['status'],count=3,truth_opened=False)))


def score(args):
    if not args.authorize_exposed_truth:raise PermissionError('Parent authorization for these3 sealed exposed-development endpoints required')
    outputs=[args.output.with_suffix('.'+ext) for ext in ('json','csv','md')]
    if any(p.exists() for p in outputs):raise FileExistsError('Preserve previous addenda; supply a new versioned output')
    closed=read(args.seal)
    if closed['status']!='THREE_EXPOSED_DEVELOPMENT_ADDENDA_SEALED' or len(closed['entries'])!=3:
        raise RuntimeError('Complete3-output seal required')
    evidence=[closed['old13_table'],closed['original_confirmation_table']]+closed['sources']
    evidence += [item for row in closed['entries'] for item in (row['receipt'],row['prediction'],row['public'])]
    if not all(unchanged(item) for item in evidence):raise RuntimeError('Sealed evidence changed')
    opened=utc();rows=[]
    old13=read(ROOT/'results/DEVELOPMENT_SCORES.json')
    confirmation=read(ROOT/'results/CONFIRMATION_SCORES.json')
    for spec,entry in zip(SPECS,closed['entries']):
        method,ident,folder,classification,public_path,truth_path=spec
        with np.load(truth_path,allow_pickle=False) as data:
            truth=data['m'].copy();support=data['mask'].ravel().astype(bool)
        rec=read(entry['receipt']['path']);public=public_metrics(public_path)
        m,yp=prediction(entry['prediction']['path'],(24,24))
        original=next(row for row in (old13['rows'] if ident==CASE else confirmation['rows'])
                      if row.get('arm',row.get('method'))==('onehop_absorption' if ident==CASE else 'onehop')
                      and (ident==CASE or row['observation_id']==ident))
        row=dict(method=method,observation_id=ident,classification=classification,
            truth_source=str(truth_path),**metrics(m,yp,public,truth,support,(24,24)))
        row.update(original_onehop_reference={key:original[key] for key in ('amplitude_normalized_nrmse','otsu_iou','relative_residual','physical_objective')},
            delta_nrmse_from_original_onehop=row['amplitude_normalized_nrmse']-original['amplitude_normalized_nrmse'],
            delta_iou_from_original_onehop=row['otsu_iou']-original['otsu_iou'],
            delta_relative_residual_from_original_onehop=row['relative_residual']-original['relative_residual'],
            delta_Phi_from_original_onehop=row['physical_objective']-original['physical_objective'],
            physical_objective_receipt_delta=row['physical_objective']-rec['physical_objective'],
            relative_residual_receipt_delta=row['relative_residual']-rec['physical_relative_residual'],
            physical_model_cost=rec['physical_model_cost'],wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],
            caller_wall_seconds=rec.get('caller_wall_seconds'),solver_import_charge=rec.get('solver_import_charge'),
            auxiliary_operation_counts=rec.get('auxiliary_operation_counts'),
            corrector_stop=rec.get('corrector_stop',rec.get('corrector',{}).get('stop')),
            cost_scope=rec['cost_scope'],failed_attempts=rec.get('failed_attempts'),
            receipt=entry['receipt']['path'],prediction=entry['prediction']['path'])
        rows.append(row)
    result=dict(status='EXPOSED_DEVELOPMENT_ADDENDA_SCORED',seal=str(args.seal),
        classification='SEPARATE_POST_CONFIRMATION_AND_CROSS_PLATFORM_DEVELOPMENT_NOT_BLIND',
        first_exposed_truth_opened_this_pass_utc=opened,rows=rows,contract=frozen_contract(),
        old13_table_preserved=unchanged(closed['old13_table']),original18_table_preserved=unchanged(closed['original_confirmation_table']),
        original_confirmation_joint_counts={row['method']:row['joint_pass'] for row in confirmation['summaries']},
        original_confirmation_denominator=6,original_confirmation_tasks=18,append_to_confirmation=False,
        fixed_case_selection='42cb was selected from public unqualified objects by minimum opaque ID, not GT material ranking',
        scope='Different mechanism revision and disclosed cross-platform replay; no independent confirmation claim or inherited6-object success-rate update',
        new_physical_LU=0,new_complex_RHS=0,hash_checks=False)
    write_new(outputs[0],result);write_csv_new(outputs[1],rows)
    lines=['# Separate post-confirmation development addenda V1','',
        '**Not blind; not part of the frozen confirmation denominator.** Both partition-average runs are POST_CONFIRMATION_DEVELOPMENT. Windows onehop is a disclosed original-R2 cross-platform replay. All3 predictions were sealed before this pass reopened their already exposed GT. Original13-row development and18-row confirmation tables remain unchanged; all original method/secondary counts remain4/6.','',
        '| Endpoint | Object | RMSE/12 | Otsu256 IoU | Relative residual | White residual | Original Phi | Joint |',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for row in rows:
        lines.append(f"|{row['method']}|{row['observation_id']}|{row['amplitude_normalized_nrmse']:.9f}|{row['otsu_iou']:.9f}|{row['relative_residual']:.12f}|{row['whitened_residual']:.6f}|{row['physical_objective']:.12g}|{'PASS' if row['joint_gate'] else 'FAIL'}|")
    lines+=['','| Endpoint | Delta NRMSE from original onehop | Delta IoU | Delta relative residual |',
        '|---|---:|---:|---:|']
    for row in rows:lines.append(f"|{row['method']}|{row['delta_nrmse_from_original_onehop']:.9f}|{row['delta_iou_from_original_onehop']:.9f}|{row['delta_relative_residual_from_original_onehop']:.12f}|")
    lines+=['','| Endpoint | Physical LU | Physical complex RHS | Receipt wall s | Receipt CPU s |',
        '|---|---:|---:|---:|---:|']
    for row in rows:lines.append(f"|{row['method']}|{row['physical_model_cost']['LU']}|{row['physical_model_cost']['actual_complex_RHS']}|{row['wall_seconds']:.6f}|{row['cpu_seconds']:.6f}|")
    lines+=['',
        'Each original run top-level physical total is charged once; no nested stage cost is added again. Local process timers and Windows adapter timers have different scopes, retained in JSON; these rows support no cross-platform speed claim. Windows import charge/caller wall are separately recorded. All correctors retain their1000-iteration cap stop. Scoring adds0 LU/0 RHS.','',
        'The new partition-average version does not replace frozen onehop, any of the original18 predictions, or the secondary choices. The42cb input was selected by public data failure and smallest opaque ID. Future claims about this revision need independent validation; this disclosed subset supplies no new blind rate. Further4000-step budget endpoints are not included here. Final scientific interpretation belongs to parent Codex.']
    with outputs[2].open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],rows=[{key:row[key]for key in ('method','observation_id','amplitude_normalized_nrmse','otsu_iou','relative_residual','whitened_residual','joint_gate','delta_nrmse_from_original_onehop')}for row in rows],old13_preserved=result['old13_table_preserved'],old18_preserved=result['original18_table_preserved'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['seal','score'])
    parser.add_argument('--seal',type=Path,default=ROOT/'audit/POST_CONFIRMATION_ADDENDA_V1_SEAL.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/DEVELOPMENT_POST_CONFIRMATION_ADDENDA_V1')
    parser.add_argument('--authorize-exposed-truth',action='store_true')
    args=parser.parse_args()
    if args.operation=='seal':seal(args)
    else:score(args)
