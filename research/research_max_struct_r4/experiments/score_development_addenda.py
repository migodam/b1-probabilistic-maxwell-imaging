"""Write-once separate development addenda. Never touches the old13 score table."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import json
import numpy as np
from offline_score_common import (read,write_new,stamp,unchanged,public_metrics,
    prediction,physical,metrics,write_csv_new,frozen_contract)

ROOT=Path(__file__).resolve().parents[1]
CASE='obs_b20d5a55677d4eca'


def utc():return datetime.now(timezone.utc).isoformat()


def specs(target):
    if target=='phase_pi':
        directory=ROOT/'results/onehop_phase_pi'
        return dict(shape=(24,24),classification='POST_SEAL_PUBLIC_DEVELOPMENT_NOT_BLIND',
            case_id=CASE,public=ROOT.parent/'research_max_struct_r2/data_public'/(CASE+'.npz'),
            truth=ROOT.parent/'research_max_struct_r2/sealed_truth'/(CASE+'.npz'),
            entries=[('onehop_phase_pi',directory/'RECEIPT.json',directory/'prediction.npz')])
    directory=ROOT/'grid32_diagnostic'
    return dict(shape=(32,32),classification='DISCLOSED_DEVELOPMENT_GRID_DIAGNOSTIC_NOT_BLIND',
        case_id='grid32_curved_bar',public=directory/'data_public/grid32_curved_bar.npz',
        truth=directory/'development_truth/grid32_curved_bar_truth.npz',
        entries=[(method,directory/'results'/method/'RECEIPT.json',directory/'results'/method/'prediction_and_trajectory.npz') for method in ('onehop','global')])


def seal(args):
    spec=specs(args.target);public=public_metrics(spec['public']);entries=[]
    for method,receipt_path,prediction_path in spec['entries']:
        rec=read(receipt_path)
        if rec['status'] not in ('SEALED_DEVELOPMENT_PREDICTION','SEALED_DISCLOSED_DEVELOPMENT_PREDICTION'):
            raise RuntimeError('Development prediction is not closed: '+method)
        m,yp=prediction(prediction_path,spec['shape'])
        entries.append(dict(method=method,receipt=stamp(receipt_path),prediction=stamp(prediction_path),
            status=rec['status'],original_Phi_pretruth=physical(m,yp,public,spec['shape'])['physical_objective']))
    result=dict(status='DEVELOPMENT_ADDENDUM_SEALED',target=args.target,
        sealed_utc=utc(),classification=spec['classification'],case_id=spec['case_id'],
        entries=entries,public=stamp(spec['public']),truth_opened=False,
        old13_table=stamp(ROOT/'results/DEVELOPMENT_SCORES.json'),old13_preserved=True,
        new_physical_LU=0,new_complex_RHS=0,hash_checks=False)
    write_new(args.seal,result)
    print(json.dumps(dict(status=result['status'],target=args.target,count=len(entries),truth_opened=False)))


def score(args):
    if not args.authorize_development_truth:
        raise PermissionError('Parent must explicitly authorize this disclosed development addendum')
    outputs=[args.output.with_suffix('.'+ext) for ext in ('json','csv','md')]
    if any(p.exists() for p in outputs):raise FileExistsError('Write-once addendum; supply a new versioned prefix')
    closed=read(args.seal);spec=specs(args.target)
    if closed['status']!='DEVELOPMENT_ADDENDUM_SEALED' or closed['target']!=args.target:
        raise RuntimeError('Matching pre-truth development seal required')
    evidence=[closed['public'],closed['old13_table']]+[item for row in closed['entries'] for item in (row['receipt'],row['prediction'])]
    if not all(unchanged(item) for item in evidence):raise RuntimeError('Evidence changed after seal')
    # Only this target's authorized disclosed-development truth; no confirmation path.
    opened=utc()
    with np.load(spec['truth'],allow_pickle=False) as data:
        truth=data['m'].copy();support=data['mask'].ravel().astype(bool)
    public=public_metrics(spec['public']);rows=[]
    for entry in closed['entries']:
        rec=read(entry['receipt']['path'])
        m,yp=prediction(entry['prediction']['path'],spec['shape'])
        row=dict(case_id=spec['case_id'],method=entry['method'],classification=spec['classification'],
            **metrics(m,yp,public,truth,support,spec['shape']))
        row.update(physical_objective_receipt_delta=row['physical_objective']-rec['physical_objective'],
            relative_residual_receipt_delta=row['relative_residual']-rec['physical_relative_residual'],
            physical_model_cost=rec['physical_model_cost'],wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],
            caller_wall_seconds=rec.get('caller_wall_seconds'),caller_cpu_seconds=rec.get('caller_cpu_seconds'),
            corrector_stop=rec.get('corrector_stop',rec.get('corrector',{}).get('stop')),
            failed_attempts=rec.get('failed_attempts'),receipt=entry['receipt']['path'],prediction=entry['prediction']['path'])
        rows.append(row)
    result=dict(status='DEVELOPMENT_ADDENDUM_SCORED',target=args.target,
        classification=spec['classification'],seal=str(args.seal),first_development_truth_opened_utc=opened,
        old13_table_preserved=unchanged(closed['old13_table']),contract=frozen_contract(),rows=rows,
        shared_generation_receipt=str(ROOT/'grid32_diagnostic/audit/generation/GENERATION_RECEIPT.json') if args.target=='grid32' else None,
        scope='Disclosed development geometry/phase control; not fresh blind sample or PDE convergence study',
        new_physical_LU=0,new_complex_RHS=0,hash_checks=False,confirmation_GT_access=False)
    write_new(outputs[0],result);write_csv_new(outputs[1],rows)
    lines=[f'# Development addendum: {args.target}','',spec['classification'],
        '', 'Separate versioned addendum; the previous13-row table remains unchanged. Original frozen RMSE/12, Otsu256, public99% CN radius, prior and material box. No new physics or GT selection; no confirmation truth access.','',
        '| Method | RMSE/12 | IoU | Relative residual | White residual | Joint | LU | Physical RHS | Wall s |',
        '|---|---:|---:|---:|---:|---|---:|---:|---:|']
    for r in rows:
        lines.append(f"|{r['method']}|{r['amplitude_normalized_nrmse']:.9f}|{r['otsu_iou']:.9f}|{r['relative_residual']:.12f}|{r['whitened_residual']:.6f}|{'PASS' if r['joint_gate'] else 'FAIL'}|{r['physical_model_cost'].get('LU')}|{r['physical_model_cost'].get('actual_complex_RHS')}|{r['wall_seconds']:.6f}|")
    lines+=['', 'All displayed physical costs are original run top-level totals; do not add nested stage costs again. Solver wall/caller wall scopes remain explicit in JSON. Grid32 shared operator generation is a separate one-time receipt, not included in each solve. Corrector stop states are retained in JSON; no convergence claim. Scientific judgment stays with parent Codex.']
    with outputs[2].open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],target=args.target,rows=rows,old13_preserved=result['old13_table_preserved'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['seal','score'])
    parser.add_argument('--target',required=True,choices=['phase_pi','grid32'])
    parser.add_argument('--seal',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--authorize-development-truth',action='store_true')
    args=parser.parse_args()
    args.seal=args.seal or ROOT/'audit'/('DEVELOPMENT_'+args.target.upper()+'_SEALED_V1.json')
    args.output=args.output or ROOT/'results'/('DEVELOPMENT_'+args.target.upper()+'_SCORES_V1')
    if args.operation=='seal':seal(args)
    else:score(args)
