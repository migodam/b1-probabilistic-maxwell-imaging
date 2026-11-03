"""Post-seal score only: published R2 truth, unchanged metrics, no physics.

Parent explicitly authorized this independent scorer after the four R4 outputs
were closed. Online modules do not import it. No GT-ranked candidate selection.
"""
from pathlib import Path
from datetime import datetime, timezone
import csv
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
R2 = ROOT.parent/'research_max_struct_r2'
R3 = ROOT.parent/'research_max_struct_r3'
sys.path.insert(0,str(REPO))
# These exact frozen functions import no Maxwell operator or inverse runner.
from research.research_max_struct_r2.experiments.score_blind import otsu, iou, prior

CASE='obs_b20d5a55677d4eca'
SPECIFICATIONS=(
    ('cooperative_clipped','cooperative_paths','cooperative_prediction.npz','cooperative'),
    ('weak_gn_clipped','cooperative_paths','weak_gn_prediction.npz','weak_gn'),
    ('box_cooperative','box_cooperative_paths','prediction.npz',None),
    ('tmatrix_completion','tmatrix_control','prediction.npz',None))


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_new(path,contents):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x')as stream:
        json.dump(contents,stream,indent=2,allow_nan=False)


def score():
    started,cpu=time.perf_counter(),time.process_time()
    destinations=[ROOT/'results'/f'DEVELOPMENT_SCORES.{ext}'for ext in ('json','csv','md')]
    if any(path.exists()for path in destinations):
        raise FileExistsError('Preserve existing post-seal development score; no metric retuning')
    entries=[];loaded=[]
    # All closure and prediction checks precede any published-truth array access.
    for arm,directory,filename,subarm in SPECIFICATIONS:
        receipt_path=ROOT/'results'/directory/'RECEIPT.json'
        prediction_path=receipt_path.parent/filename
        rec=json.loads(receipt_path.read_text())
        if rec['status']not in ('SEALED_DEVELOPMENT_PREDICTION','SEALED_DEVELOPMENT_PREDICTIONS'):
            raise RuntimeError('Unclosed R4 output: '+arm)
        with np.load(prediction_path,allow_pickle=False)as data:
            m=data['m'].copy();yp=data['Y_prediction'].copy()
        if m.shape!=(576,)or yp.shape!=(48,12)or not np.isfinite(m).all()or not np.isfinite(yp).all():
            raise ValueError('Invalid sealed prediction: '+arm)
        entries.append(dict(arm=arm,receipt=str(receipt_path.relative_to(ROOT)),
                       prediction=str(prediction_path.relative_to(ROOT)),status=rec['status'],
                       prediction_bytes=prediction_path.stat().st_size,prediction_mtime_ns=prediction_path.stat().st_mtime_ns,
                       receipt_bytes=receipt_path.stat().st_size,receipt_mtime_ns=receipt_path.stat().st_mtime_ns))
        detail=rec['arms'][subarm]if subarm else rec
        loaded.append((arm,rec,detail,m,yp,entries[-1]))
    seal=dict(status='ALL_FOUR_R4_PREDICTIONS_CLOSED_BEFORE_SCORING',checked_utc=utc(),
              classification='POST_SEAL_PUBLIC_DEVELOPMENT',case_id=CASE,entries=entries,
              authorization='Parent explicitly authorized already-published R2 truth for independent scoring only',
              truth_opened_by_this_scorer=False,truth_already_public=True,hash_checks=False,
              new_physical_LU=0,new_complex_RHS=0)
    write_new(ROOT/'audit/DEVELOPMENT_PREDICTIONS_SEALED.json',seal)
    truth_opened_utc=utc()
    with np.load(R2/'sealed_truth'/f'{CASE}.npz',allow_pickle=False)as data:
        truth=data['m'].copy();support=data['mask'].ravel().astype(bool)
    with np.load(R2/'data_public'/f'{CASE}.npz',allow_pickle=False)as data:
        Y=data['Y'].copy();sigma=float(data['sigma']);radius=float(data['noise_radius'])
    published_r2=json.loads((R2/'results/BLIND_SCORES.json').read_text())
    native=next(row for row in published_r2['rows']if row['observation_id']==CASE and row['arm']=='native')
    published_r3=json.loads((R3/'results/DEVELOPMENT_SCORES.json').read_text())
    anchor=next(row for row in published_r3['rows']if row['arm']=='scaled_original_restart')
    references=dict(
        R2_native=dict(nrmse=native['amplitude_normalized_nrmse'],otsu_iou=native['otsu_iou'],
                       relative_residual=native['relative_residual'],source='research/research_max_struct_r2/results/BLIND_SCORES.json'),
        R3_scaled_original_restart=dict(nrmse=anchor['normalized_material_error'],otsu_iou=anchor['otsu_iou'],
                       relative_residual=anchor['original_noisy_relative_residual'],source='research/research_max_struct_r3/results/DEVELOPMENT_SCORES.json'))
    rows=[]
    for arm,rec,detail,m,yp,entry in loaded:
        mask,threshold=otsu(m)
        nrmse=float(np.linalg.norm(m-truth)/np.sqrt(truth.size)/12.)
        overlap=iou(mask,support)
        residual=float(np.linalg.norm(yp-Y));relative=residual/float(np.linalg.norm(Y));white=residual/sigma
        phi=.5*relative**2+prior(m)
        violation=float(max(0.,.02-m.min(),m.max()-20.));legal=violation<=1e-12
        material=nrmse<=.15;image=material and overlap>=.7;noise=white<=radius
        expected_phi=detail.get('physical_objective',detail.get('final_Phi'))
        expected_residual=detail.get('physical_relative_residual',detail.get('physical_residual'))
        base=references['R3_scaled_original_restart']
        rows.append(dict(arm=arm,case_id=CASE,classification='POST_SEAL_PUBLIC_DEVELOPMENT',
                    amplitude_normalized_nrmse=nrmse,RMSE=12*nrmse,otsu_iou=overlap,otsu_threshold=threshold,
                    relative_residual=relative,whitened_residual=white,noise_radius99=radius,
                    physical_objective=phi,material_gate=bool(material),image_gate=bool(image),noise_gate=bool(noise),
                    material_and_data_gate=bool(material and noise and legal),joint_gate=bool(image and noise and legal),
                    material_legal=bool(legal),material_constraint_violation=violation,
                    material_min=float(m.min()),material_max=float(m.max()),
                    physical_objective_receipt_delta=phi-expected_phi,
                    relative_residual_receipt_delta=relative-expected_residual,
                    delta_nrmse_from_R3_scaled=nrmse-base['nrmse'],delta_iou_from_R3_scaled=overlap-base['otsu_iou'],
                    delta_residual_from_R3_scaled=relative-base['relative_residual'],
                    material_nrmse_improved_from_R3_scaled=bool(nrmse<base['nrmse']),
                    iou_improved_from_R3_scaled=bool(overlap>base['otsu_iou']),
                    residual_improved_from_R3_scaled=bool(relative<base['relative_residual']),
                    prediction=entry['prediction'],receipt=entry['receipt'],
                    projected_gradient_inf=rec.get('projected_gradient_inf'),corrector_stop=detail['corrector_stop']))
    tmatrix=next(row for row in rows if row['arm']=='tmatrix_completion')
    assessment=dict(
        comparison='sealed T-matrix+corrector endpoint versus published R3 scaled-original endpoint',
        residual_decreased=tmatrix['residual_improved_from_R3_scaled'],
        nrmse_decreased=tmatrix['material_nrmse_improved_from_R3_scaled'],
        iou_increased=tmatrix['iou_improved_from_R3_scaled'],
        joint_gate=tmatrix['joint_gate'],
        scope='One disclosed development object and finite budgets; no fresh success rate, no GT selection, no algorithm change')
    result=dict(status='DEVELOPMENT_SCORED',classification='POST_SEAL_PUBLIC_DEVELOPMENT',
                scored_utc=utc(),first_truth_opened_utc=truth_opened_utc,seal='audit/DEVELOPMENT_PREDICTIONS_SEALED.json',
                inherited_contract=dict(source='research/research_max_struct_r2/experiments/score_blind.py',
                    NRMSE_denominator=12,NRMSE_gate=.15,Otsu_bins=256,Otsu_IoU_gate=.7,
                    data_gate='unchanged public R2 99% CN whitened radius',sigma=sigma,noise_radius99=radius,
                    material_box=[.02,20.],legality_tolerance=1e-12),
                references=references,rows=rows,tmatrix_material_assessment=assessment,
                new_physical_LU=0,new_complex_RHS=0,hash_checks=False,
                online_algorithms_changed=False,candidate_reselection=False,
                wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
                timer_scope='closure/read checks,published truth reads,metrics; result JSON/CSV/MD serialization excluded')
    write_new(destinations[0],result)
    with destinations[1].open('x',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    lines=['# R4 independent development scores','',
           '**POST_SEAL_PUBLIC_DEVELOPMENT**: all4 outputs were closed before this scoring pass. Parent authorized already-published R2 truth for scoring only. No new physics, GT-ranked selection, or algorithm change.',
           '',f'Unchanged R2 gates: RMSE/12≤.15; prediction-only256-bin Otsu IoU≥.7; original noisy whitened residual≤{radius:.9f}; public material box.',
           '', '| Endpoint | NRMSE/12 | Otsu IoU | Relative residual | White residual | Joint |',
           '|---|---:|---:|---:|---:|---|']
    for name,ref in references.items():
        lines.append(f"|{name} (published reference)|{ref['nrmse']:.9f}|{ref['otsu_iou']:.9f}|{ref['relative_residual']:.9f}|in inherited report|FAIL|")
    for row in rows:
        lines.append(f"|{row['arm']}|{row['amplitude_normalized_nrmse']:.9f}|{row['otsu_iou']:.9f}|{row['relative_residual']:.9f}|{row['whitened_residual']:.3f}|{'PASS'if row['joint_gate']else'FAIL'}|")
    base=references['R3_scaled_original_restart']
    lines+=['',f"T-matrix endpoint versus R3 scaled-original: residual {base['relative_residual']:.9f} → {tmatrix['relative_residual']:.9f}; NRMSE {base['nrmse']:.9f} → {tmatrix['amplitude_normalized_nrmse']:.9f}; IoU {base['otsu_iou']:.9f} → {tmatrix['otsu_iou']:.9f}.",
            f"Residual decreased: {assessment['residual_decreased']}; material NRMSE decreased: {assessment['nrmse_decreased']}; IoU increased: {assessment['iou_increased']}. Original joint gate: {'PASS'if assessment['joint_gate']else'FAIL'}.",
            '', 'These metrics describe one disclosed development endpoint. They do not establish a new-blind recovery, superiority over converged CIE, or a converged optimum. Original R2 blind results remain unchanged. Scoring adds0 LU and0 RHS.']
    with destinations[2].open('x')as stream:stream.write('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],classification=result['classification'],rows=[{k:r[k]for k in
               ('arm','amplitude_normalized_nrmse','otsu_iou','relative_residual','joint_gate')}for r in rows],assessment=assessment)),flush=True)


if __name__=='__main__':
    score()
