"""Score completed development predictions against already-public R2 truth.
No forward model imported. Clean-target and original-noisy gates stay separate.
"""
from pathlib import Path
from datetime import datetime,timezone
import json,csv,time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];R2=ROOT.parent/'research_max_struct_r2'
CASE='obs_b20d5a55677d4eca'
ARMS=['noisy_strong_release','clean_native','passive_path','extra_original','globaldual_restart','scaled_original_restart',
      'onecell_roots','neighbor_swaps','all_pair_swaps','paired_neighbor_rays','dressed_neighbor_rays']

def utc():return datetime.now(timezone.utc).isoformat()
def otsu(m):
    if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
    hist,edges=np.histogram(m,bins=256);c=(edges[:-1]+edges[1:])/2
    w=np.cumsum(hist);q=np.cumsum(hist*c);den=w*(w[-1]-w);z=np.full(256,-np.inf);good=den>0
    z[good]=(q[-1]*w[good]-q[good]*w[-1])**2/den[good]
    threshold=float(c[np.argmax(z)]);return m>threshold,threshold
def iou(a,b):
    u=np.count_nonzero(a|b);return float(np.count_nonzero(a&b)/u)if u else 1.
def prior(m):
    m=m.reshape(24,24);d=np.r_[np.diff(m,axis=0).ravel(),np.diff(m,axis=1).ravel()]/19.5
    return float(1e-4*np.mean(d*d)+1e-4*np.mean(np.sqrt(d*d+.01**2)-.01))

def main():
    wall,cpu=time.perf_counter(),time.process_time();records=[]
    for arm in ARMS:
        receipt=ROOT/f'results/development/{arm}/RECEIPT.json';rec=json.loads(receipt.read_text())
        if rec['status']!='SEALED_DEVELOPMENT_PREDICTION':raise RuntimeError('Unclosed '+arm)
        pred=ROOT/rec['prediction']if'prediction'in rec else receipt.parent/'prediction.npz'
        with np.load(pred,allow_pickle=False)as d:
            if d['m'].shape!=(576,)or not all(np.isfinite(d[k]).all()for k in d.files):raise ValueError('Invalid '+arm)
        records.append((arm,receipt,rec,pred))
    seal=dict(status='ALL_DEVELOPMENT_PREDICTIONS_CLOSED',count=len(records),checked_utc=utc(),classification='PUBLISHED_R2_DEVELOPMENT_NOT_BLIND',
        truth_already_public=True,entries=[dict(arm=a,receipt=str(p.relative_to(ROOT)),prediction=str(pred.relative_to(ROOT)),
        finished_utc=r['finished_utc'],status=r['status'],prediction_bytes=pred.stat().st_size)for a,p,r,pred in records],
        no_hash_checks=True,new_forward_LU=0,new_forward_RHS=0)
    (ROOT/'audit/DEVELOPMENT_SEAL.json').write_text(json.dumps(seal,indent=2))
    with np.load(R2/f'sealed_truth/{CASE}.npz',allow_pickle=False)as d:truth=d['m'].copy();support=d['mask'].ravel().astype(bool);clean=d['Y_clean'].copy()
    with np.load(R2/f'data_public/{CASE}.npz',allow_pickle=False)as d:Y=d['Y'].copy();sigma=float(d['sigma']);radius=float(d['noise_radius'])
    rows=[]
    # Published baseline is a reference only: no additional run/charge.
    old=json.loads((R2/'results/BLIND_SCORES.json').read_text())
    ref=next(r for r in old['rows']if r['observation_id']==CASE and r['arm']=='native')
    for arm,p,rec,pred in records:
        with np.load(pred,allow_pickle=False)as d:m=d['m'].copy();yp=d['Y_prediction'].copy()
        mask,threshold=otsu(m);N=float(np.linalg.norm(m-truth)/np.sqrt(576)/12);I=iou(mask,support)
        noisy_rel=float(np.linalg.norm(yp-Y)/np.linalg.norm(Y));clean_rel=float(np.linalg.norm(yp-clean)/np.linalg.norm(clean))
        white=float(np.linalg.norm(yp-Y)/sigma);data=white<=radius;legal=bool(m.min()>=.02-1e-12 and m.max()<=20+1e-12)
        clean_arm=arm=='clean_native';target=clean if clean_arm else Y
        target_relative=clean_rel if clean_arm else noisy_rel
        phi=.5*target_relative**2+prior(m)
        row=dict(arm=arm,case_id=CASE,family='curved_bar',classification='POST_SEAL_PUBLISHED_R2_DEVELOPMENT',
          normalized_material_error=N,otsu_iou=I,otsu_threshold=threshold,fixed_midpoint_iou=iou(m>6.1,support),
          material_gate=bool(N<=.15),image_gate=bool(N<=.15 and I>=.7),material_legal=legal,
          target='PUBLISHED_CLEAN_Y_COUNTERFACTUAL'if clean_arm else'ORIGINAL_NOISY_Y',target_relative_residual=target_relative,
          target_noise_gate=None if clean_arm else bool(data),target_noise_gate_null_reason='Clean target has no CN noise; no invented clean discrepancy threshold'if clean_arm else None,
          original_noisy_relative_residual=noisy_rel,original_noisy_whitened_residual=white,original_noise_radius99=radius,
          original_noisy_data_gate=bool(data),original_noisy_joint_gate=bool(data and N<=.15 and I>=.7 and legal),
          residual_to_cleanY=clean_rel,physical_objective_to_run_target=phi,
          physical_objective_receipt_delta=phi-float(rec['physical_objective']),
          material_min=float(m.min()),material_max=float(m.max()),prediction=str(pred.relative_to(ROOT)),receipt=str(p.relative_to(ROOT)),
          wall_seconds=rec['wall_seconds'],cpu_seconds=rec['cpu_seconds'],
          projected_gradient_inf=rec.get('projected_gradient_inf'),note='Clean data change is an information counterfactual, not same-data algorithm improvement'if clean_arm else'')
        rows.append(row)
    out=dict(status='DEVELOPMENT_SCORED',scored_utc=utc(),scope=f'{len(rows)} completed controls on already-public R2 case; no fresh-blind success claims',
      inherited_metrics=dict(NRMSE_denominator=12,material_gate=.15,OtsuIoU=.7,noisy_gate='Inherited R2 99% CN whitened radius; no clean-target noise gate'),
      R2_native_reference=dict(nrmse=ref['amplitude_normalized_nrmse'],iou=ref['otsu_iou'],relative_residual=ref['relative_residual'],
        reused=True,new_cost=0,new_cost_reason='Immutable published result reference, not an R3 run'),rows=rows,new_forward_LU=0,new_forward_RHS=0,
      wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,timer_scope='Receipt/prediction/public-truth read and metrics; own JSON/CSV output excluded',
      no_hash_checks=True)
    (ROOT/'results/DEVELOPMENT_SCORES.json').write_text(json.dumps(out,indent=2,allow_nan=False))
    with(ROOT/'results/DEVELOPMENT_SCORES.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    lines=['# R3 published-case development scores','',
      f'All {len(rows)} predictions were closed before this scoring pass. The material truth is already public R2 development data. No new physical solve. N=RMSE/12; IoU is reconstruction-only256-bin Otsu. Clean-target discrepancy has no Gaussian gate.','',
      '|Control|N|IoU|Target relative residual|Original noisy white residual|Original-noisy joint|',
      '|---|---:|---:|---:|---:|---|']
    for r in rows:lines.append(f"|{r['arm']}|{r['normalized_material_error']:.6f}|{r['otsu_iou']:.6f}|{r['target_relative_residual']:.8f}|{r['original_noisy_whitened_residual']:.3f}|{'PASS'if r['original_noisy_joint_gate']else'FAIL'}|")
    lines+=['',f"Published native reference: N={ref['amplitude_normalized_nrmse']:.6f}, IoU={ref['otsu_iou']:.6f}, original-noisy residual={ref['relative_residual']:.8f}.",
      '', 'Lower data residual alone is not a material improvement. These development controls cannot modify the original R2 blind table or establish a new R3 blind rate.']
    (ROOT/'results/DEVELOPMENT_SCORES.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps([{k:r[k]for k in ['arm','normalized_material_error','otsu_iou','target_relative_residual','original_noisy_joint_gate','physical_objective_receipt_delta']}for r in rows]))

if __name__=='__main__':main()
