"""Separate post-hoc-development capacity score. No main blind result changes."""
from pathlib import Path
import json,time
import numpy as np
from score_blind import ROOT,utc,write_new,otsu,iou,prior

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    p=next((ROOT/'results/capacity_control').glob('*/RECEIPT.json'));r=json.loads(p.read_text())
    if r['status']!='SEALED_POST_HOC_PREDICTION':raise RuntimeError('Post-hoc prediction not closed')
    prediction=ROOT/r['prediction'];ident=r['case_id']
    with np.load(prediction,allow_pickle=False)as d:
        m=d['m'].copy();yp=d['Y_prediction'].copy()
    write_new(ROOT/'audit/POSTHOC_CAPACITY_SEALED.json',dict(status='POST_HOC_PREDICTION_SEALED',sealed_utc=utc(),
        original_truth_already_released=True,classification='POST_HOC_DEVELOPMENT_NOT_FRESH_BLIND',
        runner_truth_access=r['truth_access'],prediction=r['prediction'],prediction_bytes=prediction.stat().st_size,
        hashes=None,hash_reason='No hash checks per user',new_forward_LU=0,new_forward_RHS=0))
    with np.load(ROOT/f'sealed_truth/{ident}.npz',allow_pickle=False)as d:truth=d['m'];support=d['mask'].ravel().astype(bool)
    with np.load(ROOT/f'data_public/{ident}.npz',allow_pickle=False)as d:Y=d['Y'];sigma=float(d['sigma']);radius=float(d['noise_radius'])
    mask,threshold=otsu(m);N=float(np.linalg.norm(m-truth)/np.sqrt(576)/12);white=float(np.linalg.norm(yp-Y)/sigma)
    old=next(x for x in json.loads((ROOT/'results/BLIND_SCORES.json').read_text())['rows']if x['observation_id']==ident and x['arm']=='native')
    value=dict(status='POST_HOC_DEVELOPMENT_SCORED',classification='POST_HOC_DEVELOPMENT_AFTER_GT_RELEASE',
        case_id=ident,family='curved_bar',prediction=r['prediction'],source='experiments/score_capacity.py',
        amplitude_normalized_nrmse=N,otsu_iou=iou(mask,support),otsu_threshold=threshold,
        fixed_midpoint_iou=iou(m>6.1,support),relative_residual=float(np.linalg.norm(yp-Y)/np.linalg.norm(Y)),
        whitened_residual=white,noise_radius99=radius,data_gate=bool(white<=radius),material_gate=bool(N<=.15),
        joint_gate=bool(N<=.15 and iou(mask,support)>=.7 and white<=radius),
        physical_objective_receipt_delta=.5*float(np.linalg.norm(yp-Y)/np.linalg.norm(Y))**2+prior(m)-r['physical_objective'],
        baseline_native=dict(nrmse=old['amplitude_normalized_nrmse'],iou=old['otsu_iou'],relative_residual=old['relative_residual']),
        independent_claim='Residual reduction does not rescue material recovery; no change to original3-case blind denominator',
        new_forward_LU=0,new_forward_RHS=0,wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
        timer_scope='Saved post-hoc prediction,truth/public reads and metrics; own JSON write excluded')
    write_new(ROOT/'results/POSTHOC_CAPACITY_SCORE.json',value)
    print(json.dumps(value))

if __name__=='__main__':main()
