#!/usr/bin/env python3
"""Select one ordinary smoothing/shrink action on exposed validation128 only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r2.core import OUT,Book,write,save
from p3r2.classical_risk import apply_action


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation',type=Path,default=OUT/'A/anchor_data/validation.npz')
    parser.add_argument('--output',type=Path,default=OUT/'D/classical')
    parser.add_argument('--job',default='validation_classical_smoothing')
    args=parser.parse_args()
    if 'locked' in str(args.validation).lower() or 'validation' not in str(args.validation).lower():
        raise ValueError('only declared exposed validation labels allowed')
    if (args.output/'FROZEN.json').exists():raise FileExistsError('classical action already selected; no implicit reselection')
    book=Book('D',args.job,max_wall=600,analysis_only=True);status='COMPLETE';error=None
    try:
        with np.load(args.validation,allow_pickle=False) as data:
            anchor=data['chi0'];truth=data['chi_GT'];parents=data['parent'] if 'parent' in data else data['parents']
        if anchor.shape!=(128,576) or truth.shape!=anchor.shape or len(set(parents.astype(str)))!=128:
            raise ValueError('exact128 independent validation parents and material arrays required')
        if not np.isfinite(anchor).all() or not np.isfinite(truth).all():raise ValueError('finite validation arrays required')
        cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        if float(cfg['field_mean'])!=0:raise ValueError('exact symmetric-prior mean contract required')
        points=np.asarray(json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'])
        roi=np.flatnonzero(np.linalg.norm(points,axis=1)<=cfg['roi_radius'])
        candidates=[]
        for sigma in (0.,.35,.7,1.,1.5):
            for blend in (.25,.5,1.):
                for shrink in (0.,.05,.1,.2,.4):
                    prediction=apply_action(anchor,cfg,roi,sigma,blend,shrink)
                    risk=float(np.mean(np.abs(prediction-truth)**2))
                    candidates.append(dict(sigma=sigma,blend=blend,prior_shrink=shrink,validation_fullgrid_complex_MSE=risk))
        selected=min(candidates,key=lambda item:(item['validation_fullgrid_complex_MSE'],item['sigma'],item['prior_shrink'],item['blend']))
        prediction=apply_action(anchor,cfg,roi,selected['sigma'],selected['blend'],selected['prior_shrink'])
        write(args.output/'FROZEN.json',dict(selected=selected,candidates=candidates,
            baseline_validation_fullgrid_complex_MSE=float(np.mean(np.abs(anchor-truth)**2)),
            calibration_parents=parents.astype(str),validation_count=128,
            selection_scope='exposed validation calibration only; not independent performance evidence',
            action='parentwise reflected Gaussian blur -> convex blend with physics point -> shrink to exact prior material mean; outsideROI exact host',
            neural=False,density_or_posterior_claim=False,freeze_before_new_locked_observations=True))
        save(args.output/'VALIDATION_PREDICTIONS.npz',parent=parents,chi=prediction,
             parent_risk=np.mean(np.abs(prediction-truth)**2,axis=1))
        print(json.dumps(selected),flush=True)
    except BaseException as exception:status='ERROR';error=repr(exception);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
