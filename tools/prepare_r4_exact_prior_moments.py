#!/usr/bin/env python3
"""Exact original P1 unconditional prior moment action; no data or physics."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import json
from pathlib import Path
import resource
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from p3r4.core import Book,OUT
from p3r4.probability_eval import decoder_from_geometry
from p3r4.field_moments import FORMAT,load_field_moments
from p3r4.evaluation import write_new_json
from p3l.moments import gaussian_saturation_moments


def main():
    torch.set_num_threads(1)
    destination=OUT/'models/risk/prior_exact'
    if destination.exists():raise FileExistsError('exact prior artifact already exists; no implicit overwrite/replay')
    job='exact_p1_prior_moments_'+str(time.time_ns())
    audit=OUT/'AUDIT'/job;audit.mkdir(parents=True,exist_ok=False)
    write_new_json(audit/'FROZEN_CONTRACT.json',dict(status='UNCONDITIONAL_PRIOR_ONLY',
        source='unchanged p3l.moments.gaussian_saturation_moments and original p1_material.MaterialPrior.A',
        new_data=0,training_updates=0,physical_actions=0,GT_read=False,
        independence='u1,u2 mutually independent and independent of Gaussian field under ORIGINAL PRIOR ONLY; never assume posterior factorization',
        likelihood_used=False,baseline='fixed prior mean action; no observation dependence',max_wall_seconds=120))
    book=Book('I',job,analysis_only=True,max_wall=120)
    wall,cpu=time.monotonic(),time.process_time();status='COMPLETE';error=None
    try:
        cfg=json.loads((OUT/'configs/TARGET.json').read_text())
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
        geometry=decoder_from_geometry(cfg,points);A=np.asarray(geometry['A']);roi=np.asarray(geometry['roi'])
        variance_f=np.sum(A*A,axis=1);mean_f=np.full(576,cfg['field_mean']);sd_f=np.sqrt(variance_f)
        saturation=gaussian_saturation_moments(mean_f,sd_f,halfwidth=cfg['phase_transition_halfwidth'])
        mask=np.zeros(576);mask[roi]=1
        eh=saturation.mean*mask;eh2=saturation.second_moment*mask
        a,b=cfg['inclusion_increment_real_interval'];c,d=cfg['inclusion_increment_imag_interval']
        er,ei=(a+b)/2,(c+d)/2
        er2,ei2=(a*a+a*b+b*b)/3,(c*c+c*d+d*d)/3
        eri=er*ei # valid only because original prior u1,u2 are independent.
        host=complex(*cfg['host_chi']);delta_mean=complex(er,ei)
        mean=host+delta_mean*eh
        variance=(er2+ei2)*eh2-abs(delta_mean*eh)**2
        second=np.stack((host.real**2+2*host.real*er*eh+er2*eh2,
                         host.imag**2+2*host.imag*ei*eh+ei2*eh2,
                         host.real*host.imag+(host.real*ei+host.imag*er)*eh+eri*eh2))
        components=np.stack((mean.real,mean.imag,second[0],second[1],second[2]))[:,roi]
        contract=dict(format=FORMAT,kind='prior_exact',architecture='train_mean',second_moments=True,
            target_id=cfg['target_id'],prior_id=cfg['prior_id'],train_size=0,validation_size=0,
            information_set='fixed EXACT original P1 unconditional prior action; ignores observations and reads no GT or joint data',
            formula='Echi=host+Edelta*Eh; Var_complex=E|delta|²*E[h²]-|Edelta*Eh|²',
            independence_scope='ONLY ORIGINAL UNCONDITIONAL PRIOR; invalid posterior factorization is not used',
            moment_source='unchanged gaussian_saturation_moments(field_mean,sqrt(sum A_i²),halfwidth=original_tau)',
            numerical_exactness='analytic prior identity evaluated with existing double-precision recurrence/scaled quadrature; not an interval certificate',
            spatial_scope='fullfield mean and pixel-marginal joint Re/Im moments only; no full spatial covariance or joint sampler',
            likelihood_used=False,field_points=576,ROI_points=len(roi),parameter_count=0,
            physics_calls=0,optimizer_updates=0,selected_update=0,GT_read=False)
        warning=dict(prior_only='known exact unconditional prior moments; not posterior moments or posterior uncertainty given y',
            factorization='Edelta*Eh is justified only by original prior independence; posterior properties and phase must remain jointly modeled',
            scope='no pixel Gaussian density, spatial joint posterior, material sampler, or credible interval qualification',
            numerical='reuses old proven Gaussian saturation moments; floating quadrature/roundoff remains explicitly numerical')
        destination.mkdir(parents=True,exist_ok=False)
        payload=dict(format=FORMAT,kind='prior_exact',architecture='train_mean',seed=None,second_moments=True,
            moment_source='EXACT_ORIGINAL_P1_UNCONDITIONAL_PRIOR',geometry=geometry,
            y_mean=np.zeros(512),y_scale=np.ones(512),target_mean=components,target_scale=np.ones((5,1)),
            state_dict=None,selected_update=0,contract=contract,warnings=warning,
            input_statistics='unused fixed placeholders: this constant prior action never normalizes or conditions on y')
        torch.save(payload,destination/'BEST.pt')
        np.savez_compressed(destination/'PRIOR_MOMENTS.npz',material_mean=mean,complex_variance=variance,
            raw_second_moments=second,fraction_mean=eh,fraction_second_moment=eh2,
            field_mean=mean_f,field_variance=variance_f,ROI=roi,
            E_delta=np.array(delta_mean),E_abs_delta_squared=np.array(er2+ei2))
        predictor=load_field_moments(destination/'BEST.pt','cpu')
        one=predictor.predict_moments(np.zeros((1,512)))
        two=predictor.predict_moments(np.ones((1,512)))
        assert np.array_equal(one['material_mean'],two['material_mean'])
        mean_error=float(np.max(abs(one['material_mean'][0]-mean)))
        variance_error=float(np.max(abs(one['raw_complex_variance'][0]-variance)))
        assert mean_error<1e-13 and variance_error<1e-13
        symmetry=float(np.max(abs(saturation.mean-.5))) if cfg['field_mean']==0 else None
        if symmetry is not None:assert symmetry<1e-12
        assert np.all(variance>=-1e-14) and np.all(mean[mask==0]==host) and np.all(variance[mask==0]==0)
        write_new_json(destination/'PRIOR_CONTRACT.json',contract)
        result=dict(status='EXACT_ORIGINAL_PRIOR_BASELINE_READY',artifact=str(destination),checkpoint=str(destination/'BEST.pt'),
            checkpoint_format=FORMAT,loader='field',recommended_method=dict(id='prior_exact',loader='field',family='prior_exact',seed=0,checkpoint=str(destination/'BEST.pt')),
            geometry_target_id=cfg['target_id'],prior_id=cfg['prior_id'],
            material_mean_ROI_real_min=float(mean[roi].real.min()),material_mean_ROI_real_max=float(mean[roi].real.max()),
            material_mean_ROI_imag_min=float(mean[roi].imag.min()),material_mean_ROI_imag_max=float(mean[roi].imag.max()),
            mean_fullgrid_prior_risk=float(variance.mean()),
            fraction_symmetry_error=symmetry,field_variance_min=float(variance_f.min()),field_variance_max=float(variance_f.max()),
            mean_checkpoint_replay_error=mean_error,variance_checkpoint_replay_error=variance_error,
            y_independence_bitwise=True,variance_status=one['second_moment_status'],
            analytic_status='unconditional original prior identity; numerical implementation, no interval certificate',
            posterior_qualification='NOT_CLAIMED',physics_calls=0,optimizer_updates=0,new_data=0,GT_read=False,
            wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            book_receipt=str(book.jobdir/'CLOSED.json'))
        write_new_json(audit/'RESULT.json',result)
        print(json.dumps(result,indent=2),flush=True)
    except BaseException as exc:status='FAILED';error=repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()
