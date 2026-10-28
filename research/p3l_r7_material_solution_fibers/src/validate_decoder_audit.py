"""Independent, read-only E0/E1 validation after source relocation.

Original archives are never rewritten. No geometry ranking helper is reused:
SciPy gesvd checks spectra, and the original P1 decoder checks finite paths.
"""
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import sys,json,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import svd,svdvals
from scipy.special import ndtr
from scipy.optimize import brentq
import decoder_geometry as dg

P=Path(__file__).resolve().parents[1]

def inverse(prior,p):
    return np.array([prior.temperature*brentq(lambda u:float(dg.saturation(u))-v,-1,1,xtol=1e-14) for v in p])

def local_interval(prior,x,chi):
    q=x.copy();q[64:]=ndtr(q[64:]);_,a=prior.decode(q)
    T=prior.mask&(abs(a['field'])<prior.temperature);H=prior.mask&(a['field']<-prior.temperature)
    pure=prior.mask&(a['field']>=prior.temperature)
    A=prior.A[T];s=svdvals(A);tol=np.finfo(float).eps*max(A.shape)*s[0]
    rank=int(np.sum(s>tol));fullrow=rank==A.shape[0]
    no_seams=not np.any(prior.mask&(abs(a['field'])==prior.temperature))
    d=a['inclusion_increment'];fraction=a['fraction'][T]
    low=max(float(prior.real_interval[0]/d.real),float(prior.imag_interval[0]/d.imag),float(fraction.max()),0.)
    high=min(float(prior.real_interval[1]/d.real),float(prior.imag_interval[1]/d.imag))
    rec=dict(no_pure_cells=not pure.any(),no_seams=no_seams,transition_cells=int(T.sum()),A_transition_machine_rank=rank,A_transition_full_row_rank=fullrow,
       A_transition_singular_min=float(s[-1]),A_transition_singular_max=float(s[0]),
       strict_host_margin_min=float(np.min(-prior.temperature-a['field'][H])),
       analytic_property_and_fraction_alpha_open_interval=[low,high],
       theorem_applicable=bool(not pure.any() and no_seams and fullrow and low<1<high))
    if not rec['theorem_applicable']:return rec
    R=np.linalg.lstsq(A,np.eye(A.shape[0]),rcond=None)[0]
    rec['right_inverse_identity_residual']=float(np.linalg.norm(A@R-np.eye(A.shape[0])))
    margins=-prior.temperature-a['field'][H]
    AR=prior.A[H]@R
    eta=min(.25,.5*(1-low),.5*(high-1))
    for _ in range(64):
        lo,hi=1-eta,1+eta
        flo=inverse(prior,fraction/hi);fhi=inverse(prior,fraction/lo)
        slo=dg.saturation_prime(flo/prior.temperature);shi=dg.saturation_prime(fhi/prior.temperature)
        # S' is unimodal on (-1,1); its minimum on each inverse interval
        # occurs at an endpoint. This bounds |d f_T(alpha)/d alpha|.
        L=prior.temperature*fraction/(lo**2*np.minimum(slo,shi))
        bound=eta*(np.abs(AR)@L)
        ratio=float(np.max(bound/margins)) if H.any() else 0.
        if ratio<.5:break
        eta/=2
    rec.update(float64_sufficient_local_alpha_closed_interval=[1-eta,1+eta],
               interval_host_margin_consumption_bound_ratio=ratio,
               rigorous_interval_arithmetic_certificate=False,
               interval_qualification='analytical sufficient-bound formula evaluated in float64; theorem is symbolic, numerical endpoints are diagnostics')
    return rec

def main():
    start=time.perf_counter();prior,cfg,geom=dg.load_prior()
    with np.load(P/'raw/decoder/FULL66_DECODER_ARRAYS.npz',allow_pickle=False) as f:a={k:f[k] for k in f.files}
    old=json.loads((P/'audit/E1_RANK_STRATA.json').read_text())
    frozen=json.loads((P/'configs/E1_PRIOR_MANIFEST.json').read_text())
    assert np.array_equal(a['x'][:256],np.array(frozen['x']))
    ranks=[];spectral_error=[];derror=[];material_error=[];inactive=[]
    directions=json.loads((P/'raw/decoder/DIRECTIONS_AND_PROJECTORS.json').read_text())['records']
    finite=json.loads((P/'raw/decoder/FINITE_STEP_COMPARISONS.json').read_text())['records']
    for i,x in enumerate(a['x']):
        q=x.copy();q[64:]=ndtr(x[64:]);chi,aux=prior.decode(q)
        chain=np.r_[np.ones(64),np.exp(-x[64:]**2/2)/np.sqrt(2*np.pi)]
        D=prior.derivative(q)*chain
        C=np.vstack([D.real,D.imag])/24
        _,s,_=svd(C,full_matrices=False,lapack_driver='gesvd')
        rr=[int(np.sum(s>eps*s[0])) if s[0]>0 else 0 for eps in dg.RELS]
        ranks.append(rr);spectral_error.append(float(np.max(abs(s-a['singular_values'][i]))))
        derror.append(float(np.max(abs(D-a['jacobian'][i]))));material_error.append(float(np.max(abs(chi-a['material'][i]))))
        if 'saturation_inactive_exact_local_subspace' in directions[i]['direction_labels']:
            idx=directions[i]['direction_labels'].index('saturation_inactive_exact_local_subspace')
            v=np.array(directions[i]['vectors'][idx]);T=prior.mask&(abs(aux['field'])<prior.temperature)
            row=dict(parent_id=str(a['parent_id'][i]),A_T_v_residual=float(np.linalg.norm(prior.A[T]@v[:64])),property_step_zero=bool(np.all(v[64:]==0)))
            safe=old['records'][i]['inactive_first_direction_symmetric_plateau_interval']
            tested=[v for v in finite if v['parent_id']==str(a['parent_id'][i]) and v['direction']=='saturation_inactive_exact_local_subspace' and (safe is None or v['t']<safe)]
            row.update(symmetric_plateau_interval=safe,tested_steps_within_interval=[v['t'] for v in tested],max_material_RMS_inside_interval=max((v['material_RMS'] for v in tested),default=0.))
            inactive.append(row)
    exact_rank_matches=[rr==old['records'][i]['rank_by_relative'] for i,rr in enumerate(ranks)]
    with np.load(P/'raw/decoder/NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.npz',allow_pickle=False) as f:b={k:f[k] for k in f.files}
    curve=[];intervals=[]
    for parent in ['R7_E1_092','R7_E1_201','R7_E1_233','R7_E1_251']:
        idx=int(parent.rsplit('_',1)[1]);interval=local_interval(prior,a['x'][idx],a['material'][idx]);interval['parent_id']=parent;intervals.append(interval)
        inds=np.flatnonzero(b['parent_id']==parent);assert len(inds)==4
        for j in inds:
            xp=b['x'][j];q=xp.copy();q[64:]=ndtr(xp[64:]);chi,aux=prior.decode(q)
            curve.append(dict(parent_id=parent,t=float(b['t'][j]),native_properties=q[64:].tolist(),legal=prior.legal(q),properties_strictly_interior=bool(np.all((q[64:]>0)&(q[64:]<1))),
              saved_material_redecode_max_absolute_error=float(np.max(abs(chi-b['material'][j]))),material_RMS_to_parent=prior.material_rms(chi,a['material'][idx]),
              interpretation='numerical same-material witness; not a global exact floating-point identity'))
    assert all(exact_rank_matches) and max(derror)==0 and max(material_error)==0
    assert all(v['legal'] and v['properties_strictly_interior'] for v in curve)
    assert max(v['material_RMS_to_parent'] for v in curve)<1e-14
    assert max(v['max_material_RMS_inside_interval'] for v in inactive)<1e-14
    result=dict(status='PASS_INDEPENDENT_VALIDATION',new_Maxwell_RHS=0,new_forward=0,new_LU=0,NN_training=0,
       relocated_repo_root=str(dg.ROOT),relocated_output_root=str(dg.OUT),prior_source=str(sys.modules[prior.__class__.__module__].__file__),
       rank_comparison_algorithm='SciPy gesvd vs original NumPy gesdd',all_rank_cutoffs_match=all(exact_rank_matches),
       max_singular_absolute_discrepancy=max(spectral_error),max_saved_derivative_discrepancy=max(derror),max_saved_material_discrepancy=max(material_error),
       rank_mismatch_parent_ids=[str(a['parent_id'][i]) for i,v in enumerate(exact_rank_matches) if not v],
       frozen_coordinates_exactly_unchanged=True,inactive_direction_validations=inactive,radial_curve_validations=curve,
       local_IFT_interval_diagnostics=intervals,all_four_curves_four_steps_legal=all(v['legal'] and v['properties_strictly_interior'] for v in curve),
       maximum_corrected_material_RMS=max(v['material_RMS_to_parent'] for v in curve),
       original_scientific_outputs_overwritten=False,wall_seconds=time.perf_counter()-start)
    dg.write(P/'audit/E1_INDEPENDENT_VALIDATION.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ['inactive_direction_validations','radial_curve_validations','local_IFT_interval_diagnostics']},indent=2))

if __name__=='__main__':main()
