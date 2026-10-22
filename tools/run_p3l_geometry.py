"""Paid P1 projector diagnostic recipe. Importing this module performs no work.

The parent must explicitly start main(); this file does not schedule or launch
itself. It never reads GT, chains or teachers and never changes the target.
"""
from pathlib import Path
import json
import os
import sys
import traceback

for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_name,'1')

import numpy as np
import scipy
from scipy.linalg import eigh

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save,plain
from p3l.conditional import make_frame

PARENTS=('DEV_000','DEV_028')
STEPS=(.01,.005,.0025)
DIMENSION=8
JOB='geometry_projector_DEV000_DEV028'
RELATIVE_GAP_FLAG=1e-6


def _ratio(numerator,denominator):
    return float(numerator/denominator) if denominator > 0 else None


def _spectrum(q,J,D):
    """Same posterior-material-risk span as make_frame; only its projector varies."""
    frame=make_frame(q,J,DIMENSION,kind='posterior_material_risk',material_derivative=D)
    V=frame['Vw'];P=V@V.T
    Hdata=J[:,:64].T@J[:,:64];Hdata=(Hdata+Hdata.T)/2
    ev,U=eigh(Hdata)
    Croot=(U*(1/np.sqrt(1+np.maximum(ev,0))))@U.T
    M=np.real(D[:,:64].conj().T@((1/D.shape[0])*D[:,:64]));M=(M+M.T)/2
    risk=Croot@M@Croot;risk=(risk+risk.T)/2
    lam,W=eigh(risk);descending=lam[::-1]
    selected=Croot@W[:,-DIMENSION:][:,::-1]
    selected_sv=np.linalg.svd(selected,compute_uv=False)
    gap=float(descending[DIMENSION-1]-descending[DIMENSION])
    boundary_scale=max(abs(float(descending[DIMENSION-1])),abs(float(descending[DIMENSION])),np.finfo(float).tiny)
    risk_scale=max(float(np.max(abs(lam))),np.finfo(float).tiny)
    relative_gap=gap/boundary_scale
    roundoff_scale=64*np.finfo(float).eps*risk_scale
    close=bool(relative_gap<=RELATIVE_GAP_FLAG or gap<=roundoff_scale)
    singular_J=np.linalg.svd(J[:,:64],compute_uv=False)
    meta=dict(selection='posterior_material_risk; top8 of Croot M Croot mapped by Croot, Euclidean QR',
        weak_dimension=DIMENSION,gaussian_dimension=64,material_weight=1/D.shape[0],
        selection_boundary_gap=gap,gap_relative_to_boundary=float(relative_gap),
        seed_eigenvalue_gap_first_second=float(descending[0]-descending[1]),
        seed_eigenvalue_gap_second_third=float(descending[1]-descending[2]),
        gap_relative_to_risk_operator=float(gap/risk_scale),
        arithmetic_gap_diagnostic_scale=float(roundoff_scale),
        selection_gap_near_degenerate=close,
        selection_gap_status='SMALL_GAP_DERIVATIVE_INTERPRETATION_UNRESOLVED' if close else 'GAP_DIAGNOSTIC_RESOLVED_NOT_A_CERTIFICATE',
        clipped_negative_Hdata_eigenvalues=int(np.sum(ev<0)),minimum_raw_Hdata_eigenvalue=float(ev[0]),
        Hdata_negative_relative_to_spectral_scale=_ratio(max(0.,-float(ev[0])),max(float(abs(ev).max()),np.finfo(float).tiny)),
        prior_plus_GN_condition=float((1+max(0.,float(ev[-1])))/(1+max(0.,float(ev[0])))),
        selected_mapped_span_condition=_ratio(float(selected_sv[0]),float(selected_sv[-1])),
        observation_J_condition=_ratio(float(singular_J[0]),float(singular_J[-1])),
        projector_idempotence_fro=float(np.linalg.norm(P@P-P)),
        projector_symmetry_fro=float(np.linalg.norm(P-P.T)))
    arrays=dict(projector=P,Vw=V,Hdata_eigenvalues=ev,risk_eigenvalues_descending=descending,
                risk_matrix=risk,prior_GN_covariance_root=Croot,material_Gram=M,
                mapped_selected_span=selected,mapped_selected_span_singular_values=selected_sv,
                observation_J_singular_values=singular_J)
    return P,V,meta,arrays


def _rotation(P,V,Pbase,Vbase):
    delta=P-Pbase
    sv=np.linalg.svd(Vbase.T@V,compute_uv=False)
    angle=np.arccos(np.clip(sv,0,1))
    return dict(projector_rotation_operator=float(np.linalg.norm(delta,ord=2)),
                projector_rotation_frobenius=float(np.linalg.norm(delta)),
                largest_principal_angle_radians=float(np.max(angle)))


def _frozen_J(T,Cchi,q):
    """Fixed Maxwell material-response factor, exact changing P1 decoder: 0 RHS."""
    D=T.prior.derivative(q)
    complex_J=(Cchi@D).reshape(32,8,66)
    J=np.column_stack([T.white(complex_J[:,:,k]) for k in range(66)])
    return J,D


def _stratum(field,halfwidth):
    return np.where(field<=-halfwidth,-1,np.where(field>=halfwidth,1,0))


def _directional_result(projectors,h,seeds,Pbase):
    dpi=(projectors[(0,1)]-projectors[(0,-1)])/(2*h)
    dpj=(projectors[(1,1)]-projectors[(1,-1)])/(2*h)
    bracket=dpi@seeds[:,1]-dpj@seeds[:,0]
    normal=(np.eye(64)-Pbase)@bracket
    tangent=Pbase@bracket
    return dict(DP_vi=dpi,DP_vj=dpj,bracket=bracket,normal_bracket=normal,
                tangent_bracket=tangent)


def _norms(result):
    return dict(DP_vi_frobenius=float(np.linalg.norm(result['DP_vi'])),
                DP_vj_frobenius=float(np.linalg.norm(result['DP_vj'])),
                bracket_norm=float(np.linalg.norm(result['bracket'])),
                normal_bracket_norm=float(np.linalg.norm(result['normal_bracket'])),
                tangent_bracket_norm=float(np.linalg.norm(result['tangent_bracket'])),
                normal_to_bracket_ratio=_ratio(np.linalg.norm(result['normal_bracket']),np.linalg.norm(result['bracket'])))


def _convergence(results):
    pairs=[]
    for coarse,fine in zip(results[:-1],results[1:]):
        pairs.append(dict(coarse_h=coarse[0],fine_h=fine[0],
            DP_vi_difference_frobenius=float(np.linalg.norm(coarse[1]['DP_vi']-fine[1]['DP_vi'])),
            DP_vj_difference_frobenius=float(np.linalg.norm(coarse[1]['DP_vj']-fine[1]['DP_vj'])),
            bracket_difference_norm=float(np.linalg.norm(coarse[1]['bracket']-fine[1]['bracket'])),
            normal_bracket_difference_norm=float(np.linalg.norm(coarse[1]['normal_bracket']-fine[1]['normal_bracket']))))
    ratios={}
    for key in ['DP_vi_difference_frobenius','DP_vj_difference_frobenius',
                'bracket_difference_norm','normal_bracket_difference_norm']:
        ratio=_ratio(pairs[0][key],pairs[1][key])
        ratios[key]=dict(coarse_difference_over_fine_difference=ratio,
                        observed_order_log2=float(np.log2(ratio)) if ratio is not None and ratio>0 else None)
    return dict(successive_differences=pairs,descriptive_orders=ratios,
        status='THREE_STEP_DIFFERENCE_DIAGNOSTIC_NOT_AN_ERROR_CERTIFICATE',
        caution='Second-order behavior requires additional smoothness; C2 decoder seams and small selection gaps can reduce convergence. No extrapolated integrability verdict.')


def _one_parent(T,parent,jobdir):
    y=T.observation(parent)
    if y.shape!=(32,8) or not np.isfinite(y).all():
        raise ValueError('frozen online observation shape required')
    anchor_path=OUT/f'A/anchors/{parent}.npz'
    base=T.load_cache(anchor_path,owner='P3L A paid anchor; geometry diagnostic reuses only y-derived inference state')
    dest=jobdir/parent;dest.mkdir(parents=True,exist_ok=True)
    # Save raw paid derivative before frame analysis or the next physical action.
    J0=T.jacobian(base);D0=T.prior.derivative(base.q)
    save(dest/'BASE_RAW.npz',q=base.q,J=J0,material_derivative=D0,
         prediction=base.state.data,chi=base.state.chi,receiver_transfer=base.T,
         field=base.state.field,current=base.state.current,endpoint_key=np.array(base.key))
    P0,V0,meta0,arrays0=_spectrum(base.q,J0,D0)
    Cchi=np.einsum('qp,ps->qsp',base.T,base.state.field).reshape(256,576)
    Jf0,_=_frozen_J(T,Cchi,base.q)
    Pf0,Vf0,metaf0,arraysf0=_spectrum(base.q,Jf0,D0)
    seeds=V0[:,:2].copy()
    save(dest/'BASE_GEOMETRY.npz',**arrays0,direction_seeds=seeds,
         frozen_Maxwell_material_factor=Cchi,frozen_J=Jf0,
         frozen_projector=Pf0,frozen_risk_eigenvalues=arraysf0['risk_eigenvalues_descending'],
         frozen_Hdata_eigenvalues=arraysf0['Hdata_eigenvalues'])
    write(dest/'BASE.json',dict(parent=parent,anchor_path=anchor_path,anchor_owner='P3L A paid anchor',
         anchor_endpoint=base.key,properties_fixed=base.q[64:],Phi=T.phi(base,y),
         actual_full_J_geometry=meta0,frozen_geometry=metaf0,
         frozen_vs_full_J_relative=_ratio(np.linalg.norm(Jf0-J0),np.linalg.norm(J0)),
         frozen_vs_full_projector_operator=float(np.linalg.norm(Pf0-P0,ord=2)),
         pair_definition='first two base Vw columns only; basis signs do not define endpoint projectors',
         exact_fiber_claim=False,posterior_reference_claim=False))
    _,base_aux=T.prior.decode(base.q)
    roi=T.prior.mask;halfwidth=T.prior.temperature
    base_strata=_stratum(base_aux['field'],halfwidth)
    rows=[];actual_results=[];frozen_results=[]
    for h in STEPS:
        projectors={};frozen_projectors={};endpoint_rows=[]
        near_gap=meta0['selection_gap_near_degenerate']
        for direction in range(2):
            for sign in (1,-1):
                q=base.q.copy();q[:64]+=sign*h*seeds[:,direction]
                if not np.array_equal(q[64:],base.q[64:]):raise RuntimeError('properties moved')
                label=f'h{h:g}_v{direction}_{"plus" if sign==1 else "minus"}'
                counters_before=T.book.counts.copy()
                p=T.evaluate(q,purpose=f'geometry_{parent}_{label}')
                J=T.jacobian(p);D=T.prior.derivative(q)
                save(dest/f'{label}_RAW.npz',q=q,J=J,material_derivative=D,
                     prediction=p.state.data,chi=p.state.chi,endpoint_key=np.array(p.key))
                states_delta=int(T.book.counts.get('full_forward_calls',0)-counters_before.get('full_forward_calls',0))
                rhs_delta=int(T.book.counts.get('linear_solve_rhs',0)-counters_before.get('linear_solve_rhs',0))
                if states_delta!=1 or rhs_delta!=40:
                    write(dest/f'{label}_COST_DEVIATION.json',dict(new_states=states_delta,new_solver_RHS=rhs_delta,
                          expected_states=1,expected_RHS=40,continuation_stopped=True))
                    raise RuntimeError('geometry endpoint cost differs from frozen 1 state/40 RHS recipe')
                P,V,meta,arrays=_spectrum(q,J,D)
                Jf,_=_frozen_J(T,Cchi,q)
                Pf,Vf,metaf,arraysf=_spectrum(q,Jf,D)
                save(dest/f'{label}_GEOMETRY.npz',**arrays,
                     frozen_J=Jf,frozen_projector=Pf,frozen_Vw=Vf,
                     frozen_risk_eigenvalues_descending=arraysf['risk_eigenvalues_descending'],
                     frozen_Hdata_eigenvalues=arraysf['Hdata_eigenvalues'])
                _,aux=T.prior.decode(q)
                stratum=_stratum(aux['field'],halfwidth)
                row=dict(parent=parent,h=h,direction_index=direction,sign=sign,endpoint=p.key,
                    raw_path=dest/f'{label}_RAW.npz',geometry_path=dest/f'{label}_GEOMETRY.npz',
                    Gaussian_step_norm=float(np.linalg.norm(q[:64]-base.q[:64])),properties_fixed=True,
                    material_RMS=float(np.linalg.norm(p.state.chi-base.state.chi)/24),
                    full_prediction_distance_complex=float(np.linalg.norm(p.state.data-base.state.data)),
                    full_prediction_distance_whitened=float(np.linalg.norm(T.white(p.state.data-base.state.data))),
                    Phi=T.phi(p,y),backward_residual=float(p.state.residual),
                    ROI_saturation_strata_changed_from_base=int(np.sum((stratum!=base_strata)&roi)),
                    minimum_ROI_field_distance_to_saturation_seam=float(np.min(abs(abs(aux['field'][roi])-halfwidth))),
                    actual_geometry=meta,frozen_geometry=metaf,
                    actual_rotation=_rotation(P,V,P0,V0),frozen_rotation=_rotation(Pf,Vf,Pf0,Vf0),
                    frozen_vs_full_J_relative=_ratio(np.linalg.norm(Jf-J),np.linalg.norm(J)),
                    frozen_vs_full_projector_operator=float(np.linalg.norm(Pf-P,ord=2)),
                    new_solver_RHS=int(T.book.counts.get('linear_solve_rhs',0)-counters_before.get('linear_solve_rhs',0)),
                    new_states=int(T.book.counts.get('full_forward_calls',0)-counters_before.get('full_forward_calls',0)),
                    frozen_comparison_extra_RHS=0)
                write(dest/f'{label}.json',row)
                endpoint_rows.append(row);projectors[(direction,sign)]=P;frozen_projectors[(direction,sign)]=Pf
                near_gap=near_gap or meta['selection_gap_near_degenerate'] or metaf['selection_gap_near_degenerate']
        actual=_directional_result(projectors,h,seeds,P0)
        # Same base normal projection and same directions isolate the changed
        # Maxwell response from changed decoder, rather than changing the test.
        frozen=_directional_result(frozen_projectors,h,seeds,P0)
        delta={key:actual[key]-frozen[key] for key in actual}
        save(dest/f'h{h:g}_DIRECTIONAL_PROJECTOR.npz',
             **{'actual_'+k:v for k,v in actual.items()},
             **{'frozen_'+k:v for k,v in frozen.items()},
             **{'actual_minus_frozen_'+k:v for k,v in delta.items()})
        row=dict(parent=parent,h=h,actual=_norms(actual),frozen_Maxwell_changing_decoder=_norms(frozen),
                 actual_minus_frozen=_norms(delta),normal_projection='I-P_actual_base for both comparisons',
                 minimum_boundary_gap_actual=min([meta0['selection_boundary_gap']]+[r['actual_geometry']['selection_boundary_gap'] for r in endpoint_rows]),
                 any_small_selection_gap=bool(near_gap),endpoints=endpoint_rows,
                 interpretation='FINITE_DIFFERENCE_LOCAL_PROJECTOR_DIAGNOSTIC_ONLY',
                 integrability_verdict='NOT_ASSIGNED',exact_fiber_claim=False,posterior_reference_claim=False)
        rows.append(row);actual_results.append((h,actual));frozen_results.append((h,frozen))
        write(dest/f'h{h:g}_SUMMARY.json',row)
        write(dest/'PARTIAL_RESULTS.json',dict(parent=parent,rows=rows,status='IN_PROGRESS'))
        print(json.dumps(plain(dict(parent=parent,h=h,normal_bracket_norm=row['actual']['normal_bracket_norm'],
                                   small_selection_gap=near_gap))),flush=True)
    result=dict(parent=parent,status='DIAGNOSTIC_COMPLETE_NOT_AN_INTEGRABILITY_CERTIFICATE',rows=rows,
                actual_step_convergence=_convergence(actual_results),frozen_step_convergence=_convergence(frozen_results),
                caution='One pair at one anchor cannot establish integrability; nonzero normal leakage requires gap, truncation and floating-point/convergence review. No exact-fiber or posterior-reference claim.')
    write(dest/'RESULTS.json',result)
    return result


def main():
    # Explicit paid entry point: the parent starts this job, never import time.
    book=Book('C',JOB,max_wall=3600)
    status='COMPLETE';error=None
    try:
        protocol=dict(job=JOB,stage='C',parents=PARENTS,steps=STEPS,weak_dimension=DIMENSION,
            projector_kind='posterior_material_risk',direction_pairs='first two base Euclidean Vw columns',
            properties='fixed native anchor u, no probit transformation',
            formula='DP[vi]=(P(q+hvi)-P(q-hvi))/(2h); B=DP[vi]vj-DP[vj]vi; N=(I-P0)B',
            full_forward_states_planned=24,incremental_forward_and_J_RHS_planned=960,
            additional_base_J_RHS_maximum_if_no_cached_transfers=64,max_wall_seconds=3600,
            relative_selection_gap_flag=RELATIVE_GAP_FLAG,GT_used=False,chains_used=False,teachers_used=False,
            target_changed=False,flow_loop_executed=False,
            comparison='same endpoint changing decoder with one base frozen Maxwell material factor; zero additional RHS',
            expected_rows=6,status='PAID_JOB_STARTED_BY_PARENT')
        write(book.jobdir/'PROTOCOL.json',protocol)
        sources=['tools/run_p3l_geometry.py','src/p3l/conditional.py','src/p3l/core.py','src/p1_material/prior.py']
        for name in sources:
            destination=book.jobdir/'source'/name;destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes((ROOT/name).read_bytes())
        write(book.jobdir/'SOURCE_PROVENANCE.json',dict(sources=sources,byte_copies_saved=True,
             anchors=[str(OUT/f'A/anchors/{parent}.npz') for parent in PARENTS],
             permitted_measurements=[str(OUT/f'inputs/online/{parent}.npz') for parent in PARENTS],
             no_hash_check=True,no_GT_or_reference_samples=True,python=sys.version,
             numpy=np.__version__,scipy=scipy.__version__,
             thread_environment={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}))
        T=Target(book,keep_fields=True)
        write(book.jobdir/'TARGET_CONFIG_SNAPSHOT.json',T.cfg)
        results=[]
        for parent in PARENTS:
            results.append(_one_parent(T,parent,book.jobdir))
            write(book.jobdir/'RESULTS.json',dict(protocol=protocol,parents=results,
                  status='DIAGNOSTIC_COMPLETE' if len(results)==len(PARENTS) else 'IN_PROGRESS',
                  solver_RHS_so_far=int(book.counts.get('linear_solve_rhs',0)),
                  fresh_states_so_far=int(book.counts.get('full_forward_calls',0)),
                  no_integrability_or_posterior_claim=True))
        if int(book.counts.get('full_forward_calls',0))!=24:
            raise RuntimeError('geometry total fresh-state count differs from frozen recipe')
    except BaseException as exc:
        status='ERROR';error=repr(exc)
        write(book.jobdir/'ERROR.json',dict(status=status,error=error,partial_outputs_preserved=True))
        traceback.print_exc();raise
    finally:
        book.close(status,error)


if __name__=='__main__':
    main()
