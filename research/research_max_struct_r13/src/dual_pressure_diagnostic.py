"""Exposed-case directional dual-pressure correction, raw receiver metric W=I.

Implements parent THEORY_DERIVATIONS section4.1. This measures a specific
resolvent-amplified AP defect; it neither estimates an operator norm nor
certifies future DR, material recovery, or theory validity.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
from pathlib import Path
from collections import Counter
import argparse,json,time
import numpy as np
from scipy.linalg import lu_solve
from scipy.stats import chi2
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r13.runtime.field_current import DenseFieldCurrentProjectors
from research.research_max_struct_r13.src.material_events import norm_pair,finite_threshold


def _ratio(a,b):
    return float(a/b) if b>np.finfo(float).tiny else None


def _delta(after,before):
    return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)}


def _json_default(value):
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,np.ndarray):return value.tolist()
    raise TypeError(type(value).__name__)


def operator_cost_summary(physical_cost,projector_cost,rows,n,ns):
    """Count dense actions omitted by inherited model/projector setup ledgers.

    FrozenModel.forward explicitly evaluates G@J and S@J once per FULL_FORWARD;
    projector construction forms B^H B once and T@H^-1T^H once per unique port.
    Matrix scaling/conjugate-transpose materialization is not an application.
    """
    extra=Counter()
    for row in rows:extra.update(row['explicit_diagnostic_operator_cost'])
    forwards=int(physical_cost.get('FULL_FORWARD',0))
    ports=int(projector_cost['port_eigh'])
    return dict(
        Green_forward_calls=projector_cost['G_action_calls']+forwards+extra['Green_forward_calls'],
        Green_forward_complex_RHS=projector_cost['G_action_rhs']+forwards*ns+extra['Green_forward_complex_RHS'],
        Green_adjoint_calls=projector_cost['G_adjoint_action_calls']+1+extra['Green_adjoint_calls'],
        Green_adjoint_complex_RHS=projector_cost['G_adjoint_action_rhs']+n+extra['Green_adjoint_complex_RHS'],
        receiver_forward_calls=projector_cost['receiver_action_calls']+ports+forwards+extra['receiver_forward_calls'],
        receiver_forward_complex_RHS=projector_cost['receiver_action_rhs']+projector_cost['H_setup_rhs']+forwards*ns+extra['receiver_forward_complex_RHS'],
        receiver_adjoint_calls=projector_cost['receiver_adjoint_action_calls']+extra['receiver_adjoint_calls'],
        receiver_adjoint_complex_RHS=projector_cost['receiver_adjoint_action_rhs']+extra['receiver_adjoint_complex_RHS'],
        scope='projector iterations + actual constructor B^H B/port Gram + FrozenModel.forward G/S actions + explicit diagnostic actions; materialization/scaling excluded',
        disjoint_from_physical_LU_RHS='These are dense operator applications, not extra LU solves; diagnostic adjoint solves already appear in physical_model_cost')


def diagnose_pair(model,P,aX,aQ,route):
    started,wall=time.process_time(),time.perf_counter()
    physics_before=dict(model.cost);projector_before=dict(P.costs);cost=Counter()
    cX,cQ,m=P.material(aX,aQ)
    dX,dQ=aX-cX,aQ-cQ
    primeX,primeQ=P.affine_ball(cX,cQ)
    # Capture the actual just-completed PA(c) metadata before another PA call.
    info=dict(P.last_A_info)
    if 'multiplier' not in info or 'pseudoinverse_limit' not in info:
        raise RuntimeError('Projector did not expose its actual PA ball multiplier convention')
    multiplier=float(info['multiplier'])
    eX,eQ=primeX-aX,primeQ-aQ
    residual=P.scale*(model.S@primeX)-np.column_stack(P.targets)
    cost['receiver_forward_calls']+=1;cost['receiver_forward_complex_RHS']+=model.ns
    receive_rhs=model.S.conj().T@residual
    cost['receiver_adjoint_calls']+=1;cost['receiver_adjoint_complex_RHS']+=model.ns
    gh=model.G.conj().T@np.concatenate([eQ,dQ+eQ],axis=1)
    cost['Green_adjoint_calls']+=1;cost['Green_adjoint_complex_RHS']+=2*model.ns
    correction_rhs=eX+P.scale*gh[:,:model.ns]
    # One physical factorization at m_c, retained by MaxwellModel. Both
    # adjoint blocks use this same LU in one multi-RHS triangular solve.
    state=model.state(m)
    with model._physics.operation('DUAL_PRESSURE_DIAGNOSTIC_ADJOINT',rhs=2*model.ns,
                                 adjoint_RHS=2*model.ns):
        solved=lu_solve(state['lu'],np.concatenate([receive_rhs,correction_rhs],axis=1),
                        trans=2,check_finite=False)
    p,correction=solved[:,:model.ns],-solved[:,model.ns:]
    cost['batched_retained_LU_adjoint_calls']+=1
    cost['batched_retained_LU_adjoint_complex_RHS']+=2*model.ns
    compatible_eQ=P.scale*(model.G@eX)
    cost['Green_forward_calls']+=1;cost['Green_forward_complex_RHS']+=model.ns
    gap=norm_pair(dX,dQ);pair_norm=norm_pair(eX,eQ)
    correction_norm=float(np.linalg.norm(correction));dx_norm=float(np.linalg.norm(dX))
    rhs_norm=float(np.linalg.norm(correction_rhs))
    penergy=np.sum(abs(p)**2,axis=1);cenergy=np.sum(abs(correction)**2,axis=1)
    threshold=finite_threshold(cX,cQ,dX,dQ,m,P.scale,P.lower,P.upper,P.components)
    first=int(threshold['first_cell'])
    finite_multiplier=bool(np.isfinite(multiplier))
    slater_margin=float(P.delta-P.minimum_reachable_residual)
    slater=bool(P.delta>0 and slater_margin>0)
    valid_multiplier=bool(finite_multiplier and not info['pseudoinverse_limit'] and slater)
    receiver_pressure=-multiplier*P.scale*p if finite_multiplier else None
    rnorm=float(np.linalg.norm(residual))
    kkt=None
    if finite_multiplier:
        stationarity=dX+eX+P.scale*gh[:,model.ns:]+multiplier*P.scale*receive_rhs
        kkt=dict(current_stationarity_norm=float(np.linalg.norm(stationarity)),
            current_stationarity_relative=_ratio(float(np.linalg.norm(stationarity)),
                norm_pair(dX+eX,P.scale*gh[:,model.ns:])+float(np.linalg.norm(multiplier*P.scale*receive_rhs))),
            complementarity_absolute=abs(multiplier*(rnorm*rnorm-P.delta*P.delta)),
            reconstruction_discrepancy_relative=_ratio(float(np.linalg.norm(dX-receiver_pressure-correction)),dx_norm),
            diagnostic_scope='qualification residuals only; no machine-precision identity-test headline')
    pressure_rank=1+int(np.count_nonzero(penergy>penergy[first]))
    first_interior=bool(threshold['interior'][first])
    pressure_only_t=None
    if valid_multiplier and multiplier>0 and first_interior:
        pressure_only_t=_ratio(float(np.linalg.norm(cQ[first])),
            multiplier*P.scale*float(np.linalg.norm(p[first])))
    row=dict(route=route,material_source='PC of the saved public-data auxiliary A state',
        relative_pair_defect=_ratio(pair_norm,gap),pair_defect_norm=pair_norm,normal_gap=gap,
        dX_norm=dx_norm,directional_correction_norm=correction_norm,
        directional_correction_relative_to_dX=_ratio(correction_norm,dx_norm),
        correction_amplification_per_pair_defect=_ratio(correction_norm,pair_norm),
        pre_resolvent_rhs_amplification_per_pair_defect=_ratio(rhs_norm,pair_norm),
        actual_resolvent_gain_on_this_rhs=_ratio(correction_norm,rhs_norm),
        pressure_direction_cosine=_ratio(float(np.vdot(dX,receiver_pressure).real),dx_norm*float(np.linalg.norm(receiver_pressure))) if finite_multiplier else None,
        geometric_first_event_t=float(threshold['first_t']) if np.isfinite(threshold['first_t']) else None,
        first_event_cell=first,first_event_cell_material=float(m[first]),first_event_cell_interior=first_interior,
        first_event_cell_p_energy=float(penergy[first]),first_event_cell_p_energy_descending_rank=pressure_rank,
        first_event_cell_p_energy_share=_ratio(float(penergy[first]),float(penergy.sum())),
        first_event_cell_correction_energy_share=_ratio(float(cenergy[first]),float(cenergy.sum())),
        first_event_cell_pressure_only_interior_t=pressure_only_t,
        pressure_only_t_scope='uncorrected exact-pair formula; not a valid forecast when directional correction is non-negligible',
        p_energy_total=float(penergy.sum()),p_energy_top10_cells=np.argsort(-penergy)[:10].tolist(),
        actual_PA_multiplier=multiplier if finite_multiplier else None,
        multiplier_convention='P.last_A_info.multiplier=mu for half-squared raw ball constraint; T=sS and W=I; no factor-two conversion',
        projector_last_A_info={k:(None if isinstance(v,float) and not np.isfinite(v) else v) for k,v in info.items()},
        KKT_qualification=dict(raw_receiver_metric='W=I, no masks/whitening',Slater=slater,
            Slater_margin=slater_margin,delta=P.delta,minimum_reachable_residual=P.minimum_reachable_residual,
            finite_multiplier=finite_multiplier,pseudoinverse_limit=bool(info['pseudoinverse_limit']),
            finite_multiplier_interpretation_eligible=valid_multiplier,ball_residual_over_delta=_ratio(rnorm,P.delta),
            ball_active=bool(finite_multiplier and multiplier>0),
            material_amplitude_normal_defect=float(np.linalg.norm(dQ+(m/P.scale)[:,None]*dX)),
            propagation_pair_defect_compatibility=float(np.linalg.norm(eQ-compatible_eQ)),
            physical_material_state_residual=float(state['state_relative_residual']),
            numerical_KKT=kkt,invertibility_scope='existing retained Maxwell LU and finite solves; no global condition-number or pole certificate'),
        explicit_diagnostic_operator_cost=dict(cost),physical_model_delta=_delta(model.cost,physics_before),
        projector_delta=_delta(P.costs,projector_before),CPU_seconds=time.process_time()-started,
        wall_seconds=time.perf_counter()-wall,truth_access='NONE',
        interpretation='actual directional correction omitted by a scalar pair-defect ratio; not a global amplification/recovery/future-trajectory certificate')
    arrays=dict(aX=aX,aQ=aQ,cX=cX,cQ=cQ,dX=dX,dQ=dQ,eX=eX,eQ=eQ,
        a_prime_X=primeX,a_prime_Q=primeQ,material=m,p=p,correction=correction,
        receiver_rhs=receive_rhs,correction_rhs=correction_rhs,p_energy=penergy,
        correction_energy=cenergy,geometric_thresholds=threshold['thresholds'])
    return row,arrays


def run_public_case(source,out):
    wall,cpu=time.perf_counter(),time.process_time();out.mkdir(parents=True,exist_ok=False)
    with np.load(source,allow_pickle=False) as z:
        if any(any(word in k.lower() for word in ('truth','m_true','ground','label')) for k in z.files):
            raise ValueError('Truth-like field in exposed public input')
        public={k:z[k].copy() for k in z.files}
    m0=public.pop('m_init');Y=public['Y'];model=MaxwellModel(public)
    delta=float(public.get('delta',float(public['sigma'])*np.sqrt(chi2.ppf(.99,2*Y.size)/2)))
    initial=model.state(m0)
    P=DenseFieldCurrentProjectors(model.G,model.E,model.S,Y,delta,scale=6.,lower=.02,upper=20.,grouped_ports=True)
    initialX,initialQ=initial['j']/P.scale,initial['Q']
    rows=[];route_cost=[]
    for route in ('AP128','DR128_single_A_shadow'):
        before=dict(P.costs);r_cpu,r_wall=time.process_time(),time.perf_counter()
        X,Q=initialX.copy(),initialQ.copy()
        for _ in range(128):
            cX,cQ,_=P.material(X,Q)
            if route=='AP128':X,Q=P.affine_ball(cX,cQ)
            else:
                aX,aQ=P.affine_ball(2*cX-X,2*cQ-Q)
                X,Q=X+aX-cX,Q+aQ-cQ
        if route!='AP128':
            cX,cQ,_=P.material(X,Q);X,Q=P.affine_ball(cX,cQ)
        route_cost.append(dict(route=route,projector_delta=_delta(P.costs,before),
            CPU_seconds=time.process_time()-r_cpu,wall_seconds=time.perf_counter()-r_wall))
        row,arrays=diagnose_pair(model,P,X,Q,route)
        np.savez_compressed(out/(route+'_FIELDS.npz'),**arrays)
        rows.append(row)
    report=dict(status='EXPOSED_PUBLIC_DIRECTIONAL_DUAL_PRESSURE_DIAGNOSTIC',truth_access='NONE',
        source=str(source),theory='parent THEORY_DERIVATIONS.md section4.1; mathematical validity remains parent-owned',
        formula='dX + lambda*s*p = correction = -A_m^{-H}(eX+s*G^H eQ), p=A_m^{-H} S^H r_prime',
        raw_receiver_metric='W=I',projector_setup_CPU_seconds=P.setup_cpu_seconds,
        projector_setup_wall_seconds=P.setup_wall_seconds,routes=rows,preparation_cost=route_cost,
        physical_model_cost=dict(model.cost),projector_cost=dict(P.costs),
        total_dense_operator_cost=operator_cost_summary(model.cost,P.costs,rows,model.n,model.ns),
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        timer_scope='input/model/projector setup, both128step routes, diagnostics and field serialization; final report serialization excluded',
        scope='four previously exposed public cases only; one fixed projector per case; no GT/literature search/fresh or remote work')
    (out/'DUAL_PRESSURE.json').write_text(json.dumps(report,indent=2,allow_nan=False,default=_json_default)+'\n')
    return report


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if len(a.input)>4:raise ValueError('Bounded diagnostic accepts at most four exposed public cases')
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        for source in a.input:
            r=run_public_case(source.resolve(),a.output/source.stem)
            print(json.dumps(dict(case=source.stem,CPU_seconds=r['CPU_seconds'],
                rows=[{k:v for k,v in x.items() if k in ['route','relative_pair_defect','directional_correction_relative_to_dX',
                    'correction_amplification_per_pair_defect','geometric_first_event_t','first_event_cell_p_energy_descending_rank']}
                    for x in r['routes']])),flush=True)


if __name__=='__main__':main()
