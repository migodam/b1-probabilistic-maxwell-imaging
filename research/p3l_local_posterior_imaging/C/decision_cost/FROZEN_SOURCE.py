"""Frozen DEV_000 conditional-u MH decision/cost control; explicit paid job.

Authoring this file does not execute it. No truth, training or network access.
All candidate decisions are frozen before exact T-matrix/Full predictions.
This evaluates independent anchor-based decisions, never a posterior chain.
"""
from pathlib import Path
import json
import shutil
import sys
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save
from p3l.decision_bounds import (BOUND_LABEL,stability_from_saved_inverse,
                                  dual_corrected_prediction,same_uniform_decision)

PARENT='DEV_000'
JOB='decision_cost_DEV000'
STEPS=(.01,.02,.04,.08)
PROPOSAL_SEEDS=(2026102290,2026102291,2026102292,2026102293)
UNIFORM_SEED=2026102297


def _relative(value,reference):
    return float(np.linalg.norm(value-reference)/max(1e-30,np.linalg.norm(reference)))


def _without_arrays(record):
    """Large arrays live in the companion NPZ, not a JSON list expansion."""
    return {key:value for key,value in record.items() if not isinstance(value,np.ndarray)}


def _unwhite(value,sigma):
    value=np.asarray(value).reshape(8,64)*sigma/np.sqrt(2)
    return value[:,:32].T+1j*value[:,32:].T


def _project(L,rhs,S,V,j0,increment_rhs):
    start=time.perf_counter()
    LV=L@V
    reduced=V.conj().T@LV
    common_wall=time.perf_counter()-start
    start=time.perf_counter()
    projected_rhs=V.conj().T@rhs
    coefficients=np.linalg.solve(reduced,projected_rhs)
    j=V@coefficients
    direct_prediction=S@j
    direct_wall=time.perf_counter()-start
    start=time.perf_counter()
    # Same reduced operator and basis, but a separately counted RHS solve.
    increment_coefficients=np.linalg.solve(reduced,V.conj().T@increment_rhs)
    j_a17=j0+V@increment_coefficients
    a17_prediction=S@j_a17
    a17_wall=time.perf_counter()-start
    counts=dict(projected_factorizations=1,projected_linear_solves=2,
        projected_rhs_columns=16,full_physical_RHS=0,
        matrix_actions=[dict(label='L_V',rhs_columns=V.shape[1]),
          dict(label='Vstar_LV',rhs_columns=V.shape[1]),
          dict(label='Vstar_rhs',rhs_columns=8),dict(label='V_coefficients',rhs_columns=8),
          dict(label='S_j',rhs_columns=8),dict(label='Vstar_increment_rhs',rhs_columns=8),
          dict(label='V_increment_coefficients',rhs_columns=8),dict(label='S_j_A17',rhs_columns=8)],
        shared_operator=True,a17_not_independent_increment=True)
    # np.linalg.solve refactorizes on each call: count both, do not claim one.
    counts['projected_factorizations']=2
    return dict(prediction=direct_prediction,a17_prediction=a17_prediction,
        current=j,a17_current=j_a17,reduced_operator=reduced,
        reduced_rhs=projected_rhs,coefficients=coefficients,
        increment_coefficients=increment_coefficients,
        common_operator_wall_seconds=common_wall,
        direct_rhs_solve_output_wall_seconds=direct_wall,
        a17_rhs_solve_output_wall_seconds=a17_wall,
        direct_with_common_operator_wall_seconds=common_wall+direct_wall,
        a17_with_common_operator_wall_seconds=common_wall+a17_wall,
        direct_A17_relative_difference=_relative(j_a17,j),counts=counts)


def _krylov(book,L,rhs,S,R0,index):
    """Eight separately paid Krylov RHS; saved inverse has no new LU charge."""
    started=time.perf_counter()
    book._update('GMRES_FULL_ATTEMPT',states=1,rhs=8,lu=0,resolution=24,
        proposal=index,role='decision_GMRES_R0_preconditioned')
    actions=dict(operator_matvecs=0,preconditioner_matvecs=0)
    def operator(v):
        actions['operator_matvecs']+=1
        return L@v
    def preconditioner(v):
        actions['preconditioner_matvecs']+=1
        return R0@v
    Lop=LinearOperator(L.shape,matvec=operator,dtype=np.complex128)
    Mop=LinearOperator(L.shape,matvec=preconditioner,dtype=np.complex128)
    current=np.zeros_like(rhs);columns=[]
    with book.physical('decision_GMRES_R0_preconditioned_8RHS'):
        for column in range(rhs.shape[1]):
            residual_history=[]
            solution,info=gmres(Lop,rhs[:,column],M=Mop,rtol=1e-10,atol=0.,
                restart=40,maxiter=10,callback=lambda value:residual_history.append(float(value)),
                callback_type='pr_norm')
            current[:,column]=solution
            residual=L@solution-rhs[:,column]
            backward=float(np.linalg.norm(residual)/max(1e-30,np.linalg.norm(rhs[:,column])))
            columns.append(dict(source=column,info=int(info),iterations=len(residual_history),
                relative_rhs_residual=backward,preconditioned_residual_history=residual_history,
                converged=bool(info==0 and np.isfinite(backward) and backward<=1e-10)))
        prediction=S@current
    return dict(status='CONVERGED' if all(v['converged'] for v in columns) else 'NONCONVERGED_REQUIRES_FULL_FALLBACK',
        current=current,prediction=prediction,columns=columns,counts=actions,
        RHS=8,states=1,LU=0,rtol=1e-10,atol=0,restart=40,maxiter=10,
        diagnostic_residual_matvecs=8,output_matrix_rhs_columns=8,
        wall_seconds=time.perf_counter()-started,saved_inverse_setup_shared=True)


def main():
    destination=OUT/'C/decision_cost'
    if destination.exists():
        raise FileExistsError('decision-cost output exists; no implicit replay/overwrite')
    anchor_path=OUT/f'A/anchors/{PARENT}.npz'
    # This preflight reads only an already-paid cache, never a truth file.
    with np.load(anchor_path,allow_pickle=False) as data:
        if 'receiver_T' not in data.files or data['receiver_T'].shape!=(32,576):
            raise RuntimeError('already-paid anchor receiver_T is required; no fill-in solve')
        q0=data['q'].copy()
    uniform_rng=np.random.default_rng(UNIFORM_SEED)
    proposals=[]
    for step in STEPS:
        for seed in PROPOSAL_SEEDS:
            q=q0.copy();q[:64]+=step*np.random.default_rng(seed).standard_normal(64)
            proposals.append(dict(index=len(proposals),step=step,seed=seed,q=q,
                logu=float(np.log(uniform_rng.random())),
                logprior_proposal_ratio=float(-.5*(q[:64]@q[:64]-q0[:64]@q0[:64]))))
    destination.mkdir(parents=True)
    protocol=dict(status='FROZEN_BEFORE_PAID_SETUP_AND_FULL_CANDIDATE_REFERENCE',
        parent=PARENT,job=JOB,steps=STEPS,proposal_seeds=PROPOSAL_SEEDS,uniform_seed=UNIFORM_SEED,
        proposal_count=16,proposals=proposals,
        generator='fixed symmetric native z64 Gaussian random walk from same anchor; native u2 fixed',
        seed_reuse='same four fixed Gaussian directions across steps; no distance matching or selection',
        q_proposal_log_ratio=0,prior='original Gaussian64; properties conditioned at anchor',
        posterior_chain=False,GT_used=False,bound_label=BOUND_LABEL,
        decision_rule='same log uniform; freeze decision before exact T-matrix and Full candidate reference',
        bases=[8,16],full_reference='one paid Target.evaluate per candidate; shared validation/fallback',
        fallback_reuses_validation_query=True,no_free_oracle=True,
        GMRES=dict(status='FROZEN_TO_RUN',rtol=1e-10,atol=0,restart=40,maxiter=10,
            callback_type='pr_norm',preconditioner='already-paid saved R0',RHS_per_proposal=8,
            state_attempts_per_proposal=1,new_LU_per_proposal=0),
        warm_array_replays=dict(status='NOT_RUN',reason='excluded from fixed first decision-cost control'),
        interpretation='conditional floating-point algebra only; unavailable Neumann bound does not imply singularity')
    write(destination/'FROZEN_CONFIG.json',protocol)
    shutil.copyfile(Path(__file__),destination/'FROZEN_SOURCE.py')
    book=Book('C',JOB,max_wall=3600)
    status='COMPLETE';error=None;rows=[];setup_started=time.perf_counter()
    try:
        target=Target(book,keep_fields=True)
        anchor=target.load_cache(anchor_path,owner='paid DEV_000 decision-control anchor; no GT')
        if anchor.T is None:raise RuntimeError('receiver_T unexpectedly absent after paid-cache load')
        if not np.array_equal(q0,anchor.q):raise RuntimeError('anchor changed after proposal freeze')
        y=target.observation(PARENT)
        model=target.backend()
        L0=anchor.state.matrix
        S,G,incident=model.S,model.G,model.incident
        n=L0.shape[0]
        if n!=576 or anchor.state.current.shape!=(576,8):raise RuntimeError('TM24 layout required')
        inverse_started=time.perf_counter()
        with book.physical('decision_stability_inverse_576_RHS'):
            R0=model.solve_anchor(np.eye(n,dtype=complex),anchor.state,role='decision_stability_inverse')
        inverse_wall=time.perf_counter()-inverse_started
        Z0=anchor.T.conj().T
        V8=np.linalg.qr(anchor.state.current,mode='reduced')[0]
        Jwhite=target.jacobian(anchor) # cached T only; no new receiver RHS
        phi0=target.phi(anchor,y)
        setup_wall=time.perf_counter()-setup_started
        save(destination/'SETUP.npz',q0=q0,chi0=anchor.state.chi,current0=anchor.state.current,
            field0=anchor.state.field,prediction0=anchor.state.data,L0=L0,R0=R0,
            S=S,G=G,incident=incident,Z0=Z0,V8=V8,Jwhite=Jwhite,observation=y)
        write(destination/'SETUP.json',dict(status='PAID_SETUP_COMPLETE',
            inverse_RHS=576,inverse_new_LU=0,inverse_wall_seconds=inverse_wall,
            cold_setup_wall_seconds=setup_wall,anchor_phi=phi0,sigma=target.sigma,
            receiver_T_prepaid=True,new_receiver_RHS=0,
            setup_counts=dict(book.counts),bound_label=BOUND_LABEL))
        for proposal in proposals:
            index=proposal['index'];folder=destination/f'proposal_{index:02d}'
            folder.mkdir();proposal_started=time.perf_counter();start=proposal_started;q=np.asarray(proposal['q'])
            book.context=dict(parent=PARENT,proposal=index,step=proposal['step'],seed=proposal['seed'])
            chi,_=target.prior.decode(q)
            dchi=chi-anchor.state.chi
            # All candidate matrix setup uses exact signed complex delta chi.
            Lnew=np.eye(n,dtype=complex)-chi[:,None]*G
            rhs=chi[:,None]*incident
            increment_rhs=dchi[:,None]*anchor.state.field
            matrix_setup_wall=time.perf_counter()-start
            start=time.perf_counter()
            U=R0@increment_rhs
            V16=np.linalg.qr(np.column_stack((anchor.state.current,U)),mode='reduced')[0]
            enrichment_wall=time.perf_counter()-start
            start=time.perf_counter()
            stability=stability_from_saved_inverse(L0,R0,Lnew)
            stability_wall=time.perf_counter()-start
            results={};bounds={};decisions={};bound_walls={}
            for dimension,V in ((8,V8),(16,V16)):
                name=str(dimension)
                results[name]=_project(Lnew,rhs,S,V,anchor.state.current,increment_rhs)
                start=time.perf_counter()
                bounds[name]=dual_corrected_prediction(Lnew,rhs,S,V,Z0,stability,
                    observation=y,sigma=target.sigma)
                bound_walls[name]=time.perf_counter()-start
                decisions[name]=same_uniform_decision((phi0,phi0),bounds[name].get('likelihood_interval'),
                    proposal['logprior_proposal_ratio'],proposal['logu'])
            # This irreversible-on-disk ordering prevents Full reference from
            # selecting the basis, fallback action, proposal or shared uniform.
            frozen=dict(proposal=proposal,stability=stability,decisions=decisions,
                bounds={k:_without_arrays(v) for k,v in bounds.items()},
                full_candidate_prediction_available=False,exact_Tmatrix_prediction_available=False,
                bound_label=BOUND_LABEL,GT_used=False)
            write(folder/'PRE_REFERENCE_DECISIONS.json',frozen)
            pre_reference_wall=time.perf_counter()-proposal_started
            # Frozen linearized comparators are computed independently of Full.
            start=time.perf_counter()
            frozen_latent=anchor.state.data+_unwhite(Jwhite@(q-q0),target.sigma)
            frozen_material=anchor.state.data+anchor.T@increment_rhs
            frozen_linear_wall=time.perf_counter()-start
            # Exact full-dimensional T matrix gets a full paid state reservation
            # before matrix setup/LU/solve; no callback credit is created.
            start=time.perf_counter()
            book._update('EXACT_TMATRIX_FULL_ATTEMPT',states=1,lu=1,rhs=8,
                resolution=24,proposal=index,role='decision_exact_Tmatrix_576')
            with book.physical('decision_exact_Tmatrix_576_LU8'):
                RdeltaG=R0@(dchi[:,None]*G)
                operator=np.eye(n,dtype=complex)-RdeltaG
                t_rhs=R0@increment_rhs
                factor=lu_factor(operator)
                increment=lu_solve(factor,t_rhs)
                t_current=anchor.state.current+increment
                t_prediction=S@t_current
            t_wall=time.perf_counter()-start
            krylov=_krylov(book,Lnew,rhs,S,R0,index)
            # Reference occurs only after both basis decisions were frozen.
            start=time.perf_counter()
            full=target.evaluate(q,'decision_cost_full_validation_shared_fallback')
            full_wall=time.perf_counter()-start
            phi=target.phi(full,y)
            full_log_ratio=proposal['logprior_proposal_ratio']+phi0-phi
            full_action='ACCEPT' if proposal['logu']<min(0.,full_log_ratio) else 'REJECT'
            candidate_predictions=dict(full=full.state.data,exact_Tmatrix=t_prediction,
                frozen_latent_J=frozen_latent,frozen_material_linear=frozen_material,
                GMRES=krylov['prediction'])
            basis_records={}
            raw=dict(q=q,chi=chi,delta_chi=dchi,Lnew=Lnew,rhs=rhs,increment_rhs=increment_rhs,
                U=U,V8=V8,V16=V16,t_matrix_operator=operator,t_matrix_rhs=t_rhs,
                t_matrix_LU=factor[0],t_matrix_pivots=factor[1],t_matrix_increment=increment,
                t_matrix_current=t_current,full_current=full.state.current,
                GMRES_current=krylov['current'])
            for name,result in results.items():
                candidate_predictions['ROM'+name]=result['prediction']
                candidate_predictions['A17'+name]=result['a17_prediction']
                bound=bounds[name];interval=bound.get('likelihood_interval')
                action=decisions[name]['action']
                final_action=full_action if action=='REFINE' else action
                basis_records[name]=dict(ROM_A17=_without_arrays(result),
                    bound=_without_arrays(bound),decision=decisions[name],
                    final_action=final_action,shared_paid_full_fallback=action=='REFINE',
                    decision_mismatch=final_action!=full_action,
                    nonfallback_decision_mismatch=(action!='REFINE' and action!=full_action),
                    full_phi_enclosed=None if interval is None else bool(interval[0]<=phi<=interval[1]),
                    bound_wall_seconds=bound_walls[name],
                    projected_cost_note='ordinary ROM/A17 separate from bound helper; helper counts include another projected solve')
                for key,value in result.items():
                    if isinstance(value,np.ndarray):raw[f'basis{name}_{key}']=value
                for key,value in bound.items():
                    if isinstance(value,np.ndarray):raw[f'bound{name}_{key}']=value
                if bound.get('prediction') is not None:
                    candidate_predictions['dual_corrected'+name]=bound['prediction']
            errors={}
            for name,value in candidate_predictions.items():
                comparator_phi=float(np.linalg.norm(value-y)**2/target.sigma**2)
                comparator_action='ACCEPT' if proposal['logu']<min(0.,proposal['logprior_proposal_ratio']+phi0-comparator_phi) else 'REJECT'
                errors[name]=dict(relative_complex_prediction_error=_relative(value,full.state.data),
                    absolute_complex_prediction_error=float(np.linalg.norm(value-full.state.data)),
                    Phi=comparator_phi,same_uniform_uncertified_action=comparator_action,
                    same_uniform_uncertified_decision_mismatch=comparator_action!=full_action)
            raw.update({f'prediction_{k}':v for k,v in candidate_predictions.items()})
            save(folder/'ARRAYS.npz',**raw)
            row=dict(proposal=proposal,full_phi=phi,anchor_phi=phi0,full_action=full_action,
                full_log_ratio=full_log_ratio,bases=basis_records,errors=errors,
                timings=dict(matrix_setup_wall_seconds=matrix_setup_wall,
                    enrichment_wall_seconds=enrichment_wall,stability_wall_seconds=stability_wall,
                    frozen_linear_comparators_wall_seconds=frozen_linear_wall,
                    exact_Tmatrix_wall_seconds=t_wall,full_reference_wall_seconds=full_wall,
                    pre_reference_bounds_and_ROM_wall_seconds=pre_reference_wall,
                    GMRES_wall_seconds=krylov['wall_seconds'],
                    whole_proposal_including_array_persistence_wall_seconds=time.perf_counter()-proposal_started),
                setup_shared_over_16=setup_wall/16,
                enrichment_counts=dict(matrix_actions=[dict(label='R0_deltaChiE0',rhs_columns=8)],
                    QR_columns=16,new_physical_RHS=0),
                exact_Tmatrix_cost=dict(RHS=8,states=1,LU=1,
                    matrix_actions=[dict(label='R0_deltaChiG',rhs_columns=576),
                      dict(label='R0_deltaChiE0',rhs_columns=8),dict(label='S_j',rhs_columns=8)],
                    matrix_setup_included_in_wall=True),
                full_reference_cost=dict(RHS=8,states=1,LU=1,shared_validation_fallback=True),
                GMRES=dict(**_without_arrays(krylov),
                    shared_paid_full_fallback=krylov['status']!='CONVERGED',
                    final_action=full_action if krylov['status']!='CONVERGED' else errors['GMRES']['same_uniform_uncertified_action']),
                warm_array_replays=protocol['warm_array_replays'],
                bound_label=BOUND_LABEL,GT_used=False,posterior_chain=False)
            write(folder/'RESULT.json',row);rows.append(row)
            model._transfers.clear()
            del full,raw,candidate_predictions,results,bounds,factor,operator,krylov
        summary=dict(status='COMPLETE_CONDITIONAL_U_DECISION_COST_CONTROL',
            parent=PARENT,proposals=len(rows),posterior_chain=False,GT_used=False,
            bound_label=BOUND_LABEL,setup_cold_wall_seconds=setup_wall,
            setup_inverse_RHS=576,candidate_exact_Tmatrix_RHS=16*8,full_reference_RHS=16*8,
            candidate_GMRES_RHS=16*8,total_new_physical_RHS=576+16*24,
            total_new_state_attempts=48,total_new_LU_attempts=32,
            by_basis={name:dict(available_intervals=sum(row['bases'][name]['full_phi_enclosed'] is not None for row in rows),
                observed_enclosure_failures=sum(row['bases'][name]['full_phi_enclosed'] is False for row in rows),
                nonfallback_decisions=sum(row['bases'][name]['decision']['action']!='REFINE' for row in rows),
                shared_paid_full_fallbacks=sum(row['bases'][name]['shared_paid_full_fallback'] for row in rows),
                decision_mismatches=sum(row['bases'][name]['decision_mismatch'] for row in rows)) for name in ('8','16')},
            GMRES=dict(status='RUN',converged_proposals=sum(row['GMRES']['status']=='CONVERGED' for row in rows),
                full_fallbacks=sum(row['GMRES']['shared_paid_full_fallback'] for row in rows)),
            warm_array_replays=protocol['warm_array_replays'],
            interpretation='fixed proposals and shared uniforms; conditional floating-point bounds, not interval-certified; validation/fallback Full query shared and paid; no posterior mixing claim')
        write(destination/'SUMMARY.json',summary)
    except BaseException as exc:
        status='ERROR';error=repr(exc)
        write(destination/'FAILURE.json',dict(error=error,completed_proposals=len(rows),paid_costs_retained=True,
            automatic_replay=False,bound_label=BOUND_LABEL))
        raise
    finally:
        book.close(status,error)


if __name__=='__main__':main()
