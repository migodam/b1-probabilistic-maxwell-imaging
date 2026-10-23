"""Single-frequency finite-state Green reuse engineering control; NOT posterior.

Only main opens Book('P1B') after the qualified probability-reference gate.
No observations/truth are read. No physics, budget or Family I/O at import.
Uses explicit private FEM interfaces inspected in wq/fem.py; never substitutes
FE point receivers for the required free-Green volume receiver functional.
"""
from __future__ import annotations

import argparse
import csv
import gc
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.linalg import lu_factor,lu_solve,qr
from scipy.sparse.linalg import LinearOperator,gmres,norm as sparse_norm
from skfem import asm
from wq.fem import FEM,PMLConfig,_forcing,_contrast_mass,_sparse_bytes

from .accounting import ROOT,Book,BudgetStop,write,encode
from .native_world import LinearWorld
from .physical_witness import mesh_for
from .probability_family import Family


IDENTITY='P1B_SINGLE_FREQUENCY_ENGINEERING_CONTROL_NOT_FIVE_FREQUENCY_POSTERIOR'
SCOPE='green_reuse_single_frequency_v1'
CENTER=(-.8,2.15)
RADIUS=.12
SIGMA=.02
CONTRASTS=np.linspace(2.94,3.06,16)


def _ledger(book):
    if not book.path.exists():return []
    return [json.loads(line) for line in book.path.read_text().splitlines() if line.strip()]


class _Allowance:
    """Persistent <=4000 RHS sublimit inside existing P1B, serial queue only."""
    def __init__(self,book,out,cap=4000):
        self.book,self.out,self.cap=book,out,int(cap)
        self.cpu0,self.wall0=time.process_time(),time.perf_counter()
        self.first=len(_ledger(book))
    def __getattr__(self,name):return getattr(self.book,name)
    def used(self):
        return sum(r.get('rhs',0) for r in _ledger(self.book)
                   if r.get('stage')=='P1B' and r.get('scope')==SCOPE and r.get('event')=='reserve')
    def cost(self):
        rows=[r for r in _ledger(self.book)[self.first:] if r.get('scope')==SCOPE]
        rr=[r for r in rows if r.get('event')=='reserve']
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return dict(stage='P1B',scope=SCOPE,scope_rhs_cap=self.cap,
                    scope_rhs_all_attempts=self.used(),new_rhs=sum(r.get('rhs',0) for r in rr),
                    new_material_frequency_LU=sum(r.get('factorizations',0) for r in rr),
                    reserve_ids=[r['id'] for r in rr],
                    cpu_seconds=time.process_time()-self.cpu0,wall_seconds=time.perf_counter()-self.wall0,
                    process_peak_RSS_bytes=int(rss if sys.platform=='darwin' else 1024*rss),
                    cpu_accounting='measured only; existing launcher charges CPU once')
    def event(self,**row):
        self.book.event(scope=SCOPE,**row)
        with (self.out/'EVENTS.jsonl').open('a') as stream:
            stream.write(json.dumps(row,default=encode,allow_nan=False)+'\n')
        write(self.out/'LIVE_COST.json',self.cost())
    def reserve(self,phase,label,rhs=0,factorizations=0,**metadata):
        if self.used()+rhs>self.cap:
            self.event(event='engineering_scope_stop',phase=phase,label=label,requested_rhs=rhs)
            raise BudgetStop('green reuse persistent 4000 RHS sublimit')
        result=self.book.reserve(phase,label,rhs,factorizations,scope=SCOPE,**metadata)
        write(self.out/'LIVE_COST.json',self.cost())
        return result


def _timed(book,label,operation,**metadata):
    cpu,wall=time.process_time(),time.perf_counter()
    result=operation()
    book.event(event='engineering_operation',label=label,
               operation_cpu_seconds=time.process_time()-cpu,
               operation_wall_seconds=time.perf_counter()-wall,**metadata)
    return result


def _rhs(p,loads):
    return np.asfortranarray(np.column_stack([
        asm(_forcing,p.material_basis,load=loads[:,j].reshape(p._material_shape))[p._free_dofs]
        for j in range(loads.shape[1])]))


def _action(p,s0,rhs,label,trans='N'):
    rhs=np.asfortranarray(rhs)
    if rhs.ndim==1:rhs=rhs[:,None]
    p.book.reserve(p.phase,label,rhs=rhs.shape[1],factorizations=0,transpose=trans,
                   base_LU_reused=True,freq=1.)
    result=_timed(p.book,label,lambda:s0.factor.solve(rhs,trans=trans),
                  RHS=rhs.shape[1],transpose=trans)
    matrix=s0._matrix if trans=='N' else s0._matrix.conjugate().T
    residual=float(np.max(np.linalg.norm(matrix@result-rhs,axis=0)
                         /np.maximum(np.linalg.norm(rhs,axis=0),np.finfo(float).tiny)))
    p.book.event(event='base_LU_action_residual',label=label,RHS=rhs.shape[1],residual=residual)
    if not np.isfinite(residual) or residual>p.residual_tolerance:
        raise RuntimeError('base-LU action residual exceeds backend tolerance')
    return result


def _response(p,s0,epsilon1,du):
    dofs=np.zeros((p.ndofs,du.shape[1]),complex)
    dofs[p._free_dofs]=du
    sampled=np.column_stack([np.asarray(p.material_basis.interpolate(dofs[:,j])).ravel()
                             for j in range(du.shape[1])])
    total=s0.field+sampled
    current=(epsilon1-1.)[:,None]*total
    # REQUIRED volume functional; the point-probe receiver is never used here.
    return s0.direct+(p.k**2*p._receiver_green@(p.material_weights[:,None]*current)).T


def _orthogonal_columns(a):
    a=np.array(a,complex,copy=True)
    norms=np.linalg.norm(a,axis=0)
    a=a[:,norms>0]/norms[norms>0]
    q,r,pivot=qr(a,mode='economic',pivoting=True)
    diagonal=abs(np.diag(r))
    rank=int(np.sum(diagonal>max(a.shape)*np.finfo(float).eps*max(diagonal.max(initial=0.),1.)))
    if rank==0:raise RuntimeError('empty numerical reduced space')
    return q[:,:rank],dict(input_columns=a.shape[1],rank=rank,pivot_order=pivot.tolist(),
                          rank_rule='machine-precision pivoted QR; approximation ROM only, never exact Woodbury')


def _dense(book,label,a,b):
    def run():
        factor=lu_factor(a,check_finite=True)
        answer=lu_solve(factor,b,check_finite=True)
        if not np.all(np.isfinite(answer)):
            raise RuntimeError('nonfinite dense reduced/core solution; no least-squares fallback')
        return answer
    return _timed(book,label,run,dense_core_dimension=a.shape[0],dense_core_LU=1,
                  counted_as_material_frequency_LU=False)


def _snapshot(book):
    cost=book.cost()
    return {key:cost[key] for key in ('new_rhs','new_material_frequency_LU','cpu_seconds','wall_seconds')}


def _delta(after,before):return {key:after[key]-before[key] for key in before}


def _fem(seed,mesh,tx,rx,book,phase):
    return FEM(1.,tx,rx,book,h=.13,order=3,mesh=mesh,pml=PMLConfig(strength=6),
               phase=phase,receiver_mode='volume',residual_tolerance=1e-10,
               max_dofs=110000,max_factor_bytes=3500000000)


def main(gate_path='results/P1A/reference/GATE.json'):
    """No posterior calculation; probability gate is a prerequisite only."""
    gate_file=Path(gate_path)
    if not gate_file.is_absolute():gate_file=ROOT/gate_file
    gate=json.loads(gate_file.read_text())
    if gate.get('probability_reference_qualified') is not True:
        raise RuntimeError('P1B engineering control blocked: probability-reference gate not qualified')
    family=Family.load()
    q=family.walls[1].copy()
    seed=int(family.config['wall_support']['seed'])
    tx=np.asarray(family.config['acquisition']['tx'])
    rx=np.asarray(family.config['acquisition']['rx'])
    if tx.shape!=(4,2) or rx.shape!=(12,2):raise RuntimeError('requires original 4x12 acquisition')
    out=ROOT/'results/P1B/green_reuse_control'
    out.mkdir(parents=True,exist_ok=False)
    book=_Allowance(Book('P1B'),out,4000)
    rows,qualification=[],[]
    failure=None
    setup={}
    exact_status='NOT_PREPARED'
    max_memory_estimate=0
    p=s0=None
    write(out/'PLAN.json',dict(identity=IDENTITY,gate=str(gate_file),gate_snapshot=gate,
          public_family=family.config['family_id'],public_wall_branch=1,
          q_source=family.config['wall_support']['wall_parameter_file'],
          contrasts=CONTRASTS,target_center=CENTER,target_radius=RADIUS,target_sigma=SIGMA,
          freq=1.,gamma=1.,sigma_real_quadrature=.01,
          geometry=dict(h=.13,tol=2.5e-5,res=4,circle_segments=64,
                        fine_circle_segments=512,qualification_white_tolerance=.05),
          RHS_sublimit=4000,independent_engineering_comparison=True,
          replaces_five_frequency_posterior=False,
          methods=['OPTIMIZED_FULL','EXACT_WOODBURY','FIELD_GOAL_ROM_V',
                   'FIELD_GOAL_ROM_ENRICHED_QR_V_R0V','PRECONDITIONED_KRYLOV',
                   'CURRENT_GALERKIN_V','A17_ALGEBRAIC_CURRENT_FEEDBACK_V'],
          current_equation='z + deltaA R0 z = rhs_delta; du=R0 z',
          A17_claim='algebraic finite-state control only; not reproduction/certification of new A17 theory',
          trial_space_warning='same V in current equation is NOT same physical field trial space as field Galerkin(V)',
          enrichment_rights='paid R0V shared with ordinary field goal-ROM enriched QR[V,R0V]',
          exact_resource_gate=dict(active_rows_max=2048,R0U_bytes_max=3000000000,
                                   no_truncation_or_pseudoinverse=True),
          observation_reads=0,delta_log_likelihood='NOT_COMPUTED; response error only',probabilities=None))

    def persist():
        write(out/'RESULTS.json',rows)
        if rows:
            with (out/'RESULTS.csv').open('w',newline='') as stream:
                fields=list(dict.fromkeys(key for row in rows for key in row if key!='cost'))
                writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
                writer.writeheader();writer.writerows(rows)
        write(out/'SETUP.json',setup)

    try:
        # Sequential endpoint qualification: each response is a real 4-RHS Full.
        # No factor/State from these checks is retained into the comparison.
        for contrast in (float(CONTRASTS[0]),float(CONTRASTS[-1])):
            pair=[]
            for segments in (64,512):
                probe=(CENTER,RADIUS,contrast)
                geometry=mesh_for(seed,probe,h=.13,tol=2.5e-5,res=4,n=segments)
                pf=_fem(seed,geometry,tx,rx,book,f'P1B_green_qualification_n{segments}')
                sf=pf.solve(LinearWorld(seed,probe).coefficient(q,1.))
                pair.append(sf.data.copy())
                np.savez_compressed(out/f'QUALIFICATION_c{contrast:.2f}_n{segments}.npz',
                                    data=sf.data,residuals=sf.residuals,contrast=contrast,q=q)
                del pf,sf,geometry;gc.collect()
            error=float(np.linalg.norm(pair[0]-pair[1])/.01)
            qualification.append(dict(contrast=contrast,circle64_vs512_white=error,pass_=bool(error<=.05)))
            write(out/'GEOMETRY_QUALIFICATION.json',qualification)
            if error>.05:raise RuntimeError('circle64 failed .05-white fine512 endpoint qualification; branch stopped')
        baseline_probe=(CENTER,RADIUS,3.)
        geometry=mesh_for(seed,baseline_probe,h=.13,tol=2.5e-5,res=4,n=64)
        baseline_start=_snapshot(book)
        p=_fem(seed,geometry,tx,rx,book,'P1B_green_comparison')
        s0=p.solve(LinearWorld(seed,baseline_probe).coefficient(q,1.))
        A0=s0._matrix
        nfree=len(p._free_dofs)
        u0=s0.scattered_dofs[p._free_dofs].copy()
        mask=(np.linalg.norm(p.material_points-np.asarray(CENTER),axis=1)<RADIUS).astype(float)
        # A1=A0+(contrast-3)*D; rhs_delta=(contrast-3)*b_unit uses E0 TOTAL.
        D=_timed(book,'assemble_unit_deltaA',lambda:(-p.k**2*asm(
            _contrast_mass,p.material_basis,contrast=mask.reshape(p._material_shape))
            [p._free_dofs][:,p._free_dofs]).tocsc(),assembly=True)
        D.eliminate_zeros()
        bunit=_timed(book,'assemble_finite_update_rhs_unit',
                     lambda:_rhs(p,p.k**2*mask[:,None]*s0.field),assembly=True)
        binc=_timed(book,'assemble_full_rhs_contrast_increment',
                    lambda:_rhs(p,p.k**2*mask[:,None]*s0.incident),assembly=True)
        rhs0=_timed(book,'assemble_full_baseline_rhs',
                    lambda:_rhs(p,p.k**2*(s0.epsilon-1.)[:,None]*s0.incident),assembly=True)
        np.savez_compressed(out/'BASELINE.npz',data=s0.data,q=q,rhs=rhs0,
                            field=s0.field,epsilon=s0.epsilon,points=p.material_points,
                            weights=p.material_weights,residuals=s0.residuals)
        setup['baseline']=dict(cost=_delta(_snapshot(book),baseline_start),
                                free_dofs=nfree,actual_factor_bytes=p._factor_storage_bytes,
                                factor_heuristic_bytes=p.estimated_factor_bytes)

        # Ordinary goal-ROM: four paid primal R0 actions, twelve paid Hermitian
        # duals for ORIGINAL VOLUME output, not twelve point-receiver probes.
        rom_start=_snapshot(book)
        part_start=_snapshot(book)
        primal=_action(p,s0,bunit,'ROM_primal_R0_4source')
        setup['primal_R0_rhs_unit']=dict(cost=_delta(_snapshot(book),part_start),shared_by=['ROM','WOODBURY'])
        goal_rhs=_timed(book,'assemble_original_volume_goal_duals',
                        lambda:_rhs(p,p.k**2*(s0.epsilon-1.)[:,None]
                                    *p._receiver_green.T).conjugate(),assembly=True)
        part_start=_snapshot(book)
        dual=_action(p,s0,goal_rhs,'ROM_volume_goal_dual_12receiver',trans='H')
        setup['original_volume_duals']=dict(cost=_delta(_snapshot(book),part_start))
        V,qr_info=_timed(book,'ROM_normalize_pivoted_QR',
                        lambda:_orthogonal_columns(np.column_stack((primal,dual))))
        part_start=_snapshot(book)
        RV=_action(p,s0,V,'shared_R0V_current_and_enriched_field_ROM')
        setup['shared_R0V']=dict(cost=_delta(_snapshot(book),part_start),enrichment_rights='ordinary field enriched ROM and both current-space methods')
        W,w_info=_timed(book,'ordinary_enriched_QR_V_R0V',
                       lambda:_orthogonal_columns(np.column_stack((V,RV))))
        def projected(basis,label):
            a0v=_timed(book,label+'_A0_MVP',lambda:A0@basis,MVP_vectors=basis.shape[1])
            dv=_timed(book,label+'_D_MVP',lambda:D@basis,MVP_vectors=basis.shape[1])
            return basis.conj().T@a0v,basis.conj().T@dv,basis.conj().T@bunit
        G0,GD,gb=projected(V,'field_V')
        H0,HD,hb=projected(W,'field_enriched_W')
        DRV=_timed(book,'current_D_R0V_MVP',lambda:D@RV,MVP_vectors=RV.shape[1])
        feedback=V.conj().T@DRV
        current_rhs=V.conj().T@bunit
        current_gram=V.conj().T@V
        setup['shared_ROM']=dict(cost=_delta(_snapshot(book),rom_start),V=qr_info,W=w_info,
                                  shared_R0V_RHS=V.shape[1],goal_operator='baseline original volume functional',
                                  paid_enrichment_available_to_ordinary_ROM=True)
        np.savez_compressed(out/'ROM_SETUP.npz',V=V,R0V=RV,W=W,primal=primal,dual=dual,
                            rhs_unit=bunit,goal_rhs=goal_rhs)
        del dual,goal_rhs;gc.collect()

        # Exact selector update: ALL structurally active rows, no singular-value
        # truncation. C may be singular; (I+C U^T R0U) needs no C inverse.
        nonzero_rows,nonzero_columns=D.nonzero()
        idx=np.unique(nonzero_rows)
        outside_columns=np.setdiff1d(np.unique(nonzero_columns),idx)
        setup['exact_selector_support']=dict(
            active_row_count=len(idx),nonzero_column_count=len(np.unique(nonzero_columns)),
            columns_outside_active_rows=outside_columns.tolist(),
            selector_factorization_valid=not bool(len(outside_columns)))
        if len(outside_columns):
            exact_status='NOT_RUN_INVALID_SELECTOR_SUPPORT'
            raise RuntimeError('D has nonzero columns outside active row selector; U C U.T is not exact')
        r=len(idx)
        if r==0:raise RuntimeError('empty structural target update; no exact Woodbury claim')
        ru_bytes=int(nfree*r*np.dtype(complex).itemsize)
        dense_peak=int(6*r*r*np.dtype(complex).itemsize)
        max_memory_estimate=int(ru_bytes+dense_peak+2*p._factor_storage_bytes+
                                V.nbytes+RV.nbytes+W.nbytes+primal.nbytes+
                                4*64*nfree*np.dtype(complex).itemsize+
                                p._receiver_green.nbytes+_sparse_bytes(A0)+_sparse_bytes(D))
        setup['exact_resource']=dict(active_rows=r,estimated_R0U_bytes=ru_bytes,
                                     estimated_dense_workspace_peak_bytes=dense_peak,
                                     estimated_combined_major_array_peak_bytes=max_memory_estimate,
                                     estimate_excludes_opaque_basis_and_LU_workspaces=True)
        R0U=Cunit=None
        if r>2048 or ru_bytes>3000000000:
            exact_status='NOT_RUN_RESOURCE_SCALE'
        elif book.used()+r>book.cap:
            exact_status='NOT_RUN_INSUFFICIENT_REMAINING_SCOPE_RHS'
        else:
            exact_start=_snapshot(book)
            R0U=np.empty((nfree,r),complex,order='F')
            for first in range(0,r,64):
                stop=min(first+64,r)
                selector=np.zeros((nfree,stop-first),complex,order='F')
                selector[idx[first:stop],np.arange(stop-first)]=1.
                R0U[:,first:stop]=_action(p,s0,selector,'EXACT_WOODBURY_R0U_chunk64')
                del selector
            Cunit=D[idx][:,idx].toarray()
            exact_status='READY_EXACT_ALL_ACTIVE_ROWS'
            setup['exact_R0U']=dict(cost=_delta(_snapshot(book),exact_start),rank_bound=r,
                                     exact_all_active_rows=True,no_truncation=True)
        setup['exact_status']=exact_status
        setup['cold_cost_convention']=dict(
            geometry_qualification='all four paid Full endpoints are benchmark validation, never hidden as free',
            optimized_full='fresh candidate LU + four RHS; baseline/ROM setup is not required for standalone deployment',
            woodbury='baseline four RHS/LU + paid four primal R0(rhs_unit) + all r R0U RHS; no goal/ROM setup required',
            field_V='baseline four RHS/LU + four primal R0 RHS + twelve original-volume dual RHS + QR/projections',
            field_enriched_and_current='same field_V setup + rank(V) paid shared R0V RHS + own QR/projections',
            Krylov='baseline four RHS/LU preconditioner + every actual base-LU preconditioner apply paid warm',
            standalone_setup_wall='use separately measured component events; pooled setup is not credited free to any arm',
            no_teacher_or_dense_core_cost_exclusion=True,
            warm_measurement='deployment snapshot immediately after real response extraction; Full reuses sf.data',
            audit_measurement='added residual/error/formula checks separately measured after deployment; retained in audited_total_cost and full stage CPU',
            backend_required_internal_checks='existing FEM solve/base-action invariants retained in their actual operation cost')
        if R0U is None:
            # The exact-resource *potential* estimate is retained above, while
            # the enabled-method estimate excludes an unallocated giant R0U.
            max_memory_estimate=int(2*p._factor_storage_bytes+V.nbytes+RV.nbytes+
                                    W.nbytes+primal.nbytes+p._receiver_green.nbytes+
                                    _sparse_bytes(A0)+_sparse_bytes(D))
        persist()

        for number,contrast in enumerate(CONTRASTS):
            delta=float(contrast-3.)
            folder=out/f'candidate_{number:02d}';folder.mkdir()
            rhs=delta*bunit
            full_rhs=rhs0+delta*binc
            np.savez_compressed(folder/'RHS.npz',finite_update_rhs=rhs,full_rhs=full_rhs,
                                contrast=contrast,delta_contrast=delta)
            method_start=_snapshot(book)
            sf=p.solve(LinearWorld(seed,(CENTER,RADIUS,float(contrast))).coefficient(q,1.))
            # The backend already produced the required full volume response.
            # Freeze deployment cost BEFORE all finite-update audit bookkeeping.
            full_deployment_end=_snapshot(book)
            A1=sf._matrix
            owned_state_bytes=sum(getattr(s0,key).nbytes for key in
                                  ('epsilon','material_points','material_weights','incident','field','current','scattered_dofs'))
            active_extra=(R0U.nbytes+dense_peak if R0U is not None else 0)
            candidate_major_estimate=int(s0.metadata['actual_factor_storage_bytes']+
                sf.metadata['actual_factor_storage_bytes']+active_extra+2*owned_state_bytes+
                V.nbytes+RV.nbytes+W.nbytes+primal.nbytes+p._receiver_green.nbytes+
                _sparse_bytes(A0)+_sparse_bytes(A1)+_sparse_bytes(D)+12*nfree*16)
            max_memory_estimate=max(max_memory_estimate,candidate_major_estimate)
            truth=sf.data.copy()
            direct_du=sf.scattered_dofs[p._free_dofs]-u0
            matrix_identity_error=float(sparse_norm(A1-(A0+delta*D))/max(sparse_norm(A1),1e-300))
            np.savez_compressed(folder/'OPTIMIZED_FULL_REFERENCE.npz',data=truth,delta_u=direct_du,
                                full_rhs_residuals=sf.residuals,epsilon=sf.epsilon)

            def record(method,du,start,*,response=None,deployment_end=None,additional_validation=None,**metadata):
                if response is None:
                    response=_timed(book,method+'_volume_extraction',
                                    lambda:_response(p,s0,sf.epsilon,du))
                    # Every reused method freezes deployment after its actual
                    # response extraction, before added benchmark audits.
                    deployment_end=_snapshot(book)
                elif deployment_end is None:
                    raise RuntimeError('precomputed response requires an explicit deployment-end snapshot')
                if not np.all(np.isfinite(du)):
                    raise RuntimeError('nonfinite finite-update field; no response claimed')
                defect=_timed(book,method+'_full_matrix_residual_MVP',lambda:A1@du-rhs,
                               MVP_vectors=du.shape[1],cost_role='separate_validation')
                relative=np.linalg.norm(defect,axis=0)/np.maximum(np.linalg.norm(rhs,axis=0),1e-300)
                full_defect=_timed(book,method+'_full_rhs_residual_MVP',
                                    lambda:A1@(u0+du)-full_rhs,MVP_vectors=du.shape[1],
                                    cost_role='separate_validation')
                full_relative=np.linalg.norm(full_defect,axis=0)/np.maximum(np.linalg.norm(full_rhs,axis=0),1e-300)
                if not np.all(np.isfinite(response)) or not np.all(np.isfinite(relative)) or not np.all(np.isfinite(full_relative)):
                    raise RuntimeError('nonfinite response/residual; no engineering result claimed')
                # Validate the finite-increment volume formula independently on
                # the FIRST Full candidate only; never redo its deployment output.
                if method=='OPTIMIZED_FULL' and number==0:
                    extracted=_timed(book,'first_candidate_finite_volume_formula_audit',
                                     lambda:_response(p,s0,sf.epsilon,du),
                                     cost_role='separate_validation')
                    identity_error=float(np.linalg.norm(extracted-sf.data)/.01)
                    write(folder/'FINITE_VOLUME_IDENTITY_AUDIT.json',dict(
                        independent_formula_vs_backend_volume_error_white=identity_error,
                        cost_role='separate_validation',no_additional_PDE_RHS=True))
                    np.savez_compressed(folder/'FINITE_VOLUME_IDENTITY_AUDIT.npz',
                                        independently_extracted=extracted,backend_volume_data=sf.data)
                    if identity_error>1e-7:
                        raise RuntimeError('finite-update volume formula differs from backend original volume output')
                if additional_validation is not None:
                    metadata.update(_timed(book,method+'_additional_identity_audit',
                                           lambda:additional_validation(response),
                                           cost_role='separate_validation'))
                warm_deployment=_delta(deployment_end,start)
                row=dict(candidate=number,contrast=float(contrast),method=method,
                         finite_update_matrix_relative_residual_max=float(relative.max()),
                         full_rhs_relative_residual_max=float(full_relative.max()),
                         response_error_white=float(np.linalg.norm(response-truth)/.01),
                         response_error_relative=float(np.linalg.norm(response-truth)/max(np.linalg.norm(truth),1e-300)),
                         finite_response_change_white=float(np.linalg.norm(truth-s0.data)/.01),
                         matrix_assembly_identity_relative_error=matrix_identity_error,
                         warm_deployment_cost=warm_deployment,
                         warm_cost=warm_deployment,warm_rhs=warm_deployment['new_rhs'],
                         warm_material_frequency_LU=warm_deployment['new_material_frequency_LU'],
                         warm_cpu_seconds=warm_deployment['cpu_seconds'],
                         warm_wall_seconds=warm_deployment['wall_seconds'],
                         max_major_array_memory_estimate_bytes=max_memory_estimate,
                         cold_setup_cost_in_SETUP=True,
                         cost_partition='deployment through response; subsequent benchmark audits/output separate; all work remains in stage total',
                         **metadata)
                np.savez_compressed(folder/f'{method}.npz',data=response,delta_u=du,
                                    finite_rhs_relative_residuals=relative,
                                    full_rhs_relative_residuals=full_relative)
                audited_end=_snapshot(book)
                row['separate_validation_cost']=_delta(audited_end,deployment_end)
                row['audited_total_cost']=_delta(audited_end,start)
                row['post_row_persistence_accounted_in_stage_total']=True
                rows.append(row);persist()
                return response

            full_extraction=record('OPTIMIZED_FULL',direct_du,method_start,
                                   response=sf.data,deployment_end=full_deployment_end,
                                   status='COMPLETE_DIRECT_FULL_REFERENCE',setup_shared=True,
                                   original_backend_volume_output_reused=True)
            if R0U is not None:
                start=_snapshot(book)
                z=delta*primal # already-paid R0(rhs_unit), exactly reused by linearity
                C=delta*Cunit
                core=_timed(book,'exact_Woodbury_dense_core_assembly',
                            lambda:np.eye(r,dtype=complex)+C@R0U[idx],
                            dense_core_dimension=r,dense_matrix_multiplications=1)
                reduced=_dense(book,'exact_Woodbury_dense_core',core,C@z[idx])
                wb=_timed(book,'exact_Woodbury_reconstruct',lambda:z-R0U@reduced,
                          dense_matrix_shape=[nfree,r],source_columns=4)
                record('EXACT_WOODBURY',wb,start,status='EXACT_ALL_ACTIVE_ROWS_NO_TRUNCATION',
                       shared_paid_primal_R0_rhs=True,exact_setup_rhs_included_in_SETUP=True)
                del C,core,reduced,wb,z
            else:
                rows.append(dict(candidate=number,contrast=float(contrast),method='EXACT_WOODBURY',
                                 status=exact_status,active_rows=r,estimated_R0U_bytes=ru_bytes));persist()
            for name,basis,ga,gd,br in [('FIELD_GOAL_ROM_V',V,G0,GD,gb),
                                       ('FIELD_GOAL_ROM_ENRICHED_QR_V_R0V',W,H0,HD,hb)]:
                start=_snapshot(book)
                coefficients=_dense(book,name+'_dense_core',ga+delta*gd,delta*br)
                du=basis@coefficients
                record(name,du,start,status='APPROXIMATE_FIELD_GALERKIN',
                       reduced_rank=basis.shape[1],same_paid_R0V_enrichment_rights=True)
                del du,coefficients
            start=_snapshot(book)
            cgal=_dense(book,'current_Galerkin_dense_core',current_gram+delta*feedback,delta*current_rhs)
            du_g=RV@cgal
            response_g=record('CURRENT_GALERKIN_V',du_g,start,
                              status='APPROXIMATE_CURRENT_LOAD_GALERKIN',current_space_rank=V.shape[1])
            start=_snapshot(book)
            ca17=_dense(book,'A17_algebraic_feedback_dense_core',
                         np.eye(V.shape[1],dtype=complex)+delta*feedback,delta*current_rhs)
            du_a=RV@ca17
            def same_current_space_audit(response):
                same_space_error=float(np.linalg.norm(du_a-du_g)/max(np.linalg.norm(du_g),1e-300))
                response_error=float(np.linalg.norm(response-response_g)/.01)
                if same_space_error>1e-9 or response_error>1e-8:
                    raise RuntimeError('same-current-space A17/Galerkin algebraic identity check failed')
                return dict(same_current_space_relative_difference=same_space_error,
                            same_current_space_response_difference_white=response_error)
            response_a=record('A17_ALGEBRAIC_CURRENT_FEEDBACK_V',du_a,start,
                              additional_validation=same_current_space_audit,
                              status='SAME_CURRENT_SPACE_ALGEBRAIC_CONTROL_NOT_THEORY_REPRODUCTION',
                              independent_algorithmic_gain_claim=False)
            del cgal,ca17,du_g,du_a,response_g,response_a

            start=_snapshot(book)
            counts=dict(matrix_vector_products=0,matrix_vector_wall_seconds=0.,preconditioner_RHS=0)
            convergence=[]
            def mvp(v):
                tic=time.perf_counter();answer=A1@v
                counts['matrix_vector_products']+=1
                counts['matrix_vector_wall_seconds']+=time.perf_counter()-tic
                return answer
            def precondition(v):
                counts['preconditioner_RHS']+=1
                return _action(p,s0,np.asarray(v)[:,None],'KRYLOV_base_LU_apply_1RHS')[:,0]
            op=LinearOperator(A1.shape,matvec=mvp,dtype=complex)
            pre=LinearOperator(A1.shape,matvec=precondition,dtype=complex)
            krylov=np.empty_like(rhs)
            infos=[]
            for source in range(rhs.shape[1]):
                def callback(value):
                    convergence.append(dict(source=source,preconditioned_relative_residual=float(value)))
                    write(folder/'KRYLOV_PROGRESS.json',convergence)
                krylov[:,source],info=gmres(op,rhs[:,source],M=pre,rtol=1e-10,atol=0.,
                                           restart=30,maxiter=20,callback=callback,callback_type='pr_norm')
                infos.append(int(info))
            book.event(event='Krylov_MVP_and_base_action_cost',candidate=number,**counts,GMRES_info=infos)
            record('PRECONDITIONED_KRYLOV',krylov,start,
                   status='CONVERGED_GMRES' if all(i==0 for i in infos) else 'GMRES_NOT_CONVERGED',
                   GMRES_info=infos,**counts)
            del sf,krylov,A1,truth,direct_du,full_extraction
            gc.collect() # only baseline State/LU and the backend's current factor remain
    except Exception as exc:
        failure=dict(error_type=type(exc).__name__,error=str(exc),
                     status='INCOMPLETE_BUDGET_STOP' if isinstance(exc,BudgetStop) else 'INCOMPLETE_ENGINEERING_BRANCH',
                     no_recovery_solve=True)
        write(out/'FAILURE.json',failure)
    finally:
        persist()
        summary=dict(identity=IDENTITY,status='ENGINEERING_CONTROL_COMPLETE' if failure is None
                     else 'INCOMPLETE_PRESERVED_RESULTS',failure=failure,
                     geometry_qualification=qualification,exact_status=exact_status,
                     completed_method_rows=len(rows),cost=book.cost(),setup=setup,
                     max_major_array_memory_estimate_bytes=max_memory_estimate,
                     no_posterior_or_deltaPhi_claim=True,observations_read=False,
                     five_frequency_probability_result_replaced=False,
                     A17_scope='finite-state current-feedback algebra only; no novelty/certification claim',
                     ordinary_goal_ROM_has_equal_paid_R0V_enrichment_rights=True,
                     prior_minimum_feature_status='UNRESOLVED_CONDITIONAL_FIXED_SUPPORT',probabilities=None)
        write(out/'SUMMARY.json',summary)
        write(out/('COMPLETE.json' if failure is None else 'INCOMPLETE.json'),summary)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate',default='results/P1A/reference/GATE.json')
    args=parser.parse_args()
    result=main(args.gate)
    if result['failure'] is not None:raise SystemExit(2)
