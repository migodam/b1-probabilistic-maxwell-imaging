"""Saved-current invariants and data redundancy; no forward or inverse solve.

Parent-supplied identities are diagnostic hypotheses here. No exact physical
states are generated to recheck them. Corrections minimize physical current
Frobenius norm in each stored stage's free current span, with numerical Gram
rank and unresolved residuals exposed. No material candidate is constructed.
"""
from pathlib import Path
from collections import Counter,defaultdict
from contextlib import contextmanager
from datetime import datetime,timezone
import json,sys,time,traceback
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
R2=ROOT.parent/'research_max_struct_r2'
R3=ROOT.parent/'research_max_struct_r3'
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import fourier_apcs

CASE='obs_b20d5a55677d4eca'
EPS=np.finfo(float).eps


def utc():return datetime.now(timezone.utc).isoformat()


def dump(path,obj):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def ratio(numerator,denominator):
    return float(la.norm(numerator)/max(float(la.norm(denominator)),np.finfo(float).tiny))


class Ledger:
    def __init__(self):
        self.counts=Counter();self.wall=defaultdict(float);self.cpu=defaultdict(float)

    @contextmanager
    def operation(self,name,**counts):
        self.counts[name]+=1;self.counts.update(counts)
        begun,cpu=time.perf_counter(),time.process_time()
        try:yield
        except Exception:
            self.counts['failed_auxiliary_operations']+=1
            raise
        finally:
            self.wall[name]+=time.perf_counter()-begun;self.cpu[name]+=time.process_time()-cpu


def rowspace_distance(M,Vd,ledger):
    with ledger.operation('ROWSPACE_DIAGNOSTIC',small_row_SVD=2,projected_row_products=4):
        residual=M-(M@Vd)@Vd.conj().T
        _,s,vh=la.svd(M,full_matrices=False,check_finite=False)
        tolerance=max(M.shape)*EPS*float(s[0])
        rank=int(np.count_nonzero(s>tolerance))
        Qrow=vh[:rank]
        transverse=Qrow-(Qrow@Vd)@Vd.conj().T
        sins=la.svdvals(transverse,check_finite=False)
        info=dict(raw_relative_Frobenius_distance=ratio(residual,M),row_rank=rank,
                  row_rank_tolerance=tolerance,orthonormal_row_relative_Frobenius_distance=ratio(transverse,Qrow),
                  principal_angle_sines=sins.tolist(),largest_principal_angle_sine=float(sins.max(initial=0)))
    return info,residual,transverse


def defects(E,Gamma,J,ledger):
    with ledger.operation('SAVED_CURRENT_INVARIANTS',source_current_products=2,
                          Gamma_current_products=1,current_Gram_products=1):
        C=E.T@J;reciprocity=C-C.T
        EJ=E.conj().T@J
        lhs=(EJ-EJ.conj().T)/(2j)
        rhs=J.conj().T@(Gamma@J)
        optical=lhs-rhs
        optical_den=float(la.norm(lhs)+la.norm(rhs))
        info=dict(current_fro=float(la.norm(J)),reciprocity_defect_fro=float(la.norm(reciprocity)),
                  reciprocity_relative_to_ETJ=ratio(reciprocity,C),
                  optical_defect_fro=float(la.norm(optical)),
                  optical_relative_to_sum_sides=float(la.norm(optical)/max(optical_den,np.finfo(float).tiny)),
                  optical_lhs_fro=float(la.norm(lhs)),optical_rhs_fro=float(la.norm(rhs)))
    return info,reciprocity,optical


class MinimumReciprocityCorrection:
    def __init__(self,E,U,family,ledger):
        self.E,self.U,self.family,self.ledger=E,U,family,ledger
        p=E.shape[1]
        self.pairs=np.asarray([(i,j)for i in range(p)for j in range(i+1,p)],int)
        with ledger.operation('RECIPROCITY_CONSTRAINT_SETUP',constraint_Gram_builds=1,
                              constraint_Gram_SVD=1,constraint_dimension=len(self.pairs),
                              E_basis_products=1,source_Gram_products=1):
            self.B=E.T@U
            K=self.B@self.B.conj().T
            ii,jj=self.pairs.T
            # M(DeltaA)_ij=(B DeltaA)_ij-(B DeltaA)_ji, i<j.
            Gram=(K[ii[:,None],ii[None,:]]*(jj[:,None]==jj[None,:])
                  -K[ii[:,None],jj[None,:]]*(jj[:,None]==ii[None,:])
                  -K[jj[:,None],ii[None,:]]*(ii[:,None]==jj[None,:])
                  +K[jj[:,None],jj[None,:]]*(ii[:,None]==ii[None,:]))
            self.Gram=Gram
            self.left,self.s,self.vh=la.svd(Gram,full_matrices=False,check_finite=False)
            self.tolerance=len(self.pairs)*EPS*float(self.s[0])
            self.keep=self.s>self.tolerance
            self.info=dict(family=family,orthonormal_spatial_rank=U.shape[1],
                           real_current_dof=2*U.shape[1]*p,complex_constraints=len(self.pairs),
                           constraint_Gram_numerical_rank=int(self.keep.sum()),
                           constraint_Gram_singular_values=self.s.tolist(),
                           constraint_Gram_rank_tolerance=self.tolerance,
                           Gram_rank_rule='66*machine_epsilon*largest_Gram_singular_value; fixed numerical rule, no sweep',
                           source_basis_action_fro=float(la.norm(self.B)),
                           source_basis_action_relative_to_ET=float(la.norm(self.B)/la.norm(E.T)),
                           orthonormal_basis_defect=float(la.norm(U.conj().T@U-np.eye(U.shape[1]))),
                           minimum_scope='minimum current Frobenius norm in retained numerical constraint-Gram range; no exact-rank certificate')

    def correction(self,J,defect):
        with self.ledger.operation('RECIPROCITY_MINIMUM_CORRECTION',constraint_small_RHS=1,
                                   constraint_adjoint_actions=1,physical_current_basis_products=1):
            ii,jj=self.pairs.T
            target=-defect[ii,jj]
            lam=self.vh[self.keep].conj().T@((self.left[:,self.keep].conj().T@target)/self.s[self.keep])
            deltaA=np.zeros((self.U.shape[1],self.E.shape[1]),complex)
            for (i,j),weight in zip(self.pairs,lam):
                deltaA[:,j]+=self.B[i].conj()*weight
                deltaA[:,i]-=self.B[j].conj()*weight
            deltaJ=self.U@deltaA
            remaining=target-self.Gram@lam
            info=dict(correction_current_fro=float(la.norm(deltaJ)),
                      correction_relative_to_stored_current=ratio(deltaJ,J),
                      correction_coefficient_fro=float(la.norm(deltaA)),
                      unresolved_upper_constraint_fro=float(la.norm(remaining)),
                      unresolved_upper_constraint_relative=ratio(remaining,target),
                      target_upper_constraint_fro=float(la.norm(target)))
        return deltaJ,deltaA,info


def run():
    started,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/current_invariants';out.mkdir(parents=True,exist_ok=True)
    path=out/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve single diagnostic disposition; no tuning')
    ledger=Ledger();rows=[];arrays={}
    receipt=dict(status='RUNNING',started_utc=utc(),case_id=CASE,truth_access='NONE',
                 classification='SAVED_AUXILIARY_CURRENT_DIAGNOSTIC_NO_PHYSICS',
                 parent_supplied_identities=dict(reciprocity='E.T J is symmetric',
                    optical='Im_H(E.conj().T J)=J.conj().T Gamma J; Gamma=(G-G.conj().T)/(2j)'),
                 configuration=dict(dpcs_rtol=1e-3,Fourier_stages=[3,6,10],
                    constraint_Gram_shape=[66,66],numerical_rank='dimension*eps*smax',
                    objective='minimum Frobenius norm of physical-coordinate deltaJ within current family',
                    source_fields=['G','E0','S','Y','grid_shape'],stored_fields=['J'],threads=1),
                 new_physical_LU=0,new_complex_RHS=0,inverse_runs=0,hash_checks=False,
                 limits='No success tolerance, optical correction, material readout, physical solve, or theorem validity decision')
    dump(path,receipt)
    try:
        with np.load(R2/'data_public'/f'{CASE}.npz',allow_pickle=False)as f:
            G,E,S,Y,shape=[f[k].copy()for k in ('G','E0','S','Y','grid_shape')]
        n,p=E.shape
        with ledger.operation('RECEIVER_SVD_SETUP',receiver_SVD=1,DPCS_current_products=2):
            Us,s,Vh=la.svd(S,full_matrices=False,check_finite=False)
            rank=int(np.count_nonzero(s>1e-3*s[0]));Vd=Vh[:rank].conj().T
            JD=Vd@((Us[:,:rank].conj().T@Y)/s[:rank,None])
            arrays.update(receiver_singular_values=s,Vd=Vd,JD=JD)
        rT,tres,torth=rowspace_distance(E.T,Vd,ledger)
        rH,hres,horth=rowspace_distance(E.conj().T,Vd,ledger)
        arrays.update(ET_transverse=tres,EH_transverse=hres,ET_orthogonal_row_transverse=torth,EH_orthogonal_row_transverse=horth)
        with ledger.operation('GAMMA_DATA_REDUNDANCY',Gamma_projection_products=4,Gamma_left_projection_products=2):
            Gamma=(G-G.conj().T)/(2j)
            projected=Vd@(Vd.conj().T@Gamma@Vd)@Vd.conj().T
            left=Gamma-Vd@(Vd.conj().T@Gamma)
            gamma_info=dict(full_fro=float(la.norm(Gamma)),projected_fro=float(la.norm(projected)),
                            leakage_relative_to_full=ratio(Gamma-projected,Gamma),
                            one_sided_leakage_relative_to_full=ratio(left,Gamma))
            arrays.update(Gamma=Gamma,Gamma_data_projection=projected)
        receipt['data_redundancy']=dict(receiver_retained_rank=rank,receiver_singular_cutoff=float(1e-3*s[0]),
                                      ET_rowspace=rT,EH_rowspace=rH,Gamma=gamma_info,
                                      interpretation='Distances only; small values do not supply a success threshold or certify exact rowspace inclusion')
        bases={};solvers={};basis_info={}
        for mf in [3,6,10]:
            with ledger.operation('NATIVE_APCS_BASIS',orthonormal_APCS_SVD=1):
                Z,meta=fourier_apcs(shape,mf,Vd,ledger.counts)
                u,sz,_=la.svd(Z,full_matrices=False,check_finite=False)
                tol=max(Z.shape)*EPS*max(float(sz[0]),1.)
                U=u[:,sz>tol]
                family=f'native_mf{mf}'
                bases[family]=U;arrays[family+'_U']=U
                basis_info[family]=dict(**meta,orthonormal_rank=U.shape[1],rank_tolerance=tol,
                                       DPCS_overlap_fro=float(la.norm(Vd.conj().T@U)))
            solvers[family]=MinimumReciprocityCorrection(E,U,family,ledger)
        with ledger.operation('RELEASED_CURRENT_BASIS',released_basis_SVD=1):
            Z=np.column_stack((Vd,bases['native_mf10']))
            u,sr,_=la.svd(Z,full_matrices=False,check_finite=False)
            tol=max(Z.shape)*EPS*float(sr[0]);U=u[:,sr>tol]
            family='released_DPCS_mf10';bases[family]=U;arrays[family+'_U']=U
            basis_info[family]=dict(orthonormal_rank=U.shape[1],rank_tolerance=tol,
                                   data_locked=False,scope='actual R3 final-stage free Vd plus APCS family')
        solvers[family]=MinimumReciprocityCorrection(E,U,family,ledger)
        receipt['basis_info']=basis_info
        receipt['constraint_families']={key:solver.info for key,solver in solvers.items()}
        for key,solver in solvers.items():
            arrays[key+'_constraint_Gram']=solver.Gram
            arrays[key+'_constraint_singular_values']=solver.s
        streams=[('R2_native',R2/'results/blind'/CASE/'native'),
                 ('R2_global',R2/'results/blind'/CASE/'global'),
                 ('R3_noisy_strong_release',R3/'results/development/noisy_strong_release')]
        for stream,directory in streams:
            for stage,mf in enumerate([3,6,10]):
                source=directory/f'stage_{stage}.npz'
                with np.load(source,allow_pickle=False)as f:J=f['J'].copy()
                if J.shape!=(n,p)or not np.isfinite(J).all():raise ValueError('Invalid saved current '+str(source))
                family='released_DPCS_mf10'if stream=='R3_noisy_strong_release'and stage==2 else f'native_mf{mf}'
                solver=solvers[family];U=bases[family]
                before,reciprocity,optical=defects(E,Gamma,J,ledger)
                with ledger.operation('STORED_CURRENT_FAMILY_AUDIT',current_projection_products=4,data_current_products=1):
                    centered=J-JD
                    membership=centered-U@(U.conj().T@centered)
                    dpcs_change=Vd.conj().T@centered
                    data_residual=S@J-Y
                dJ,dA,corr=solver.correction(J,reciprocity)
                after,r_after,o_after=defects(E,Gamma,J+dJ,ledger)
                with ledger.operation('CORRECTION_VISIBILITY_AUDIT',correction_data_products=1,DPCS_correction_products=1):
                    data_change=S@dJ
                    locked_change=Vd.conj().T@dJ
                    correction_visible=dict(full_data_change_relative_to_Y=ratio(data_change,Y),
                         corrected_auxiliary_data_relative_residual=ratio(data_residual+data_change,Y),
                         DPCS_correction_fro=float(la.norm(locked_change)),
                         DPCS_correction_relative_to_current=ratio(locked_change,J))
                key=f'{stream}_stage{stage}'
                arrays.update({key+'_deltaJ':dJ,key+'_deltaA':dA,key+'_reciprocity_before':reciprocity,
                               key+'_reciprocity_after':r_after,key+'_optical_before':optical,key+'_optical_after':o_after})
                row=dict(stream=stream,stage=stage,mf=mf,family=family,source=str(source.relative_to(REPO)),
                         data_locked=family!='released_DPCS_mf10',before=before,after_reciprocity_correction=after,
                         stored_affine_family_relative_defect=ratio(membership,J),
                         stored_DPCS_offset_relative_defect=ratio(dpcs_change,J),
                         stored_auxiliary_data_relative_residual=ratio(data_residual,Y),
                         minimum_reciprocity_correction=corr,correction_visibility=correction_visible)
                rows.append(row)
                receipt.update(completed_currents=len(rows),rows=rows,auxiliary_counts=dict(ledger.counts))
                dump(path,receipt)
                print(json.dumps(dict(current=key,reciprocity=before['reciprocity_relative_to_ETJ'],
                     optical=before['optical_relative_to_sum_sides'],
                     correction_ratio=corr['correction_relative_to_stored_current'],
                     remaining=corr['unresolved_upper_constraint_relative'])),flush=True)
        np.savez_compressed(out/'ARRAYS.npz',**arrays)
        receipt.update(status='COMPLETED_SAVED_CURRENT_DIAGNOSTIC',finished_utc=utc(),rows=rows,
                       completed_currents=len(rows),arrays='ARRAYS.npz',
                       retained_Gram_rank_caveat='Gram squares constraint conditioning; pseudoinverse discards singular values below recorded floating cutoff. Remaining defect reports the unresolved component; no numerical-rank sweep.')
        lines=['# Saved-current collective invariants','',
               'No GT, inverse run, or new physical state. Parent-supplied reciprocal/real-material identities are evaluated only on saved auxiliary currents. No success threshold or theorem verification is introduced.',
               '',f'Retained receiver rank at relative1e-3: {rank}.',
               f"E.T rowspace raw relative distance={rT['raw_relative_Frobenius_distance']:.9g}; orthonormal-row distance={rT['orthonormal_row_relative_Frobenius_distance']:.9g}; largest principal-angle sine={rT['largest_principal_angle_sine']:.9g}.",
               f"E.conj().T rowspace raw relative distance={rH['raw_relative_Frobenius_distance']:.9g}; orthonormal-row distance={rH['orthonormal_row_relative_Frobenius_distance']:.9g}; largest principal-angle sine={rH['largest_principal_angle_sine']:.9g}.",
               f"Gamma leakage outside P Gamma P, relative Frobenius={gamma_info['leakage_relative_to_full']:.9g}.",
               '', '| Stored current | Free family | Reciprocity defect / ||E.T J|| | Optical defect / sum sides | Minimum ||deltaJ||/||J|| | Remaining reciprocity constraint fraction | Data change / ||Y|| |',
               '|---|---|---:|---:|---:|---:|---:|']
        for row in rows:
            b=row['before'];c=row['minimum_reciprocity_correction'];v=row['correction_visibility']
            lines.append(f"|{row['stream']} stage{row['stage']}|{row['family']}|{b['reciprocity_relative_to_ETJ']:.6g}|{b['optical_relative_to_sum_sides']:.6g}|{c['correction_relative_to_stored_current']:.6g}|{c['unresolved_upper_constraint_relative']:.6g}|{v['full_data_change_relative_to_Y']:.6g}|")
        lines+=['', 'Each correction uses an orthonormal physical-current basis and the66-complex-constraint Gram pseudoinverse. Its reported minimum applies to the retained numerical Gram range; discarded constraint components and actual post-correction defects are preserved. Basis orthonormalization changes only this diagnostic norm representation, not the historical CIE optimization.',
                '', 'Native APCS corrections preserve the locked DPCS coordinates to floating arithmetic, while weak/unretained receiver channels can change. R3 final released-stage corrections use its actual free[Vd,APCS] family and therefore can alter locked data components. These controls are not interchangeable.',
                '', 'Small rowspace/Gamma leakage means these quantities are largely visible in the retained measurement subspace for this acquisition; the measured minimum correction describes the cost of using the remaining freedom. Neither fact alone establishes feasibility of the full real-material state equation or an improved reconstruction. Parent owns the algorithm decision.',
                '', 'All detailed norms, optical changes after reciprocity correction, numerical ranks, SVDs, Gram solves and timings are in RECEIPT.json; matrices and corrections are in ARRAYS.npz. New physical LU/RHS and inverse runs are all zero. No source/current coefficient tuning or material scoring was performed.']
        (out/'NOTE.md').write_text('\n'.join(lines)+'\n')
    except Exception as exc:
        receipt.update(status='FAILED_PAID_DIAGNOSTIC',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        receipt.update(auxiliary_counts=dict(ledger.counts),operation_wall_seconds=dict(ledger.wall),
                       operation_cpu_seconds=dict(ledger.cpu),wall_seconds=time.perf_counter()-started,
                       cpu_seconds=time.process_time()-cpu,
                       timer_scope='public/storedcurrent reads,basis/SVD/Gram/constraint work,arrays/note writes;imports/finalreceiptwrite excluded')
        dump(path,receipt)


if __name__=='__main__':
    with threadpool_limits(limits=1):run()
