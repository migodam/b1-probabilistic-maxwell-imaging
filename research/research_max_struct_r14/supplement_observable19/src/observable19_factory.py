"""Fixed rank19 background bases; inherited R13 PR and R10 full derivative.

Only U and the original LowRankDenominator projection mode vary. No material
DOF, source metric, current family, optimizer or physical forward is replaced.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import importlib.util,json,sys,time
import numpy as np
from scipy.linalg import svd,lu_solve
_frozen_path=Path(__file__).resolve().parents[1]/'runtime/frozen_full_current.py'
_spec=importlib.util.spec_from_file_location('observable19_frozen_full_current',_frozen_path)
original=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=original;_spec.loader.exec_module(original)
from research.research_max_struct_r11.src.material_readouts import prepare_full as original_prepare_full,full_material_readout as original_readout
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
from research.research_max_struct_r10.src import lowrank_denominator as original_denominator

ARMS=('g19_internal_galerkin','g19_observation_galerkin','g19_internal_right','g19_observation_right','g19_anchor_observation_galerkin','g19_anchor_observation_right')
RANK=19
SNAPSHOT_CALLS=(1,128,512,1500)


def write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def svd_basis_headers(G,S):
    _,sg,vg=svd(G,full_matrices=False,check_finite=False)
    _,so,vo=svd(S,full_matrices=False,check_finite=False)
    tol_g=np.finfo(float).eps*max(G.shape)*sg[0]
    tol_o=np.finfo(float).eps*max(S.shape)*so[0]
    return dict(G=vg[:RANK].conj().T,O=vo[:RANK].conj().T),dict(
        singular_G=sg,singular_S=so,
        numerical_G_rank=int(np.count_nonzero(sg>tol_g)),numerical_S_rank=int(np.count_nonzero(so>tol_o)),
        numerical_G_threshold=float(tol_g),numerical_S_threshold=float(tol_o),
        original_G_threshold_0p1_rank=int(np.count_nonzero(sg>=.1*sg[0])))


def spectrum_metadata(values,not_computed):
    if values is None:return dict(singular_values=None,sigma19=None,sigma20=None,relative_gap=None,null_reason=not_computed)
    values=np.asarray(values,float);s19=float(values[18]) if len(values)>=19 else None;s20=float(values[19]) if len(values)>=20 else None
    reason='fewer than20 singular values' if len(values)<20 else 'sigma1 is zero' if values[0]==0 else None
    return dict(singular_values=values.tolist(),sigma19=s19,sigma20=s20,relative_gap=(s19-s20)/float(values[0]) if reason is None else None,null_reason=reason,gap_scope='descriptive only; fixed19 never changed by gap')


def bind_runner(arm,anchor_m0=None):
    if arm not in ARMS:raise ValueError('Unregistered fixed arm: '+str(arm))
    projection='right' if arm.endswith('_right') else 'galerkin'
    selection='anchor_observation' if '_anchor_observation_' in arm else 'observation' if '_observation_' in arm else 'internal'
    context=dict(arm=arm,setup={},readout_calls=0,snapshots=[])

    def prepare(model,config):
        started=time.process_time();out=Path(config['receipt_path']).parent;cost=Counter()
        context['setup']=dict(status='RUNNING',arm=arm,rank=RANK,selection=selection,projection=projection)
        try:
            if model.G.shape!=(576,576) or tuple(model.shape)!=(24,24):raise ValueError('This fixed followup retains the original 576 real material DOF')
            if arm=='g19_internal_galerkin':
                # Exact inherited prepare path, including its original 3 SVDs.
                spectra=[]
                def observed_svd(*args,**kwargs):
                    answer=original_denominator.svd(*args,**kwargs);spectra.append(answer[1].copy());return answer
                bases_observed=original._local_function(original_denominator.prepare_bases,svd=observed_svd)
                prepare_observed=original._local_function(original_denominator._prepare,prepare_bases=bases_observed)
                green_observed=original._local_function(original_denominator.prepare_green_strong,_prepare=prepare_observed)
                full_observed=original._local_function(original_prepare_full,prepare_green_strong=green_observed)
                supplied=full_observed(model,config)
                operator=supplied['matrix_denominator'];old=json.loads((out/'DENOMINATOR_SETUP.json').read_text())
                cost.update(dict(G_SVD=int(old['setup_G_SVD']),S_SVD=int(old['setup_S_SVD']),projected_basis_SVD=int(old['setup_projected_basis_SVD']),Green_basis_and_Galerkin_setup=1,dense_Green_complex_RHS=operator.rank))
                if operator.rank!=RANK:raise ValueError(f'Original G strong(.1) rank is {operator.rank}, not fixed19; no rank adjustment is permitted')
                headers=dict(singular_G=spectra[0],singular_S=spectra[1])
                info=dict(old,original_prepare_full_code_unchanged=True,private_SVD_observers_only=True,original_G_threshold_0p1_rank=operator.rank,basis_rule='Exact original prepare_full G-strong basis; required to equal first19 G right singular vectors')
            else:
                bases,headers=svd_basis_headers(model.G,model.S)
                cost.update(dict(G_SVD=1,S_SVD=1,projected_basis_SVD=0))
                if selection=='anchor_observation':
                    state=getattr(model,'_current',None);physical=getattr(model,'_physics',None)
                    expected=np.full(model.n,.1) if anchor_m0 is None else np.asarray(anchor_m0,float)
                    if not np.array_equal(expected,np.full(model.n,.1)):raise ValueError('Public anchor m0 must be the common original uniform .1')
                    if state is None or physical is None or state is not physical.state or 'lu' not in state or state['lu'] is not physical.factor or not np.array_equal(state['m'],expected):
                        raise RuntimeError('Matching already-paid initial state/LU unavailable; no extra forward or replacement anchor permitted')
                    before_physics=dict(model.cost);nr=model.S.shape[0]
                    with physical.operation('OBSERVABLE19_ANCHOR_TRANSPORT',rhs=nr,adjoint_RHS=nr,observable19_anchor_adjoint_RHS=nr):
                        O0=lu_solve(state['lu'],model.S.conj().T,trans=2,check_finite=False).conj().T
                    _,sa,va=svd(O0,full_matrices=False,check_finite=False)
                    threshold=np.finfo(float).eps*max(O0.shape)*sa[0]
                    bases['A']=va[:RANK].conj().T;headers['numerical_anchor_observation_rank']=int(np.count_nonzero(sa>threshold));headers['numerical_anchor_observation_threshold']=float(threshold);headers['singular_O0']=sa
                    cost['anchor_observation_SVD']+=1
                    delta={k:model.cost.get(k,0)-before_physics.get(k,0) for k in set(model.cost)|set(before_physics)}
                    headers['anchor_physical_cost_delta']=delta
                    headers['anchor_definition']='O0=S(I-D_m0 G)^-1; current constitutive defect to ports, NOT material Jacobian'
                    headers['anchor_cached_initial_LU_reused']=True;headers['anchor_new_LU']=delta.get('LU',0);headers['anchor_adjoint_complex_RHS']=nr
                    np.savez_compressed(out/'ANCHOR_OBSERVATION_OPERATOR.npz',O0=O0,m0=expected,singular_O0=sa)
                    chosen='A';rank_key='numerical_anchor_observation_rank'
                else:
                    chosen='O' if selection=='observation' else 'G';rank_key='numerical_S_rank' if chosen=='O' else 'numerical_G_rank'
                if headers[rank_key]<RANK or bases[chosen].shape[1]!=RANK:raise ValueError(f'{selection} public operator supports rank {headers[rank_key]}, below fixed19; no rank lowering')
                operator=LowRankDenominator(model.G,bases[chosen],projection)
                cost.update(dict(Green_basis_and_Galerkin_setup=int(projection=='galerkin'),Green_right_projection_setup=int(projection=='right'),dense_Green_complex_RHS=RANK))
                supplied=dict(config,matrix_denominator=operator,matrix_inner_maxiter=8)
                info={k:v for k,v in headers.items() if not isinstance(v,np.ndarray)}
                info.update(original_prepare_full_direct=False,basis_rule='first19 right singular vectors of unscaled O0 at public m0=.1' if chosen=='A' else 'first19 right singular vectors of unscaled S; same right subspace as common-scalar W S' if chosen=='O' else 'first19 right singular vectors of G')

            info.update(status='FIXED19_PUBLIC_SETUP_COMPLETE',truth_access='NONE',study='EXPOSED24_FOLLOWUP',arm=arm,
                rank=operator.rank,material_real_DOF=model.n,projection=operator.projection,
                K='P G P = U (U* G U) U*' if projection=='galerkin' else 'G P = (G U) U*',
                U_source='G right singular vectors' if selection=='internal' else 'unscaled O0 right singular vectors; same common-scalar whitened subspace' if selection=='anchor_observation' else 'unscaled S right singular vectors; same common-scalar whitened subspace',
                whitening='Original common scalar CN row whitening leaves right subspace unchanged. Basis SVD uses unscaled S/O0; no inverse-singular-value weighting and no noise sigma input to factory.',
                matrix_inner_maxiter=8,material_derivative='Original moving-denominator matrix_material_readout full analytic derivative; original old/seed comparison',
                initialize_physical_current=True,current_control='restricted nativeFourier/DPCS, original real PR+',
                setup_operation_counts=dict(cost),setup_complete_cpu_seconds=time.process_time()-started,
                full_U_shape=list(operator.U.shape),full_B_shape=list(operator.B.shape),U_B_file='DENOMINATOR_BASIS.npz',
                header_rank_policy='Fixed19; deficient requested operator is reported as failure, never lowered',
                geometry_scope='Auxiliary target metric only; physical model remains original full G')
            info['singular_spectra']={name:spectrum_metadata(headers.get(key),not_computed='O0 is not formed for this arm' if name=='O0' else 'not available') for name,key in [('G','singular_G'),('S','singular_S'),('O0','singular_O0')]}
            np.savez_compressed(out/'DENOMINATOR_SPECTRA.npz',**{key:headers[key] for key in ('singular_G','singular_S','singular_O0') if key in headers})
            if arm!='g19_internal_galerkin':np.savez_compressed(out/'DENOMINATOR_BASIS.npz',U=operator.U,B=operator.B)
            info['setup_complete_cpu_seconds']=time.process_time()-started
            write(out/'DENOMINATOR_SETUP.json',info);context['setup']=info
            return supplied
        except BaseException as exc:
            info=dict(context['setup'],status='SETUP_FAILED_OR_INTERRUPTED',error=repr(exc),setup_operation_counts=dict(cost),setup_complete_cpu_seconds=time.process_time()-started,truth_access='NONE',rank_lowering=False)
            write(out/'DENOMINATOR_SETUP_FAILED.json',info);context['setup']=info
            raise

    def readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
        context['readout_calls']+=1;call=context['readout_calls'];saved=None
        if call in SNAPSHOT_CALLS:
            op=config['matrix_denominator']
            saved=dict(J=J.copy(),Q=Q.copy(),m_before=m.copy(),beta=np.asarray(beta),rho=np.asarray(rho),Os=Os.copy(),Ls=Ls.copy(),U=op.U.copy(),B=op.B.copy(),readout_call=np.asarray(call),arm=np.asarray(arm),projection=np.asarray(op.projection),snapshot_kind=np.asarray('MATERIAL_READOUT_INPUT_OUTPUT_NOT_FULL_PR_CHECKPOINT'))
        try:
            result=original_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost)
        except BaseException as exc:
            if saved is not None:save_readout(saved,'INTERRUPTED',config,cost,error=repr(exc))
            raise
        if saved is not None:
            saved['m_after']=np.asarray(result[0]).copy()
            save_readout(saved,'RETURNED',config,cost,readout_metadata=result[1])
        return result

    def save_readout(arrays,status,config,cost,**extra):
        started=time.process_time();folder=Path(config['receipt_path']).parent/'MATERIAL_READOUT_SNAPSHOTS';folder.mkdir(exist_ok=True)
        call=int(arrays['readout_call']);stem=folder/f'readout_{call:05d}'
        arrays['status']=np.asarray(status);np.savez_compressed(str(stem)+'.npz',**arrays)
        item=dict(readout_call=call,status=status,arm=arm,array_file=str(stem.name)+'.npz',snapshot_kind='Readout inputs and returned material only; NOT a complete PR state or optimizer checkpoint',contains_returned_material='m_after' in arrays,truth_access='NONE',**extra)
        write(str(stem)+'.json',item);context['snapshots'].append(item)
        cost['observable19_readout_snapshot_files']+=1;cost['observable19_readout_snapshot_cpu_microseconds']+=int(1e6*(time.process_time()-started))

    bound=original._local_function(original.run_full_current,prepare_full=prepare,full_material_readout=readout)

    def run(model,Y,m0,config,*args,**kwargs):
        cfg=dict(config);caller=cfg.get('cost_callback')
        def observe_receipt(receipt):
            # R13 hardcodes three preparation-SVD labels; retain that record,
            # then report actual factory operations, with no optimizer change.
            receipt['inherited_denominator_setup_cost_before_factory_accounting']=receipt.get('denominator_setup_cost')
            receipt['denominator_setup_cost']=context['setup'].get('setup_operation_counts',{})
            receipt['observable19_setup']=context['setup'];receipt['observable19_arm']=arm
            receipt['observable19_readout_calls']=context['readout_calls'];receipt['observable19_readout_snapshots']=context['snapshots']
            receipt['observable19_source_scope']='Only private prepare/readout observers; original R13 run_full_current code, R8 PR, R10 full moving derivative, C2 physical corrector unchanged'
            if cfg.get('receipt_path'):write(cfg['receipt_path'],receipt)
            if caller is not None:caller(receipt)
        cfg['cost_callback']=observe_receipt
        return bound(model,Y,m0,cfg,*args,**kwargs)
    run.prepare_for_validation=prepare
    return run,context
