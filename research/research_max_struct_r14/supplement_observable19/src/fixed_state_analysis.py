"""Parent-sealed fixed-state OFFLINE_DIAGNOSTIC; no GT and no online changes.

Implemented now, to run only after all144 main predictions are terminal and a
parent-authorized observable19 seal exists. Seven cases x saved calls1/128/512,
using only actual RETURNED internal-Galerkin inputs. Missing slots stay missing.
Each valid common state pays one original Maxwell gradient/state and one nr-RHS
adjoint transport; six original LBFGS8 readouts are diagnostic proposals only.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
import argparse,csv,json,math,sys,time,traceback,warnings
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True

FIXED_CASES=tuple(f'r14_new_{i:04d}' for i in (0,2,9,12,15,20,21))
FIXED_CALLS=(1,128,512)
BASELINE='g19_internal_galerkin'
PUBLIC_KEYS={'G','E0','S','Y','sigma','delta','m_init','lower','upper','grid_shape','coords','centers','sources','receivers','k','radius','cell_radius','components','inverse_side','forward_side'}


def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    if hasattr(value,'tolist'):return clean(value.tolist())
    if hasattr(value,'item') and not hasattr(value,'shape') or getattr(value,'shape',None)==():
        try:return clean(value.item())
        except (ValueError,AttributeError):pass
    if isinstance(value,float) and not math.isfinite(value):return None
    if isinstance(value,Path):return str(value)
    return value


def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(clean(value),f,indent=2,allow_nan=False,ensure_ascii=False);f.write('\n')


def table(out,name,rows):
    rows=clean(rows);write(out/(name+'.json'),rows);keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+'.csv')).open('x',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()})


def stamp(path):
    p=Path(path).resolve();return dict(path=str(p),exists=p.is_file(),bytes=p.stat().st_size if p.is_file() else None)


def delta(after,before):
    return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before) if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float))}


def finite_number(v):
    try:x=float(v);return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None


def norms(x,np):
    if x is None:return None
    return finite_number(np.linalg.norm(x))


def ratio(numerator,denominator):
    if numerator is None or denominator is None:return None,'nonfinite or unavailable operand'
    if denominator==0:return None,'zero denominator'
    result=finite_number(numerator/denominator)
    return result,None if result is not None else 'nonfinite ratio'


def direction(x,y,np):
    nx,ny=norms(x,np),norms(y,np)
    if nx is None or ny is None:return dict(angle_degrees=None,cosine=None,null_reason='unavailable/nonfinite gradient or step')
    if nx==0 or ny==0:return dict(angle_degrees=None,cosine=None,null_reason='zero norm')
    cosine=float(np.clip(np.dot(np.asarray(x).ravel(),np.asarray(y).ravel())/(nx*ny),-1.,1.))
    return dict(angle_degrees=float(np.degrees(np.arccos(cosine))),cosine=cosine,null_reason=None)


def scalar_snapshot(z,key):
    return z[key].item()


def load_returned_snapshot(path,call,np):
    if not path.is_file():return None,'MISSING_BASELINE_READOUT_SNAPSHOT'
    with np.load(path,allow_pickle=False) as z:
        status=str(scalar_snapshot(z,'status')) if 'status' in z else None
        if status!='RETURNED':return None,'BASELINE_SNAPSHOT_NOT_RETURNED:'+str(status)
        required=('J','Q','m_before','m_after','beta','rho','Os','Ls','U','B','readout_call')
        missing=[k for k in required if k not in z]
        if missing:return None,'MISSING_SNAPSHOT_FIELDS:'+','.join(missing)
        state={k:z[k].copy() for k in required}
    if int(state['readout_call'])!=call:raise ValueError('Snapshot cumulative call identity differs')
    m=state['m_before'];J,Q=state['J'],state['Q']
    if np.iscomplexobj(m):raise ValueError('Original material inputs must be real')
    if m.shape!=(576,) or J.shape!=Q.shape or J.ndim!=2 or J.shape[0]!=576:raise ValueError('Unqualified material/current snapshot dimensions')
    if state['Os'].shape!=(J.shape[1],J.shape[1]) or state['Ls'].shape!=state['Os'].shape:raise ValueError('Source metric dimensions differ')
    if not all(np.isfinite(state[k]).all() for k in ('m_before','J','Q','Os','Ls','U','B')):raise ValueError('Nonfinite saved source input')
    return state,None


def missing_arm_row(cid,call,arm,reason,source):
    return dict(case_id=cid,readout_call=call,arm=arm,status='MISSING_OR_UNUSABLE_SOURCE',diagnostic_kind='OFFLINE_DIAGNOSTIC',source=str(source),reason=reason,
        raw_projected_current_energy_ratio=None,transported_current_omission_relative_norm=None,full_gradient_L2=None,
        diagnostic_proposal_Phi=None,diagnostic_Phi_change=None,new_physical_cost=None)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--seal',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start,wall=time.process_time(),time.perf_counter();root=a.root.resolve();legacy=a.legacy_repo.resolve();out=a.out.resolve()
    sys.path[:0]=[str(root/'src'),str(root/'runtime'),str(root.parent/'src'),str(legacy),str(legacy/'vendor')]
    # Stdlib-only authorized seal, exact frozen configuration, complete terminal
    # queue and result inventory checks precede even output directory creation.
    import analyze_results as authorized
    seal,cfg,tasks,account,jobs=authorized.authorize(root,a.seal.resolve())
    if cfg.get('exposure_status')!='EXPOSED24_FOLLOWUP' or cfg.get('basis_rank')!=19:raise ValueError('Frozen exposed fixed19 cohort required')
    if out.exists() or not (out.is_relative_to(root/'audit') or out.is_relative_to(root/'analysis')):raise ValueError('New output directory under supplement audit/ or analysis/ required')
    original={ (t['case_id'],t['method']):t for t in tasks };arms=tuple(authorized.METHODS)
    for cid in FIXED_CASES:
        if any((cid,arm) not in original for arm in arms):raise ValueError('Fixed case/arm absent from authorized manifest')
    prior=cfg['protocol'].get('prior',{})
    if any(float(prior.get(k,float('nan')))!=v for k,v in [('smooth',1e-4),('tv',1e-4),('epsilon',.01)]):raise ValueError('Original C2 prior contract changed')
    if list(cfg['protocol'].get('material_box',[]))!=[.02,20.]:raise ValueError('Original real material box changed')
    out.mkdir(parents=True,exist_ok=False)
    plan=[]
    for cid in FIXED_CASES:
        base_folder=root/'RESULTS'/authorized.CAMPAIGN/original[(cid,BASELINE)]['task_id']
        for call in FIXED_CALLS:
            plan.append(dict(case_id=cid,readout_call=call,source_snapshot=stamp(base_folder/'MATERIAL_READOUT_SNAPSHOTS'/f'readout_{call:05d}.npz'),DPCS_RULE=stamp(base_folder/'DPCS_RULE.json'),DPCS_basis=stamp(base_folder/'RECEIVER_DPCS.npz'),
                frozen_arm_bases=[dict(arm=arm,file=stamp(root/'RESULTS'/authorized.CAMPAIGN/original[(cid,arm)]['task_id']/'DENOMINATOR_BASIS.npz')) for arm in arms]))
    write(out/'AUTHORIZED_FIXED_INPUT_PLAN.json',dict(parent_seal=str(a.seal.resolve()),schema=seal['schema_version'],parent_authorized=seal['authorized'],manifest_snapshot=cfg,fixed_cases=list(FIXED_CASES),fixed_readout_calls=list(FIXED_CALLS),expected_source_slots=21,expected_arm_rows=126,plans=plan,truth_access='NONE',no_missing_state_reruns=True,main_source_or_predictions_modified=False,
        meanings=dict(port_ratio='||O_m(I-P)e||/||O_m e||, raw Frobenius norms; may exceed1, not material information percentage',O_m='S(I-D_mG)^-1 from one retained original Maxwell LU; common scalar CN whitening cancels in relative norms; not material Jacobian',proposals='Original LBFGS8 moving-denominator helper outputs; OFFLINE_DIAGNOSTIC, never main candidates')))
    import numpy as np
    from scipy.linalg import lu_solve,svdvals
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator,matrix_material_readout
    from research.research_max_struct_r11.src.material_readouts import gradient_parts

    class SavedDenominator(LowRankDenominator):
        """No constructor setup: original methods on the sealed U/B factors."""
        def __init__(self,U,B,projection,cost):
            self.U=np.asarray(U,complex).copy();self.B=np.asarray(B,complex).copy();self.n,self.rank=self.U.shape;self.projection=projection;self.cost=cost
            if self.U.shape!=(576,19) or self.B.shape!=self.U.shape:raise ValueError('Saved denominator must be576x19 U/B')
            if not np.isfinite(self.U).all() or not np.isfinite(self.B).all():raise ValueError('Nonfinite sealed U/B')

    basis={};basis_errors={};publics={};public_errors={};dpcs={};dpcs_errors={};sources=[]
    for cid in FIXED_CASES:
        task=original[(cid,BASELINE)];input_path=(root/task['input_file']).resolve()
        if input_path!=root.parent/'data_public/independent24'/(cid+'.npz'):raise ValueError('Private/unexpected input path')
        sources.append(stamp(input_path))
        try:
            with np.load(input_path,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files if k in PUBLIC_KEYS}
            if public['G'].shape!=(576,576) or public['S'].shape[1]!=576 or tuple(public['grid_shape'])!=(24,24):raise ValueError('Original public grid changed')
            if not np.array_equal(public['m_init'],np.full(576,.1)):raise ValueError('Original common public initial material changed')
            if float(public.get('lower',.02))!=.02 or float(public.get('upper',20.))!=20.:raise ValueError('Public box differs from frozen contract')
            publics[cid]=public
        except Exception as exc:public_errors[cid]=repr(exc)
        folder=root/'RESULTS'/authorized.CAMPAIGN/task['task_id']
        try:
            rule=json.loads((folder/'DPCS_RULE.json').read_text());rank=rule['rank']
            if not isinstance(rank,int) or rank<0:raise ValueError('No valid saved DPCS rule rank')
            Vd=None;reason=None
            try:
                with np.load(folder/'RECEIVER_DPCS.npz',allow_pickle=False) as z:Vd=z['Vd'].copy()
                if Vd.shape!=(576,rank):raise ValueError('Saved Vd shape differs from rank')
            except Exception as exc:Vd=None;reason=repr(exc)
            dpcs[cid]=dict(rank=rank,Vd=Vd,basis_null_reason=reason,rule=rule)
        except Exception as exc:dpcs_errors[cid]=repr(exc)
        for arm in arms:
            folder=root/'RESULTS'/authorized.CAMPAIGN/original[(cid,arm)]['task_id'];file=folder/'DENOMINATOR_BASIS.npz';sources.append(stamp(file))
            try:
                setup=json.loads((folder/'DENOMINATOR_SETUP.json').read_text());projection='right' if arm.endswith('_right') else 'galerkin'
                if setup.get('rank')!=19 or setup.get('projection')!=projection:raise ValueError('Saved denominator arm/rank/projection metadata differs')
                with np.load(file,allow_pickle=False) as z:U=z['U'].copy();B=z['B'].copy()
                if U.shape!=(576,19) or B.shape!=U.shape:raise ValueError('Fixed19 factors missing/incompatible')
                if not np.isfinite(U).all() or not np.isfinite(B).all():raise ValueError('Nonfinite sealed U/B')
                basis[(cid,arm)]=dict(U=U,B=B,projection=projection,setup=setup)
            except Exception as exc:basis_errors[(cid,arm)]=repr(exc)

    arm_rows=[];state_rows=[];events=[];failures=[];new_physics=Counter();aux_total=Counter();completed_states=0;interrupted=False;active_event=None
    def event(label,cid,call,arm,before,model,counter,cpu0,wall0,status,error=None):
        after=dict(model.cost) if model is not None else {};physical=delta(after,before)
        new_physics.update(physical);aux_total.update(counter)
        item=dict(label=label,case_id=cid,readout_call=call,arm=arm,status=status,error=error,CPU_seconds=time.process_time()-cpu0,wall_seconds=time.perf_counter()-wall0,new_physical_cost=physical,new_auxiliary_cost=dict(counter))
        events.append(item)
        with (out/'COST_EVENTS.jsonl').open('a',encoding='utf-8') as f:f.write(json.dumps(clean(item),allow_nan=False)+'\n')
        return physical

    try:
        for item in plan:
            cid,call=item['case_id'],item['readout_call'];source=Path(item['source_snapshot']['path']);state_id=f'{cid}__readout_{call:05d}';directory=out/state_id;directory.mkdir()
            try:common,missing=load_returned_snapshot(source,call,np)
            except Exception as exc:common=None;missing='INVALID_SAVED_BASELINE_INPUT:'+repr(exc)
            if cid in public_errors:missing='PUBLIC_INPUT_UNAVAILABLE:'+public_errors[cid];common=None
            if cid in dpcs_errors:missing='SAVED_DPCS_RULE_UNAVAILABLE:'+dpcs_errors[cid];common=None
            if common is None:
                state_rows.append(dict(case_id=cid,readout_call=call,status='MISSING_OR_UNUSABLE_SOURCE',reason=missing,source=str(source),physical_state_evaluated=False))
                arm_rows.extend(missing_arm_row(cid,call,arm,missing,source) for arm in arms);write(directory/'MISSING_SOURCE.json',state_rows[-1]);continue
            if (cid,BASELINE) in basis:
                saved=basis[(cid,BASELINE)]
                if not np.array_equal(common['U'],saved['U']) or not np.array_equal(common['B'],saved['B']):
                    reason='BASELINE_SNAPSHOT_U_B_DIFFERS_FROM_SEALED_DENOMINATOR_FACTORS';state_rows.append(dict(case_id=cid,readout_call=call,status='SOURCE_LINEAGE_MISMATCH',reason=reason));arm_rows.extend(missing_arm_row(cid,call,arm,reason,source) for arm in arms);continue
            public=publics[cid];m=np.asarray(common['m_before'],float);J,Q=common['J'],common['Q'];beta=float(common['beta']);rho=float(common['rho']);Os,Ls=common['Os'],common['Ls'];ns=J.shape[1]
            if public['Y'].shape[1]!=ns:
                reason='SAVED_SOURCE_COUNT_DIFFERS_FROM_PUBLIC_DATA';state_rows.append(dict(case_id=cid,readout_call=call,status='SOURCE_LINEAGE_MISMATCH',reason=reason));arm_rows.extend(missing_arm_row(cid,call,arm,reason,source) for arm in arms);continue
            model=None;state=None;Phi0=None;gphysical=None;Om=None;Oe=None;e=J-m[:,None]*Q;common_arrays=dict(m=m,J=J,Q=Q,e=e,beta=np.asarray(beta),rho=np.asarray(rho),Os=Os,Ls=Ls)
            common_cost=Counter();cpu0,wall0=time.process_time(),time.perf_counter();before={};common_errors=[]
            active_event=dict(label='common state interrupted/unhandled',cid=cid,call=call,arm=None,before=before,model=None,counter=common_cost,cpu0=cpu0,wall0=wall0)
            try:
                model=MaxwellModel(public);active_event['model']=model;common_cost['diagnostic_original_Phi_gradient_attempts']+=1;Phi0,gphysical=model.data_value_gradient(m,public['Y']);state=model._current;common_cost['diagnostic_original_Phi_gradient_completed']+=1;common_cost['diagnostic_original_physical_prior_evaluations']+=1
                if state is None or not np.array_equal(state['m'],m) or state['lu'] is not model._physics.factor:raise RuntimeError('Original common-state LU was not retained')
                common_arrays.update(physical_current=state['j'],physical_Q=state['Q'],physical_Y=state['Y'],physical_Phi_gradient=gphysical,physical_LU_factors=state['lu'][0],physical_LU_pivots=state['lu'][1])
            except Exception as exc:common_errors.append(dict(component='common_original_Phi_gradient',error=repr(exc),traceback=traceback.format_exc()))
            if state is not None:
                try:
                    nr=model.S.shape[0]
                    with model._physics.operation('FIXED_STATE_PORT_TRANSPORT',rhs=nr,adjoint_RHS=nr,diagnostic_port_transport_adjoint_RHS=nr):
                        Om=lu_solve(state['lu'],model.S.conj().T,trans=2,check_finite=False).conj().T
                    Oe=Om@e;common_cost['diagnostic_Om_dense_actions']+=1;common_cost['diagnostic_Om_dense_complex_RHS']+=ns;common_arrays.update(O_m=Om,Oe=Oe)
                except Exception as exc:common_errors.append(dict(component='common_Om',error=repr(exc),traceback=traceback.format_exc()))
            Ge=public['G']@e;common_cost['diagnostic_dense_G_actions']+=1;common_cost['diagnostic_dense_G_complex_RHS']+=ns;common_arrays['G_e']=Ge
            np.savez_compressed(directory/'COMMON_STATE.npz',**common_arrays)
            common_physical=event('common original Maxwell state/gradient and retainedLU port transport',cid,call,None,before,model,common_cost,cpu0,wall0,'PARTIAL_FAILURE' if common_errors else 'COMPLETE',common_errors or None);active_event=None
            state_row=dict(case_id=cid,readout_call=call,status='COMMON_STATE_READY' if not common_errors else 'COMMON_PHYSICS_PARTIAL_FAILURE',source=str(source),beta=beta,rho=rho,sources=ns,physical_state_evaluated=state is not None,physical_Phi_before=finite_number(Phi0),physical_gradient_L2=norms(gphysical,np),raw_auxiliary_current_defect_L2=norms(e,np),O_m_available=Om is not None,physical_cost=common_physical,common_errors=common_errors,DPCS_rank=dpcs[cid]['rank'],DPCS_basis_null_reason=dpcs[cid]['basis_null_reason'])
            state_rows.append(state_row);failures.extend(dict(case_id=cid,readout_call=call,**error) for error in common_errors)
            if model is None:
                arm_rows.extend(missing_arm_row(cid,call,arm,'ORIGINAL_MODEL_SETUP_FAILED',source) for arm in arms);continue
            # Recreate the original C2 material-readout config without any SVD.
            try:
                readout_cfg=solver.arm_config('a_plus',dict(host=.1,material_box=[.02,20.],dpcs_rank=dpcs[cid]['rank']),public['Y'],model.E,[dict(mf=12,maxiter=1500)])
                readout_cfg.update(matrix_inner_maxiter=8,aplus_regularizer=model.reg,initialize_physical_current=True,physical_sigma=float(np.linalg.norm(public['Y'])),fft_geometry=dict(coords=public['coords'],radius=float(public.get('cell_radius',public['radius'])),k=float(public['k'])))
                plain_cfg={k:v for k,v in readout_cfg.items() if not callable(v) and k not in ('aplus_regularizer','matrix_denominator','_counter','_aplus_prior_context')}
                write(directory/'REBUILT_READOUT_CONFIG.json',dict(source='Original R8 benchmark_solver.arm_config(a_plus) and original C2 runtime flags; no current solver invoked',saved_DPCS_rank=dpcs[cid]['rank'],actual_material_beta_from_snapshot=beta,actual_rho_from_snapshot=rho,Os_Ls='saved common snapshot, not recomputed',config=plain_cfg))
            except Exception as exc:
                reason='ORIGINAL_READOUT_CONFIG_REBUILD_FAILED:'+repr(exc);failures.append(dict(case_id=cid,readout_call=call,component='readout_config',error=repr(exc)));arm_rows.extend(missing_arm_row(cid,call,arm,reason,source) for arm in arms);continue
            lo=np.full(576,.02);hi=np.full(576,20.);paid_physical_states=[]
            if state is not None:paid_physical_states.append((m.copy(),state,finite_number(Phi0)))
            for arm in arms:
                if (cid,arm) not in basis:
                    arm_rows.append(missing_arm_row(cid,call,arm,'SEALED_U_B_UNAVAILABLE:'+basis_errors.get((cid,arm),'UNKNOWN'),source));continue
                entry=basis[(cid,arm)];cost=Counter();operator=SavedDenominator(entry['U'],entry['B'],entry['projection'],cost);arrays={};warnings_seen=[];errors=[]
                row=dict(case_id=cid,readout_call=call,arm=arm,status='COMPLETE',diagnostic_kind='OFFLINE_DIAGNOSTIC',source=str(source),rank=19,projection=entry['projection'],beta=beta,rho=rho,
                    physical_Phi_before=finite_number(Phi0),raw_e_L2=norms(e,np),original_DPCS_rank=dpcs[cid]['rank'],T19_min_singular=None,T19_condition=None,
                    material_solver='Original R10 matrix_material_readout moving derivative; original scalar seed/old comparison; maxiter8',main_candidate_pool_modified=False,truth_access='NONE')
                folder=directory/arm;folder.mkdir();cpu0,wall0=time.process_time(),time.perf_counter();before=dict(model.cost)
                active_event=dict(label='arm interrupted/unhandled',cid=cid,call=call,arm=arm,before=before,model=model,counter=cost,cpu0=cpu0,wall0=wall0)
                def reg(z):cost['diagnostic_original_regularizer_evaluations']+=1;return model.reg.evaluate(z)[:2]
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        warnings.simplefilter('always')
                        U=operator.U;Pe=U@(U.conj().T@e);omitted=e-Pe;GPe=public['G']@Pe;Gmissing=public['G']@omitted;output_trim=GPe-U@(U.conj().T@GPe)
                        cost['diagnostic_orthogonal_projection_actions']+=2;cost['diagnostic_orthogonal_projection_complex_RHS']+=2*ns;cost['diagnostic_dense_G_actions']+=2;cost['diagnostic_dense_G_complex_RHS']+=2*ns
                        arrays.update(U=U,B=operator.B,projected_e=Pe,omitted_e=omitted,G_projected_e=GPe,G_input_omitted=Gmissing,PGP_extra_output_trim=output_trim,K_e=operator.action(e))
                        ne,np_e=norms(e,np),norms(Pe,np);ratio1,why=ratio(np_e*np_e if np_e is not None else None,ne*ne if ne is not None else None);row.update(raw_projected_current_energy_ratio=ratio1,raw_projected_energy_null_reason=why)
                        for key,num,den in [('Green_input_omission_relative_norm',norms(Gmissing,np),norms(Ge,np)),('PGP_extra_output_trim_relative_to_GPe',norms(output_trim,np),norms(GPe,np))]:row[key],row[key+'_null_reason']=ratio(num,den)
                        row['Green_input_omission_absolute_norm']=norms(Gmissing,np);row['PGP_extra_output_trim_absolute_norm']=norms(output_trim,np);row['PGP_output_trim_applied_to_this_arm']=entry['projection']=='galerkin'
                        if Om is not None:
                            transported=Om@omitted;cost['diagnostic_Om_dense_actions']+=1;cost['diagnostic_Om_dense_complex_RHS']+=ns;arrays['O_m_omitted_e']=transported
                            row['transported_current_omission_relative_norm'],row['transported_current_omission_null_reason']=ratio(norms(transported,np),norms(Oe,np))
                            row['transported_omission_absolute_norm']=norms(transported,np);row['transported_raw_defect_absolute_norm']=norms(Oe,np)
                        else:row.update(transported_current_omission_relative_norm=None,transported_current_omission_null_reason='common O_m unavailable')
                        Vd=dpcs[cid]['Vd']
                        if Vd is not None:
                            cross=U.conj().T@Vd;overlap=float(np.vdot(cross,cross).real);PDe=Vd@(Vd.conj().T@e);PDPe=Vd@(Vd.conj().T@Pe);PDmissing=Vd@(Vd.conj().T@omitted)
                            arrays.update(DPCS_Vd=Vd,U_adjoint_Vd=cross,DPCS_e=PDe,DPCS_projected_e=PDPe,DPCS_omitted_e=PDmissing)
                            cost['diagnostic_DPCS_projection_actions']+=3;cost['diagnostic_DPCS_projection_complex_RHS']+=3*ns;cost['diagnostic_basis_DPCS_overlap_products']+=1
                            row.update(basis_DPCS_trace_overlap=overlap,basis_fraction_in_DPCS=overlap/19,DPCS_fraction_in_basis=overlap/Vd.shape[1] if Vd.shape[1] else None,DPCS_overlap_null_reason='DPCS rank0' if Vd.shape[1]==0 else None)
                            row['DPCS_omitted_e_relative_norm'],row['DPCS_omitted_e_null_reason']=ratio(norms(PDmissing,np),norms(PDe,np))
                        else:row.update(basis_DPCS_trace_overlap=None,basis_fraction_in_DPCS=None,DPCS_fraction_in_basis=None,DPCS_overlap_null_reason=dpcs[cid]['basis_null_reason'])
                        # Build only the inherited19x19 Woodbury T; this SVD is
                        # diagnostic and never controls whether readout is run.
                        try:
                            z=m/(1+beta*m);T=np.eye(19)-U.conj().T@(z[:,None]*operator.B);cost['diagnostic_T19_assemblies']+=1;cost['diagnostic_T19_SVD_attempts']+=1;singular=svdvals(T,check_finite=False);cost['diagnostic_T19_SVD']+=1
                            arrays.update(T19=T,T19_singular_values=singular);small=float(singular[-1]);row.update(T19_min_singular=small,T19_condition=finite_number(singular[0]/small) if small>0 else None,T19_condition_null_reason='zero minimum singular value' if small==0 else 'nonfinite condition number' if finite_number(singular[0]/small) is None else None)
                        except Exception as exc:errors.append(dict(component='T19_diagnostic_only',error=repr(exc)));row['T19_condition_null_reason']=repr(exc)
                        try:
                            cost['diagnostic_gradient_parts_attempts']+=1;parts=gradient_parts(operator,J,Q,m,beta,Os,rho,reg);cost['diagnostic_gradient_parts_evaluations']+=1;cost['diagnostic_gradient_parts_source_metric_products']+=1
                            for name,key in [('full','full_gradient'),('const','constitutive_gradient'),('metric','metric_gradient'),('prior','prior_gradient')]:
                                g=np.asarray(parts[key]);arrays['g_'+name]=g;row[name+'_gradient_L2']=norms(g,np);row[name+'_gradient_Linf']=finite_number(np.max(np.abs(g)))
                                d=direction(g,gphysical,np) if gphysical is not None else direction(None,None,np)
                                row[name+'_gradient_to_physical_gradient_angle_degrees']=d['angle_degrees'];row[name+'_gradient_to_physical_gradient_cosine']=d['cosine'];row[name+'_gradient_direction_null_reason']=d['null_reason']
                            arrays['transported_auxiliary_defect']=parts['transported_defect'];row['auxiliary_material_value_at_common_m']=finite_number(parts['value'])
                        except Exception as exc:errors.append(dict(component='original_gradient_parts',error=repr(exc),traceback=traceback.format_exc()))
                        try:
                            local=dict(readout_cfg,matrix_denominator=operator,regularizer=reg);cost['diagnostic_readout_source_whitening_dense_products']+=2;cost['diagnostic_LBFGS8_readout_attempts']+=1
                            proposal,readout_meta=matrix_material_readout(None,J.copy(),Q.copy(),m.copy(),beta,Os.copy(),Ls.copy(),rho,(lo.copy(),hi.copy()),reg,local,cost)
                            cost['diagnostic_LBFGS8_readout_completed']+=1;proposal=np.asarray(proposal,float);arrays['diagnostic_proposal_m']=proposal;step=proposal-m;arrays['diagnostic_material_step']=step;row.update(diagnostic_material_step_L2=norms(step,np),readout_metadata=readout_meta)
                            if arm=='g19_internal_galerkin':
                                row['repeated_internal_readout_saved_m_after_Linf_difference']=float(np.max(np.abs(proposal-common['m_after'])))
                                row['repeated_internal_readout_saved_m_after_L2_difference']=float(np.linalg.norm(proposal-common['m_after']))
                                row['repeat_difference_is_diagnostic_not_gate']=True
                            d=direction(step,-gphysical,np) if gphysical is not None else direction(None,None,np);row['diagnostic_step_to_negative_physical_gradient_angle_degrees']=d['angle_degrees'];row['diagnostic_step_direction_null_reason']=d['null_reason']
                            for name in ('full','const','metric','prior'):
                                if 'g_'+name in arrays:row[name+'_gradient_work_on_diagnostic_step']=finite_number(np.dot(arrays['g_'+name],step))
                            try:
                                reused=next(((st,val) for old,st,val in paid_physical_states if np.array_equal(old,proposal)),None)
                                if reused is not None:post,val=reused;cost['diagnostic_reused_paid_physical_states']+=1
                                else:
                                    post=model.state(proposal);residual=post['Y']-public['Y'];cost['diagnostic_original_regularizer_evaluations']+=1;val=.5*float(np.vdot(residual,residual).real/np.vdot(public['Y'],public['Y']).real)+float(model.reg.evaluate(proposal)[0]);paid_physical_states.append((proposal.copy(),post,float(val)))
                                arrays.update(diagnostic_physical_current=post['j'],diagnostic_physical_Q=post['Q'],diagnostic_physical_Y=post['Y'])
                                row['diagnostic_proposal_Phi']=finite_number(val);row['diagnostic_Phi_change']=finite_number(val-Phi0) if val is not None and Phi0 is not None else None
                                row['diagnostic_physical_relative_residual']=float(np.linalg.norm(post['Y']-public['Y'])/np.linalg.norm(public['Y']))
                                row['diagnostic_proposal_legal']=bool(np.isfinite(proposal).all() and proposal.min()>=.02-1e-12 and proposal.max()<=20+1e-12)
                            except Exception as exc:errors.append(dict(component='diagnostic_proposal_original_Phi',error=repr(exc),traceback=traceback.format_exc()))
                        except Exception as exc:errors.append(dict(component='original_LBFGS8_readout',error=repr(exc),traceback=traceback.format_exc()))
                        warnings_seen=[dict(category=w.category.__name__,message=str(w.message)) for w in caught]
                except Exception as exc:errors.append(dict(component='space_diagnostic',error=repr(exc),traceback=traceback.format_exc()))
                except BaseException as exc:
                    errors.append(dict(component='interrupted_arm',error=repr(exc),traceback=traceback.format_exc()));raise
                finally:
                    np.savez_compressed(folder/'DIAGNOSTIC_ARRAYS.npz',**arrays)
                    row['warnings']=warnings_seen;row['errors']=errors;row['status']='INTERRUPTED' if any(e['component']=='interrupted_arm' for e in errors) else 'PARTIAL_DIAGNOSTIC_FAILURE' if errors else 'COMPLETE';row['new_auxiliary_cost']=dict(cost)
                    row['new_physical_cost']=event('same-input space/gradient/original LBFGS8 and proposal Phi',cid,call,arm,before,model,cost,cpu0,wall0,row['status'],errors or None);active_event=None
                    arm_rows.append(row);write(folder/'DIAGNOSTIC.json',row);failures.extend(dict(case_id=cid,readout_call=call,arm=arm,**error) for error in errors)
            completed_states+=1
    except BaseException as exc:
        interrupted=True;failures.append(dict(component='unhandled_or_interrupted_main',error=repr(exc),traceback=traceback.format_exc()))
        raise
    finally:
        if active_event is not None:
            event(**active_event,status='INTERRUPTED_OR_UNHANDLED',error='Partial active action retained; see main failure record')
            active_event=None
        table(out,'COMMON_STATES',state_rows);table(out,'SAME_STATE_SIX_ARM_DIAGNOSTICS',arm_rows);table(out,'FAILURES',failures);table(out,'VALIDATION_COST_EVENTS',events);write(out/'INPUT_SOURCES.json',sources)
        receipt=dict(status='INTERRUPTED_OR_FAILED_RETAINED' if interrupted else 'COMPLETE_WITH_RETAINED_MISSING_OR_FAILURES' if failures or any(r['status']!='COMPLETE' for r in arm_rows) else 'COMPLETE',diagnostic_kind='OFFLINE_DIAGNOSTIC',parent_authorized_seal=str(a.seal.resolve()),same_frozen_config=True,truth_access='NONE',fixed_cases=list(FIXED_CASES),fixed_readout_calls=list(FIXED_CALLS),expected_source_slots=21,actual_common_source_rows=len(state_rows),expected_arm_rows=126,actual_arm_rows=len(arm_rows),arm_status_counts=dict(Counter(r['status'] for r in arm_rows)),common_states_with_physics=sum(r.get('physical_state_evaluated') is True for r in state_rows),new_original_physical_cost=dict(new_physics),new_auxiliary_cost=dict(aux_total),missing_or_failed_components=len(failures),new_basis_SVD=0,new_G_times_U_setup=0,T19_SVD_is_paid_diagnostic_only=True,original_LBFGS8_is_offline_only=True,main_predictions_or_candidate_pool_modified=False,CPU_seconds=time.process_time()-start,whole_process_CPU_seconds_at_receipt=time.process_time(),wall_seconds=time.perf_counter()-wall,CPU_scope='Imports, source reads, original physical state/gradient+nr transport, all lowrank diagnostics/readouts and proposal forwards, copies and output; eventCPU subscopes are not added again. Final receipt/seal serialization follows timestamp.',interpretation='Current defect/projector diagnostics, not material information percentage, material-Jacobian identifiability, recovery predictor, new algorithm or novelty judgment',finished_utc=datetime.now(timezone.utc).isoformat())
        write(out/'RECEIPT.json',receipt)
        write(out/'COUNTER_SCOPE.json',dict(physical='Original model aggregate actual_complex_RHS includes adjoint roles; do not add role counters again. Each common m state/gradient is evaluated once; its retained LU supports exactly one nr-adjoint transport solve. Proposal forward states may reuse already-paid identical m states.',small_system='denominator_small_LU/triangular RHS and diagnostic_T19_SVD are lowrank operations, never added as576-dimensional physicalLU/RHS.',prior='diagnostic_original_regularizer_evaluations records actual reg wrapper calls; inherited material_prior_evaluations are logical helper counts and must not be added as disjoint reg evaluations.',setup='Sealed U/B loaded by thin subclass, no constructor G@U, no basis or DPCS SVD. Existing DPCS_RULE/Vd reused; missing data remain null.',new_proposals='Original full moving denominator LBFGS8, seed/old comparison unchanged; every proposal remains OFFLINE_DIAGNOSTIC and never enters main selection. Condition numbers do not gate proposals.',failures='Inherited small-system counters may count completed calls only; failed attempts retain elapsedCPU/physical counters, warnings and error context. No zero-cost claim for an interrupted operation.'))
        inventory=[stamp(p) for p in sorted(out.rglob('*')) if p.is_file()]
        write(out/'DIAGNOSTIC_OUTPUT_SEAL.json',dict(schema_version='observable19_fixed_state_outputs_v1',source_parent_seal=str(a.seal.resolve()),diagnostic_only=True,authorized_scoring=False,truth_access='NONE',files=inventory,provenance='Existence/bytes only; no hash; same-size replacement undetectable',main_predictions_changed=False))
        text=['# 固定同状态解释性诊断','', '固定7案例、internal_galerkin readout1/128/512实际RETURNED输入；缺失/中断不补跑。所有结果OFFLINE_DIAGNOSTIC，不进入主候选池，不读取GT。', '', '原完整144终态父授权seal与逐项冻结manifest/inventory先核验，之后才创建此新输出。六个K直接从原档U/B载入，无新basis/DPCS SVD或G@U setup。每共同m一次原物理state/gradient，保留LU后一次nr-adjoint RHS构造O_m；各原LBFGS8提案的物理Phi另外付费，不能拿这些validation费用当在线算法成本。', '', '||O_m(I-P)e||/||O_m e||可能大于1，包含相位抵消，不能标为材料信息覆盖百分比。正交current能量、Green输入遗漏、PGP额外输出裁剪与DPCS重合均是current空间描述；没有D_Q材料切向就不作材料Jacobian/可辨识性结论。', '', '小19阶T最小奇异值/条件数只是离线诊断，不据此改变rank、停止或拒绝原readout。所有失败/警告/零分母原因与费用保留。GT若以后用于解释，需父线程另在本诊断输出封存后授权附加；此脚本没有GT读取路径。','',json.dumps(clean(receipt),ensure_ascii=False,indent=2)]
        (out/'REPORT_ZH.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
        print(json.dumps(clean(receipt),indent=2))
    return 0

if __name__=='__main__':raise SystemExit(main())
