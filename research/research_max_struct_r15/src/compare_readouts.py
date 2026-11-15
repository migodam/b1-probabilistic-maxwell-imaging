"""Fixed public-state material-readout counterfactual; implement/run separately.

Four public IDs x four original LBFGS8 readouts. No current/PR updates, main
candidate export, corrector, GT, online selection, rank search or new SVD.
"""
from __future__ import annotations
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,csv,importlib.util,json,math,sys,time,traceback
from pathlib import Path
from collections import Counter
CASES=tuple(f'r14_new_{i:04d}' for i in (9,20,21,23))
_ACTIVE_OUT=None
_START_WALL=time.perf_counter()
ARMS=(('G19_nohop','basis_G19',False),('S19_nohop','basis_S19',False),('U11_nohop','U11',False),('U11_onehop','U11',True))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path);p.add_argument('--out',type=Path);p.add_argument('--geometry-directions',action='store_true');a=p.parse_args()
    root=a.root.resolve();legacy=(a.legacy_repo or root.parent.parent).resolve();out=(a.out or root/'audit/readout_counterfactual_v1').resolve();wall=time.perf_counter();cpu=time.process_time()
    if out.exists() or not out.is_relative_to(root/'audit'):raise ValueError('New R15 audit directory required')
    global _ACTIVE_OUT
    out.mkdir(parents=True);_ACTIVE_OUT=out;sys.path[:0]=[str(root/'src'),str(root/'runtime'),str(legacy),str(legacy/'vendor')]
    import numpy as np
    import dualmetric_core as codec
    from queue_interventions import source_checkpoint,geometry_map
    from twofold_denominator import TwofoldDenominator
    from research.research_max_struct_r10.src.lowrank_denominator import matrix_material_readout
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    def encode(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    def write(path,value):
        with Path(path).open('x') as f:json.dump(value,f,default=encode,indent=2,allow_nan=False);f.write('\n')
    def changes(after,before):return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before) if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float))}
    def norm(x):return float(np.linalg.norm(x))
    def projected(g,m,lo,hi):
        residual=m-np.clip(m-g,lo,hi);kkt=g.copy();blocked=((m<=lo+1e-8)&(g>0))|((m>=hi-1e-8)&(g<0));kkt[blocked]=0
        return dict(projected_descent_residual_L2=norm(residual),projected_descent_residual_Linf=float(np.max(abs(residual))),
            KKT_projected_gradient_L2=norm(kkt),KKT_projected_gradient_Linf=float(np.max(abs(kkt))),blocked_coordinates=int(np.count_nonzero(blocked))),residual,kkt
    def box(m,lo,hi):return dict(min=float(m.min()),max=float(m.max()),lower_active=int(np.count_nonzero(m<=lo+1e-8)),upper_active=int(np.count_nonzero(m>=hi-1e-8)),legal=bool(np.isfinite(m).all() and np.all(m>=lo-1e-12) and np.all(m<=hi+1e-12)))
    rows=[];cases=[];issues=[]
    try:geometry=geometry_map(root)
    except Exception as exc:geometry={};issues.append(dict(scope='optional_real_geometry_metadata',error=repr(exc)))
    for cid in CASES:
        c=time.process_time();w=time.perf_counter();folder=out/cid;folder.mkdir();source=None;state=None;model=None;origin=None;basis=None;initial_phi=None;grad=None;case_error=None;case_initial_cost=None
        try:
            source=source_checkpoint(root,cid)
            if source['path'] is None:raise FileNotFoundError('No complete source checkpoint under frozen512 availability policy')
            state=codec._load_state_npz(source['path']);m=np.asarray(state['m'],float).copy();J=state['J'].copy();Q=state['Q'].copy();Os=state['Os'].copy();Ls=state['Ls'].copy()
            beta=float(state['beta1']);rho=float(state['frozen_config'].get('rho',1.));cfg=dict(state['frozen_config'])
            fixed=Path(source['path']).parent/state['r15']['fixed_cache_relative']
            with np.load(fixed,allow_pickle=False) as z:basis={key:z[key].copy() for key in ('basis_G19','basis_S19','U11')}
            public_path=root/'data_public/independent24'/(cid+'.npz')
            with np.load(public_path,allow_pickle=False) as z:public={key:z[key].copy() for key in z.files}
            model=MaxwellModel(public);Y=public['Y'];lo=np.broadcast_to(public['lower'],m.shape).copy();hi=np.broadcast_to(public['upper'],m.shape).copy()
            # Exactly one baseline physical forward/retained LU and one adjoint.
            # The unchanged API includes the original physical material prior.
            initial_phi,grad=model.data_value_gradient(m,Y);origin=model.state(m);case_initial_cost=dict(model.cost)
            if not np.isfinite(initial_phi) or not np.isfinite(grad).all():raise FloatingPointError('Nonfinite shared physical Phi/gradient')
            physical_pg,pgr,kg=projected(grad,m,lo,hi)
            np.savez_compressed(folder/'SOURCE_INPUT.npz',m=m,J=J,Q=Q,Os=Os,Ls=Ls,beta=np.asarray(beta),rho=np.asarray(rho),physical_gradient=grad,
                physical_projected_descent_residual=pgr,physical_KKT_projected_gradient=kg,initial_physical_Phi=np.asarray(initial_phi))
            setup=dict(case_id=cid,status='READY',source_selection=source,source_fixed_cache=str(fixed),source_stage=state['stage'],source_stage_iteration=state['iteration'],source_global_iteration=state['iterations'],
                beta_current=state['beta0'],beta_material=beta,rho=rho,initial_physical_Phi=initial_phi,initial_box=box(m,lo,hi),physical_gradient_norm=norm(grad),physical_projected_gradient=physical_pg,
                initial_physical_cost=case_initial_cost,checkpoint_inherited_auxiliary=state['operation_counts'],checkpoint_inherited_projection=state['r15']['projector_cost'],
                regularizer_source='Unchanged MaxwellModel(public).reg; same original callback as main runner',actual_prior={k:getattr(model.reg,k,None) for k in ('smooth','tv','epsilon')},
                geometry_metadata=geometry.get(cid),truth_access='NONE')
            write(folder/'SOURCE_RECEIPT.json',setup)
        except Exception as exc:case_error=repr(exc);issues.append(dict(case_id=cid,scope='case_initial_state_or_gradient',error=case_error))
        for arm,key,one_hop in ARMS:
            ac=time.process_time();aw=time.perf_counter();counter=Counter();operator=None;arrays_out={};physical_before=dict(model.cost) if model else None
            row=dict(case_id=cid,readout=arm,basis_key=key,one_hop=one_hop,status='RUNNING',source_selection=source,
                truth_access='NONE',online_selection=False,current_or_PR_changed=False,original_LBFGS_maxiter=8,no_new_SVD=True,
                inherited_auxiliary=state['operation_counts'] if state else None,inherited_projection=state['r15']['projector_cost'] if state else None)
            try:
                if case_error:raise RuntimeError('Case prerequisite failed: '+case_error)
                operator=TwofoldDenominator(model.G,basis[key],tuple(model.shape),one_hop=one_hop,cost=counter)
                reg=lambda x:model.reg.evaluate(x)[:2]
                local=dict(cfg,matrix_denominator=operator,matrix_inner_maxiter=8,regularizer=reg,aplus_regularizer=model.reg)
                # Diagnostic objective/gradients are separately paid; original
                # seed, old/seed acceptance, LBFGS8, box and gtol are untouched.
                aux_old,g_old=operator.value_gradient(J,Q,m,beta,Os,rho,reg)
                if not np.isfinite(aux_old) or not np.isfinite(g_old).all():raise FloatingPointError('Nonfinite old auxiliary objective/gradient')
                returned,metadata=matrix_material_readout(operator.blocks,J.copy(),Q.copy(),m.copy(),beta,Os,Ls,rho,(lo,hi),reg,local,counter)
                newm=np.asarray(returned,float);aux_new,g_new=operator.value_gradient(J,Q,newm,beta,Os,rho,reg)
                if not np.isfinite(aux_new) or not np.isfinite(g_new).all():raise FloatingPointError('Nonfinite new auxiliary objective/gradient')
                if not box(newm,lo,hi)['legal']:raise ValueError('Readout returned invalid/out-of-box material')
                dm=newm-m;linear=float(np.dot(grad,dm));dot_negative=-linear;denom=norm(dm)*norm(grad)
                cosine=float(np.clip(dot_negative/denom,-1.,1.)) if denom>0 else None
                old_pg,old_pgr,old_kg=projected(g_old,m,lo,hi);new_pg,new_pgr,new_kg=projected(g_new,newm,lo,hi)
                arrays_out=dict(m_new=newm,delta_m=dm,auxiliary_gradient_old=g_old,auxiliary_gradient_new=g_new,
                    auxiliary_projected_descent_residual_old=old_pgr,auxiliary_projected_descent_residual_new=new_pgr,auxiliary_KKT_projected_gradient_old=old_kg,auxiliary_KKT_projected_gradient_new=new_kg)
                # New material: forward and prior ONLY; no physical gradient,
                # correction, current optimization or candidate-pool insertion.
                newstate=model.state(newm);res=newstate['Y']-Y;D=.5*float(np.vdot(res,res).real)/float(np.vdot(Y,Y).real);R=float(model.reg.evaluate(newm)[0]);actual_phi=D+R
                if not np.isfinite(actual_phi):raise FloatingPointError('Nonfinite new physical Phi')
                row.update(status='COMPLETE',source_stage=state['stage'],source_stage_iteration=state['iteration'],source_global_iteration=state['iterations'],beta_current=state['beta0'],beta_material=beta,
                    selected_basis_rank=operator.rank,operator_setup=operator.setup_receipt,initial_physical_Phi=initial_phi,new_physical_Phi=actual_phi,new_physical_D=D,new_physical_R=R,
                    delta_m_norm=norm(dm),delta_dot_negative_physical_gradient=dot_negative,cosine_to_negative_physical_gradient=cosine,
                    angle_to_negative_physical_gradient_degrees=float(np.degrees(np.arccos(cosine))) if cosine is not None else None,
                    angle_null_reason=None if cosine is not None else 'ZERO_MATERIAL_STEP_OR_PHYSICAL_GRADIENT',linear_predicted_delta_Phi=linear,linear_predicted_Phi=initial_phi+linear,
                    actual_delta_Phi=actual_phi-initial_phi,actual_minus_linear_delta_Phi=actual_phi-initial_phi-linear,old_auxiliary_objective=aux_old,new_auxiliary_objective=aux_new,
                    auxiliary_delta=aux_new-aux_old,old_auxiliary_projected_gradient=old_pg,new_auxiliary_projected_gradient=new_pg,old_box=box(m,lo,hi),new_box=box(newm,lo,hi),original_readout_metadata=metadata,
                    diagnostic_auxiliary_value_gradient_evaluations=2,diagnostic_scope='Two explicit objective/gradient diagnostics outside original readout, charged separately in same fresh counter; no physical gradient at new material')
                arrays_out.update(physical_j_new=newstate['j'],physical_Q_new=newstate['Q'],physical_Y_prediction_new=newstate['Y'])
                if a.geometry_directions and cid in geometry:
                    try:
                        with np.load(geometry[cid]['path'],allow_pickle=False) as z:frames={k:z[k].copy() for k in ('U11','U10','U01','U00')}
                        e=float(np.vdot(dm,dm).real);row['real_material_direction_fractions']={k:float(np.linalg.norm(U.conj().T@dm)**2)/e if e>0 else None for k,U in frames.items()};row['real_direction_null_reason']=None if e>0 else 'ZERO_MATERIAL_STEP'
                    except Exception as exc:row['real_direction_null_reason']=repr(exc);issues.append(dict(case_id=cid,readout=arm,scope='optional_real_geometry_direction',error=repr(exc)))
            except Exception as exc:
                row.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc());issues.append(dict(case_id=cid,readout=arm,error=repr(exc)))
            finally:
                row.update(fresh_auxiliary_cost=dict(counter),fresh_physical_cost=changes(dict(model.cost),physical_before) if model and physical_before is not None else None,
                    CPU_seconds=time.process_time()-ac,wall_seconds=time.perf_counter()-aw,cost_scope='Fresh zero-start readout setup/inner/diagnostic counters; checkpoint inherited counts separate, never added. Physical per-arm model delta only; shared baseline gradient billed once per case.')
                if arrays_out:np.savez_compressed(folder/(arm+'.npz'),**arrays_out)
                write(folder/(arm+'_RECEIPT.json'),row);rows.append(row)
        case=dict(case_id=cid,status='FAILED_INITIAL_STATE' if case_error else 'COMPLETE',initial_error=case_error,source_selection=source,
            initial_physical_cost=case_initial_cost,case_physical_total=dict(model.cost) if model else None,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,
            arm_rows=[r['readout'] for r in rows if r['case_id']==cid],scope='Case CPU includes source load, one baseline physical gradient, all four readouts and writes; nested arm CPU is not added again',truth_access='NONE')
        write(folder/'CASE_RECEIPT.json',case);cases.append(case)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/'COMPARISON.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows({k:json.dumps(v,default=encode) if isinstance(v,(list,dict,np.ndarray)) else v for k,v in r.items()} for r in rows)
    write(out/'COMPARISON.json',rows);write(out/'CASES.json',cases);write(out/'ISSUES.json',issues)
    receipt=dict(status='COMPLETE_WITH_FAILURES' if any(r['status']=='FAILED' for r in rows) else 'COMPLETE',expected_cases=4,expected_readouts=16,completed_readouts=sum(r['status']=='COMPLETE' for r in rows),
        source_policy='queue_interventions.source_checkpoint: requested512, latest complete numeric<=512 else earliest complete; metadata-only availability selection',
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,whole_process_CPU_seconds_at_receipt=time.process_time(),truth_access='NONE',current_PR_updates=0,main_candidate_pool_writes=0,corrector_calls=0,online_selection=False,
        physical_LU_total=sum(c['case_physical_total'].get('LU',0) for c in cases if c['case_physical_total']),actual_complex_RHS_total=sum(c['case_physical_total'].get('actual_complex_RHS',0) for c in cases if c['case_physical_total']),
        cost_scope='One case physical gradient (retained initial LU+adjoint); four new-material forwards+prior only (exact cache hits retained). Case/arm/internal timers nested, do not sum; failed work charged. Auxiliary setup does not perform a new SVD.',
        projected_gradient_definitions='descent residual=m-clip(m-g,lo,hi); KKT projected gradient removes blocked descent at box-active tolerance1e-8. These are fixed diagnostics, not changed optimizer stopping rules.',issues=len(issues))
    write(out/'RECEIPT.json',receipt)
    (out/'README.md').write_text('# Fixed-state material-readout counterfactual\n\n'+json.dumps(receipt,indent=2)+'\n\nEach arm starts from the same frozen full m,J,Q,beta,Os,Ls and original prior. Only one inherited matrix_material_readout is applied per arm, with original seed/bounds/LBFGS8/acceptance. No current/PR state is changed or continued. No new result is inserted into a main candidate pool. No GT is read.\n\nPhysical-gradient alignment and linear prediction use the single original full physical Phi gradient at the shared source. New-material Phi uses only a full physical forward and prior, no adjoint/corrector. Nonlinear finite changes are not guarantees or online selectors. Parent owns interpretation.\n\nCheckpoint inherited auxiliary/projector counters are preserved as provenance; every counterfactual operator starts a fresh Counter. Initial physical-gradient costs are shared once per case, per-arm physics deltas and full case totals are distinct; do not add these nested views. Optional real material sector directions use only saved geometry, no SVD. Failed case prerequisites still produce four failed readout rows and receipts.\n')
    print(json.dumps(receipt))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        if _ACTIVE_OUT is not None:
            target=_ACTIVE_OUT/'FAILED_RECEIPT.json'
            if not target.exists():target.write_text(json.dumps(dict(status='COUNTERFACTUAL_FAILED',error=repr(exc),whole_process_CPU_seconds=time.process_time(),wall_since_script_import=time.perf_counter()-_START_WALL,truth_access='NONE'),indent=2)+'\n')
        raise
