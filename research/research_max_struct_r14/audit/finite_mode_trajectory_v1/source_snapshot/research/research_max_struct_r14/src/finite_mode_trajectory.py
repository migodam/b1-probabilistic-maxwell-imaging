"""Declared leading-real-mode finite-program diagnostic; not an inverse method.

The artificial root and saved ideal01 map are the only material reference.
One fixed epsilon, four pairs, two signs, two explicit PR-memory policies.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime,timezone
import json
from pathlib import Path
import shutil
import sys
import time
import traceback


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args();root=args.root.resolve();repo=root.parents[1]
    out=(args.out or root/'audit/finite_mode_trajectory_v1').resolve();out.mkdir(parents=True,exist_ok=False)
    sys.path[:0]=[str(root/'src'),str(repo)]
    cpu,wall=time.process_time(),time.perf_counter()
    import numpy as np
    from dualmetric_core import StatefulDualMetric,_load_state_npz,_save_state_npz
    from run_vector_factorial import dump
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r13.src.score_results import table
    source=root/'audit/local_closed_loop_v1'
    costs=Counter();events=[];failures=[];rows=[];physical_rows=[];trajectory_summaries=[];first_steps={}
    model=None;status='RUNNING';selection={};state_source=source/'ROOT_FULL_STATE.npz'
    def norm(value):return float(np.linalg.norm(value))
    def delta(after,before):return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before)
        if isinstance(after.get(k,0),(int,float)) and isinstance(before.get(k,0),(int,float)) and after.get(k,0)!=before.get(k,0)}
    def charged(label,function,solver=None,metadata=None):
        start_cpu,start_wall=time.process_time(),time.perf_counter()
        before=dict(model.cost) if model is not None else {};aux=dict(solver.cost) if solver is not None else {}
        event=dict(label=label,metadata=metadata)
        try:
            result=function();event['status']='COMPLETE';return result
        except BaseException as exc:
            event.update(status='FAILED',error=repr(exc));raise
        finally:
            event.update(CPU_seconds=time.process_time()-start_cpu,wall_seconds=time.perf_counter()-start_wall,
                physical_operation_delta=delta(dict(model.cost) if model is not None else {},before),
                auxiliary_operation_delta=delta(dict(solver.cost) if solver is not None else {},aux))
            events.append(event)
    try:
        saved=charged('read complete artificial-root checkpoint',lambda:_load_state_npz(state_source))
        def read_public():
            with np.load(source/'ROOT_PUBLIC_AND_PHYSICAL.npz',allow_pickle=False) as z:
                return {key:z[key].copy() for key in z.files}
        public=charged('read declared artificial-root public operators and saved physical fields',read_public)
        maps={}
        for pair in ('00','01','10','11'):
            def read_map():
                with np.load(source/f'IDEAL_{pair}.npz',allow_pickle=False) as z:return z['T'].copy()
            maps[pair]=charged('read saved ideal material map',read_map,metadata=dict(pair=pair))
        eigenvalues,eigenvectors=charged('leading ideal01 eigendecomposition',lambda:np.linalg.eig(maps['01']))
        costs['leading_ideal01_eigendecompositions']+=1
        index=int(np.argmax(abs(eigenvalues)));value=eigenvalues[index];vector=eigenvectors[:,index]
        tolerance=1e-10;imaginary_relative=norm(vector.imag)/max(norm(vector),1e-300)
        selection=dict(selection='largest modulus eigenvalue of saved IDEAL_01, no replacement or direction search',
            leading_index=index,eigenvalue_real=float(value.real),eigenvalue_imag=float(value.imag),
            eigenvector_relative_imaginary_norm=imaginary_relative,
            real_tolerance=tolerance,expected_parent_reported_rho=1.0027124491177037,
            real_mode_eligible=bool(abs(value.imag)<=tolerance*(1+abs(value.real)) and imaginary_relative<=tolerance))
        dump(out/'MODE_SELECTION.json',selection)
        if not selection['real_mode_eligible']:
            status='NOT_RUN_COMPLEX_LEADING_MODE';return 0
        h=vector.real.copy();h/=norm(h);pivot=int(np.argmax(abs(h)))
        if h[pivot]<0:h=-h
        rho=float(value.real);mroot=public['m_anchor'].copy();Jroot=public['J_anchor'].copy();Qroot=public['Q_anchor'].copy();Y=public['Y_anchor'].copy()
        epsilon=.001
        if np.any(mroot-epsilon*abs(h)<.02) or np.any(mroot+epsilon*abs(h)>20):
            raise ValueError('Declared fixed epsilon leaves original material box; no replacement amplitude permitted')
        original_flags=deepcopy(saved['frozen_config']);extended=deepcopy(saved)
        # The sole numerical schedule extension is explicit and frozen.
        extended['frozen_config']['stages'][0]['maxiter']=32
        extended['done']=False
        if (extended['stage'],extended['iteration'],extended['iterations'])!=(0,0,0):
            raise ValueError('Input is not the prescribed untouched root iteration0 state')
        specification=dict(status='FROZEN_BEFORE_FIRST_TRAJECTORY',utc=datetime.now(timezone.utc).isoformat(),
            source_checkpoint=str(state_source),source_public=str(source/'ROOT_PUBLIC_AND_PHYSICAL.npz'),
            ideal01=str(source/'IDEAL_01.npz'),scope='offline finite directional trajectory at declared artificial root, not inverse recovery',
            epsilon=epsilon,steps=32,pairs=['00','01','10','11'],signs=[1,-1],
            memory_policies=['native_PR','reset_before_every_update'],fixed_snapshots=[1,8,32],
            total_native_updates=512,root_mode_selection=selection,mode_norm=norm(h),fixed_sign_pivot=pivot,
            mode_sign_at_pivot=float(h[pivot]),root_reference='declared artificial m_anchor only; no GT/test arrays',
            original_frozen_flags=original_flags,extended_frozen_flags=extended['frozen_config'],
            explicit_changes=['stage0.maxiter:1->32','done=False; rootstage/iteration/globaliteration stay0',
                'initial m=mroot+sign*.001*h only; rootA/J/Q/DPCS/weights unchanged',
                'reset policy clears old_grad/old_direction BEFORE every update; nativepolicy preserves original PR'],
            physical_audit='force new original64x64 LU for each1/8/32 snapshot,6 sourceRHS; physical audits never feed auxiliary update',
            no_amplitude_search=True,no_alternative_mode=True,no_new_optimizer=True,
            prohibited_claims=['FrechetJacobian existence','infiniteLyapunov exponent','global recovery','novel theorem','continuum/PDE certificate'])
        dump(out/'FROZEN_CONFIG.json',specification)
        snapshot=out/'source_snapshot';snapshot.mkdir();files=[]
        sources=[Path(__file__).resolve(),root/'src/dualmetric_core.py',root/'src/local_closed_loop.py',
            repo/'research/research_max_struct_r7/src/aplus_readout.py',
            repo/'research/research_max_struct_r4/src/cie_onehop.py',
            repo/'research/research_fbg_r3_experiments/src/cie_tsom_baseline.py',
            state_source,source/'ROOT_PUBLIC_AND_PHYSICAL.npz',*[source/f'IDEAL_{pair}.npz' for pair in maps]]
        def freeze():
            for path in sources:
                destination=snapshot/path.relative_to(repo);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
                files.append(dict(source=str(path),snapshot=str(destination),bytes=destination.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
            dump(snapshot/'PROVENANCE.json',dict(status='FROZEN_BEFORE_TRIALS',files=files,GT_access='no test/GT arrays; only declared artificialroot'))
        charged('freeze exact source and artificial-root inputs before trials',freeze)
        def setup_model():
            nonlocal model
            model=MaxwellModel(public);model.reg.smooth=0.;model.reg.tv=0.
        charged('existing64-cell original Maxwell model and zero prior setup',setup_model)
        cfg=deepcopy(extended['frozen_config']);cfg.update(regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg)
        expected={}
        for pair,T in maps.items():
            vectors=[h.copy()]
            for step in range(32):vectors.append(T@vectors[-1]);costs['ideal_repeated_material_map_vectors']+=1
            expected[pair]=np.asarray(vectors)
        np.savez_compressed(out/'MODE_AND_IDEAL_PREDICTIONS.npz',h=h,eigenvalue=value,
            eigen_residual=maps['01']@h-rho*h,mode_sign_pivot=pivot,
            pair_ids=np.asarray(list(maps)),ideal_T_h=np.asarray([maps[pair]@h for pair in maps]),
            ideal_material_traces=np.asarray([expected[pair] for pair in maps]),
            rho01_powers=rho**np.arange(33),epsilon=epsilon,mroot=mroot,Jroot=Jroot,Qroot=Qroot,Yroot=Y)
        def diagnostics(sol):
            port=sol.S@sol.J-Y;sol.cost['S_columns']+=sol.ns
            defect=sol.J-sol.m[:,None]*sol.Q
            current=sol._P(sol.m,defect,sol.pair[0]);material=sol._P(sol.m,defect,sol.pair[1])
            sol.cost['source_metric_products']+=3
            lower=sol.m<=sol.lo+1e-8;upper=sol.m>=sol.hi-1e-8
            return dict(auxiliary_port_relative=norm(port)/norm(Y),constitutive_defect_relative_to_J=norm(defect)/max(norm(sol.J),1e-300),
                constitutive_defect_norm=norm(defect),auxiliary_data_loss=.5*float(np.vdot(port,port@sol.Od).real),
                auxiliary_current_state_loss=.5*sol.rho*float(np.vdot(current,current@sol.Os).real),
                auxiliary_material_state_loss=.5*sol.rho*float(np.vdot(material,material@sol.Os).real),
                lower_active=int(lower.sum()),upper_active=int(upper.sum())),dict(port=port,defect=defect,lower_active_mask=lower,upper_active_mask=upper)
        def physical_audit(sol,folder,step,label):
            # Explicitly pay a fresh audit even when the material repeats.
            model._current=None
            state=model.state(sol.m);residual=state['Y']-Y
            D=.5*norm(residual)**2/norm(Y)**2;R=float(model.reg.evaluate(sol.m)[0])
            np.savez_compressed(folder/f'PHYSICAL_{step:02d}.npz',m=sol.m,physical_J=state['j'],physical_Q=state['Q'],
                physical_Y=state['Y'],physical_residual=residual,auxiliary_J=sol.J,auxiliary_Q=sol.Q)
            return dict(label=label,step=step,Phi=D+R,D=D,R=R,relative_residual=norm(residual)/norm(Y),
                physical_constitutive_relative=float(state['state_relative_residual']),physical_cost=dict(model.cost),fresh_LU_dimension=64,source_RHS=6)
        for pair in maps:
            for sign in (1,-1):
                for policy in ('native_PR','reset_before_every_update'):
                    label=f'pair{pair}__sign{sign:+d}__{policy}';folder=out/'trajectories'/label;folder.mkdir(parents=True)
                    trial_cpu,trial_wall=time.process_time(),time.perf_counter();before_physical=dict(model.cost)
                    trial_rows=[];trace_m=[];trace_J=[];trace_Q=[];trial_counts=Counter();trial_status='RUNNING';sol=None
                    try:
                        sol=charged(label+' complete-state restore',lambda:StatefulDualMetric.from_state_dict(model,Y,cfg,extended,pair=pair))
                        sol.m=mroot+sign*epsilon*h
                        charged(label+' fullstart export',lambda:sol.checkpoint(folder/'START_STATE.npz'),solver=sol)
                        trace_m.append(sol.m.copy());trace_J.append(sol.J.copy());trace_Q.append(sol.Q.copy())
                        for step in range(1,33):
                            previous_m,previous_J,previous_Q=sol.m.copy(),sol.J.copy(),sol.Q.copy()
                            if policy=='reset_before_every_update':sol.old_grad=sol.old_direction=None
                            before=dict(sol.cost)
                            row=charged(label+' native completed update',sol.step,solver=sol,metadata=dict(step=step))
                            trial_counts.update(delta(dict(sol.cost),before))
                            diag,arrays=charged(label+' cached-state observational diagnostics',lambda:diagnostics(sol),solver=sol,metadata=dict(step=step))
                            dm=sol.m-mroot
                            result=dict(label=label,pair=pair,sign=sign,memory_policy=policy,step=step,
                                material_gain=norm(dm)/epsilon,signed_mode_coefficient=float(np.dot(h,dm)/epsilon),
                                sign_aligned_mode_gain=float(sign*np.dot(h,dm)/epsilon),
                                material_step_over_epsilon=norm(sol.m-previous_m)/epsilon,
                                current_root_relative=norm(sol.J-Jroot)/norm(Jroot),field_root_relative=norm(sol.Q-Qroot)/norm(Qroot),
                                current_relative_step=norm(sol.J-previous_J)/max(norm(previous_J),1e-300),
                                field_relative_step=norm(sol.Q-previous_Q)/max(norm(previous_Q),1e-300),
                                ideal_repeated_map_gain=norm(expected[pair][step]),
                                ideal_signed_mode_coefficient=float(sign*np.dot(h,expected[pair][step])),
                                ideal01_rho_power=float(rho**step),ideal01_signed_rho_prediction=float(sign*rho**step),
                                alpha=row['line_alpha'],PR_plus=row['PR_plus'],direction_restarted=row['direction_restarted'],
                                material_readout_nit=row.get('material_inner_nit'),material_fallback=row.get('material_fallback'),
                                material_scalar_seed_preferred=None if 'material_old_block_loss' not in row else row['material_scalar_seed_block_loss']<row['material_old_block_loss'],
                                Arc_guard=row.get('aplus_inherited_gtol'),Arc_projected_gradient_inf=row.get('aplus_projected_gradient_inf'),
                                coefficient_relative_change_squared=row['coefficient_relative_change_squared'],stage_end=row['stage_end'],
                                **diag)
                            rows.append(result);trial_rows.append(result);trace_m.append(sol.m.copy());trace_J.append(sol.J.copy());trace_Q.append(sol.Q.copy())
                            if step==1:first_steps[(pair,sign,policy)]={key:deepcopy(getattr(sol,key)) for key in ('m','A','J','Q','old_grad','old_direction')}
                            if step in (1,8,32):
                                charged(label+' full completed checkpoint',lambda:sol.checkpoint(folder/f'STATE_{step:02d}.npz'),solver=sol,metadata=dict(step=step))
                                np.savez_compressed(folder/f'DIAGNOSTICS_{step:02d}.npz',**arrays,ideal_T_h=maps[pair]@h,ideal_Tpow_h=expected[pair][step])
                                paid=charged(label+' original physical audit',lambda:physical_audit(sol,folder,step,label),solver=sol,metadata=dict(step=step,LU_dimension=64,source_RHS=6))
                                physical_rows.append(paid);result.update(physical_Phi=paid['Phi'],physical_relative_residual=paid['relative_residual'])
                        trial_status='COMPLETE'
                    except Exception as exc:
                        trial_status='FAILED';failure=dict(label=label,error=repr(exc),traceback=traceback.format_exc());failures.append(failure)
                        (folder/'FAILURE.txt').write_text(failure['traceback'])
                    finally:
                        np.savez_compressed(folder/'TRACE.npz',m=np.asarray(trace_m),J=np.asarray(trace_J),Q=np.asarray(trace_Q))
                        dump(folder/'HISTORY.json',trial_rows)
                        trajectory_summaries.append(dict(label=label,pair=pair,sign=sign,memory_policy=policy,status=trial_status,
                            completed_updates=len(trial_rows),CPU_seconds=time.process_time()-trial_cpu,wall_seconds=time.perf_counter()-trial_wall,
                            native_step_operation_counts=dict(trial_counts),physical_delta=delta(dict(model.cost),before_physical),
                            final_material_gain=trial_rows[-1]['material_gain'] if trial_rows else None,
                            final_sign_aligned_mode_gain=trial_rows[-1]['sign_aligned_mode_gain'] if trial_rows else None,
                            final_ideal01_rho_prediction=float(rho**32)))
                        dump(folder/'RECEIPT.json',trajectory_summaries[-1])
        equality=[]
        for pair in maps:
            for sign in (1,-1):
                a=first_steps.get((pair,sign,'native_PR'));b=first_steps.get((pair,sign,'reset_before_every_update'))
                arrays=None if a is None or b is None else {key:bool(np.array_equal(a[key],b[key])) for key in a}
                equality.append(dict(pair=pair,sign=sign,arrays_equal=arrays,all_first_step_arrays_equal=None if arrays is None else all(arrays.values()),
                    use='actual observed comparison, never a stopping/admission gate'))
        dump(out/'FIRST_STEP_POLICY_EQUALITY.json',equality)
        table(out,'TRAJECTORIES',trajectory_summaries);table(out,'FINITE_PROGRAM_ROWS',rows);table(out,'PAID_PHYSICAL_AUDITS',physical_rows)
        (out/'INTERPRETATION_BOUNDARIES.md').write_text('''# Finite-program directional diagnostic\n\nThe h direction is the largest-modulus real eigenmode of saved IDEAL01, normalized once with fixed sign; it is not the previous largest-difference h1. Epsilon is exactly.001 with no amplitude/direction search. The sole stage extension is maxiter1→32, preserving fullmf4/native64coefficients, beta0=1/beta1=.2, sourceweights, zero prior, rootA/J/Q and flags.\n\nThese16 trajectories compare actual32-step inherited PR/Arc readout arithmetic with saved ideal exact-block-map T^n h and the ideal01 rho^n counterfactual. NativePR and clearing oldgradient/direction before everyupdate are distinct declared policies; all8 first-step pairs are compared as observations, not gates.\n\nAt this finite epsilon, floating floors, rootseed reconstruction gaps, stationary-root choices and Arc projected-gradient guards may affect directional response. No FrechetJacobian existence, infinite-timeLyapunov exponent, globalrecovery, newtheorem or continuum/PDE certificate follows. Only the declared artificialroot is used; no GT/test arrays are inputs. Parent owns mathematical/scientific interpretation.\n''')
        status='COMPLETE' if not failures else 'PARTIAL_WITH_FAILURES'
    except Exception as exc:
        failures.append(dict(label='setup or summary',error=repr(exc),traceback=traceback.format_exc()));status='FAILED';traceback.print_exc()
    finally:
        dump(out/'FAILURES.json',failures);dump(out/'COST_LEDGER.json',events)
        receipt=dict(status=status,trajectories=len(trajectory_summaries),completed_trajectories=sum(row['status']=='COMPLETE' for row in trajectory_summaries),
            native_completed_updates=len(rows),physical_audits=len(physical_rows),physical_cost={} if model is None else dict(model.cost),
            setup_ideal_counts=dict(costs),CPU_seconds=time.process_time()-cpu,whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
            mode_selection=selection,failed_attempts=len(failures),CUDA_seconds=0,
            source_and_config_frozen_before_trials=True,amplitude_search=False,GT_access='declared artificialroot only; no GT/test arrays',
            preflight_note='prior read-only source/root/eigen inspection occurred before this timed runner; that inspection CPU was not separately instrumented',
            scientific_acceptance='parent controlled finite-scale interpretation only; no derivative/Lyapunov/recovery theorem claim')
        dump(out/'RECEIPT.json',receipt);print(json.dumps(dict(status=status,trajectories=receipt['completed_trajectories'],updates=len(rows),physical_audits=len(physical_rows),physical_cost=receipt['physical_cost'],CPU_seconds=receipt['CPU_seconds'])))
    return 0 if status in ('COMPLETE','NOT_RUN_COMPLEX_LEADING_MODE') else 1


if __name__=='__main__':sys.exit(main())
