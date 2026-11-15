#!/usr/bin/env python3
"""Authorized OFFLINE_DIAGNOSTIC of old exposed R14 G19/DR endpoints only.

No optimization, candidate selection, new SVD, or new R15 prediction access.
Saved real/complex conditional geometry is loaded through EXECUTION maps.
"""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import argparse,csv,json,shutil,sys,time,traceback
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
MODULE_CPU=time.process_time();MODULE_WALL=time.perf_counter()
import numpy as np
from scipy.linalg import lu_solve
ROOT=Path(__file__).resolve().parents[1];WORKTREE=ROOT.parents[1]
sys.path.insert(0,str(WORKTREE))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

BANDS=('U11','U10','U01','U00')
SCOPE='OFFLINE_DIAGNOSTIC_OLD_PUBLIC_PREDICTIONS_ONLY'


def clean(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,Path):return str(v)
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v


def dump(path,value):path.write_text(json.dumps(clean(value),indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf8')


def table(path,rows):
    if not rows:path.write_text('status\nNO_ROWS\n');return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in rows:writer.writerow({k:json.dumps(clean(v),ensure_ascii=False) if isinstance(v,(dict,list,tuple)) else v for k,v in r.items()})


def sq(x):return float(np.vdot(x,x).real)


def ratio(n,d):
    if not np.isfinite(n) or not np.isfinite(d):return None,'nonfinite numerator or denominator'
    if d<=0:return None,'zero denominator'
    return float(n/d),None


def fraction_fields(n,d,prefix='fraction'):
    r,why=ratio(n,d);return {prefix:r,prefix+'_null_reason':why}


def project(U,x):return U@(U.conj().T@x)


def load_geometry(root):
    maps={};cache={};sources=[]
    for domain in ('complex','real'):
        maps[domain]={}
        for p in sorted((root/domain).glob('EXECUTION_*.json')):
            info=json.loads(p.read_text())
            if info.get('status')!='COMPLETE':continue
            name=Path(info['public_path']).name
            group=info['operator_group'];key=(domain,group)
            if key not in cache:
                folder=p.parent/group
                with np.load(folder/'SPACES.npz',allow_pickle=False) as z:arrays={k:z[k].copy() for k in z.files}
                metadata=json.loads((folder/'METADATA.json').read_text())
                if metadata['scalar_domain']!=('real_material' if domain=='real' else 'current_complex'):raise ValueError('Geometry scalar domain mismatch')
                if metadata['tau']!=.1 or metadata['numerical_rank_rtol']!=1e-12:raise ValueError('Geometry thresholds differ from authorization')
                cache[key]=(arrays,metadata)
                sources.append(dict(domain=domain,group=group,SPACES=str(folder/'SPACES.npz'),METADATA=str(folder/'METADATA.json')))
            if name in maps[domain]:raise ValueError('Duplicate public basename in geometry map: '+name)
            maps[domain][name]=(cache[key][0],cache[key][1],str(p))
    return maps,sources


def spectral_bands(arrays,operator):
    V=arrays[operator+'_right_frame'];s=arrays[operator+'_singular_values']
    full=np.zeros(V.shape[1]);full[:len(s)]=s
    maximum=float(s[0]) if len(s) else 0.
    if maximum<=0:raise ValueError('Zero operator spectrum not in this frozen cohort')
    rel=full/maximum
    return dict(strong=V[:,rel>=.1],weak_observable=V[:,(rel<.1)&(rel>1e-12)],numerical_null=V[:,rel<=1e-12])


def measured_contributions(model,state,components,Qbase,label):
    n,ns=Qbase.shape;k=len(components)
    rhs=(components[:,:,None]*Qbase[None]).transpose(1,0,2).reshape(n,k*ns)
    with model._physics.operation(label,rhs=k*ns,**{label+'_complex_RHS':k*ns}):
        result=lu_solve(state['lu'],rhs,check_finite=False)
    return result.reshape(n,k,ns).transpose(1,0,2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--authorization',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();authorization=args.authorization.resolve();out=args.out.resolve()
    if out.exists():raise FileExistsError('Retain original diagnostics; choose a fresh output: '+str(out))
    auth=json.loads(authorization.read_text())
    if auth.get('authorized') is not True or auth.get('expected')!=48:raise ValueError('Explicit authorization for48 old endpoints required')
    items=auth['items'];keys=[(r['case_id'],r['method']) for r in items]
    cases=sorted(set(c for c,m in keys))
    if len(items)!=48 or len(set(keys))!=48 or len(cases)!=24 or any({m for c,m in keys if c==case}!={'g19','dr'} for case in cases):raise ValueError('Authorization must contain precisely g19/dr ×24 old cases')
    oldroot=(WORKTREE/'research/research_max_struct_r14/RESULTS/independent24_v1').resolve()
    for item in items:
        if item['scope']!=SCOPE or item.get('new_online_predictions_allowed') is not False:raise ValueError('Wrong authorization scope')
        prediction=Path(item['prediction']).resolve()
        if not prediction.is_relative_to(oldroot) or prediction.name!='primary.npz':raise ValueError('Only old R14 primary predictions allowed')
        if item['truth_key']!='material_coarse':raise ValueError('Only authorized coarse truth key allowed')
    out.mkdir(parents=True);(out/'cases').mkdir();(out/'public_inputs').mkdir()
    shutil.copyfile(authorization,out/'AUTHORIZATION_AT_ACTIONS.json')
    receipt=dict(status='RUNNING',scope=SCOPE,new_online_predictions_allowed=False,new_candidate_selection=False,
                 authorization=str(authorization),expected=48,new_SVD=0,physical_model_source=str(Path(sys.modules[MaxwellModel.__module__].__file__).resolve()),
                 thresholds=dict(tau=.1,numerical_rtol=1e-12),BLAS_threads=1,errors=[],closure_role='Wiring metadata only, not a new theorem/proof claim')
    endpoints=[];bands=[];background=[];cross=[];cost_events=[];failures=[];total=Counter();case_stats=[]
    try:
        maps,cache_sources=load_geometry(Path(auth['geometry_cache']))
        receipt['geometry_sources']=cache_sources
        for case in cases:
            case_cpu,case_wall=time.process_time(),time.perf_counter();model=None
            selected=[r for r in items if r['case_id']==case]
            first=selected[0];truth_state=None;case_error=None
            try:
                if any(r['public_input']!=first['public_input'] or r['truth']!=first['truth'] for r in selected):raise ValueError('Case input/truth inconsistent across two authorized methods')
                public_path=Path(first['public_input']);truth_path=Path(first['truth'])
                complex_arrays,complex_meta,complex_execution=maps['complex'][public_path.name]
                real_arrays,real_meta,real_execution=maps['real'][public_path.name]
                with np.load(public_path,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
                with np.load(truth_path,allow_pickle=False) as z:mgt=np.asarray(z[first['truth_key']],float).reshape(-1).copy()
                shutil.copyfile(public_path,out/'public_inputs'/public_path.name)
                model=MaxwellModel(public)
                # Restore the original public R14 smooth/TV prior through inherited defaults.
                shape=tuple(public['grid_shape']);sigma=float(public['sigma']);Y=public['Y']
                if sigma<=0:raise ValueError('Positive public sigma required for requested whitened norm')
                gt_before=Counter(model.cost);c,w=time.process_time(),time.perf_counter()
                truth_state=model.state(mgt)
                gt_cost=dict(Counter(model.cost)-gt_before)
                cost_events.append(dict(case_id=case,method='SHARED_GT_REFERENCE',status='COMPLETE',CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,cost=gt_cost))
                folder=out/'cases'/case;folder.mkdir()
                np.savez_compressed(folder/'COARSE_GT_REFERENCE.npz',m_gt=mgt,J_gt=truth_state['j'],Q_gt=truth_state['Q'],Y_gt_coarse=truth_state['Y'])
                metadata=dict(case_id=case,public_input=str(public_path),truth=str(truth_path),truth_key=first['truth_key'],truth_use=SCOPE,
                              complex_geometry_execution=complex_execution,real_geometry_execution=real_execution,
                              complex_ranks=complex_meta['ranks'],real_ranks=real_meta['ranks'],sigma=sigma,sources=model.ns,receivers=model.nr,material_dof=model.n,
                              public_inverse_side=int(public.get('inverse_side',shape[0])),public_forward_side=int(public.get('forward_side',shape[0])),
                              GT_reference_scope='Same coarse inverse-model Maxwell reference only; cross-grid acquisition truth forward is not replaced or certified',
                              coarse_GT_vs_observation_raw_L2=sq(truth_state['Y']-Y)**.5,
                              coarse_GT_vs_observation_whitened_energy=sq((truth_state['Y']-Y)/sigma),
                              no_noise_certificate=True,GT_reference_cost=gt_cost)
                dump(folder/'CASE_METADATA.json',metadata)
                sbands=spectral_bands(real_arrays,'S');gbands=spectral_bands(real_arrays,'G')
            except BaseException as exc:
                case_error=repr(exc)
                failures.append(dict(case_id=case,method='CASE_SETUP',error=case_error,traceback=traceback.format_exc()))
            for item in selected:
                method=item['method'];row=dict(case_id=case,method=method,status='RUNNING',scope=SCOPE)
                before=Counter(model.cost) if model else Counter();c,w=time.process_time(),time.perf_counter()
                # Stage all rows per endpoint so a failure never leaves false COMPLETE summaries.
                staged_bands=[];staged_background=[];staged_cross=[]
                try:
                    if case_error:raise RuntimeError('Shared case setup failed: '+case_error)
                    prediction=Path(item['prediction'])
                    if prediction.stat().st_size!=item['prediction_bytes']:raise ValueError('Authorized prediction byte-size changed')
                    with np.load(prediction,allow_pickle=False) as z:mhat=np.asarray(z['m'],float).reshape(model.n).copy()
                    state=model.state(mhat);h=mhat-mgt;htotal=sq(h)
                    components=np.array([project(real_arrays[b],h) for b in BANDS],float)
                    chord=measured_contributions(model,state,components,truth_state['Q'],'OFFLINE_FINITE_CHORD')
                    local=measured_contributions(model,state,components,state['Q'],'OFFLINE_LOCAL_JACOBIAN')
                    phi,gphi=model.data_value_gradient(mhat,Y)
                    prior,gprior,_=model.reg.evaluate(mhat)
                    data_phi=phi-prior;gdata=gphi-gprior
                    dj=state['j']-truth_state['j'];dy=state['Y']-truth_state['Y'];dq=state['Q']-truth_state['Q']
                    chordG=np.einsum('ij,bjs->bis',model.G,chord);chordS=np.einsum('ri,bis->brs',model.S,chord)
                    localG=np.einsum('ij,bjs->bis',model.G,local);localS=np.einsum('ri,bis->brs',model.S,local)
                    model.cost['OFFLINE_dense_G_output_matrix_products']+=2
                    model.cost['OFFLINE_dense_G_output_complex_columns']+=8*model.ns
                    model.cost['OFFLINE_dense_S_output_matrix_products']+=2
                    model.cost['OFFLINE_dense_S_output_complex_columns']+=8*model.ns
                    material_fractions={};current_fractions={}
                    for i,b in enumerate(BANDS):
                        em=sq(components[i]);ej=sq(project(complex_arrays[b],dj))
                        material_fractions[b]=fraction_fields(em,htotal)
                        current_fractions[b]=fraction_fields(ej,sq(dj))
                        staged_bands.append(dict(case_id=case,method=method,band=b,scalar_domain='real material; complex current separate',
                            material_rank=real_meta['ranks'][b],current_rank=complex_meta['ranks'][b],material_squared_error=em,
                            material_error_fraction=material_fractions[b]['fraction'],material_error_fraction_null_reason=material_fractions[b]['fraction_null_reason'],
                            complex_current_squared_error=ej,complex_current_error_fraction=current_fractions[b]['fraction'],complex_current_fraction_null_reason=current_fractions[b]['fraction_null_reason'],
                            chord_current_energy=sq(chord[i]),chord_G_field_energy=sq(chordG[i]),chord_S_port_energy=sq(chordS[i]),chord_whitened_port_energy=sq(chordS[i]/sigma),
                            local_current_action_energy=sq(local[i]),local_G_field_action_energy=sq(localG[i]),local_S_port_action_energy=sq(localS[i]),local_whitened_port_action_energy=sq(localS[i]/sigma),
                            local_whitened_action_per_material_energy=ratio(sq(localS[i]/sigma),em)[0],local_gain_null_reason=ratio(sq(localS[i]/sigma),em)[1],
                            physical_Phi_gradient_dot_h_band=float(gphi@components[i]),data_Phi_gradient_dot_h_band=float(gdata@components[i])))
                    saved_background={}
                    for operator,groups in [('S_real',sbands),('G_real',gbands)]:
                        for name,U in groups.items():
                            hb=project(U,h);en=sq(hb);saved_background[operator+'_'+name]=hb
                            staged_background.append(dict(case_id=case,method=method,operator=operator,band=name,rank=U.shape[1],material_squared_error=en,**fraction_fields(en,htotal)))
                    hnull=saved_background['S_real_numerical_null']
                    null_chord=measured_contributions(model,state,hnull[None],truth_state['Q'],'OFFLINE_EXACT_SNULL_FINITE_CHORD')[0]
                    null_local=measured_contributions(model,state,hnull[None],state['Q'],'OFFLINE_EXACT_SNULL_LOCAL_JACOBIAN')[0]
                    null_S_chord=model.S@null_chord;null_S_local=model.S@null_local
                    model.cost['OFFLINE_dense_S_output_matrix_products']+=2
                    model.cost['OFFLINE_dense_S_output_complex_columns']+=2*model.ns
                    for kind,total_array,components_out in [('current',dj,chord),('G_field',dq,chordG),('S_port',dy,chordS),('whitened_S_port',dy/sigma,chordS/sigma)]:
                        isolated=[sq(v) for v in components_out];total_energy=sq(total_array)
                        staged_cross.append(dict(case_id=case,method=method,quantity=kind,total_actual_chord_energy=total_energy,sum_isolated_component_energies=sum(isolated),
                            signed_cross_energy=total_energy-sum(isolated),isolated_band_energies=dict(zip(BANDS,isolated)),
                            cross_over_total=ratio(total_energy-sum(isolated),total_energy)[0],cross_ratio_null_reason=ratio(total_energy-sum(isolated),total_energy)[1],
                            scope='Isolated propagated components are not an additive percentage allocation'))
                    local_total=np.sum(local,axis=0)
                    local_S_total=np.sum(localS,axis=0)
                    for kind,total_array,components_out in [('local_current',local_total,local),('local_G_field',np.sum(localG,axis=0),localG),('local_S_port',local_S_total,localS),('local_whitened_S_port',local_S_total/sigma,localS/sigma)]:
                        isolated=[sq(v) for v in components_out];total_energy=sq(total_array)
                        staged_cross.append(dict(case_id=case,method=method,quantity=kind,total_actual_chord_energy=None,total_local_action_energy=total_energy,sum_isolated_component_energies=sum(isolated),signed_cross_energy=total_energy-sum(isolated),isolated_band_energies=dict(zip(BANDS,isolated)),
                            cross_over_total=ratio(total_energy-sum(isolated),total_energy)[0],cross_ratio_null_reason=ratio(total_energy-sum(isolated),total_energy)[1],scope='Local derivative actions, not finite-chord state difference'))
                    model.cost['OFFLINE_dense_S_material_null_matrix_products']+=1
                    model.cost['OFFLINE_dense_S_material_null_columns']+=1
                    model.cost['OFFLINE_saved_basis_orthogonal_projection_actions']+=14
                    row.update(status='COMPLETE',prediction=str(prediction),sources=model.ns,receivers=model.nr,material_dof=model.n,material_squared_error=htotal,material_L2_error=htotal**.5,
                        physical_Phi=phi,data_Phi=data_phi,prior_value=prior,physical_relative_residual=ratio(sq(state['Y']-Y),sq(Y))[0]**.5 if sq(Y)>0 else None,
                        physical_model_residual=state['state_relative_residual'],coarse_GT_reference_residual=truth_state['state_relative_residual'],
                        true_current_error_energy=sq(dj),true_G_field_error_energy=sq(dq),true_port_error_energy=sq(dy),true_whitened_port_error_energy=sq(dy/sigma),
                        material_conditional_fractions=material_fractions,complex_current_conditional_fractions=current_fractions,
                        real_S_null_error_fraction=ratio(sq(hnull),htotal)[0],real_S_null_fraction_null_reason=ratio(sq(hnull),htotal)[1],
                        exact_real_Snull_background_S_action_energy=sq(model.S@hnull),
                        exact_real_Snull_chord_port_energy=sq(null_S_chord),exact_real_Snull_chord_whitened_port_energy=sq(null_S_chord/sigma),
                        exact_real_Snull_local_port_energy=sq(null_S_local),exact_real_Snull_local_whitened_port_energy=sq(null_S_local/sigma),
                        exact_real_Snull_local_gain_per_material_energy=ratio(sq(null_S_local/sigma),sq(hnull))[0],exact_real_Snull_local_gain_null_reason=ratio(sq(null_S_local/sigma),sq(hnull))[1],
                        physical_Phi_gradient_dot_exact_Snull=float(gphi@hnull),
                        wiring_material_component_closure_relative=ratio(sq(components.sum(axis=0)-h),htotal)[0],
                        wiring_current_chord_closure_relative=ratio(sq(chord.sum(axis=0)-dj),sq(dj))[0],
                        wiring_S_port_chord_closure_relative=ratio(sq(chordS.sum(axis=0)-dy),sq(dy))[0],
                        public_forward_side=int(public.get('forward_side',shape[0])),public_inverse_side=int(public.get('inverse_side',shape[0])),
                        reference_scope='Coarse inverse-model truth reference; cross-grid observations not noise-certified',noise_certificate_claim=False)
                    folder=out/'cases'/case
                    np.savez_compressed(folder/(method+'_ENDPOINT_ARRAYS.npz'),m_hat=mhat,m_gt=mgt,h=h,h_conditional=components,
                        J_hat=state['j'],Q_hat=state['Q'],Y_hat=state['Y'],true_current_difference=dj,true_G_field_difference=dq,true_port_difference=dy,
                        finite_chord_current_contributions=chord,finite_chord_G_field_contributions=chordG,finite_chord_S_port_contributions=chordS,
                        local_current_actions=local,local_G_field_actions=localG,local_S_port_actions=localS,
                        physical_Phi_gradient=gphi,data_Phi_gradient=gdata,exact_Snull_chord_current=null_chord,exact_Snull_local_current=null_local,
                        exact_Snull_chord_S_port=null_S_chord,exact_Snull_local_S_port=null_S_local,
                        **{'h_'+k:v for k,v in saved_background.items()})
                    bands+=staged_bands;background+=staged_background;cross+=staged_cross
                except BaseException as exc:
                    row.update(status='FAILED',error=repr(exc),traceback=traceback.format_exc())
                    failures.append(dict(case_id=case,method=method,error=repr(exc),traceback=traceback.format_exc()))
                finally:
                    now=Counter(model.cost) if model else Counter();delta={k:now[k]-before[k] for k in now if now[k]!=before[k]}
                    row.update(CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,physical_cost=delta)
                    cost_events.append(dict(case_id=case,method=method,status=row['status'],CPU_seconds=row['CPU_seconds'],wall_seconds=row['wall_seconds'],cost=delta))
                    endpoints.append(row)
                    if (out/'cases'/case).is_dir():dump(out/'cases'/case/(method+'_RECEIPT.json'),row)
            if model:
                total.update(model.cost)
                case_stats.append(dict(case_id=case,status='COMPLETE' if case_error is None else 'FAILED',physical_cost=dict(model.cost),physical_wall_by_operation=dict(model._physics.timings),physical_CPU_by_operation=dict(model._physics.cpu_timings),CPU_seconds=time.process_time()-case_cpu,wall_seconds=time.perf_counter()-case_wall))
            dump(out/'PROGRESS.json',dict(completed_rows=len(endpoints),failed_rows=sum(r['status']!='COMPLETE' for r in endpoints),last_case=case,whole_CPU_seconds=time.process_time()-MODULE_CPU,whole_wall_seconds=time.perf_counter()-MODULE_WALL))
        table(out/'ENDPOINTS.csv',endpoints);table(out/'CONDITIONAL_BANDS.csv',bands);table(out/'BACKGROUND_MATERIAL_BANDS.csv',background);table(out/'PROPAGATED_CROSS_ENERGIES.csv',cross);table(out/'COST_EVENTS.csv',cost_events);table(out/'FAILURES.csv',failures)
        dump(out/'ENDPOINTS.json',endpoints);dump(out/'FAILURES.json',failures);dump(out/'CASE_COSTS.json',case_stats)
        good=[r for r in endpoints if r['status']=='COMPLETE']
        receipt.update(status='COMPLETE_WITH_FAILURES' if failures else 'COMPLETE',completed_endpoints=len(good),failed_endpoints=48-len(good),physical_cost=dict(total),GT_coarse_reference_cases=sum(x['method']=='SHARED_GT_REFERENCE' and x['status']=='COMPLETE' for x in cost_events),
            physics_counter_scope='Each case model summed once. ENDPOINT cost deltas exclude shared GT; CASE_COSTS includes GT and endpoints. Operation-specific roles overlap aggregate actual_complex_RHS: never add both.',
            exact_Snull_extra_scope='Parent-approved separate ns finite chord plus ns local RHS per endpoint; not extra conditional band or additive share.',
            diagonal_sigma_scope='Scalar public sigma used only to report whitened sensitivity; no noise or GT certificate.',
            errors=failures)
        lines=['# Old exposed R14 endpoint geometry — offline diagnostic','',
               'This is posthoc analysis of48 already-public G19/DR primary endpoints (24 cases). It is not a new confirmation cohort, predictor, candidate selector or new R15 online result. Authorized material_coarse is used only here.',
               '',f"Completed endpoints: {len(good)}/48. Failures: {48-len(good)}. Shared coarse GT forward states: {receipt['GT_coarse_reference_cases']}/24.",
               '', 'Saved real-domain conditional projections describe material error. Saved complex-domain projections describe the actual current difference. They are not interchangeable. All four conditional groups, real background S strong/weak-observable/numerical-null and global G groups retain separate labels, ranks and denominators.',
               '', 'Finite chord contributions use Lhat inverse on h_band Qgt with one grouped4ns solve; local derivative actions use Qhat in a distinct4ns solve. Their G/S isolated energies and signed cross terms are recorded separately. No isolated-energy percentages are forced to sum to100%. Closure defects are wiring metadata only, not a new proof result.',
               '', 'The exact numerical real-S null component has separately charged chord/local ns solves. A small raw background S h_null does not imply its total-field/resolvent-transported local or finite-chord port action vanishes. Actual per-endpoint energies are in ENDPOINTS.csv; no material-identifiability or recovery conclusion is inferred from the fixed basis alone.',
               '', 'Coarse GT Maxwell states are inverse-model references. Cases generated on a finer grid retain their acquisition mismatch; coarse GT versus measured observations does not produce a noise certificate. Physical Phi and its gradient come from the original model and include its unchanged regularizer; data-only value/gradient are also separated.',
               '', '| Method | endpoints | pooled real S-null material-error fraction | pooled whitened S-null local-action energy / material energy |',
               '|---|---:|---:|---:|']
        for method in ('g19','dr'):
            rr=[r for r in good if r['method']==method]
            null_energy=sum(r['material_squared_error']*(r['real_S_null_error_fraction'] or 0) for r in rr)
            total_energy=sum(r['material_squared_error'] for r in rr)
            fr=ratio(null_energy,total_energy)[0];gain=ratio(sum(r['exact_real_Snull_local_whitened_port_energy'] for r in rr),null_energy)[0]
            lines.append(f'| {method} | {len(rr)} | {fr if fr is not None else "undefined"} | {gain if gain is not None else "undefined"} |')
        lines+=['','Ratios are pooled squared-norm summaries, not averages of percentages and not algorithm qualification. All endpoint values and null reasons remain available.','',
                'Artifacts: ENDPOINTS.csv/json; CONDITIONAL_BANDS.csv; BACKGROUND_MATERIAL_BANDS.csv; PROPAGATED_CROSS_ENERGIES.csv; CASE_COSTS.json; COST_EVENTS.csv; per-case public inputs, GT reference, endpoint arrays and metadata. These files contain offline authorized truth-derived material/error fields and must not become inputs to online inverse code.','']
        (out/'REPORT.md').write_text('\n'.join(lines),encoding='utf8')
    except BaseException as exc:
        receipt.update(status='FAILED',fatal_error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        receipt.update(CPU_seconds=time.process_time()-MODULE_CPU,wall_seconds=time.perf_counter()-MODULE_WALL,whole_process_CPU_seconds_at_receipt=time.process_time(),completed_utc=datetime.now(timezone.utc).isoformat())
        dump(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k] for k in ['status','completed_endpoints','failed_endpoints','physical_cost','CPU_seconds','wall_seconds']},indent=2))


if __name__=='__main__':main()
