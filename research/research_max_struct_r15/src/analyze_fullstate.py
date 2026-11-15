"""No-GT saved-state dynamics for R15; never invokes a solver or forward model.

Reads actual complete-state codec and public-derived checkpoints only. Default
output is a NEW analysis/fullstate_v1. Missing/failed expected entries remain.
Finite differences are not local spectra, Lyapunov exponents or chaos tests.
"""
from __future__ import annotations
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,csv,importlib.util,json,math,sys,time
from pathlib import Path
from datetime import datetime,timezone
METHODS=('twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset')
CASES=tuple(f'r14_new_{i:04d}' for i in range(24))
MAIN='twofold_exposed24_v1';FORK='fullstate_interventions_v1'
_ACTIVE_OUT=None
_START_WALL=time.perf_counter()


def finite(x):return isinstance(x,(int,float)) and math.isfinite(x)
def fresh(cumulative,inherited):
    if not isinstance(cumulative,dict) or not isinstance(inherited,dict):return None
    return {k:cumulative.get(k,0)-inherited.get(k,0) for k in set(cumulative)|set(inherited) if isinstance(cumulative.get(k,0),(int,float)) and isinstance(inherited.get(k,0),(int,float))}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path);p.add_argument('--out',type=Path);a=p.parse_args();root=a.root.resolve();out=(a.out or root/'analysis/fullstate_v1').resolve();start=time.perf_counter();cpu=time.process_time()
    if out.exists() or not out.is_relative_to(root/'analysis'):raise ValueError('New R15 analysis directory required')
    global _ACTIVE_OUT
    out.mkdir(parents=True);_ACTIVE_OUT=out;legacy=(a.legacy_repo or root.parent.parent).resolve();sys.path[:0]=[str(root/'runtime'),str(root/'src'),str(legacy),str(legacy/'vendor')]
    import numpy as np
    # Load the real inherited safe codec, no checkpoint field-name emulation,
    # optimizer/model construction, or numerical state reinitialization.
    source=root/'runtime/dualmetric_core.py';spec=importlib.util.spec_from_file_location('r15_analysis_state_codec',source);codec=importlib.util.module_from_spec(spec);spec.loader.exec_module(codec)
    issues=[];read_paths=set();cache={};basis_cache={}
    def read(path,default=None):
        path=Path(path);read_paths.add(str(path))
        if not path.is_file():return default
        try:return json.loads(path.read_text())
        except Exception as exc:issues.append(dict(path=str(path),issue='JSON_READ_ERROR',error=repr(exc)));return default
    def arrays(path):
        path=Path(path);read_paths.add(str(path))
        if not path.is_file():return None
        try:
            with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}
        except Exception as exc:issues.append(dict(path=str(path),issue='ARRAY_READ_ERROR',error=repr(exc)));return None
    def state(path):
        path=Path(path);read_paths.add(str(path))
        if not path.is_file():return None
        try:return codec._load_state_npz(path)
        except Exception as exc:issues.append(dict(path=str(path),issue='COMPLETE_STATE_READ_ERROR',error=repr(exc)));return None
    def nrm(x):return float(np.linalg.norm(x))
    def scalar(x):return x.item() if isinstance(x,np.ndarray) and x.ndim==0 else x
    def view(s,full=False):
        if s is None:return None
        try:
            if full:
                result=dict(m=s['m'],J=s['J'],Q=s['Q'],A=s['A'],memory_X=s['r15']['memory_X'],memory_Q=s['r15']['memory_Q'],
                    stage=s['stage'],iteration=s['iteration'],beta=[s['beta0'],s['beta1']],old_grad=s['old_grad'],old_direction=s['old_direction'])
            else:result={k:scalar(s[k]) for k in ('m','J','Q','A','memory_X','memory_Q','stage','iteration','beta')}
            if not all(np.isfinite(result[k]).all() for k in ('m','J','Q','A','memory_X','memory_Q','beta')):raise ValueError('Nonfinite saved state')
            return result
        except Exception as exc:issues.append(dict(issue='SAVED_STATE_SCHEMA_OR_FINITE_ERROR',error=repr(exc)));return None
    def phase(s):return dict(stage=int(s['stage']),iteration=int(s['iteration']),beta=np.asarray(s['beta'],float).tolist()) if s else dict(stage=None,iteration=None,beta=None)
    def compare(x,y,eps=None):
        if x is None or y is None:return dict(available=False,reason='MISSING_COMPARISON_STATE')
        if any(x[k].shape!=y[k].shape for k in ('m','J','Q','memory_X','memory_Q')):
            return dict(available=False,reason='SAVED_STATE_SHAPE_MISMATCH',before_stage=phase(y),after_stage=phase(x))
        result=dict(available=True,**{'delta_'+k+'_norm':nrm(x[k]-y[k]) if x[k].shape==y[k].shape else None for k in ('m','J','Q','A','memory_X','memory_Q')})
        result['coefficient_comparison_null_reason']=None if x['A'].shape==y['A'].shape else 'DIFFERENT_NATIVE_COEFFICIENT_MASK_SHAPE'
        for key in ('old_grad','old_direction'):
            if key in x and key in y:
                result[key+'_both_none']=x[key] is None and y[key] is None
                result[key+'_distance']=nrm(x[key]-y[key]) if x[key] is not None and y[key] is not None and x[key].shape==y[key].shape else None
        result['scaled_field_state_distance']=float(np.hypot(nrm((x['J']-y['J'])/6.),nrm(x['Q']-y['Q']))) if x['J'].shape==y['J'].shape else None
        result['normal_memory_distance']=float(np.hypot(nrm(x['memory_X']-y['memory_X']),nrm(x['memory_Q']-y['memory_Q'])))
        result.update(before_stage=phase(y),after_stage=phase(x),same_stage_and_beta=bool(x['stage']==y['stage'] and np.array_equal(x['beta'],y['beta'])))
        if eps is not None:
            result.update(realized_perturbation_norm=eps,material_finite_gain=result['delta_m_norm']/eps if eps>0 else None,
                scaled_field_finite_gain=result['scaled_field_state_distance']/eps if eps>0 and result['scaled_field_state_distance'] is not None else None,
                normal_memory_finite_gain=result['normal_memory_distance']/eps if eps>0 else None,
                coefficient_finite_gain=result['delta_A_norm']/eps if eps>0 and result['delta_A_norm'] is not None else None,
                gain_null_reason=None if eps>0 else 'ZERO_REALIZED_PERTURBATION',interpretation='Finite material-only perturbation response per realized input norm; mixed material/field physical units are not identified; never a local spectral radius')
        return result
    geometry={}
    for path in sorted((root/'audit/rank_all24/real').glob('EXECUTION_*.json')):
        row=read(path,{})
        if row.get('status')!='COMPLETE':continue
        cid=str(row.get('public_path','')).replace('\\','/').rsplit('/',1)[-1].removesuffix('.npz');group=row.get('operator_group')
        if cid in CASES and group:
            if Path(group).name!=group or group in ('.','..'):raise ValueError('Unsafe saved geometry group')
            geometry[cid]=dict(path=root/'audit/rank_all24/real'/group/'SPACES.npz',execution=str(path),group=group)
    def basis(cid):
        g=geometry.get(cid)
        if g is None:return None
        key=str(g['path'])
        if key not in basis_cache:
            metadata=read(g['path'].parent/'METADATA.json',{})
            if metadata.get('scalar_domain')!='real_material':issues.append(dict(case_id=cid,issue='NOT_REAL_MATERIAL_BASIS'));basis_cache[key]=None
            else:
                try:
                    with np.load(g['path'],allow_pickle=False) as z:basis_cache[key]={k:z[k].copy() for k in ('U11','U10','U01','U00')}
                    read_paths.add(key)
                except Exception as exc:issues.append(dict(path=key,issue='BASIS_READ_ERROR',error=repr(exc)));basis_cache[key]=None
        return basis_cache[key]
    def sectors(dm,cid):
        b=basis(cid);energy=float(np.vdot(dm,dm).real)
        if b is None:return dict(sector_available=False,sector_null_reason='MISSING_REAL_GEOMETRY')
        d=dict(sector_available=True,material_change_norm=float(np.sqrt(energy)),real_sector_ranks={k:v.shape[1] for k,v in b.items()})
        for k,U in b.items():
            e=float(np.linalg.norm(U.conj().T@dm)**2);d[k+'_material_energy']=e;d[k+'_material_fraction']=e/energy if energy>0 else None
        d['material_outside_real_U01_fraction']=1-d['U01_material_fraction'] if energy>0 else None
        d['zero_material_change']=energy==0;d['fraction_null_reason']='ZERO_MATERIAL_CHANGE' if energy==0 else None
        d['geometry_scope']='Real material four-sector diagnostic; current routing H is not a material-step constraint'
        return d
    # 480 expected before/after pairs; no object/cycle dropped for missing files.
    cycles=[];maincost=[];mainq=read(root/'queues'/MAIN/'RECEIPT.json',{});mainjobs={x['task_id']:x for x in mainq.get('jobs',[])}
    for cid in CASES:
        for method in METHODS:
            task=cid+'__'+method;folder=root/'RESULTS'/MAIN/task;r=read(folder/'RECEIPT.json',{});job=mainjobs.get(task,{})
            maincost.append(dict(case_id=cid,method=method,status=r.get('status','MISSING'),external_child_CPU_seconds=(job.get('external') or {}).get('child_process_CPU_seconds'),
                internal_CPU_seconds=r.get('whole_process_CPU_seconds',r.get('CPU_seconds')),physical_cost=r.get('physical_cost'),auxiliary_cost=r.get('auxiliary_cost'),projector_cost=r.get('fixed_projection_cost'),projector_source_field='fixed_projection_cost',scope='Reused original full-child costs, not added to cycle-window differences'))
            for index in range(1,5):
                bp=folder/f'CYCLE_{index:03d}_BEFORE.npz';ap=folder/f'CYCLE_{index:03d}_AFTER.npz';before=state(bp);after=state(ap)
                row=dict(case_id=cid,method=method,cycle=index,before_path=str(bp),after_path=str(ap),available=before is not None and after is not None,run_status=r.get('status','MISSING'),geometry_source=geometry.get(cid,{}).get('execution'))
                if not row['available']:row['reason']='MISSING_OR_INVALID_BEFORE_AFTER_COMPLETE_STATE';cycles.append(row);continue
                try:
                    x,y=view(after,True),view(before,True);event=next((e for e in reversed(after['r15']['cycle_history']) if e.get('cycle')==index),{})
                    row.update(compare(x,y),**sectors(x['m']-y['m'],cid),before_phase=phase(y),after_phase=phase(x),
                        PR_reset_observed=after['old_grad'] is None and after['old_direction'] is None,
                        PR_present_before=before['old_grad'] is not None,before_normal_memory_norm=float(np.hypot(nrm(y['memory_X']),nrm(y['memory_Q']))),after_normal_memory_norm=float(np.hypot(nrm(x['memory_X']),nrm(x['memory_Q']))),
                        pending_before=before['r15']['pending_cycle'] is not None,pending_after=after['r15']['pending_cycle'] is not None,
                        auxiliary_window_cost=fresh(after['operation_counts'],before['operation_counts']),projector_window_cost=fresh(after['r15']['projector_cost'],before['r15']['projector_cost']),event=event,
                        current_projection_relative_loss=event.get('current_projection_relative_loss'),full_shadow_gap=event.get('full_shadow_gap'),routed_shadow_gap=event.get('routed_shadow_gap'),
                        scope='Saved completed before/after cycle window; prior setup and outside-window/checkpoint I/O remain in child CPU, not free')
                except Exception as exc:row.update(available=False,reason='STATE_SCHEMA_OR_SHAPE_ERROR',error=repr(exc));issues.append(dict(case_id=cid,method=method,cycle=index,issue='CYCLE_ERROR',error=repr(exc)))
                cycles.append(row)
    # Fixed56 fork inventory is independent of receipt availability.
    specs=[(f'r14_new_{i:04d}',m,'base',1) for i in (9,20,21,23) for m in METHODS]
    specs += [(f'r14_new_{i:04d}',m,b,s) for i in (21,23) for m in ('twofold_hop','twofold_dr','twofold_reset') for b in ('U11','U01','U00') for s in (-1,1)]
    q=read(root/'queues'/FORK/'RECEIPT.json',{});jobs={x['task_id']:x for x in q.get('jobs',[])};forkcost=[];phis=[];native=[];points=[];gains=[];base_initial={};state16={}
    for cid in [f'r14_new_{i:04d}' for i in (9,20,21,23)]:
        for method in METHODS:
            folder=root/'RESULTS'/FORK/(cid+'__'+method+'__base');base_initial[cid,method]=view(state(folder/'FORK_INITIAL.npz'),True);state16[cid,method]=view(arrays(folder/'STATE_016.npz'))
    initial_comparisons=[]
    for (cid,method),initial in base_initial.items():
        initial_comparisons.append(dict(case_id=cid,method=method,reference='twofold_hop',comparison=compare(initial,base_initial[cid,'twofold_hop']),scope='Actual complete fork initial arrays/PR history/reflection memory; method labels intentionally differ'))
    for cid,method,band,sign in specs:
        ident=cid+'__'+method+'__'+band+(str(sign) if band!='base' else '');folder=root/'RESULTS'/FORK/ident;r=read(folder/'RECEIPT.json',{});job=jobs.get(ident,{});external=(job.get('external') or {}).get('child_process_CPU_seconds')
        cost=dict(task_id=ident,case_id=cid,method=method,band=band,sign=sign,run_status=r.get('status','MISSING'),queue_status=job.get('status','MISSING'),completed_updates=r.get('completed_updates'),
            external_child_CPU_seconds=external,external_CPU_null_reason=None if finite(external) else 'NO_EXTERNAL_WHOLE_CHILD_CPU',internal_whole_CPU_seconds=r.get('whole_process_CPU_seconds'),
            physical_cost=r.get('physics_cost'),auxiliary_inherited=r.get('auxiliary_inherited'),auxiliary_cumulative=r.get('auxiliary_cumulative'),auxiliary_new=fresh(r.get('auxiliary_cumulative'),r.get('auxiliary_inherited')),
            projector_inherited=r.get('fixed_projection_inherited'),projector_cumulative=r.get('fixed_projection_cumulative'),projector_new=fresh(r.get('fixed_projection_cumulative'),r.get('fixed_projection_inherited')),
            source_checkpoint=r.get('source_checkpoint'),queue_source_selection=job.get('source_checkpoint'),failure=r.get('error'),cost_scope='Fresh external child CPU; inherited counters subtracted, physical child-model costs fresh, background RHS separate; internal CPU nested, never summed')
        forkcost.append(cost);hist=read(folder/'HISTORY.json',[]);nh=read(folder/'NATIVE_HISTORY.json',[])
        if not hist:phis.append(dict(task_id=ident,case_id=cid,method=method,band=band,sign=sign,available=False,reason='NO_PAID_PHI_HISTORY',run_status=cost['run_status']))
        for h in hist:phis.append(dict(h,task_id=ident,case_id=cid,method=method,band=band,sign=sign,available=True,run_status=cost['run_status']))
        for h in nh:native.append(dict(h,task_id=ident,case_id=cid,method=method,band=band,sign=sign,run_status=cost['run_status']))
        if band=='base':
            initial=base_initial[cid,method]
            for index in (1,16,64):
                s=view(arrays(folder/f'STATE_{index:03d}.npz'));hop=state16.get((cid,'twofold_hop')) if index==16 else view(arrays(root/'RESULTS'/FORK/(cid+'__twofold_hop__base')/f'STATE_{index:03d}.npz'))
                points.append(dict(case_id=cid,method=method,local_update=index,run_status=cost['run_status'],distance_from_shared_initial=compare(s,initial),distance_from_hop_at_same_local_update=compare(s,hop),phase=phase(s),initial_phase=phase(initial),state_available=s is not None))
        else:
            perturb=arrays(folder/'PERTURBATION.npz');initial=view(state(folder/'FORK_INITIAL.npz'),True);s=view(arrays(folder/'STATE_016.npz'));base=state16[cid,method];basepath=root/'RESULTS'/FORK/(cid+'__'+method+'__base');br=read(basepath/'RECEIPT.json',{})
            row=dict(task_id=ident,case_id=cid,method=method,band=band,sign=sign,run_status=cost['run_status'],base_run_status=br.get('status','MISSING'),source_checkpoint_equal=bool(r.get('source_checkpoint') and r.get('source_checkpoint')==br.get('source_checkpoint')),
                before_phase=phase(initial),after_phase=phase(s),base_after_phase=phase(base),base_initial_phase=phase(base_initial[cid,method]))
            if perturb is None or not {'requested','realized','unperturbed_material'}<=set(perturb) or not all(np.isfinite(v).all() for v in perturb.values()):row.update(available=False,reason='MISSING_OR_INVALID_SAVED_REALIZED_PERTURBATION')
            else:
                eps=nrm(perturb['realized']);requested=perturb['requested'];realized=perturb['realized'];row.update(compare(s,base,eps),requested_perturbation_norm=nrm(requested),clipped_coordinates=int(np.count_nonzero(requested!=realized)),
                    initial_difference=compare(initial,base_initial[cid,method]),realized_material_perturbation_sectors=sectors(realized,cid),
                    perturbation_reconstruction_error=nrm(initial['m']-perturb['unperturbed_material']-realized) if initial else None,
                    no_local_spectrum_claim=True,stage_warning='Different stage/beta at the16-step endpoints; not a local spectrum' if s is not None and base is not None and (s['stage']!=base['stage'] or not np.array_equal(s['beta'],base['beta'])) else 'Finite response, not a local spectrum even when stage/beta match')
            base_native=read(basepath/'NATIVE_HISTORY.json',[])
            row['perturbed_stage_beta_path']=[dict(local_update=x.get('local_update'),stage=x.get('stage'),beta_current=x.get('beta_current'),beta_material=x.get('beta_material')) for x in nh]
            row['base_stage_beta_path']=[dict(local_update=x.get('local_update'),stage=x.get('stage'),beta_current=x.get('beta_current'),beta_material=x.get('beta_material')) for x in base_native if x.get('local_update',0)<=16]
            row['within_window_stage_or_beta_change']=any(len({(x.get('stage'),x.get('beta_current'),x.get('beta_material')) for x in path})>1 for path in (row['perturbed_stage_beta_path'],row['base_stage_beta_path']))
            gains.append(row)
    def encode(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        raise TypeError(type(x).__name__)
    def output(name,rows):
        (out/(name+'.json')).write_text(json.dumps(rows,default=encode,indent=2,allow_nan=False)+'\n');keys=list(dict.fromkeys(k for r in rows for k in r))
        with (out/(name+'.csv')).open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows({k:json.dumps(v,default=encode) if isinstance(v,(dict,list,np.ndarray)) else v for k,v in r.items()} for r in rows)
    for name,data in [('MAIN_CYCLES',cycles),('MAIN_REUSED_COSTS',maincost),('FORK_COSTS',forkcost),('PAID_PHI_TRAJECTORIES',phis),('FORK_NATIVE_DYNAMICS',native),('BASE_STATE_DISTANCES',points),('BASE_INITIAL_COMPARISONS',initial_comparisons),('FINITE_RESPONSES',gains)]:output(name,data)
    # Saved scientific plots only, no interpolation/forward/success labels.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figures=out/'figures';figures.mkdir()
    fig,axes=plt.subplots(4,3,figsize=(18,16))
    for ri,index in enumerate((9,20,21,23)):
        cid=f'r14_new_{index:04d}'
        for method in METHODS:
            h=[x for x in phis if x['case_id']==cid and x['method']==method and x['band']=='base' and x['available'] and finite(x.get('Phi'))]
            axes[ri,0].plot([x['CPU_seconds'] for x in h],[x['Phi'] for x in h],'.-',label=method,color=plt.get_cmap('tab10')(METHODS.index(method)))
            n=[x for x in native if x['case_id']==cid and x['method']==method and x['band']=='base' and x.get('cycle_metadata')]
            axes[ri,1].plot([x['local_update'] for x in n],[x['cycle_metadata'].get('full_shadow_gap',np.nan) for x in n],'.-',label=method+' full',color=plt.get_cmap('tab10')(METHODS.index(method)))
            axes[ri,1].plot([x['local_update'] for x in n],[x['cycle_metadata'].get('routed_shadow_gap',np.nan) for x in n],':',alpha=.7,color=plt.get_cmap('tab10')(METHODS.index(method)))
            d=[x for x in points if x['case_id']==cid and x['method']==method and x['distance_from_hop_at_same_local_update'].get('available')]
            axes[ri,2].plot([x['local_update'] for x in d],[x['distance_from_hop_at_same_local_update']['delta_m_norm'] for x in d],'.-',label=method)
        for ci,title in enumerate(('Paid original Phi vs native CPU','Full(solid)/routed(dotted) gap vs update','Material distance from same-time hop')):axes[ri,ci].set_title(cid+' '+title,fontsize=9)
    axes[0,0].legend(fontsize=7);fig.tight_layout();fig.savefig(figures/'samecheckpoint_phi_gap_material.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(18,6))
    for i,row in enumerate(gains):
        if not row.get('available'):continue
        label=row['case_id'][-4:]+'/'+row['method'].removeprefix('twofold_')+'/'+row['band']+'/'+str(row['sign']);color='tab:blue' if row.get('same_stage_and_beta') and not row.get('within_window_stage_or_beta_change') else 'tab:red'
        for ax,key in zip(axes[:2],('material_finite_gain','scaled_field_finite_gain')):
            if finite(row.get(key)):ax.scatter(i,row[key],color=color)
        axes[2].scatter(i,row.get('clipped_coordinates',0),color=color)
    for ax,title in zip(axes,('16-update material finite gain','16-update scaled field finite gain','Actually clipped coordinates')):ax.set_title(title);ax.set_xlabel('Blue=same stage/beta throughout; red=endpoint or within-window change');ax.grid(alpha=.2)
    ticklabels=[r['case_id'][-4:]+'/'+r['method'].removeprefix('twofold_')+'/'+r['band']+'/'+str(r['sign']) for r in gains]
    for ax in axes:ax.set_xticks(range(len(gains)),ticklabels,rotation=90,fontsize=6)
    fig.suptitle('Finite saved-state responses; no spectral/chaos claim');fig.tight_layout();fig.savefig(figures/'finite_response_gains_box_stage.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(18,6))
    for mi,method in enumerate(METHODS):
        group=[x for x in cycles if x['method']==method and x.get('available')]
        for x in group:
            if finite(x.get('material_outside_real_U01_fraction')):axes[0].scatter(mi,x['material_outside_real_U01_fraction'],s=12,alpha=.5)
            if finite(x.get('current_projection_relative_loss')):axes[1].scatter(mi,x['current_projection_relative_loss'],s=12,alpha=.5)
        for bi,band in enumerate(('U11','U10','U01','U00')):
            vs=[x.get(band+'_material_fraction') for x in group if finite(x.get(band+'_material_fraction'))]
            if vs:axes[2].scatter([mi+.12*(bi-1.5)]*len(vs),vs,s=10,alpha=.4,label=band if mi==0 else None,color=plt.get_cmap('tab10')(bi))
    for ax in axes:ax.set_xticks(range(5),METHODS,rotation=30,ha='right')
    axes[0].set_title('Cycle Δm energy outside real U01');axes[1].set_title('Saved current injection projection loss');axes[2].set_title('Cycle Δm energy fractions by real sector');axes[2].legend();fig.tight_layout();fig.savefig(figures/'cycle_material_bands_projection_loss.png',dpi=150);plt.close(fig)
    output('ISSUES',issues);values=[x['external_child_CPU_seconds'] for x in forkcost]
    receipt=dict(status='COMPLETE_WITH_RETAINED_MISSING_OR_ISSUES' if issues or any(not x['available'] for x in cycles) or any(not x.get('available') for x in gains) else 'COMPLETE',created_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',physics_calls=0,
        expected_main_cycle_pairs=480,available_main_cycle_pairs=sum(x['available'] for x in cycles),expected_forks=56,received_fork_receipts=sum(x['run_status']!='MISSING' for x in forkcost),expected_perturbed_responses=36,available_perturbed_responses=sum(x.get('available') is True for x in gains),
        external_fork_CPU_total=sum(values) if all(finite(v) for v in values) else None,measured_external_fork_CPU_subtotal=sum(v for v in values if finite(v)),missing_external_CPU_tasks=[x['task_id'] for x in forkcost if not finite(x['external_child_CPU_seconds'])],
        queue_wall_seconds=q.get('queue_wall_seconds'),queue_controller_CPU_seconds=q.get('controller_CPU_seconds'),queue_status=q.get('status'),
        analysis_CPU_seconds=time.process_time()-cpu,analysis_wall_seconds=time.perf_counter()-start,state_codec=str(source),read_paths=sorted(read_paths),issues=len(issues),
        scope='No-GT reused paid state analysis; fixed finite responses and real material sector labels; stages/betas retained; no cycles/chaos/spectral or recovery conclusion')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    (out/'REPORT.md').write_text('# R15 full-state saved dynamics (no GT)\n\n'+json.dumps({k:v for k,v in receipt.items() if k!='read_paths'},indent=2)+'\n\nMAIN_CYCLES keeps every24×5×4 expected pair, actual real-sector Δm fractions, PR reset, current/field/reflection changes, saved gaps/injection loss and counter windows. A material update is not constrained to current-routing H.\n\nFINITE_RESPONSES compares perturbed STATE016 to the same-method unperturbed base STATE016 and divides by saved realized clipped input norm. These are finite material/field response gains, not local eigenvalues or spectral radii; differing endpoint stage/beta are explicitly flagged. Coefficient comparisons use the saved native coordinate frame. All scalar stages/betas and initial-state comparisons remain in tables.\n\nActual external fork CPU is fresh cost; cumulative inherited checkpoint counters are subtracted. Main full-child costs are reused separately. Cycle windows exclude earlier projector setup; all setup/serialization remains charged to whole child. No child wall sums, nested CPU sums or counter double charges are claimed.\n\nMissing/failed states remain explicit. No optimizer, forward, truth array, scientific success label, limit-cycle or chaos interpretation is produced; final interpretation belongs to parent.\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='read_paths'}))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        if _ACTIVE_OUT is not None:
            failed=_ACTIVE_OUT/'FAILED_RECEIPT.json'
            if not failed.exists():failed.write_text(json.dumps(dict(status='ANALYSIS_FAILED',error=repr(exc),whole_process_CPU_seconds=time.process_time(),wall_since_script_import=time.perf_counter()-_START_WALL,truth_access='NONE',physics_calls=0),indent=2)+'\n')
        raise
