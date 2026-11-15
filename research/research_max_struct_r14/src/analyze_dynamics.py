"""No-GT mechanical aggregation of already-paid native/intervention diagnostics.

No Maxwell model, optimizer, NN, material truth, success labels, or numerical
NPZ payload is imported/read. Every run requires a new output directory.
"""
from __future__ import annotations
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import csv
import html
import json
import math
from pathlib import Path
import time
import zipfile

EPS=2.220446049250313e-16
COST_KEYS=('PR_CG_steps','material_inner_evaluations','aplus_sweeps','Green_forward_calls',
    'Green_adjoint_calls','coefficient_FFT_calls','coefficient_IFFT_calls',
    'onehop_block_matrix_apply_calls','source_metric_products')
SOURCES=[]
ISSUES=[]


def finite(v):return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)


def read(path,required=False):
    path=Path(path).resolve()
    if any(p.lower() in ('sealed_truth','truth','data_private','private_operational') for p in path.parts):
        raise ValueError('Private input forbidden: '+str(path))
    if not path.is_file():
        if required:raise FileNotFoundError(path)
        ISSUES.append(dict(path=str(path),issue='MISSING_JSON_INPUT'));return None
    before=path.stat()
    try:value=json.loads(path.read_text(encoding='utf-8'))
    except Exception as exc:
        ISSUES.append(dict(path=str(path),issue='UNREADABLE_JSON_INPUT',error=repr(exc)));return None
    after=path.stat();SOURCES.append(dict(path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,
        unchanged_during_read=(before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)))
    if SOURCES[-1]['unchanged_during_read'] is False:ISSUES.append(dict(path=str(path),issue='INPUT_CHANGED_DURING_READ'))
    return value


def write_json(path,value):
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def table(out,name,rows):
    write_json(out/(name+'.json'),rows);keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+'.csv')).open('x',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
        for r in rows:writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def add(*values):return sum(values) if all(finite(v) for v in values) else None


def subtract(a,b):return a-b if finite(a) and finite(b) else None


def scalar_range(v):
    if finite(v):return v,v
    if isinstance(v,list):
        values=[]
        for x in v:
            lo,hi=scalar_range(x)
            if lo is not None:values.extend([lo,hi])
        if values:return min(values),max(values)
    return None,None


def numerical_floor(*values):
    numbers=[abs(v) for v in values if finite(v)]
    return 128*EPS*max(numbers+[1e-300])


def summarize_steps(rows,identity):
    output=dict(identity,native_recorded_steps=len(rows),
        valid_F_steps=sum(finite(r.get('delta_F')) for r in rows),
        source_scope='Recorded completed native updates; iterations are dependent within trajectories,not independent samples')
    for field in ('delta_F','current_step_objective_change','material_current_potential_work'):
        values=[r[field] for r in rows if finite(r.get(field))]
        output[field+'_sum']=sum(values) if values else None
        output[field+'_positive_fraction']=sum(v>0 for v in values)/len(values) if values else None
        output[field+'_max_positive']=max([0.]+values) if values else None
        output[field+'_above_scale_floor_count']=sum(r[field]>r['scale_machine_floor'] for r in rows if finite(r.get(field)))
        output[field+'_positive_sum']=sum(v for v in values if v>0) if values else None
    output.update(restarts=sum(r.get('direction_restarted') is True for r in rows),
        noncontiguous_rows=sum(r.get('previous_iteration_contiguous') is False for r in rows),
        stage_entry_rows=sum(r.get('stage_entry') is True for r in rows),
        max_scale_machine_floor=max((r['scale_machine_floor'] for r in rows),default=None),
        floor_scope='128*binary64eps*max(abs(before,after,aftercurrent),tiny); descriptive roundoff scale,not error certificate or quality gate')
    return output


def native(root,campaigns):
    steps=[];trajectories=[];stage_summaries=[];configs=[];counter_rows=[];checkpoint_rows=[]
    for campaign in campaigns:
        cfg=read(root/'configs'/(campaign+'.json'),required=True)
        for task in cfg['tasks']:
            folder=root/'RESULTS'/campaign/task['task_id'];path=folder/'NATIVE_HISTORY.json'
            receipt=read(folder/'RECEIPT.json');history=read(path)
            identity=dict(campaign=campaign,task_id=task['task_id'],case_id=task['case_id'],pair=str(task.get('pair',task.get('method'))))
            actual=read(folder/'ACTUAL_CONFIG.json') if (folder/'ACTUAL_CONFIG.json').is_file() else None
            configs.append(dict(identity,source=str(folder/'ACTUAL_CONFIG.json'),actual_config=actual,
                manifest_protocol=cfg.get('protocol'),task=task))
            if not isinstance(history,list):
                trajectories.append(dict(identity,status='NO_NATIVE_HISTORY',run_status=(receipt or {}).get('status'),recorded_updates=0));continue
            previous=None;local=[];stage_groups=defaultdict(list)
            for index,row in enumerate(history):
                if 'global_iteration' not in row or 'data_before_current' not in row:continue
                previous_prior=previous.get('regularizer') if previous is not None else 0. if row.get('global_iteration')==1 else None
                before=add(row.get('data_before_current'),row.get('state_before_current'),previous_prior)
                mid=add(row.get('data_after_current'),row.get('state_after_current'),previous_prior)
                after=add(row.get('data_term'),row.get('state_term_current'),row.get('regularizer'))
                delta=subtract(after,before);work=subtract(after,mid)
                pr_lo,pr_hi=scalar_range(row.get('PR_plus'));alpha_lo,alpha_hi=scalar_range(row.get('line_alpha'))
                counter=row.get('operation_counts',{});prevcount=previous.get('operation_counts',{}) if previous else {}
                entry=previous is None or row.get('stage')!=previous.get('stage')
                r=dict(identity,source_row=index,global_iteration=row.get('global_iteration'),stage=row.get('stage'),
                    iteration=row.get('iteration'),stage_entry=entry,
                    previous_iteration_contiguous=previous is None and row.get('global_iteration')==1 or previous is not None and row['global_iteration']==previous['global_iteration']+1,
                    beta_current=row.get('beta_current'),beta_material=row.get('beta_material'),
                    stage_beta_current=row.get('stage_beta_current'),stage_beta_material=row.get('stage_beta_material'),mf=row.get('mf'),
                    previous_regularizer=previous_prior,F_before=before,F_after_current=mid,F_after=after,delta_F=delta,
                    current_step_objective_change=row.get('current_step_objective_change'),
                    material_current_potential_work=work,scale_machine_floor=numerical_floor(before,mid,after),
                    decomposition_roundoff=subtract(delta,add(row.get('current_step_objective_change'),work)),
                    before_minus_previous_after=subtract(before,previous.get('_F_after')) if previous else None,
                    stage_change_scope='Stage-entry before terms already use new metric/basis;prior carried acrossstage;not one global potential',
                    same_current_material_potential=row.get('pair',identity['pair']) in ('00','11'),
                    PR_min=pr_lo,PR_max=pr_hi,line_alpha_min=alpha_lo,line_alpha_max=alpha_hi,
                    direction_restarted=row.get('direction_restarted'),readout_material_change=row.get('readout_material_change'),
                    material_inner_nit=row.get('material_inner_nit'),material_inner_success=row.get('material_inner_success'),
                    material_fallback=row.get('material_fallback'),aplus_sweeps=row.get('aplus_sweeps'),
                    stage_end=row.get('stage_end'),stage_stop=row.get('stage_stop'),
                    native_CPU_seconds=row.get('CPU_seconds'),native_wall_seconds=row.get('wall'),
                    first_operation_count_scope='First recorded counters include setup;per-step delta unknown' if previous is None else None)
                for key in COST_KEYS:
                    r['cumulative_'+key]=counter.get(key)
                    r['delta_'+key]=subtract(counter.get(key),prevcount.get(key)) if previous else None
                steps.append(r);local.append(r);stage_groups[row.get('stage')].append(r)
                previous=dict(row,_F_after=after)
            final_counts=(receipt or {}).get('auxiliary_cost',{}) or {}
            last_counts=previous.get('operation_counts',{}) if previous else {}
            final_last_delta={k:subtract(final_counts.get(k),last_counts.get(k)) for k in set(final_counts)|set(last_counts)}
            checkpoint_files=sorted(folder.glob('CHECKPOINT_*.npz'))
            for cp in checkpoint_files:
                # List names only;no array payload/model/state contraction.
                try:
                    with zipfile.ZipFile(cp) as z:names=z.namelist()
                    checkpoint_rows.append(dict(identity,path=str(cp),bytes=cp.stat().st_size,members=names,
                        numerical_payload_read=False,scope='Public-derived checkpoint inventory only'))
                except Exception as exc:ISSUES.append(dict(path=str(cp),issue='CHECKPOINT_INVENTORY_ERROR',error=repr(exc)))
            counter_rows.append(dict(identity,last_recorded_native_counts=last_counts,actual_final_auxiliary_counts=final_counts,
                final_minus_last_recorded_counts=final_last_delta,physical_cost=(receipt or {}).get('physical_cost'),
                scope='Final receipt includes already-paid partial/overrun work;history only completed returned updates;missing counters null'))
            summary=summarize_steps(local,identity)
            summary.update(run_status=(receipt or {}).get('status'),recorded_updates=len(local),
                last_recorded_global_iteration=previous.get('global_iteration') if previous else None,
                final_receipt_PR_steps=final_counts.get('PR_CG_steps'),
                completed_but_unrecorded_PR_steps=subtract(final_counts.get('PR_CG_steps'),len(local)) if 'PR_CG_steps' in final_counts else None,
                all_paid_auxiliary_counts=final_counts,termination=(receipt or {}).get('termination_reason'),
                actual_final_CPU_seconds=(receipt or {}).get('CPU_seconds'),
                last_native_CPU_seconds=previous.get('CPU_seconds') if previous else None,
                missing_last_overrun_scope='An interrupted step may lack a returned native row;final-minus-last counters preserved;absence is not a completed-step or success claim')
            trajectories.append(summary)
            for stage,values in stage_groups.items():stage_summaries.append(summarize_steps(values,dict(identity,stage=stage)))
    return dict(NATIVE_STEPS=steps,TRAJECTORIES=trajectories,NATIVE_STAGE_SUMMARY=stage_summaries,
        SOURCE_CONFIGURATIONS=configs,ACTUAL_COST_COUNTERS=counter_rows,CHECKPOINT_INVENTORY=checkpoint_rows)


def flatten(value,prefix=''):
    result={}
    if isinstance(value,dict):
        for key,item in value.items():result.update(flatten(item,prefix+key+'.'))
    elif isinstance(value,list):result[prefix[:-1]]=value
    else:result[prefix[:-1]]=value
    return result


def interventions(root):
    snapshots=[];branch_steps=[];branches=[];responses=[];comparisons=[];costs=[];receipts=[];configs=[];repeats=[];failures=[];summaries=[]
    parent=root/'audit/interventions'
    for folder in sorted(p for p in parent.iterdir() if p.is_dir()):
        if not (folder/'CONFIGURATION.json').is_file():continue
        conf=read(folder/'CONFIGURATION.json');receipt=read(folder/'RECEIPT.json') or {}
        case_id=folder.name.split('_k')[0];identity=dict(intervention=folder.name,case_id=case_id,
            source_global_iteration=(conf or {}).get('checkpoint_global_iteration'),source_stage=(conf or {}).get('checkpoint_stage'))
        configs.append(dict(identity,configuration=conf));receipts.append(dict(identity,**receipt))
        branchlist=read(folder/'BRANCHES.json') or []
        discovered=sorted((folder/'branches').glob('*/BRANCH.json'))
        if len(discovered)!=8:ISSUES.append(dict(intervention=folder.name,issue='BRANCH_COUNT_NOT8',actual=len(discovered)))
        if len(branchlist)!=len(discovered):ISSUES.append(dict(intervention=folder.name,issue='BRANCH_INDEX_VS_FILES_COUNT_MISMATCH',index=len(branchlist),files=len(discovered)))
        for path in discovered:
            branch=read(path) or {};initial=branch.get('initial_snapshot',{});name=branch.get('name',path.parent.name)
            branchid=dict(identity,branch=name,pair=branch.get('pair'),reset_PR=branch.get('reset_PR'))
            branches.append(dict(branchid,status=branch.get('status'),completed_updates=branch.get('completed_updates'),
                declared_update_budget=branch.get('update_budget'),physical_Phi_change=branch.get('physical_Phi_change'),
                same_state_scope='Same saved full m/A/J/Q/DPCS/source weights/stage;reset only clears old gradient/direction',
                comparison_scope='Fixed1/8/32 updates,not equalCPU;all paid costs separate'))
            for step in branch.get('steps',[]):branch_steps.append(dict(branchid,**step))
            actual_files={int(p.name.split('_')[1]):p for p in path.parent.glob('step_*_physical.json')}
            for count in (0,1,8,32):
                source=actual_files.get(count)
                snap=read(source) if source else None
                if snap is None:
                    snapshots.append(dict(branchid,relative_update=count,available=False));continue
                aux=add(snap.get('auxiliary_data_source_metric_loss'),snap.get('auxiliary_state_current_metric_loss'),snap.get('R'))
                initial_aux=add(initial.get('auxiliary_data_source_metric_loss'),initial.get('auxiliary_state_current_metric_loss'),initial.get('R'))
                compact={k:v for k,v in snap.items() if k not in ('physical_before','arrays')}
                compact.update(dict(branchid,relative_update=count,available=True,source=str(source),
                    delta_physical_Phi_from_branch_start=subtract(snap.get('Phi'),initial.get('Phi')),
                    delta_physical_D_from_branch_start=subtract(snap.get('D'),initial.get('D')),
                    delta_physical_R_from_branch_start=subtract(snap.get('R'),initial.get('R')),
                    current_metric_auxiliary_potential=aux,delta_current_metric_auxiliary_potential=subtract(aux,initial_aux),
                    delta_constitutive_defect_norm=subtract(snap.get('auxiliary_constitutive_defect_norm'),initial.get('auxiliary_constitutive_defect_norm')),
                    last_material_step_delta_physical_Phi=subtract(snap.get('Phi'),(snap.get('physical_before') or {}).get('Phi')),
                    raw_constitutive_defect_vector_angle_degrees=None,
                    raw_defect_angle_reason='Not present in paid source;not recomputed',
                    angle_scope='Available material-step angle uses paid previous-m originalPhi gradient;not a defect-vector angle'))
                snapshots.append(compact)
        repeat=read(folder/'repeat_identical_step/CONTROL.json')
        if repeat is not None:repeats.append(dict(identity,**repeat))
        summary=read(folder/'perturbations/SUMMARY.json') or {}
        summaries.append(dict(identity,status=summary.get('status'),direction_info=summary.get('direction_info'),interpretation=summary.get('interpretation')))
        for row in summary.get('responses',[]):responses.append(dict(identity,**flatten(row)))
        for row in summary.get('comparisons',[]):comparisons.append(dict(identity,**flatten(row)))
        ledger=read(folder/'COST_LEDGER.json') or []
        for row in ledger:costs.append(dict(identity,**row))
        fs=read(folder/'FAILURES.json')
        if isinstance(fs,list):
            for row in fs:failures.append(dict(identity,failure=row))
        elif fs:failures.append(dict(identity,failure=fs))
    batch=read(parent/'BATCH_RECEIPT.json')
    return dict(INTERVENTION_BRANCHES=branches,INTERVENTION_SNAPSHOTS=snapshots,INTERVENTION_NATIVE_STEPS=branch_steps,
        FINITE_JVP_RESPONSES=responses,FINITE_JVP_COMPARISONS=comparisons,INTERVENTION_COST_EVENTS=costs,
        INTERVENTION_RECEIPTS=receipts,INTERVENTION_CONFIGURATIONS=configs,REPEATED_STEP_CONTROLS=repeats,
        INTERVENTION_FAILURES=failures,PERTURBATION_SCOPE=summaries,INTERVENTION_BATCH=[batch] if batch else [])


def plots(out,data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    paths=[];dest=out/'FIGURES';dest.mkdir()
    trajectories=data['TRAJECTORIES'];steps=data['NATIVE_STEPS']
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for pi,pair in enumerate(('00','01','10','11')):
        group=[r for r in trajectories if r.get('pair')==pair]
        for ax,key in zip(axes,['delta_F_max_positive','material_current_potential_work_positive_fraction','material_current_potential_work_sum']):
            values=[r[key] for r in group if finite(r.get(key))]
            ax.scatter([pi]*len(values),values,s=22,alpha=.65);ax.set_xticks(range(4),['00','01','10','11']);ax.set_ylabel(key)
    axes[0].set_yscale('symlog',linthresh=1e-13);fig.suptitle('Each dot is one recorded trajectory;event counts are not independent samples');fig.tight_layout()
    path=dest/'native_mechanical_summary.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    case_ids=sorted({r['case_id'] for r in trajectories})
    for offset in range(0,len(case_ids),4):
        subset=case_ids[offset:offset+4];fig,axes=plt.subplots(len(subset),2,figsize=(13,3*len(subset)),squeeze=False)
        for ri,cid in enumerate(subset):
            for pair in ('00','01','10','11'):
                rows=[r for r in steps if r['case_id']==cid and r['pair']==pair]
                axes[ri,0].plot([r['global_iteration'] for r in rows],[r['F_after'] for r in rows],label=pair,lw=.9)
                # Cumulative work displayed within each stage;never telescope differing stage potentials.
                acc=0.;stage=None;xs=[];ys=[]
                for r in rows:
                    if stage!=r['stage']:acc=0.;stage=r['stage']
                    if finite(r['material_current_potential_work']):acc+=r['material_current_potential_work']
                    xs.append(r['global_iteration']);ys.append(acc)
                axes[ri,1].plot(xs,ys,label=pair,lw=.9)
            axes[ri,0].set_yscale('log');axes[ri,0].set_ylabel(cid+'\nF_after');axes[ri,1].set_ylabel('stage-local cumulative material work')
            for ax in axes[ri]:ax.set_xlabel('Recorded global update');ax.grid(alpha=.2)
        axes[0,0].legend();fig.suptitle('Native traces;stage changes can change potential definition');fig.tight_layout()
        path=dest/f'native_traces_{offset//4+1:02d}.png';fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path.relative_to(out)))
    snaps=data['INTERVENTION_SNAPSHOTS'];names=sorted({r['intervention'] for r in snaps})
    for offset in range(0,len(names),4):
        subset=names[offset:offset+4];fig,axes=plt.subplots(len(subset),2,figsize=(13,3.2*len(subset)),squeeze=False)
        for ri,name in enumerate(subset):
            for branch in sorted({r['branch'] for r in snaps if r['intervention']==name}):
                group=sorted([r for r in snaps if r['intervention']==name and r['branch']==branch and r.get('available')],key=lambda r:r['relative_update'])
                for ax,key in zip(axes[ri],['delta_physical_Phi_from_branch_start','delta_current_metric_auxiliary_potential']):
                    ax.plot([r['relative_update'] for r in group],[r.get(key) for r in group],marker='o',label=branch,lw=.9);ax.axhline(0,color='k',lw=.5);ax.set_ylabel(name+'\n'+key,fontsize=8);ax.set_xlabel('Fixed native updates')
            axes[ri,0].set_xticks([0,1,8,32]);axes[ri,1].set_xticks([0,1,8,32])
        axes[0,1].legend(fontsize=6,ncol=2);fig.suptitle('Already-paid same-state interventions;not equalCPU/not GT recovery');fig.tight_layout()
        path=dest/f'intervention_deltas_{offset//4+1:02d}.png';fig.savefig(path,dpi=150);plt.close(fig);paths.append(str(path.relative_to(out)))
    comparisons=data['FINITE_JVP_COMPARISONS'];fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,kind,field in [(axes[0],'central_JVP_epsilon_convergence','relative_difference.m'),(axes[1],'directional_additivity','relative_defect.m')]:
        for pi,pair in enumerate(('01','11')):
            values=[r[field] for r in comparisons if r.get('kind')==kind and r.get('pair')==pair and finite(r.get(field))]
            ax.scatter([pi]*len(values),values,s=20,alpha=.6)
        ax.set_xticks([0,1],['01','11']);ax.set_ylabel(field);ax.set_title(kind);ax.set_yscale('symlog',linthresh=1e-12)
    fig.suptitle('Existing finite-response diagnostics only;branch flags/actual epsilon retained in tables');fig.tight_layout()
    path=dest/'finite_JVP_existing_diagnostics.png';fig.savefig(path,dpi=160);plt.close(fig);paths.append(str(path.relative_to(out)))
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--campaigns',default='dualmetric_pilot12_v1')
    args=parser.parse_args();root=args.root.resolve();out=args.out.resolve();campaigns=[c.strip() for c in args.campaigns.split(',') if c.strip()]
    if any(c not in ('dualmetric_pilot12_v1','dualmetric_remaining36_v1') for c in campaigns):raise ValueError('Only declared main-quartet native histories supported')
    if 'dualmetric_pilot12_v1' in campaigns and not (root/'audit/PILOT_DOWNLOAD_EXTRACT.json').is_file():
        raise RuntimeError('Pilot extraction marker absent;no native analysis before complete download handoff')
    if out.exists() or out==root or any(out.is_relative_to(root/p) for p in ('RESULTS','sealed_truth','data_public','src','configs','audit/interventions')):
        raise ValueError('New separate analysis output directory required')
    cpu0,wall0=time.process_time(),time.perf_counter();started=datetime.now(timezone.utc).isoformat()
    out.mkdir(parents=True,exist_ok=False)
    marker=read(root/'audit/PILOT_DOWNLOAD_EXTRACT.json') if 'dualmetric_pilot12_v1' in campaigns else None
    data=native(root,campaigns);data.update(interventions(root))
    for name,rows in data.items():table(out,name,rows)
    table(out,'ISSUES',ISSUES);table(out,'SOURCE_INPUTS',SOURCES)
    figures=plots(out,data)
    complete=sum(r.get('run_status')=='COMPLETE' for r in data['TRAJECTORIES'])
    receipt=dict(status='NO_GT_MECHANICAL_AGGREGATION_COMPLETE',started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),
        campaigns=campaigns,expected_tasks=len(data['TRAJECTORIES']),
        trajectory_rows=len(data['TRAJECTORIES']),complete_run_receipts=complete,recorded_native_steps=len(data['NATIVE_STEPS']),
        interventions=len(data['INTERVENTION_RECEIPTS']),branches=len(data['INTERVENTION_BRANCHES']),
        intervention_snapshots=len(data['INTERVENTION_SNAPSHOTS']),
        intervention_actual_paid_LU=sum(r.get('physical_operation_counts',{}).get('LU') for r in data['INTERVENTION_RECEIPTS']) if all(finite(r.get('physical_operation_counts',{}).get('LU')) for r in data['INTERVENTION_RECEIPTS']) else None,
        intervention_actual_paid_complex_RHS=sum(r.get('physical_operation_counts',{}).get('actual_complex_RHS') for r in data['INTERVENTION_RECEIPTS']) if all(finite(r.get('physical_operation_counts',{}).get('actual_complex_RHS')) for r in data['INTERVENTION_RECEIPTS']) else None,
        finite_JVP_response_rows=len(data['FINITE_JVP_RESPONSES']),
        existing_comparison_rows=len(data['FINITE_JVP_COMPARISONS']),failures=len(data['INTERVENTION_FAILURES']),issues=len(ISSUES),
        extraction_marker=marker,formula=dict(F_before='data_before_current+state_before_current+R_previous_native_row;first uniformR=0',
            F_after='data_term+state_term_current+regularizer',
            material_work='F_after−(data_after_current+state_after_current+R_previous_native_row)',
            stage='Carry prior acrossstage;new-stage beforefield already uses changed metric;no crossstage common-potential claim'),
        diagnostics='00/11 samepotential per-stage expectation;raw positive maxima and128epsscale comparison are descriptive,not acceptance gates;01/10 material work signs/counts reported without recovery inference',
        cost_scope='No new physics/actions;all original paid partial/overrun counts retained;receipt totals and nested event costs are separate views,never added twice',
        intervention_scope='Seven exposed posthoc checkpoints;8same-fullstate branches includingPR reset;fixedupdates not equalCPU,not independent material samples',
        unavailable_defect_angle='Raw defect-vector angle not in paid JSON;null,no recomputation;paid materialstep/negativephysicalgradientangle retained',
        no_GT_read=True,no_family_read=True,no_success_classification=True,new_forward=0,new_adjoint=0,new_optimizer=0,new_NN=0,
        scientific_review='PARENT_REQUIRED',CPU_seconds=time.process_time()-cpu0,wall_seconds=time.perf_counter()-wall0,
        timer_scope='Analysis entry through exports/plots;finalnotes/receipt/reportwrite andexitexcluded',figures=figures)
    write_json(out/'RECEIPT.json',receipt)
    notes=['# Mechanical dynamics notes — no GT',
        f"Analyzed {len(data['TRAJECTORIES'])} registered trajectories, {len(data['NATIVE_STEPS'])} recorded native updates, and {len(data['INTERVENTION_RECEIPTS'])} already-paid intervention checkpoints.",
        'F_before uses the previous native material prior, including across stage boundaries; the first uniform prior is0. F_after is evaluated in the current block metric after readout. Stage changes can change the potential definition, so stagewise work sums do not establish one global Lyapunov function.',
        'TRAJECTORIES/NATIVE_STAGE_SUMMARY retain raw maximum positive total/current/material changes and counts relative to a128eps-scaled arithmetic reference. That reference is neither a validated floating error bound nor a scientific quality gate. Iterations within each trajectory are dependent.',
        'INTERVENTION_SNAPSHOTS compare originalPhi, data/prior, auxiliary current-metric potential and constitutive defect after0/1/8/32 native updates. Every comparison restores the same saved full state; reset branches additionally clearPR memory. Fixed update counts do not equal CPU.',
        'Finite JVP/additivity tables copy existing paid central differences, actual feasible epsilons and branch/active-set flags. Small finite defects do not prove differentiability, contraction, a basin, or material recovery. A false recorded nonsmooth flag only means available flags did not change; unrecorded branches remain unknown.',
        'The only available angle is the paid latest material displacement versus negative originalPhi gradient at the previous material. A raw constitutive-defect vector angle is unavailable and stays null.',
        'Final counters may exceed last recorded native counters because budget interruption can occur inside a step or after its completion. Missing last rows are not added as synthetic iterations. Full receipt costs and per-event deltas are nested views, not additive ledgers.',
        'No truth, family/regime label or recovery threshold was read/applied. No Maxwell/adjoint/optimizer/NN call was made. Scientific interpretation remains with parent. Re-run with a new --out after remaining36 downloads; never overwrite this pilot artifact.']
    with (out/'MECHANICAL_NOTES.md').open('x') as f:f.write('\n\n'.join(notes)+'\n')
    body=['<!doctype html><meta charset="utf-8"><title>No-GT dynamics</title><style>body{font:16px system-ui;max-width:1400px;margin:auto;padding:30px}img{max-width:100%}</style>',
        '<h1>No-GT mechanical dynamics</h1>','<p>No new physics or material scoring. Same-stage potential changes, finite recorded response diagnostics and all prior costs remain separate. Parent scientific review required.</p>',
        '<p><a href="MECHANICAL_NOTES.md">Mechanical notes</a></p>']
    for name in data:body.append(f'<a href="{name}.csv">{name} CSV</a> · <a href="{name}.json">JSON</a><br>')
    for path in figures:body.append('<p><img src="'+html.escape(path)+'"></p>')
    with (out/'report.html').open('x') as f:f.write('\n'.join(body))
    print(json.dumps({k:receipt[k] for k in ('status','trajectory_rows','recorded_native_steps','interventions','branches','issues','CPU_seconds','wall_seconds')}));return 0


if __name__=='__main__':raise SystemExit(main())
