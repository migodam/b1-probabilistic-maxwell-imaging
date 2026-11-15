"""Post-seal pair01-versus-six paid-pool cost comparisons for independent24.

No inverse, optimization, forward/adjoint solve or source physics is imported.
All endpoint/curve candidate plans are written before importing the truth helper.
Uses the original base eligibility, min(Phi,candidate), physics counters and
material metrics. All original whole-child OS costs remain charged in full.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import html
import importlib.util
import json
from pathlib import Path
import sys
import time
import score_dualmetric as base
import score_extensions as extension

CANDIDATE='pair01'
REFERENCES=('cie','dual_cie','g19','dr','pair11','pair01_to11')
PREFIXES=(15,30,60,90,120)


def read_context(task,item,issues):
    """Same public HISTORY/receipt parsing and missing rules as base.public_plans."""
    ident=task['task_id'];history=[];receipt=None
    try:
        if item['history_file']['exists']:
            history=base.read(item['history_path'])
            if not isinstance(history,list) or any(r.get('candidate')!=i for i,r in enumerate(history)):
                raise ValueError('History/pool index schema is not consecutive')
    except Exception as exc:
        issues.append(dict(task_id=ident,issue='INVALID_PUBLIC_HISTORY',error=repr(exc)));history=[]
    try:
        if item['receipt_file']['exists'] and item.get('receipt_error') is None:
            receipt=base.read(item['receipt_path'])
    except Exception as exc:
        issues.append(dict(task_id=ident,issue='INVALID_PUBLIC_RECEIPT',error=repr(exc)))
    return dict(task=task,item=item,history=history,receipt=receipt)


def final_CPU(c):
    r=c['receipt'] or {}
    return r.get('CPU_seconds',r.get('whole_process_CPU_seconds'))


def selected_physics(history,lu,rhs):
    if lu is None or rhs is None:return None
    usable=[]
    for row in history:
        counts=base.cumulative_physics(row.get('physics_cost'))
        if (base.eligible(row) and counts['LU'] is not None and counts['complex_RHS'] is not None
                and counts['LU']<=lu and counts['complex_RHS']<=rhs):
            usable.append(row)
    return min(usable,key=lambda r:(r['Phi'],r['candidate'])) if usable else None


def public_pair_plans(caches,cases):
    lookup={(c['task']['case_id'],c['task']['method']):c for c in caches.values()}
    plans=[];limits=[]
    for cid in cases:
        for ref in REFERENCES:
            arms=[lookup[(cid,CANDIDATE)],lookup[(cid,ref)]]
            ends=[final_CPU(c) for c in arms]
            cpu_known=all(base.finite(t) and t>=0 for t in ends)
            cpu=min([120.]+ends) if cpu_known else None
            finals=[base.cumulative_physics((c['receipt'] or {}).get('physical_cost')) for c in arms]
            physics_known=all(v['LU'] is not None and v['complex_RHS'] is not None for v in finals)
            lu=min(v['LU'] for v in finals) if physics_known else None
            rhs=min(v['complex_RHS'] for v in finals) if physics_known else None
            limits.append(dict(case_id=cid,candidate=CANDIDATE,reference=ref,
                common_CPU_seconds=cpu,CPU_limit_available=cpu_known,
                candidate_final_internal_process_CPU_seconds=ends[0],reference_final_internal_process_CPU_seconds=ends[1],
                common_physical_LU=lu,common_physical_complex_RHS=rhs,physics_limits_available=physics_known,
                candidate_actual_final_physics=finals[0],reference_actual_final_physics=finals[1],
                CPU_scope='min(120,both original RECEIPT final internal processCPU); not last native proposal row or external childCPU',
                physics_scope='Componentwise min of both final paid original physical LU and aggregate complex RHS; both limits must hold',
                paid_pool_only=True,spent_cost_refunded=False))
            requests=[('CPU_endpoint','CPU',120,cpu),('physics_endpoint','physics',None,None)]
            requests.extend(('CPU_curve','CPU',p,min(float(p),cpu) if cpu is not None else None) for p in PREFIXES)
            for kind,budget,prefix,cutoff in requests:
                for side,c,fin,end in zip(('candidate','reference'),arms,finals,ends):
                    t=c['task'];selected=(base.pick(c['history'],cutoff) if cutoff is not None else None) if budget=='CPU' else selected_physics(c['history'],lu,rhs)
                    known=cpu_known if budget=='CPU' else physics_known
                    reason=('AVAILABLE' if selected else 'NO_ELIGIBLE_COMPLETED_PAID_CANDIDATE') if known else ('MISSING_FINAL_CPU_COUNTER' if budget=='CPU' else 'MISSING_FINAL_PHYSICS_COUNTER')
                    plans.append(dict(plan_id=f'{cid}__{CANDIDATE}_vs_{ref}__{kind}__{prefix}__{side}',
                        case_id=cid,comparison=CANDIDATE+'_vs_'+ref,candidate_method=CANDIDATE,reference_method=ref,
                        role=side,solver=t['method'],task_id=t['task_id'],campaign=t['campaign'],kind=kind,budget=budget,
                        nominal_CPU_prefix=prefix,cutoff_CPU_seconds=cutoff if budget=='CPU' else None,
                        common_recorded_process_CPU=cpu,common_physical_LU=lu if budget=='physics' else None,
                        common_physical_complex_RHS=rhs if budget=='physics' else None,
                        limits_available=known,public_availability_reason=reason,
                        selection_kind='pool',candidate=selected['candidate'] if selected else None,
                        metadata=selected,source_file=c['item']['candidate_pool_path'],source_exists=c['item']['candidate_pool_file']['exists'],
                        actual_final_internal_process_CPU_seconds=end,actual_final_physics_counts=fin,
                        selection_rule='original base.eligible + min(Phi,candidate), never material error or GT',
                        physics_counter_rules='base.cumulative_physics: actual_complex_RHS aggregate preferred; never add role counters twice',
                        observation_contract=t.get('observation_contract',{}),run_status=c['item']['status']))
    return plans,limits


def full_cost_ledger(root,cfg,tasks,caches):
    qpath=root/'queues'/cfg['campaign']/'RECEIPT.json'
    jobs={j['task_id']:j for j in base.read(qpath).get('jobs',[])} if qpath.is_file() else {}
    ledger=[]
    for t in tasks:
        c=caches[t['task_id']];r=c['receipt'] or {};j=jobs.get(t['task_id'],{})
        external=(j.get('external') or {}).get('child_process_CPU_seconds',j.get('child_process_CPU_seconds'))
        ledger.append(dict(task_id=t['task_id'],case_id=t['case_id'],solver=t['method'],run_status=c['item']['status'],
            used_in_these_comparisons=t['method'] in (CANDIDATE,)+REFERENCES,
            external_child_CPU_seconds=external,internal_final_process_CPU_seconds=final_CPU(c),
            external_wall_seconds=j.get('wall_seconds'),internal_wall_seconds=r.get('wall_seconds'),
            raw_physical_cost=r.get('physical_cost'),actual_final_physics=base.cumulative_physics(r.get('physical_cost')),
            auxiliary_cost=r.get('auxiliary_cost',r.get('auxiliary_receipt')),core=j.get('cpu'),core_type=j.get('core_type'),
            all_original_spent_CPU_charged=True,paid_pool_cutoffs_refund_nothing=True,
            actual_corrector_policy=cfg['protocol'].get('corrector_policy')))
    return ledger


def score_selected(plans,caches,truths,helper,np,issues):
    wanted=defaultdict(set)
    for p in plans:
        if p['candidate'] is not None:wanted[p['task_id']].add(p['candidate'])
    materials={};pool_errors={}
    for ident,indices in wanted.items():
        c=caches[ident]
        try:
            if not c['item']['candidate_pool_file']['exists']:raise FileNotFoundError(c['item']['candidate_pool_path'])
            with np.load(c['item']['candidate_pool_path'],allow_pickle=False) as z:pool=z['m']
            if len(pool)!=len(c['history']):raise ValueError('Material pool count differs from paid history')
            for i in indices:materials[(ident,i)]=pool[i].copy()
        except Exception as exc:
            pool_errors[ident]=repr(exc);issues.append(dict(task_id=ident,issue='INVALID_OR_MISSING_MATERIAL_POOL',error=repr(exc)))
    output=[]
    for p in plans:
        meta=p['metadata'] or {};m=materials.get((p['task_id'],p['candidate']));metrics=helper.material_metrics(m,truths.get(p['case_id']),.02,20.)
        r={k:v for k,v in p.items() if k not in ('metadata','source_exists','observation_contract')};r.update(metrics)
        if not p['source_exists']:r.update(available=False,availability_reason='MISSING_FIXED_PREDICTION_OR_POOL')
        elif not p['limits_available']:r.update(available=False,availability_reason=p['public_availability_reason'])
        elif p['candidate'] is None:r.update(available=False,availability_reason='NO_ELIGIBLE_COMPLETED_PAID_CANDIDATE')
        elif p['task_id'] in pool_errors:r.update(available=False,availability_reason='INVALID_OR_MISSING_MATERIAL_POOL')
        for k in ('CPU_seconds','wall_seconds','Phi','D','R','relative_residual','residual_over_delta','noise_feasible','legal','eligible','phase','deadline','lower_active','upper_active'):r[k]=meta.get(k)
        counts=base.cumulative_physics(meta.get('physics_cost'));r['selected_physical_LU']=counts['LU'];r['selected_physical_complex_RHS']=counts['complex_RHS']
        r['physical_cost_at_candidate']=helper.physical_counts(meta.get('physics_cost',{}))
        r['model_mismatch_status']=p['observation_contract'].get('model_mismatch_status');r['truth_available']=truths.get(p['case_id']) is not None
        r['accuracy_available']=base.accuracy_available(r);r['recovered_posthoc']=base.recovered(r);output.append(r)
    return output,materials


def paired_tables(rows,materials,np):
    lookup=defaultdict(dict)
    for r in rows:lookup[(r['case_id'],r['reference_method'],r['kind'],r['nominal_CPU_prefix'])][r['role']]=r
    diffs=[]
    metrics=('NRMSE','RMSE','relativeL2','IoU','Phi','D','R','residual_over_delta','CPU_seconds','selected_physical_LU','selected_physical_complex_RHS')
    for (cid,ref,kind,prefix),arms in lookup.items():
        a,b=arms['candidate'],arms['reference'];ok=base.accuracy_available(a) and base.accuracy_available(b)
        r=dict(case_id=cid,candidate=CANDIDATE,reference=ref,comparison=CANDIDATE+'_vs_'+ref,kind=kind,budget=a['budget'],nominal_CPU_prefix=prefix,
            paired_available=ok,candidate_available=a['available'],reference_available=b['available'],candidate_reason=a.get('availability_reason'),reference_reason=b.get('availability_reason'),
            candidate_task_id=a['task_id'],reference_task_id=b['task_id'],candidate_selection=a['candidate'],reference_selection=b['candidate'],
            common_CPU_seconds=a['common_recorded_process_CPU'],effective_CPU_cutoff_seconds=a['cutoff_CPU_seconds'],common_physical_LU=a['common_physical_LU'],common_physical_complex_RHS=a['common_physical_complex_RHS'],
            candidate_recovered_posthoc=base.recovered(a),reference_recovered_posthoc=base.recovered(b),
            delta_scope='candidate pair01 minus reference; fixed original posthoc thresholds are descriptive, never candidate selection')
        for k in metrics:
            av,bv=a.get(k),b.get(k);r['candidate_'+k]=av;r['reference_'+k]=bv
            r['delta_'+k]=av-bv if ok and base.finite(av) and base.finite(bv) else None
        ma=materials.get((a['task_id'],a['candidate']));mb=materials.get((b['task_id'],b['candidate']))
        if ma is not None and mb is not None and ma.shape==mb.shape:
            gap=ma-mb;rmse=float(np.linalg.norm(gap)/np.sqrt(gap.size));amp=a.get('amplitude')
            r['pair_material_RMSE']=rmse;r['pair_material_amplitude_normalized_RMSE']=rmse/amp if base.finite(amp) and amp>0 else None
            r['pair_material_Linf']=float(np.max(np.abs(gap)));r['pair_material_relative_L2_over_pair01']=float(np.linalg.norm(gap)/np.linalg.norm(ma)) if np.linalg.norm(ma)>0 else None
            r['support_disagreement_cells_at_original_chi1']=int(np.count_nonzero((ma>1)!=(mb>1)))
        diffs.append(r)
    summaries=[];gainloss=[]
    for ref in REFERENCES:
        for kind,prefix in [('CPU_endpoint',120),('physics_endpoint',None)]+[('CPU_curve',p) for p in PREFIXES]:
            g=[r for r in diffs if r['reference']==ref and r['kind']==kind and r['nominal_CPU_prefix']==prefix];valid=[r for r in g if r['paired_available']]
            out=dict(candidate=CANDIDATE,reference=ref,comparison=CANDIDATE+'_vs_'+ref,kind=kind,nominal_CPU_prefix=prefix,expected_cases=24,rows=len(g),paired_available=len(valid),missing_cases=[r['case_id'] for r in g if not r['paired_available']],
                candidate_recovered_count=sum(r['candidate_recovered_posthoc'] for r in g),reference_recovered_count=sum(r['reference_recovered_posthoc'] for r in g),recovery_denominator=24,
                scope='Fixed independent24 posthoc paired paid-pool comparison; no population/novelty certificate; missing cases retained; no cutoff refund')
            for k in ('NRMSE','IoU','Phi','residual_over_delta','pair_material_amplitude_normalized_RMSE'):
                key=k if k.startswith('pair_material') else 'delta_'+k;v=[r[key] for r in valid if base.finite(r.get(key))]
                out[key+'_mean']=float(np.mean(v)) if v else None;out[key+'_median']=float(np.median(v)) if v else None;out[key+'_n']=len(v)
            summaries.append(out)
            gainloss.append(dict(candidate=CANDIDATE,reference=ref,kind=kind,nominal_CPU_prefix=prefix,denominator=24,
                gain_cases=[r['case_id'] for r in valid if r['candidate_recovered_posthoc'] and not r['reference_recovered_posthoc']],
                loss_cases=[r['case_id'] for r in valid if r['reference_recovered_posthoc'] and not r['candidate_recovered_posthoc']],
                both_cases=[r['case_id'] for r in valid if r['candidate_recovered_posthoc'] and r['reference_recovered_posthoc']],
                neither_cases=[r['case_id'] for r in valid if not r['candidate_recovered_posthoc'] and not r['reference_recovered_posthoc']],
                unavailable_cases=[r['case_id'] for r in g if not r['paired_available']],scope='Original NRMSE<=.15 and IoU>=.70 conjunction, posthoc only; continuous fields/differences retained'))
    return diffs,summaries,gainloss


def figures(out,differences,summaries,np):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'savefig.bbox':'tight','axes.spines.top':False,'axes.spines.right':False})
    colors=dict(zip(REFERENCES,['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']));paths=[]
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ref in REFERENCES:
        g=sorted([r for r in summaries if r['kind']=='CPU_curve' and r['reference']==ref],key=lambda r:r['nominal_CPU_prefix'])
        for ax,k in zip(axes,('delta_NRMSE_mean','delta_IoU_mean')):ax.plot([r['nominal_CPU_prefix'] for r in g],[r[k] if r[k] is not None else np.nan for r in g],c=colors[ref],marker='.',label=ref)
    for ax,label in zip(axes,('Mean paired Δ NRMSE (pair01 − reference)','Mean paired Δ IoU (pair01 − reference)')):ax.set_xlabel('Nominal paid CPU prefix, clipped to each pair finalCPU');ax.set_ylabel(label);ax.set_xticks(PREFIXES);ax.axhline(0,c='gray',ls=':',lw=.7);ax.grid(alpha=.2)
    axes[0].legend(ncol=3);fig.suptitle('24 cases: separate two-arm CPU limits; public minΦ candidate plans');fig.tight_layout();fig.text(.5,-.02,'Original full OS costs remain charged. Missing paired observations stay missing; mean denominators are in tables.',ha='center',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(out/('PAIRWISE_CPU_CURVES.'+ext),dpi=160)
    plt.close(fig);paths.append('PAIRWISE_CPU_CURVES.png')
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ref in REFERENCES:
        g=[r for r in differences if r['kind']=='physics_endpoint' and r['reference']==ref and r['paired_available']]
        for ax,k in zip(axes,('delta_NRMSE','delta_IoU')):
            values=sorted(r[k] for r in g if base.finite(r.get(k)))
            if values:ax.step(values,np.arange(1,len(values)+1)/24,where='post',c=colors[ref],label=ref)
    for ax,label in zip(axes,('Paired Δ NRMSE (pair01 − reference)','Paired Δ IoU (pair01 − reference)')):ax.set_xlabel(label);ax.set_ylabel('Cumulative fraction of all 24 cases');ax.set_ylim(0,1.03);ax.axvline(0,c='gray',ls=':',lw=.7);ax.grid(alpha=.2)
    axes[0].legend(ncol=3);fig.suptitle('Endpoint comparison within each two-arm componentwise min(LU, complex RHS)');fig.tight_layout();fig.text(.5,-.02,'Paid-pool endpoint ECDF, not a physics-cost trajectory or a rerun. Missing cases leave the ECDF below 1.',ha='center',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(out/('PAIRWISE_PHYSICS_ENDPOINT_ECDF.'+ext),dpi=160)
    plt.close(fig);paths.append('PAIRWISE_PHYSICS_ENDPOINT_ECDF.png')
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--seal',type=Path,help='Default: ROOT/audit/independent24_SEAL_REMOTE_AUTHORIZED.json')
    parser.add_argument('--history-root',type=Path,help='Default: ROOT/runtime/scoring_history; reused score_results.py helper')
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    start_cpu,start_wall=time.process_time(),time.perf_counter();root=args.root.resolve();out=args.out.resolve()
    history=(args.history_root or root/'runtime/scoring_history').resolve();seal_path=(args.seal or root/'audit/independent24_SEAL_REMOTE_AUTHORIZED.json').resolve()
    seal,cfg,cases,methods,tasks,account=extension.authorize(root,seal_path)
    if cfg['campaign']!='independent24_v1' or len(cases)!=24 or len(tasks)!=264 or float(cfg['protocol']['CPU_cap'])!=120.:raise ValueError('Only the separately authorized complete independent24 campaign is supported')
    if out.exists() or out==root or out.is_relative_to(history) or any(out.is_relative_to(root/n) for n in ('RESULTS','queues','sealed_truth','data_public','src','configs','audit','runtime')):
        raise ValueError('New analysis directory outside existing results,truth,source,authorization and runtime required')
    out.mkdir(parents=True,exist_ok=False);issues=[]
    # Loading receipts for all eleven arms retains every original cost; only seven
    # arms' candidate histories/pools are needed by the fixed six comparisons.
    caches={}
    for task in tasks:
        item=account[task['task_id']]
        if task['method'] in (CANDIDATE,)+REFERENCES:c=read_context(task,item,issues)
        else:
            receipt=None
            try:
                if item['receipt_file']['exists'] and item.get('receipt_error') is None:receipt=base.read(item['receipt_path'])
            except Exception as exc:issues.append(dict(task_id=task['task_id'],issue='INVALID_PUBLIC_RECEIPT',error=repr(exc)))
            c=dict(task=task,item=item,receipt=receipt,history=[])
        caches[task['task_id']]=c
    compare_caches={k:c for k,c in caches.items() if c['task']['method'] in (CANDIDATE,)+REFERENCES}
    plans,limits=public_pair_plans(compare_caches,cases);ledger=full_cost_ledger(root,cfg,tasks,caches)
    # Every selection is final on disk before importing the truth helper or
    # reading any GT/family labels. No GT can choose a prefix or missing fallback.
    base.table(out,'PUBLIC_PAIRWISE_LIMITS_BEFORE_TRUTH',limits)
    base.table(out,'PUBLIC_PAIRWISE_PLANS_BEFORE_TRUTH',plans)
    base.table(out,'FULL_ORIGINAL_COST_LEDGER_BEFORE_TRUTH',ledger)
    base.write_json(out/'PROTOCOL_BEFORE_TRUTH.json',dict(campaign=cfg['campaign'],authorized_seal=str(seal_path),candidate=CANDIDATE,references=list(REFERENCES),cases=list(cases),endpoint_arm_plans=24*6*2*2,CPU_curve_arm_plans=24*6*2*5,total_plans=len(plans),original_methods=list(methods),original_tasks=len(tasks),eligibility='Exactly base.eligible; physical pool candidates additionally within both component limits; min(Phi,candidate).',CPU_final_definition='Original RECEIPT final CPU_seconds or whole_process_CPU_seconds; not last proposal NATIVE_HISTORY or external OS cost',CPU_prefixes=list(PREFIXES),physics_curve='Endpoint ECDF only, no invented fractional-physics budgets',recovery_rule='Inherited posthoc NRMSE<=.15 AND IoU>=.70; never selection',cost_rule='All original OS costs across 264 unique tasks retained, no cutoff refunds; per-comparison costs reuse arms and must not be summed as new jobs',truth_access_at_plan_write='NONE',new_inverse=False,new_forward=False))
    import numpy as np
    helper_path=history/'src/score_results.py';sys.path.insert(0,str(history/'src'))
    spec=importlib.util.spec_from_file_location('r14_pairwise_cost_old_metrics',helper_path);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    truths={};truth_sources=[]
    for cid in cases:
        truth_spec=dict(path=str(root/'sealed_truth/independent24'/(cid+'.npz')),key='material_coarse',support_rule='fixed_chi_gt_1',amplitude_rule='max_truth_minus_host',host=.1,object_id=cid)
        truth_sources.append(dict(case_id=cid,**truth_spec))
        try:truths[cid]=helper.load_truth(truth_spec)
        except Exception as exc:truths[cid]=None;issues.append(dict(case_id=cid,issue='AUTHORIZED_TRUTH_UNAVAILABLE',error=repr(exc)))
    rows,materials=score_selected(plans,compare_caches,truths,helper,np,issues)
    differences,summaries,gainloss=paired_tables(rows,materials,np)
    for budget,kind in [('CPU','CPU_endpoint'),('PHYSICS','physics_endpoint')]:
        base.table(out,budget+'_SELECTED_ARMS',[r for r in rows if r['kind']==kind]);base.table(out,budget+'_PAIRED_DIFFERENCES',[r for r in differences if r['kind']==kind]);base.table(out,budget+'_SUMMARY',[r for r in summaries if r['kind']==kind]);base.table(out,budget+'_POSTHOC_GAIN_LOSS',[r for r in gainloss if r['kind']==kind])
    for name,data in [('CPU_CURVE_SELECTED_ARMS',[r for r in rows if r['kind']=='CPU_curve']),('CPU_CURVE_PAIRED_DIFFERENCES',[r for r in differences if r['kind']=='CPU_curve']),('CPU_CURVE_SUMMARY',[r for r in summaries if r['kind']=='CPU_curve']),('TRUTH_SOURCES',truth_sources),('ISSUES',issues)]:base.table(out,name,data)
    figures_out=figures(out,differences,summaries,np)
    ext=[r['external_child_CPU_seconds'] for r in ledger];known=[x for x in ext if base.finite(x)]
    receipt=dict(status='OFFLINE_PAIRWISE_COST_ANALYSIS_COMPLETE',campaign=cfg['campaign'],cases=24,comparisons=6,expected_endpoint_arm_rows=576,actual_endpoint_arm_rows=sum(r['kind'].endswith('endpoint') for r in rows),expected_endpoint_case_pairs=288,actual_endpoint_case_pairs=sum(r['kind'].endswith('endpoint') for r in differences),expected_CPU_curve_arm_rows=1440,actual_CPU_curve_arm_rows=sum(r['kind']=='CPU_curve' for r in rows),source_seal=str(seal_path),metrics_helper=str(helper_path),original_unique_tasks=264,original_run_status_counts=dict(Counter(r['run_status'] or 'MISSING' for r in ledger)),full_original_OS_child_CPU_seconds=sum(known) if len(known)==len(ext) else None,measured_original_OS_child_CPU_subtotal=sum(known),missing_original_external_CPU_tasks=[r['task_id'] for r in ledger if not base.finite(r['external_child_CPU_seconds'])],all_original_costs_still_charged=True,cutoff_refund=False,curve_scope='CPU prefixes original15/30/60/90/120 clipped independently for each two-arm comparison; physics figure is endpoint-difference ECDF only',new_inverse=0,new_forward=0,new_physical_LU=0,new_complex_RHS=0,new_optimizer=0,truth_scope='Authorized independent24 only, after all public plans written',issues=len(issues),offline_CPU_seconds=time.process_time()-start_cpu,offline_wall_seconds=time.perf_counter()-start_wall,figures=figures_out,finished_utc=datetime.now(timezone.utc).isoformat(),scientific_review='Parent owns scientific conclusions; no final novelty or recovery-causality judgment')
    base.write_json(out/'RECEIPT.json',receipt)
    text=['# Independent24：pair01的六组成对成本附注','', '固定pair01分别比较cie、dual_cie、g19、dr、pair11、pair01_to11，24case全保留。CPU端点取min(120,两臂原receipt最终processCPU)，物理端点取两臂最终physical LU/aggregate complex RHS的逐分量共同上限。候选必须满足原base eligible/legal/CPU≤120规则，再按公开min(Φ,candidate)选择；缺失费用不补0，缺失候选不以GT另选。','', '所有端点及CPU曲线计划先写盘，再导入旧score_results helper并读取已授权independent24真值。CPU_CURVE沿用原15/30/60/90/120标称prefix，每case各自clip到本次两臂共同CPU；不能视作全组11方法共同上限。物理小图只画端点paired-difference ECDF，未增加新的物理预算网格。','', 'NRMSE/IoU/RMSE/relativeL2完全复用旧评分函数。CPU/PHYSICS_PAIRED_DIFFERENCES另有选中两臂材料数组间连续RMSE、truth-amplitude归一化RMSE、Linf及固定χ>1 support disagreement，避免仅看阈值过线。差值均pair01−reference。原NRMSE≤.15且IoU≥.70只作后验gain/loss描述，不参与任何选择。','', '原264个任务的全部OS子进程费用保留在FULL_ORIGINAL_COST_LEDGER_BEFORE_TRUTH；成对cutoff只是原已付费pool的离线比较，不退还后续CPU/物理费用，不是假装提前停止的inverse重跑。pair01在六比较中复用，不能把六比较费用相加当六份新实验。辅助FFT/块动作/设置/失败费用照原receipt保留，未被LU/RHS归一化为免费。','', '这份附注不更改已有score_extensions或任何原预测/评分；不作最终科学、新颖性或恢复因果解释。', '',json.dumps(receipt,ensure_ascii=False,indent=2)]
    (out/'REPORT_ZH.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    (out/'report.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>Pairwise paid-pool costs</title><style>body{font-family:system-ui;max-width:1400px;margin:25px auto;line-height:1.6}pre{white-space:pre-wrap}img{width:100%}</style><h1>Independent24：六组成对成本比较</h1><pre>'+html.escape('\n'.join(text))+'</pre>'+''.join(f'<img src="{p}">' for p in figures_out)+'</html>',encoding='utf-8')
    print(json.dumps(receipt,indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
