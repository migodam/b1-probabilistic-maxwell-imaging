"""Authorized post-seal scoring of the fixed exposed16 shared-readout control.

All view/prefix decisions are frozen from public Phi before any truth loader.
No forward, optimizer, prediction edit, or new method/threshold is used.
"""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import argparse
import csv
from datetime import datetime,timezone
import itertools
import json
from pathlib import Path
import resource
import sys
import time


def read(path):return json.loads(Path(path).read_text())


def dump(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def valid(row):
    return bool(row and row.get('eligible') and row.get('legal') and isinstance(row.get('Phi'),(int,float)))


def validate(root,manifest,seal_path,authorization_path):
    approval=read(authorization_path);seal=read(seal_path)
    if approval.get('authorized_GT_scoring') is not True or approval.get('expected_tasks')!=16:
        raise PermissionError('Explicit parent authorization for the fixed16tasks is required')
    if Path(approval['prediction_seal']).resolve()!=seal_path.resolve():
        raise ValueError('Parent authorization references another prediction seal')
    expected=sorted(task['task_id'] for task in manifest['tasks']);actual=sorted(item['task_id'] for item in seal['inventory'])
    if len(expected)!=16 or expected!=actual or seal.get('terminal_tasks')!=16:
        raise ValueError('Unexpected or nonterminal shared16task inventory')
    verified=0
    for item in seal['inventory']:
        for record in item['files']:
            stat=(root/record['path']).stat()
            if stat.st_size!=record['bytes'] or stat.st_mtime_ns!=record['mtime_ns']:
                raise ValueError('Sealed prediction artifact changed: '+record['path'])
            verified+=1
    return approval,seal,verified


def public_plans(root,manifest):
    """Complete public selection plan; never imports or reads material truth."""
    saved_views=[];plans=[];histories={};raw_audits=[]
    for task in manifest['tasks']:
        folder=root/'RESULTS/shared_readout4_v1'/task['task_id']
        history=read(folder/'HISTORY.json') if (folder/'HISTORY.json').exists() else []
        histories[task['task_id']]=history;eligible=[record for record in history if valid(record)]
        metadata={view:read(folder/(view+'.json')) if (folder/(view+'.json')).exists() else None
                  for view in ('initial','raw','corrected','primary','secondary')}
        for view in metadata:
            saved_views.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],view=view,
                exists=(folder/(view+'.npz')).exists(),source=str(folder/(view+'.npz')),
                candidate=None if metadata[view] is None else metadata[view].get('candidate'),
                Phi=None if metadata[view] is None else metadata[view].get('Phi'),
                metadata_eligible=valid(metadata[view]),truth_access='NONE'))
        raw_candidates=[record for record in eligible if record['phase'] in ('initial','proposal')]
        raw=min(raw_candidates,key=lambda record:record['Phi']) if raw_candidates else None
        raw_source='initial' if raw and raw['phase']=='initial' else 'raw'
        if raw and (metadata[raw_source] is None or metadata[raw_source]['candidate']!=raw['candidate']):raw_source='pool'
        raw_audits.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],
            archived_raw_candidate=None if metadata['raw'] is None else metadata['raw']['candidate'],
            archived_raw_Phi=None if metadata['raw'] is None else metadata['raw']['Phi'],
            initial_Phi=None if metadata['initial'] is None else metadata['initial']['Phi'],
            public_raw_candidate=None if raw is None else raw['candidate'],public_raw_Phi=None if raw is None else raw['Phi'],
            raw_initial_fallback=raw_source=='initial',raw_source=raw_source))
        for view in ('initial','raw','corrected','primary','secondary'):
            fallback=False;source_view=raw_source if view=='raw' else view
            chosen=raw if view=='raw' else metadata[view]
            if view=='secondary' and not (folder/'secondary.npz').exists():
                source_view='primary';chosen=metadata['primary'];fallback=True
            source=str(folder/(source_view+'.npz')) if source_view!='pool' else str(folder/'CANDIDATE_POOL.npz')
            available=chosen is not None and (Path(source).exists())
            if source_view=='pool':available=available and chosen['candidate']<len(history)
            expected=None
            if view=='raw':expected=raw
            elif view=='primary':expected=min(eligible,key=lambda record:record['Phi']) if eligible else None
            elif view=='corrected':
                candidates=[record for record in eligible if record['phase']=='corrector']
                expected=min(candidates,key=lambda record:record['Phi']) if candidates else None
            elif view=='secondary' and not fallback:
                candidates=[record for record in eligible if record.get('noise_feasible')]
                expected=min(candidates,key=lambda record:(record['R'],record['Phi'])) if candidates else None
            plans.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],view=view,
                source_view=source_view,source=source,available=bool(available),secondary_fallback=fallback,
                candidate=None if chosen is None else chosen['candidate'],public_metadata=chosen,
                public_selector_expected_candidate=None if expected is None else expected['candidate'],
                archived_selector_consistent=None if expected is None or chosen is None else expected['candidate']==chosen['candidate'],
                selection='archived paid originalPhi/ball-prior rule, raw min(initial,proposal); chosen before GT'))
        for point in (15,30,60,90,120):
            candidates=[record for record in eligible if record['CPU_seconds']<=point and record['CPU_seconds']<=120]
            chosen=min(candidates,key=lambda record:record['Phi']) if candidates else None
            plans.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],view='paid_prefix',CPU_prefix=point,
                source_view='pool',source=str(folder/'CANDIDATE_POOL.npz'),available=chosen is not None,
                candidate=None if chosen is None else chosen['candidate'],public_metadata=chosen,
                selection='paid legal eligible minPhi with processCPU<=prefix; no GT-best choice'))
    return saved_views,plans,histories,raw_audits


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--history-root',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--manifest',type=Path)
    parser.add_argument('--seal',type=Path);parser.add_argument('--authorization',type=Path)
    a=parser.parse_args();root=a.root.resolve();out=a.out.resolve();cpu,wall=time.process_time(),time.perf_counter()
    manifest=read(a.manifest or root/'configs/shared_readout4_v1.json')
    seal_path=(a.seal or root/'predictions/shared_readout4_v1/attempt_02/PREDICTION_SEAL.json').resolve()
    authorization=(a.authorization or root/'audit/shared_readout4_v1/PARENT_SCORING_AUTHORIZATION.json').resolve()
    approval,seal,verified=validate(root,manifest,seal_path,authorization)
    saved_views,plans,histories,raw_audits=public_plans(root,manifest)
    out.mkdir(parents=True,exist_ok=False)
    dump(out/'PUBLIC_SAVED_VIEWS_BEFORE_TRUTH.json',saved_views)
    dump(out/'PUBLIC_SELECTION_PLAN_BEFORE_TRUTH.json',plans)
    dump(out/'PUBLIC_RAW_AUDIT_BEFORE_TRUTH.json',raw_audits)
    dump(out/'AUTHORIZATION_AND_SEAL_VALIDATION.json',dict(status='VALIDATED_AND_SELECTION_PLAN_WRITTEN_BEFORE_FIRST_TRUTH_ACCESS',
        authorization=approval,seal=str(seal_path),validated_artifact_files=verified,
        selection_plans=len(plans),utc=datetime.now(timezone.utc).isoformat(),predictions_modified=False,seal_modified=False))
    sys.path[:0]=[str(root/'src'),str(root.parents[1])]
    import numpy as np
    from research.research_max_struct_r13.src.score_results import truth_spec,load_truth,material_metrics,physical_counts,table
    tasks={task['task_id']:task for task in manifest['tasks']};cases=manifest['cases'];pairs=('00','01','10','11')
    truths={};truthrecords=[]
    for case in cases:
        spec=truth_spec(a.history_root.resolve(),dict(case_id=case));stamp=datetime.now(timezone.utc).isoformat();truth=load_truth(spec)
        if truth is None:raise ValueError('Inherited truth is absent in supplied actualhistoryroot: '+case)
        truths[case]=truth;truthrecords.append(dict(case_id=case,truth_source=spec['path'],key=spec['key'],
            support_rule=spec['support_rule'],amplitude_rule=spec.get('amplitude_rule'),host=spec.get('host'),
            amplitude=truth['amplitude'],accessed_utc=stamp,scope='authorized offline only after publicselectionplan and sealvalidation'))
    pools={};predictions={};rows=[];costs=[];materials=[];truth_arrays=[];supports=[];truth_supports=[];labels=[]
    receipts={}
    for task in manifest['tasks']:
        folder=root/'RESULTS/shared_readout4_v1'/task['task_id'];receipt=read(folder/'RECEIPT.json');receipts[task['task_id']]=receipt
        with np.load(folder/'CANDIDATE_POOL.npz',allow_pickle=False) as z:pools[task['task_id']]=z['m'].copy()
        costs.append(dict(task_id=task['task_id'],case_id=task['case_id'],pair=task['pair'],status=receipt['status'],
            same_host_group='shared_readout4_v1 existingAppleCPU only',body_CPU_seconds=receipt['CPU_seconds'],
            wall_seconds=receipt['wall_seconds'],CPU_cap=120,proposal_absoluteCPU_cap=72,
            corrector_contract='originalR14 remaining absolute120 allowance, including unusedproposal time',
            auxiliary_cost=receipt['auxiliary_cost'],physical_cost=receipt['physical_cost'],**physical_counts(receipt['physical_cost']),
            corrector_termination=receipt['termination_reason'].get('corrector')))
    for plan in plans:
        task=tasks[plan['task_id']];truth=truths[plan['case_id']];m=None
        if plan['available']:
            if plan['source_view']=='pool':m=pools[plan['task_id']][plan['candidate']].copy()
            else:
                with np.load(plan['source'],allow_pickle=False) as z:m=z['m'].copy()
        metadata=plan['public_metadata'];metrics=material_metrics(m,truth,.02,20.)
        row=dict(task_id=plan['task_id'],case_id=plan['case_id'],pair=plan['pair'],method='pair'+plan['pair'],view=plan['view'],
            material_solver='shared_LBFGS60',source_view=plan['source_view'],selection_candidate=plan['candidate'],
            secondary_fallback=plan.get('secondary_fallback',False),archived_selector_consistent=plan.get('archived_selector_consistent'),
            CPU_prefix=plan.get('CPU_prefix'),run_terminal_before_prefix=bool(receipts[plan['task_id']]['CPU_seconds']<plan.get('CPU_prefix',0)),
            D=None if metadata is None else metadata['D'],R=None if metadata is None else metadata['R'],Phi=None if metadata is None else metadata['Phi'],
            actual_candidate_CPU_seconds=None if metadata is None else metadata['CPU_seconds'],
            physical_relative_residual=None if metadata is None else metadata['relative_residual'],
            noise_feasible=None if metadata is None else metadata['noise_feasible'],
            model_mismatch_status=task['observation_contract'].get('model_mismatch_status'),
            scope='four exposedcase commoninnerprogram mechanismcontrol, no newconfirmation',**metrics)
        rows.append(row)
        if plan['view']!='paid_prefix':
            predictions[(plan['case_id'],plan['pair'],plan['view'])]=m
            materials.append(np.full(truth['m'].shape,np.nan) if m is None else m)
            truth_arrays.append(truth['m']);supports.append(np.zeros(truth['m'].shape,bool) if m is None else m>1.);truth_supports.append(truth['support'])
            labels.append(plan['task_id']+'__'+plan['view'])
    views=[row for row in rows if row['view']!='paid_prefix'];prefixes=[row for row in rows if row['view']=='paid_prefix']
    summary=[]
    for pair in pairs:
        for view in ('initial','raw','corrected','primary','secondary'):
            selected=[row for row in views if row['pair']==pair and row['view']==view];available=[row for row in selected if row['available']]
            summary.append(dict(pair=pair,view=view,expected_cases=4,available=len(available),
                **{'mean_'+key:float(np.mean([row[key] for row in available])) if available else None for key in ('NRMSE','IoU','RMSE','relativeL2','D','R','Phi')},
                scope='paired fixedfour exposedcases; descriptive mechanical summary'))
    lookup={(row['case_id'],row['pair'],row['view']):row for row in views}
    differences=[];difference_arrays=[];difference_labels=[];corrector=[]
    for case in cases:
        for view in ('raw','corrected','primary','secondary'):
            for first,second in itertools.combinations(pairs,2):
                x,y=predictions.get((case,first,view)),predictions.get((case,second,view));available=x is not None and y is not None
                d=None if not available else x-y;amp=truths[case]['amplitude']
                differences.append(dict(case_id=case,view=view,pair_a=first,pair_b=second,available=available,
                    material_L2_distance=None if d is None else float(np.linalg.norm(d)),
                    material_amplitude_normalized_RMSE_distance=None if d is None else float(np.linalg.norm(d)/np.sqrt(d.size)/amp),
                    symmetric_material_relativeL2_distance=None if d is None else float(np.linalg.norm(d)/max(np.linalg.norm(x),np.linalg.norm(y),1e-300)),
                    material_Linf_distance=None if d is None else float(np.max(abs(d))),
                    support_flips=None if d is None else int(np.count_nonzero((x>1)!=(y>1))),
                    Phi_a=lookup[(case,first,view)]['Phi'],Phi_b=lookup[(case,second,view)]['Phi']))
                difference_arrays.append(np.full(truths[case]['m'].shape,np.nan) if d is None else d);difference_labels.append(case+'__'+view+'__'+first+'minus'+second)
        for pair in pairs:
            x,y=lookup[(case,pair,'raw')],lookup[(case,pair,'corrected')]
            corrector.append(dict(case_id=case,pair=pair,
                **{prefix+key:row.get(key) for prefix,row in (('raw_',x),('corrected_',y)) for key in ('NRMSE','IoU','D','R','Phi')},
                NRMSE_change=None if x.get('NRMSE') is None or y.get('NRMSE') is None else y['NRMSE']-x['NRMSE']))
    for name,data in (('SAVED_VIEWS',saved_views),('SCORES',views),('PAID_PREFIX_SCORES',prefixes),('METHOD_VIEW_SUMMARY',summary),
        ('MATERIAL_ARRAY_DIFFERENCES',differences),('CORRECTOR_DECOMPOSITION',corrector),('TASK_COSTS',costs),('TRUTH_SOURCES',truthrecords),('RAW_INITIAL_AUDIT',raw_audits)):
        table(out,name,data)
    np.savez_compressed(out/'SCORED_MATERIAL_ROWS_AND_MASKS.npz',m=np.asarray(materials),truth=np.asarray(truth_arrays),
        prediction_support=np.asarray(supports),truth_support=np.asarray(truth_supports),row_labels=np.asarray(labels),
        support_threshold=np.asarray(1.),box=np.asarray([.02,20.]))
    np.savez_compressed(out/'MATERIAL_DIFFERENCE_ARRAYS.npz',difference=np.asarray(difference_arrays),labels=np.asarray(difference_labels))
    figures=render(out,predictions,truths,lookup,cases,pairs,np)
    table(out,'FIGURE_MANIFEST',figures)
    queue=read(root/'queues/shared_readout4_v1/attempt_02/RECEIPT.json')
    totals=dict(successful_children_wholeOS_CPU_seconds=queue['solver_children_whole_OS_CPU_seconds'],
        prior_failed_bootstrap_children_OS_CPU_seconds=queue['prior_failed_attempt_OS_CPU_seconds'],
        cumulative_children_OS_CPU_seconds=queue['cumulative_solver_children_OS_CPU_seconds'],
        controller_CPU_seconds=queue['controller_CPU_seconds'],controller_wall_seconds=queue['controller_wall_seconds'],
        sum_child_body_CPU_seconds=sum(row['body_CPU_seconds'] for row in costs),
        physical_LU=sum(row['physical_LU'] for row in costs),physical_complex_RHS=sum(row['physical_complex_RHS'] for row in costs),
        failed_attempt_physics_LU=0,failed_attempt_physics_RHS=0,offline_scoring_LU=0,offline_scoring_RHS=0,CUDA_seconds=0,
        cost_counter_reconciliation=str(root/'audit/shared_readout4_v1/TASK_COSTS_AND_COUNTER_RECONCILIATION.csv'),
        normalization='bodyCPU and OSchildCPU are separate clockboundaries, never added together; hostgroup local16 only')
    dump(out/'COST_TOTALS.json',totals)
    report(out,views,summary,corrector,totals,raw_audits,cases,pairs)
    receipt=dict(status='COMPLETE',tasks=16,saved_view_score_rows=len(views),primary_rows=sum(row['view']=='primary' for row in views),
        available_saved_views=sum(row['available'] for row in views),paid_prefix_rows=len(prefixes),material_difference_rows=len(differences),
        actual_archived_saved_view_files=sum(row['exists'] for row in saved_views),
        secondary_primary_fallback_rows=sum(row['secondary_fallback'] for row in views if row['view']=='secondary'),
        available_saved_views_definition='available reported estimators including explicitly labeled secondaryfallback; see actual_archived_saved_view_files',
        public_selector_mismatches=sum(plan.get('archived_selector_consistent') is False for plan in plans),
        figure_pairs=len(figures),truth_files_accessed=len(truthrecords),public_plan_written_before_truth=True,
        GT_access='explicit parentauthorized offline after unchangedattempt02seal',predictions_modified=False,seal_modified=False,
        offline_CPU_seconds=time.process_time()-cpu,whole_offline_process_CPU_seconds=time.process_time(),offline_wall_seconds=time.perf_counter()-wall,
        offline_peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        forward_LU=0,solver_RHS=0,optimizer_calls=0,CUDA_seconds=0,
        scientific_scope='commonLBFGS16mechanism only; no finalnoveltyjudgment, no MacWindows absolute speedcomparison')
    dump(out/'SCORING_RECEIPT.json',receipt)
    print(json.dumps(dict(status='COMPLETE',views=len(views),available=receipt['available_saved_views'],prefixes=len(prefixes),figures=len(figures),CPU_seconds=receipt['offline_CPU_seconds'])))
    return 0


def render(out,predictions,truths,lookup,cases,pairs,np):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    folder=out/'FIGURES';folder.mkdir();figures=[]
    for view in ('raw','primary'):
        fig,axes=plt.subplots(4,5,figsize=(15,13),squeeze=False,constrained_layout=True)
        for i,case in enumerate(cases):
            for j,(pair,m) in enumerate([(None,truths[case]['m'])]+[(pair,predictions.get((case,pair,view))) for pair in pairs]):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                if m is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center',transform=ax.transAxes)
                else:ax.imshow(m.reshape(24,24),origin='lower',cmap='viridis',vmin=0,vmax=20,interpolation='nearest')
                if pair is None:ax.set_title(case+'\nOffline GT',fontsize=10)
                else:
                    row=lookup[(case,pair,view)];metric='UNAVAILABLE' if not row['available'] else f"NRMSE {row['NRMSE']:.5f} | IoU {row['IoU']:.4f}"
                    ax.set_title('pair'+pair+' | shared LBFGS60\n'+metric,fontsize=10)
        fig.colorbar(ScalarMappable(norm=Normalize(0,20),cmap='viridis'),ax=axes.ravel().tolist(),shrink=.8,label='shared real chi = epsilon_r - 1; absolute scale [0,20]')
        fig.suptitle('Shared material program, '+view+' | four exposed cases × four metric pairs\nOriginal24×24 material arrays; no interpolation, GT selection or new optimizer',fontsize=13)
        name='MATERIALS_'+view.upper()
        for extension in ('png','pdf'):fig.savefig(folder/(name+'.'+extension),dpi=170,bbox_inches='tight')
        plt.close(fig);figures.append(dict(view=view,png=str(folder/(name+'.png')),pdf=str(folder/(name+'.pdf')),scale=[0,20],cases=4,metric_pairs=4))
    return figures


def report(out,rows,summary,corrector,totals,raw_audits,cases,pairs):
    lookup={(row['case_id'],row['pair'],row['view']):row for row in rows}
    lines=['# Shared L-BFGS-B60：固定四对象机制控制离线评分','',
        '父线程已验收source/solver/cost并另行授权attempt02 seal后的GT评分。所有saved-view、raw min(initial,proposal)、secondaryfallback和paidCPU-prefix选择计划都在读取GT前按公开Phi/eligibility写盘。原seal/预测保持不变；没有forward/optimizer/新参数搜索。','',
        '本组为0002/0017/0032/0037四个已曝光对象×四metric。Material统一原chi坐标、同current-m起点、同box/prior、L-BFGS-B60/ftol1e-12/gtol1e-8/maxls25、完整metric导数、无scalarseed/Arc/MM，旧m非增保障。它控制materialtarget与有限innerprogram同时改变的混杂；不是新主算法、全48恢复结果或独立确认。','',
        '## 归一化与选择','',
        '复用原R13 truth_spec/load_truth/material_metrics。Amplitude NRMSE=全域RMSE/(max(coarseGT)-.1)，另报relativeL2与boxrange_NRMSE；固定prediction chi>1与truth chi>1支持。Phi=.5全Frobenius相对residual平方+原prior；GT不选择candidate、阶段、prefix或threshold。',
        '',f"保存view共80行（含16initial），其中raw/corrected/primary/secondary64行；paidprefix15/30/60/90/120共80行。Raw按公开Phi回退initial的条数为{sum(row['raw_initial_fallback'] for row in raw_audits)}。缺失corrected保留missing；缺失secondary明确fallbackprimary。",'',
        '## 每对象Primary完整四组合','']
    for case in cases:
        lines += ['### '+case,'','| pair | NRMSE | IoU | Phi |','|---|---:|---:|---:|']
        for pair in pairs:
            row=lookup[(case,pair,'primary')]
            lines.append(f"| {pair} | {row['NRMSE']:.7f} | {row['IoU']:.6f} | {row['Phi']:.9g} |")
        lines.append('')
    lines += ['## 描述性汇总与材料数组差','', '| pair | primary mean NRMSE | primary mean IoU |','|---|---:|---:|']
    for pair in pairs:
        row=next(row for row in summary if row['pair']==pair and row['view']=='primary')
        lines.append(f"| {pair} | {row['mean_NRMSE']:.7f} | {row['mean_IoU']:.6f} |")
    lines += ['', 'MATERIAL_ARRAY_DIFFERENCES保留每case、每raw/corrected/primary/secondary、全部六pair-pair比较的L2、原amplitude归一化RMS、relativeL2、Linf、supportflips及共同Phi。NPZ保存完整材料/GT/固定mask及数组差；相似NRMSE不代表相同材料数组。',
        '', f"Corrector后NRMSE增加的条数为{sum(row['NRMSE_change'] is not None and row['NRMSE_change']>0 for row in corrector)}；所有D/R/Phi和材料变化逐条保留。描述性统计不作最终新颖性判断。",'',
        '## 同host费用和失败','',
        f"成功16children完整OSCPU={totals['successful_children_wholeOS_CPU_seconds']:.6f}s；旧psutil-bootstrap失败CPU={totals['prior_failed_bootstrap_children_OS_CPU_seconds']:.6f}s，累计={totals['cumulative_children_OS_CPU_seconds']:.6f}s。Controller CPU={totals['controller_CPU_seconds']:.6f}s。旧16次setupfail均在物理前，0LU/RHS，完整日志/原失败seal保留。",
        '', f"实际physical LU={totals['physical_LU']}，complex RHS={totals['physical_complex_RHS']}。Inherited M1 whitening计数的名义offset单列于已完成costreconciliation；CPU/LU/RHS不改。Offline评分/图形CPU单列SCORING_RECEIPT，额外LU/RHS/optimizer/CUDA=0。",'',
        '## 解释边界','',
        '图MATERIALS_RAW与MATERIALS_PRIMARY均为4×(GT+4pairs)，共享绝对chi[0,20]色标，不插值。只有本地同16组内CPU具有共同host边界；不比较Mac与旧Windows绝对速度。',
        '', '本地算力使本组各臂都完成3000nativeupdates；原Arc对照在旧host/72proposal cap下的迭代能力不同。即使与原Arc结果有差异，也不能把全部变化归因于materialsolver。此次共同innerprogram只提供这四个曝光对象上的机制辨别，主综合与最终新颖性/科学解释由父线程负责。']
    (out/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':sys.exit(main())
