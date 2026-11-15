"""Parent-authorized post-seal vector4 metrics and figures; no forward/optimizer.

Inherited R13 truth, support and amplitude contracts are used unchanged.
Candidate choice uses only frozen paid eligibility/Phi, never material truth.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from collections import Counter
from datetime import datetime,timezone
import json
from pathlib import Path
import resource
import sys
import time


def read(path):return json.loads(Path(path).read_text())


def validate(args,manifest):
    seal=read(args.seal);approval=read(args.authorization)
    if approval.get('status')!='PARENT_AUTHORIZED_OFFLINE_GT_SCORING_AFTER_EXISTING_SEAL':
        raise PermissionError('Explicit reviewed parent scoring authorization is required')
    if Path(approval['reviewed_seal']).resolve()!=args.seal.resolve():
        raise ValueError('Scoring authorization references a different seal')
    if seal.get('terminal_tasks')!=24 or seal.get('successful')!=24 or seal.get('failed_or_missing')!=0:
        raise ValueError('Expected reviewed24 terminal/24 complete/0failed seal')
    expected=sorted(task['task_id'] for task in manifest['tasks'])
    actual=sorted(task['task_id'] for task in seal['inventory'])
    if expected!=actual or len(set(actual))!=24:raise ValueError('Unexpected/missing/duplicate sealed task')
    validated=0
    for task in seal['inventory']:
        for entry in task['files']:
            path=args.root/entry['path'];stat=path.stat()
            if stat.st_size!=entry['bytes'] or stat.st_mtime_ns!=entry['mtime_ns']:
                raise ValueError('Sealed artifact changed before offline truth access: '+str(path))
            validated+=1
    return seal,approval,validated


def raw_audit(args,manifest):
    rows=[]
    for task in manifest['tasks']:
        folder=args.root/'RESULTS/vector4_v1'/task['task_id']/'attempt_01'
        history=read(folder/'HISTORY.json');saved=read(folder/'raw.json')
        eligible=[row for row in history if row['eligible'] and row['legal'] and row['phase'] in ('initial','proposal')]
        paid=min(eligible,key=lambda row:row['Phi']) if eligible else None
        source='initial' if paid is not None and paid['phase']=='initial' else 'raw'
        rows.append(dict(task_id=task['task_id'],stored_raw_candidate=saved['candidate'],
            paid_initial_or_proposal_min_candidate=None if paid is None else paid['candidate'],
            stored_raw_Phi=saved['Phi'],paid_min_Phi=None if paid is None else paid['Phi'],
            same_candidate=paid is not None and saved['candidate']==paid['candidate'],
            same_Phi=paid is not None and saved['Phi']==paid['Phi'],raw_scoring_source=source))
    return dict(status='COMPLETE',tasks=len(rows),affected_predictions=sum(not row['same_candidate'] for row in rows),
        rows=rows,latent_source_issue='initial implementation chose proposal-phase raw; explicit paid-initial+proposal audit resolves canonical raw without rewriting predictions',
        prediction_files_modified=False,forward_LU=0,GT_access='NONE at raw selection audit')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--legacy-repo',type=Path,required=True)
    parser.add_argument('--history-root',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True);parser.add_argument('--seal',type=Path,required=True)
    parser.add_argument('--authorization',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.root=args.root.resolve();args.out=args.out.resolve()
    wall=time.perf_counter();cpu=time.process_time();manifest=read(args.manifest)
    seal,approval,validated=validate(args,manifest)
    raw=raw_audit(args,manifest)
    args.out.mkdir(parents=True,exist_ok=False)
    sys.path[:0]=[str(args.root/'src'),str(args.legacy_repo.resolve())]
    from run_vector_factorial import dump
    dump(args.out/'RAW_SELECTION_AUDIT.json',raw)
    dump(args.out/'AUTHORIZATION_AND_SEAL_VALIDATION.json',dict(status='VALIDATED_BEFORE_FIRST_TRUTH_ACCESS',
        validated_files=validated,approval=approval,seal=str(args.seal.resolve()),
        utc=datetime.now(timezone.utc).isoformat(),raw_audited_tasks=24,raw_affected_predictions=raw['affected_predictions']))
    import numpy as np
    from research.research_max_struct_r13.src.score_results import truth_spec,load_truth,material_metrics,physical_counts,table
    methods=manifest['methods'];cases=manifest['cases'];views=('raw','corrected','primary','secondary')
    truthdata={};truthrecords=[];pictures={};rows=[];prefixes=[];costs=[];corrections=[];incumbent_curves=[];failures=[]
    material_arrays=[];prediction_support=[];truth_arrays=[];truth_support=[]
    for case in cases:
        spec=truth_spec(args.history_root.resolve(),dict(case_id=case))
        accessed=datetime.now(timezone.utc).isoformat();truth=load_truth(spec)
        if truth is None:raise ValueError('Authorized inherited vector truth is absent: '+case)
        if truth['m'].shape!=(216,):raise ValueError('Inherited vector truth shape changed')
        truthdata[case]=truth
        truthrecords.append(dict(case_id=case,object_id=spec['object_id'],path=spec['path'],
            support_rule=spec['support_rule'],amplitude_rule=spec['amplitude_rule'],host=spec['host'],
            amplitude=truth['amplitude'],truth_accessed_at_utc=accessed,
            qualification='authorized offline old/exposed truth; no optimization/selection access'))
    raw_lookup={row['task_id']:row for row in raw['rows']}
    for task in manifest['tasks']:
        case=task['case_id'];method=task['method'];truth=truthdata[case]
        domain='crossgrid8_to6_repeated_frame' if 'crossgrid' in case else 'matched6_old_three_objects'
        object_id=case.split('_crossgrid')[0]
        folder=args.root/'RESULTS/vector4_v1'/task['task_id']/'attempt_01'
        receipt=read(folder/'RECEIPT.json');history=read(folder/'HISTORY.json')
        with np.load(folder/'CANDIDATE_POOL.npz',allow_pickle=False) as archive:pool=archive['m'].copy()
        if len(pool)!=len(history):raise ValueError('Candidate pool/history mismatch: '+task['task_id'])
        common=dict(task_id=task['task_id'],case_id=case,object_id=object_id,domain=domain,method=method,
            exposed_development=True,independent_geometries_overall=3,run_status=receipt['status'])
        cost=dict(common,body_reported_whole_process_CPU_seconds=receipt['whole_process_CPU_seconds'],
            actual_peak_RSS_bytes=receipt['actual_peak_RSS_bytes'],completed_material_cycles=receipt['completed_material_cycles'],
            saved_full_state_cycles=receipt['saved_full_state_cycles'],operation_counts_all=receipt['operation_counts_all'],
            **physical_counts(receipt['operation_counts_all']),CPU_cap=120,proposal_absolute_CPU_cap=72,corrector_CPU_allowance=48,
            corrector_termination=receipt['termination'].get('corrector'),timer_scope=receipt['timer_scope'])
        costs.append(cost)
        viewrows={}
        for view in views:
            source=raw_lookup[task['task_id']]['raw_scoring_source'] if view=='raw' else view
            path=folder/(source+'.npz');meta_path=folder/(source+'.json')
            m=metadata=None
            if path.exists() and meta_path.exists():
                with np.load(path,allow_pickle=False) as archive:m=archive['m'].copy()
                metadata=read(meta_path)
            metrics=material_metrics(m,truth,0.,8.)
            row=dict(common,view=view,source_view=source,**metrics,
                candidate=None if metadata is None else metadata['candidate'],
                D=None if metadata is None else metadata['D'],R=None if metadata is None else metadata['R'],
                Phi=None if metadata is None else metadata['Phi'],CPU_seconds=None if metadata is None else metadata['CPU_seconds'],
                relative_residual=None if metadata is None else metadata['relative_residual'],
                noise_feasible=None if metadata is None else metadata['noise_feasible'],
                noise_qualification=None if metadata is None else metadata['noise_qualification'],
                secondary_fallback_to_primary=bool(metadata.get('secondary_fallback_to_primary',False)) if metadata else None,
                selected_by='frozen original paid Phi / original ball-prior rule; no GT selection')
            rows.append(row);viewrows[view]=row
            pictures[(case,method,view)]=m
            material_arrays.append(np.full(216,np.nan) if m is None else m)
            prediction_support.append(np.zeros(216,bool) if m is None else m>1.)
            truth_arrays.append(truth['m']);truth_support.append(truth['support'])
        rawrow,corrected=viewrows['raw'],viewrows['corrected']
        corrections.append(dict(common,
            raw_Phi=rawrow['Phi'],corrected_Phi=corrected['Phi'],
            Phi_change=None if corrected['Phi'] is None or rawrow['Phi'] is None else corrected['Phi']-rawrow['Phi'],
            raw_D=rawrow['D'],corrected_D=corrected['D'],raw_R=rawrow['R'],corrected_R=corrected['R'],
            raw_NRMSE=rawrow.get('NRMSE'),corrected_NRMSE=corrected.get('NRMSE'),
            NRMSE_change=None if corrected.get('NRMSE') is None or rawrow.get('NRMSE') is None else corrected['NRMSE']-rawrow['NRMSE'],
            raw_IoU=rawrow.get('IoU'),corrected_IoU=corrected.get('IoU'),termination=receipt['termination'].get('corrector')))
        for point in (15,30,60,90,120):
            eligible=[record for record in history if record['eligible'] and record['legal']
                      and record['CPU_seconds']<=point and record['CPU_seconds']<=120]
            selected=min(eligible,key=lambda record:record['Phi']) if eligible else None
            index=None if selected is None else selected['candidate']
            metrics=material_metrics(None if index is None else pool[index],truth,0.,8.)
            prefixes.append(dict(common,view='paid_prefix',CPU_prefix=point,selection_candidate=index,**metrics,
                D=None if selected is None else selected['D'],R=None if selected is None else selected['R'],Phi=None if selected is None else selected['Phi'],
                actual_candidate_CPU_seconds=None if selected is None else selected['CPU_seconds'],
                run_already_terminal_before_prefix=receipt['whole_process_CPU_seconds']<point,
                interpretation='held incumbent after actual termination; no invented continued work or GT-best candidate'))
        best=None
        for record in history:
            if not record['eligible'] or not record['legal']:continue
            if best is None or record['Phi']<best['Phi']:
                best=record
                metrics=material_metrics(pool[record['candidate']],truth,0.,8.)
                incumbent_curves.append(dict(common,view='paid_incumbent',candidate=record['candidate'],
                    phase=record['phase'],CPU_seconds=record['CPU_seconds'],Phi=record['Phi'],D=record['D'],R=record['R'],**metrics))
    summary=[]
    for domain in ('matched6_old_three_objects','crossgrid8_to6_repeated_frame'):
        for method in methods:
            for view in views:
                subset=[row for row in rows if row['domain']==domain and row['method']==method and row['view']==view]
                available=[row for row in subset if row['available']]
                summary.append(dict(domain=domain,method=method,view=view,expected_conditions=3 if domain.startswith('matched') else 1,
                    available=len(available),independent_geometries_in_this_domain=3 if domain.startswith('matched') else 1,
                    independent_geometries_overall_not_sum_of_domains=3,
                    **{'mean_'+key:float(np.mean([row[key] for row in available])) if available else None
                       for key in ('NRMSE','IoU','RMSE','relativeL2','mean_amplitude_bias','D','R','Phi')}))
    for name,data in (('SCORES',rows),('PAID_PREFIX_SCORES',prefixes),('METHOD_VIEW_SUMMARY',summary),
        ('CORRECTOR_DECOMPOSITION',corrections),('COST_LEDGER',costs),('PAID_INCUMBENT_CURVES',incumbent_curves),
        ('TRUTH_SOURCES',truthrecords),('FAILURES',failures)):
        table(args.out,name,data)
    np.savez_compressed(args.out/'SCORED_MATERIAL_ROWS_AND_MASKS.npz',
        m=np.asarray(material_arrays),truth=np.asarray(truth_arrays),prediction_support=np.asarray(prediction_support),
        truth_support=np.asarray(truth_support),available=np.asarray([row['available'] for row in rows]),
        row_labels=np.asarray([row['task_id']+'__'+row['view'] for row in rows]),
        support_threshold=np.asarray(1.),box=np.asarray([0.,8.]))
    figure_manifest=render(args.out,rows,incumbent_curves,pictures,truthdata,costs,methods,cases,views,np)
    table(args.out,'FIGURE_MANIFEST',figure_manifest)
    queue=read(args.root/'queues/vector4_v1/attempt_01/RECEIPT.json')
    totals=dict(solver_children_whole_OS_CPU_seconds=queue['whole_children_OS_CPU_seconds'],
        all_controller_subprocess_OS_CPU_seconds=queue['all_controller_subprocess_OS_CPU_seconds'],
        controller_CPU_seconds=queue['controller_CPU_seconds'],controller_wall_seconds=queue['controller_wall_including_source_freeze_seconds'],
        sum_child_body_reported_CPU_seconds=sum(row['body_reported_whole_process_CPU_seconds'] for row in costs),
        physical_LU=sum(row['physical_LU'] for row in costs),
        physical_source_RHS=sum(row['physical_source_RHS'] for row in costs),
        physical_adjoint_RHS=sum(row['physical_adjoint_RHS'] for row in costs),
        physical_complex_RHS=sum(row['physical_complex_RHS'] for row in costs),
        max_child_actual_peak_RSS_bytes=max(row['actual_peak_RSS_bytes'] for row in costs),
        memory_scope='maximum of24 OS child high-water marks; not simultaneous aggregate resident memory',
        offline_scoring_physical_LU=0,offline_scoring_complex_RHS=0,CUDA_seconds=0)
    dump(args.out/'COST_TOTALS.json',totals)
    report(args.out,rows,summary,corrections,totals,raw,methods,cases)
    receipt=dict(status='COMPLETE',authorization=str(args.authorization.resolve()),seal=str(args.seal.resolve()),
        expected_tasks=24,scored_rows=len(rows),available_rows=sum(row['available'] for row in rows),
        paid_prefix_rows=len(prefixes),truth_files_accessed=len(truthrecords),independent_old_geometries=3,
        figure_count=len(figure_manifest),original_predictions_modified=False,original_seal_modified=False,
        GT_access='authorized offline only after validated unchangedseal/rawaudit',selection='no GT selection',
        forward_LU=0,complex_solver_RHS=0,CUDA_seconds=0,
        offline_CPU_seconds=time.process_time()-cpu,offline_whole_process_CPU_seconds=time.process_time(),
        offline_wall_seconds=time.perf_counter()-wall,
        offline_peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
        inheritance='R13 truth_spec/load_truth/material_metrics; original max_truth-.1 amplitude and chi>1 support',
        scientific_acceptance='parent scientific and figure review outstanding')
    dump(args.out/'SCORING_RECEIPT.json',receipt)
    print(json.dumps(dict(status='COMPLETE',rows=len(rows),available=receipt['available_rows'],
        figures=len(figure_manifest),offline_CPU_seconds=receipt['offline_CPU_seconds'],out=str(args.out))))
    return 0


def render(out,rows,curves,pictures,truths,costs,methods,cases,views,np):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    folder=out/'FIGURES';folder.mkdir();manifest=[]
    names={'pair00':'00 P0/P0 global','pair01':'01 P0/P1 global','pair10':'10 P1/P0 global','pair11':'11 P1/P1 global',
           'CIE_sourcewise':'CIE sourcewise','Dual_sourcewise':'Dual sourcewise'}
    colors={method:plt.get_cmap('tab10')(index) for index,method in enumerate(methods)}
    lookup={(row['case_id'],row['method'],row['view']):row for row in rows}
    def save(fig,name,scope):
        for extension in ('png','pdf'):fig.savefig(folder/(name+'.'+extension),dpi=170,bbox_inches='tight')
        plt.close(fig);manifest.append(dict(name=name,png=str(folder/(name+'.png')),pdf=str(folder/(name+'.pdf')),scope=scope))
    for case in cases:
        for view in views:
            fig,axes=plt.subplots(3,7,figsize=(20,8.6),squeeze=False,constrained_layout=True)
            values=[truths[case]['m']]+[pictures[(case,method,view)] for method in methods]
            for column,value in enumerate(values):
                metric=None if column==0 else lookup[(case,methods[column-1],view)]
                title='Offline GT' if column==0 else names[methods[column-1]]+f"\nNRMSE {metric['NRMSE']:.3f} | IoU {metric['IoU']:.3f}" if metric.get('NRMSE') is not None else names[methods[column-1]]+'\nUNAVAILABLE'
                planes=[] if value is None else [value.reshape(6,6,6)[:,:,3].T,value.reshape(6,6,6)[:,3,:].T,value.reshape(6,6,6)[3,:,:].T]
                for index,(plane_name,axes_name) in enumerate((('XY at z=+.02lambda',('x','y')),('XZ at y=+.02lambda',('x','z')),('YZ at x=+.02lambda',('y','z')))):
                    ax=axes[index,column]
                    if value is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center',transform=ax.transAxes)
                    else:ax.imshow(planes[index],origin='lower',cmap='viridis',vmin=0,vmax=8,interpolation='nearest',extent=[-.12,.12,-.12,.12])
                    ax.set_xticks([-.1,0,.1]);ax.set_yticks([-.1,0,.1]);ax.tick_params(labelsize=7)
                    ax.set_xlabel(axes_name[0]+' / lambda0',fontsize=8)
                    if column==0:ax.set_ylabel(plane_name+'\n'+axes_name[1]+' / lambda0',fontsize=8)
                    if index==0:ax.set_title(title,fontsize=9)
            fig.colorbar(ScalarMappable(norm=Normalize(0,8),cmap='viridis'),ax=axes.ravel().tolist(),shrink=.8,label='shared real chi = epsilon_r - 1; absolute scale [0,8]')
            qualification='repeated-frame crossgrid stress; not fourth object' if 'crossgrid' in case else 'exposed matched finite-DDA case'
            fig.suptitle(case+' | '+view+' | '+qualification+'\nTrue3D6^3 arrays: complete positive-central XY/XZ/YZ cell planes, index3; no interpolation',fontsize=12)
            save(fig,case+'__'+view+'__central_slices','GT and all6methods, 3complete central planes, fixed [0,8] scale; slice inspection does not show every voxel')
    fig,axes=plt.subplots(4,2,figsize=(14,15),squeeze=False,constrained_layout=True)
    for index,case in enumerate(cases):
        for method in methods:
            values=[row for row in curves if row['case_id']==case and row['method']==method]
            axes[index,0].step([row['CPU_seconds'] for row in values],[row['Phi'] for row in values],where='post',color=colors[method],label=names[method])
            axes[index,1].step([row['CPU_seconds'] for row in values],[row['NRMSE'] for row in values],where='post',color=colors[method],label=names[method])
        axes[index,0].set_yscale('log');axes[index,0].set_ylabel('Original vector Phi, log');axes[index,1].set_ylabel('Amplitude NRMSE of same paid Phi incumbent')
        for ax in axes[index]:ax.set_title(case,fontsize=10);ax.set_xlabel('actual cumulative child process CPU seconds');ax.grid(alpha=.2)
    axes[0,1].legend(fontsize=8);fig.suptitle('Paid original-Phi incumbent curves; no GT selection; trajectories end at actual convergence, not padded work',fontsize=12)
    save(fig,'PAID_COST_CURVES','four conditions separated; inherited vectorPhi and GT scored on the same frozen public-Phi incumbents')
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),constrained_layout=True)
    for method in methods:
        for case in cases:
            a=lookup[(case,method,'raw')];b=lookup[(case,method,'corrected')]
            marker='^' if 'crossgrid' in case else 'o'
            axes[0].scatter(a['Phi'],b['Phi'],color=colors[method],marker=marker,s=60,label=names[method] if case==cases[0] else None)
            axes[1].scatter(a['NRMSE'],b['NRMSE'],color=colors[method],marker=marker,s=60)
    for index,key in enumerate(('Phi','NRMSE')):
        values=[row[key] for row in rows if row['view'] in ('raw','corrected') and row.get(key) is not None]
        lower,upper=min(values),max(values);axes[index].plot([lower,upper],[lower,upper],'k--',lw=1)
        axes[index].set_xlabel('raw '+key);axes[index].set_ylabel('corrected '+key);axes[index].grid(alpha=.2)
    axes[0].set_xscale('log');axes[0].set_yscale('log');axes[0].legend(fontsize=8)
    from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
    axes[0].xaxis.set_major_locator(FixedLocator([2e-4,1e-3,5e-3]))
    axes[0].xaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.0e}'))
    axes[0].xaxis.set_minor_locator(NullLocator())
    fig.suptitle('Same original physical corrector: Phi and material changes are separate; triangles repeat the crossgrid frame')
    save(fig,'CORRECTOR_PHI_AND_MATERIAL_SCATTER','all24pairs raw->corrected; diagonal reference; lowerPhi need not lowermaterialerror')
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),constrained_layout=True)
    for method in methods:
        for case in cases:
            row=lookup[(case,method,'primary')];marker='^' if 'crossgrid' in case else 'o'
            axes[0].scatter(row['Phi'],row['NRMSE'],color=colors[method],marker=marker,s=60,label=names[method] if case==cases[0] else None)
            axes[1].scatter(row['relative_residual'],row['IoU'],color=colors[method],marker=marker,s=60)
    axes[0].set_xscale('log');axes[0].set_xlabel('original vector Phi, log');axes[0].set_ylabel('primary material amplitude NRMSE');axes[0].legend(fontsize=8)
    axes[1].set_xlabel('physical relative residual');axes[1].set_ylabel('fixed chi>1 support IoU');axes[1].set_ylim(-.03,1.03)
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Physical fit and material recovery reported separately; three independent oldobjects plus repeated-frame stress')
    save(fig,'PRIMARY_FIT_VS_MATERIAL_SCATTER','all6methods x4conditions; no pooled superiority/continuum claim')
    return manifest


def report(out,rows,summary,corrections,totals,raw,methods,cases):
    lookup={(row['case_id'],row['method'],row['view']):row for row in rows}
    lines=['# Vector4：封存后离线材料评分','',
        '24/24原任务完成后，主线程审阅固定公式/目标/费用/预测seal并另行授权GT评分。原seal和预测数组均未修改；本程序没有forward、optimizer、新调参或CUDA。',
        '', '## 合同和对象数','',
        '四组two-metric maps为已有CM+RR vector-DDA的3D几何对应，四组合都用原full analytic material derivative与L-BFGS-B inner12；不是R7 Arc或authorcode复现。强对照为原R13 CIE_sourcewise/Dual_sourcewise。每voxel一个实chi，216材料变量，共享xyz及60源；[0,8]box，原prior和per-source物理目标无1/2。',
        '', '四个条件包含三个已曝光独立几何；crossgrid8→6 frame重复旧frame，不算第四独立材料。matched/crossgrid分开汇总。strict48CPU corrector allowance和absolute120整进程cap保留；实际converged trajectories远早于很多公共prefix，后续prefix只保持已有paid incumbent，不伪造继续计算。',
        '', '## 指标与选择','',
        '完全复用R13 truth_spec/load_truth/material_metrics：amplitude NRMSE = whole-box RMSE / (max(coarse truth)-.1)，另报boxrange_NRMSE和relativeL2；固定prediction chi>1对truth chi>1支持；support内signed amplitude bias另报。没有新jointsuccess阈值。GT从未选择candidate或prefix。',
        '', f"raw incumbent审计24/24：受latent proposal-only输出问题影响的预测数为{raw['affected_predictions']}。本次全部raw都已经等于min(initial,proposal)的paid Phi incumbent；未改数组，latent来源与future一行report修复单独记录。",'',
        '## Primary：每个条件的完整六方法比较','']
    for case in cases:
        lines += ['### '+case,'','| method | NRMSE | IoU | Phi | physical residual |','|---|---:|---:|---:|---:|']
        for method in methods:
            row=lookup[(case,method,'primary')]
            lines.append(f"| {method} | {row['NRMSE']:.6f} | {row['IoU']:.6f} | {row['Phi']:.9g} | {row['relative_residual']:.6f} |")
        lines.append('')
    harmed=[row for row in corrections if row['NRMSE_change'] is not None and row['NRMSE_change']>0]
    lines += ['## Corrector与负面结果','',
        f'24条中，原物理corrector后material NRMSE增加的条数为{len(harmed)}；D/R/Phi变化与material变化逐条保留在CORRECTOR_DECOMPOSITION。不得把较低Phi直接称为材料恢复提升。各方法全部raw/corrected/primary/secondary和secondaryfallback均有完整行，未按GT挑选阶段。',
        '', '## 实际费用','',
        f"24个solver children完整OS CPU合计 {totals['solver_children_whole_OS_CPU_seconds']:.6f}秒；controller CPU {totals['controller_CPU_seconds']:.6f}秒，含sourcefreeze的controller wall {totals['controller_wall_seconds']:.6f}秒。子receipt累计CPU {totals['sum_child_body_reported_CPU_seconds']:.6f}秒是另一个body/serialization边界，不能与OS总数重复相加。",
        '', f"实际physical LU={totals['physical_LU']}，source RHS={totals['physical_source_RHS']}，adjoint RHS={totals['physical_adjoint_RHS']}，合计disjoint solver RHS={totals['physical_complex_RHS']}。basis/densecurrent/materialHop/setup/failed/output记录在原run COST_LEDGER和本表nestedcounts；不要再次加nested delta。",
        '', '离线评分/图形CPU、wall、RSS单列于SCORING_RECEIPT；离线physical LU/RHS=0。没有GPU/CUDA；旧Windows/其他物理域时间不并入本地速度结论。',
        '', '## 图形与范围','',
        '每个条件的raw/corrected/primary/secondary均显示GT+全部六方法的完整XY/XZ/YZ正中央cellplane(index3, +.02lambda)，共享绝对chi[0,8]色标与物理坐标；不插值。图像来自真正3D648-current/216-material数组，centralplane并不展示每个voxel，也不代替全volume/continuumreview。',
        '', 'PAID_COST_CURVES仅画按原Phi选择的paid incumbent；CORRECTOR_PHI_AND_MATERIAL_SCATTER分离fit与material；PRIMARY_FIT_VS_MATERIAL_SCATTER保留同一case不同物理解。跨格噪声球不覆盖模型失配；没有新Mie/PDE证明。',
        '', '最终科学解释、所有图形审阅和publication验收仍由主线程负责。']
    (out/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf8')


if __name__=='__main__':sys.exit(main())
