"""Optional saved-output probability figures: no inference, physics or truth I/O.

Read whitelist: P1A/reference/{RESULTS,GATE}.json,
P1B/{inference,point_validation}/RESULTS.json, and only completed offline score
SUMMARY/PARENT_SCORES/PAIRED_CI files. No observations, rules, parent parameters,
generation plans, TRUTH_ACCESS or OFFLINE_TRUTH file is opened here.
"""
from pathlib import Path
import json
import math
import numpy as np
import matplotlib.pyplot as plt

ARMS=('direct_point','progressive_point','direct_SMC_mean','progressive_SMC_mean')


def build_probability_figures(root,sources,errors,finish,table,read_json,atomic,write_csv):
    root=Path(root).resolve();used=[];fig_sources=[]
    def load(relative):
        path=root/relative
        if not path.exists():return None
        # All callers below use literal public-result/score names only.
        used.append(relative)
        return read_json(path,root,sources,errors)
    def records(value):
        if isinstance(value,list):return [r for r in value if isinstance(r,dict)]
        if isinstance(value,dict):
            for key in ('rows','results','observations'):
                if isinstance(value.get(key),list):return records(value[key])
            return [value]
        return []
    def num(v):return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
    def tagged(row,**context):
        result=dict(row)
        for key,value in context.items():
            if key in result and result[key]!=value:
                if key=='source':result['recorded_source']=result[key]
                else:errors.append(dict(path=inf_path,kind='DisplayContextMismatch',message=f'{key}: saved row differs from enclosing group'))
            result[key]=value
        return result
    def finish_plot(fig,name,title,caption,inputs):
        finish(fig,name,title,caption)
        fig_sources.append(dict(figure=f'figures/saved_results/{name}.png',SVG=f'figures/saved_results/{name}.svg',saved_inputs=inputs))
    def blank(name,title,inputs,note='NOT_RUN — required saved output fields are absent'):
        fig,ax=plt.subplots(figsize=(8,3.2));ax.text(.5,.5,note,ha='center',va='center',transform=ax.transAxes,wrap=True);ax.set_axis_off()
        finish_plot(fig,name,title,'缺失结果保持NOT_RUN；空白占位不表示该阶段完成，不生成概率或补全轨迹。',inputs)

    ref_path='results/P1A/reference/RESULTS.json';gate_path='results/P1A/reference/GATE.json'
    inf_path='results/P1B/inference/RESULTS.json';point_path='results/P1B/point_validation/RESULTS.json'
    references=records(load(ref_path));gate=load(gate_path);inference=records(load(inf_path));validation=records(load(point_path))
    point_lookup={r.get('observation'):r for r in validation if 'observation' in r}
    stage_rows=[dict(stage=name,status='SAVED_RECORDS' if values else 'NOT_RUN',saved_records=len(values)) for name,values in (('P1A_reference',references),('P1B_inference',inference),('P1B_point_validation',validation))]
    masses=[];order_errors=[];points=[];smc_rows=[];smc_runs=[];cdfs=[]
    for r in references:
        obs=r.get('observation');p=r.get('branch1_posterior_mass')
        if num(p):masses.append(dict(cohort='P1A',observation=obs,reference='continuous_order9',branch1_mass=p,status=r.get('reference_status')))
        order_errors.append(dict(cohort='P1A',observation=obs,event_5_9_error=r.get('quadrature_event_error'),CDF_5_9_max_error=r.get('quadrature_CDF_max_error'),status=r.get('reference_status')))
    for r in inference:
        obs=r.get('observation');continuous=r.get('continuous_ref',{});finite=r.get('finite_ref',{});ensembles=r.get('ensemble',{})
        for label,summary in (('continuous_order9',continuous),('finite_enumeration',finite)):
            mass=summary.get('branch_mass')
            if isinstance(mass,list) and len(mass)==2 and all(num(v) for v in mass):
                masses.append(dict(cohort='P1B',observation=obs,reference=label,branch1_mass=mass[1],status=summary.get('status')))
        order_errors.append(dict(cohort='P1B',observation=obs,event_5_9_error=continuous.get('event_order_error'),CDF_5_9_max_error=continuous.get('CDF_order_max_error'),status=continuous.get('status')))
        methods=r.get('methods',{});validated=point_lookup.get(obs,{}).get('methods',{})
        for name,point in methods.items():
            if not isinstance(point,dict):continue
            actual=validated.get(name,point)
            points.append(dict(observation=obs,method=name,branch=actual.get('branch'),x=actual.get('x'),Phi=actual.get('Phi'),full_validation_status=actual.get('full_validation_status','NEEDS_FULL_VALIDATION'),fallback_used=actual.get('fallback_used'),weight_semantics='POINT_DECISION_NO_POSTERIOR_MASS'))
        for schedule,summary in ensembles.items():
            if not isinstance(summary,dict):continue
            smc_rows.append(tagged({k:v for k,v in summary.items() if not isinstance(v,(dict,list))},observation=obs,schedule=schedule))
            mass=summary.get('average_branch_mass')
            if isinstance(mass,list) and len(mass)==2 and all(num(v) for v in mass):masses.append(dict(cohort='P1B',observation=obs,reference=schedule+'_SMC_ensemble',branch1_mass=mass[1],status=summary.get('probability_claim')))
        for schedule,runs in r.get('SMC_runs',{}).items():
            if isinstance(runs,list):
                for run in runs:
                    if isinstance(run,dict):smc_runs.append(tagged({k:v for k,v in run.items() if not isinstance(v,(dict,list))},observation=obs,schedule=schedule))
        for cut,values in continuous.get('CDF',{}).items():
            if not isinstance(values,dict):continue
            try:x=float(cut)
            except(ValueError,TypeError):continue
            row=dict(observation=obs,cut=x,continuous_order5=values.get('order5'),continuous_order9=values.get('order9'),finite=finite.get('CDF',{}).get(cut))
            for schedule,summary in ensembles.items():
                if isinstance(summary,dict):row[schedule]=summary.get('average_CDF',{}).get(cut)
            cdfs.append(row)

    if masses:
        cohorts=sorted({r['cohort'] for r in masses});fig,axes=plt.subplots(len(cohorts),1,figsize=(10,3.5*len(cohorts)),squeeze=False)
        for ax,cohort in zip(axes[:,0],cohorts):
            for label in sorted({r['reference'] for r in masses if r['cohort']==cohort}):
                rs=[r for r in masses if r['cohort']==cohort and r['reference']==label]
                ax.plot([r['observation'] for r in rs],[r['branch1_mass'] for r in rs],'o-',ms=3,label=label)
            ax.set_ylim(-.03,1.03);ax.set_title(cohort+' saved observations (separate cohort)');ax.set_xlabel('Saved observation ID');ax.set_ylabel('Posterior mass of branch1');ax.legend(fontsize=8)
        finish_plot(fig,'P_branch_posterior_mass','受控新prior：真实branch posterior mass','P1A与P1B观察群体分面，不把同编号误当同一数据。仅使用已保存后验质量；profile点从不赋予后验质量。资格标志原样保留在CSV。',[ref_path,inf_path])
    else:blank('P_branch_posterior_mass','Branch posterior mass',[ref_path,inf_path])
    usable=[r for r in order_errors if num(r.get('event_5_9_error')) or num(r.get('CDF_5_9_max_error'))]
    if usable:
        fig,axes=plt.subplots(1,2,figsize=(10,4))
        for ax,key in zip(axes,('event_5_9_error','CDF_5_9_max_error')):
            for cohort in sorted({r['cohort'] for r in usable}):
                rs=[r for r in usable if r['cohort']==cohort and num(r.get(key))]
                if rs:ax.plot([r['observation'] for r in rs],[r[key] for r in rs],'o-',ms=3,label=cohort)
            ax.set_title(key);ax.set_xlabel('Cohort-specific observation ID');ax.set_ylabel('Saved absolute order5/order9 difference');ax.legend()
        finish_plot(fig,'P_reference_order_errors','概率参考：5/9阶误差','只画保存的branch事件与CDF阶数变化；这是经验求积误差诊断，不称连续精确后验。未通过资格的观察不删除。',[ref_path,inf_path,gate_path])
    else:blank('P_reference_order_errors','概率参考：5/9阶误差',[ref_path,inf_path,gate_path])

    p1b_mass=[r for r in masses if r['cohort']=='P1B' and r['reference']=='finite_enumeration']
    if points and p1b_mass:
        fig,ax=plt.subplots(figsize=(10,4.5));ax.plot([r['observation'] for r in p1b_mass],[r['branch1_mass'] for r in p1b_mass],'o-',label='finite marginal branch1 mass')
        for method in sorted({r['method'] for r in points}):
            rs=[r for r in points if r['method']==method and num(r.get('branch'))]
            ax.scatter([r['observation'] for r in rs],[r['branch'] for r in rs],marker='x',s=50,label=method+' branch DECISION')
        ax.set_ylim(-.08,1.08);ax.set_xlabel('P1B observation');ax.set_ylabel('Mass (line) / selected branch label (cross)');ax.legend(fontsize=8)
        finish_plot(fig,'P_profile_vs_marginal_branch','Profile点与marginal branch分开','叉号只表示profile ML/MAP点选择了哪个离散branch，不是质量0/1的后验。可用full point validation替换近似点；fallback和资格状态保留。',[inf_path,point_path])
    else:blank('P_profile_vs_marginal_branch','Profile点与marginal branch',[inf_path,point_path])
    if cdfs:
        ids=list(dict.fromkeys(r['observation'] for r in cdfs))[:4];fig,axes=plt.subplots(1,len(ids),figsize=(4*len(ids),4),squeeze=False)
        for ax,obs in zip(axes[0],ids):
            rs=sorted((r for r in cdfs if r['observation']==obs),key=lambda r:r['cut'])
            for key in ('continuous_order5','continuous_order9','finite','direct','progressive'):
                valid=[r for r in rs if num(r.get(key))]
                if valid:ax.plot([r['cut'] for r in valid],[r[key] for r in valid],'o-',ms=3,label=key)
            for point in points:
                if point['observation']==obs and num(point.get('x')):ax.axvline(point['x'],ls=':',lw=1,label=point['method']+' point')
            ax.set_ylim(-.03,1.03);ax.set_title(f'observation {obs}');ax.set_xlabel('Declared CDF checkpoint x');ax.set_ylabel('Saved CDF');ax.legend(fontsize=6)
        finish_plot(fig,'P_CDF_and_profile_points','保存CDF与profile点位置','按源文件顺序只展示前4个观察，其余完整输出到CSV；线连接声明检查点，不生成密度或未保存CDF值。竖线是点估计，不是posterior mean材料世界。',[inf_path,point_path])
    else:blank('P_CDF_and_profile_points','保存CDF与profile点位置',[inf_path,point_path])
    if smc_rows:
        metrics=[('branch_mass_max_error_vs_finite','Branch mass error vs exact finite reference'),('CDF_max_error_vs_finite','CDF error vs exact finite reference'),('average_surviving_initial_ancestors','Average surviving INITIAL ancestors'),('branch_MISS_run_rate','Eligible-branch MISS run rate')]
        fig,axes=plt.subplots(2,2,figsize=(11,7))
        for ax,(key,title) in zip(axes.ravel(),metrics):
            count=0
            for schedule in ('direct','progressive'):
                rs=[r for r in smc_rows if r['schedule']==schedule and num(r.get(key))]
                if rs:ax.plot([r['observation'] for r in rs],[r[key] for r in rs],'o-',ms=3,label=schedule);count+=len(rs)
            if not count:ax.text(.5,.5,'NOT_RECORDED',ha='center',transform=ax.transAxes)
            ax.set_title(title,fontsize=10);ax.set_xlabel('Observation ID')
            if count:ax.legend(fontsize=8)
        finish_plot(fig,'P_SMC_accuracy_ancestry_misses','SMC direct/progressive：参考误差与退化','只比较同一观察、同一finite target的保存指标；同seed重复运行不算独立parent。漏分支的资格阈值与reference标志沿用源输出；不把ancestor存活数或SMC结束当方法成功。',[inf_path])
    else:blank('P_SMC_accuracy_ancestry_misses','SMC direct/progressive对照',[inf_path])

    # Completed scorer output only. Never follow truth provenance paths.
    offline_dirs=[]
    for base in (root/'results/P1B/offline',root/'results/P1B/offline_scores'):
        if base.exists():offline_dirs.extend([base]+[p for p in base.iterdir() if p.is_dir()])
    completed=[]
    for folder in offline_dirs:
        marker=folder/'COMPLETE.json'
        if marker.exists():
            record=load(str(marker.relative_to(root)))
            if isinstance(record,dict) and record.get('status')=='OFFLINE_SCORES_READY_FOR_PARENT_REVIEW':completed.append(folder)
    parent_scores=[];contrasts=[];offline_summaries=[]
    for folder in completed:
        prefix=str(folder.relative_to(root));summary=load(prefix+'/SUMMARY.json');parents=records(load(prefix+'/PARENT_SCORES.json'));paired=load(prefix+'/PAIRED_CI.json')
        if isinstance(summary,dict):offline_summaries.append(tagged(summary,source=prefix+'/SUMMARY.json'))
        parent_scores.extend(tagged(r,source=prefix+'/PARENT_SCORES.json') for r in parents)
        if isinstance(paired,dict):
            for family,values in paired.items():
                if isinstance(values,list):contrasts.extend(tagged(r,source=prefix+'/PAIRED_CI.json',contrast_family=family) for r in values if isinstance(r,dict))
        elif isinstance(paired,list):contrasts.extend(tagged(r,source=prefix+'/PAIRED_CI.json') for r in paired if isinstance(r,dict))
    stage_rows.append(dict(stage='P1B_offline_scoring',status='COMPLETED_SCORER_FILES_ONLY' if completed else 'NOT_RUN',saved_score_folders=len(completed)))
    risk=[r for r in parent_scores if r.get('arm') in ARMS and num(r.get('position_squared_error'))]
    if risk:
        fig,ax=plt.subplots(figsize=(10,4.5))
        for arm in ARMS:
            rs=sorted((r for r in risk if r['arm']==arm),key=lambda r:r['parent'])
            if rs:ax.plot([r['parent'] for r in rs],[r['position_squared_error'] for r in rs],'o-',ms=3,label=arm)
        ax.set_xlabel('Independent parent ID (noise replicas already averaged)');ax.set_ylabel('Saved parent mean squared position loss');ax.legend(fontsize=8)
        finish_plot(fig,'P_offline_parent_risk_four_arms','四臂离线parent风险','仅从已COMPLETE的评分PARENT_SCORES读取，绝不读取parent真值。完整parent均值与缺失覆盖按源表保留；不能把两次噪声重复当两个独立parent。',[str(p.relative_to(root))+'/PARENT_SCORES.json' for p in completed])
    else:blank('P_offline_parent_risk_four_arms','四臂离线parent风险',['completed offline/PARENT_SCORES.json'])
    ci=[r for r in contrasts if r.get('loss')=='position_squared_error' and num(r.get('gain_a_minus_b')) and isinstance(r.get('paired_percentile_CI95'),list) and len(r['paired_percentile_CI95'])==2 and all(num(v) for v in r['paired_percentile_CI95'])]
    if ci:
        fig,ax=plt.subplots(figsize=(10,max(3.5,.65*len(ci))))
        for i,r in enumerate(ci):
            lo,hi=r['paired_percentile_CI95'];ax.hlines(i,lo,hi,color='#0369a1');ax.scatter(r['gain_a_minus_b'],i,color='#0369a1')
        ax.axvline(0,color='gray',ls=':');ax.set_yticks(range(len(ci)),[f"{r.get('contrast_family','')} | {r['arm_a']} minus {r['arm_b']}" for r in ci],fontsize=8)
        ax.set_xlabel('Saved paired parent loss gain a-b (positive favors b)')
        finish_plot(fig,'P_offline_paired_CI','离线paired parent CI','直接画已完成评分提供的percentile95%区间，不重新bootstrap、不静默删除不完整parent。描述性区间不是自动GO；多假设收益与progressive schedule收益分组保留。',[str(p.relative_to(root))+'/PAIRED_CI.json' for p in completed])
    else:blank('P_offline_paired_CI','离线paired parent CI',['completed offline/PAIRED_CI.json'],'NOT_RUN / no saved complete paired-parent interval')
    outputs=dict(branch_mass=masses,reference_order_errors=order_errors,profile_points=points,CDF_checkpoints=cdfs,SMC_ensemble=smc_rows,SMC_runs=smc_runs,offline_parent_scores=parent_scores,offline_contrasts=contrasts)
    for name,data in outputs.items():write_csv(root/'tables/saved_results'/f'P_{name}.csv',data)
    snapshot=dict(status='SEALED_SAVED_RESULTS' if (root/'CLOSED.json').exists() else 'WORK_IN_PROGRESS',new_physics=0,new_inference=0,truth_access=False,read_whitelist_only=True,stage_status=stage_rows,P1A_saved_gate=gate,offline_score_summaries=offline_summaries,offline_score_folders=[str(p.relative_to(root)) for p in completed],source_files=used,scope='new controlled fixed-support prior only; not arbitrary walls; original minimum_feature UNRESOLVED; outside-family mass unquantified',**outputs)
    atomic(root/'tables/saved_results/PROBABILITY_DISPLAY_SNAPSHOT.json',json.dumps(snapshot,ensure_ascii=False,indent=2,allow_nan=False))
    atomic(root/'tables/saved_results/PROBABILITY_FIGURE_SOURCES.json',json.dumps(fig_sources,ensure_ascii=False,indent=2))
    section='<section><h2>概率参考与四臂结果 · 可选保存输出视图</h2><p>只针对受控新fixed-support prior，非任意墙体；原minimum_feature语义UNRESOLVED。profile点不等于marginal/evidence；无GT读取、无推断或PDE执行，无最终科学GO。</p>'+table(stage_rows)+'<h3>保存的参考资格状态</h3>'+table([gate] if isinstance(gate,dict) else [])+'<h3>Profile点与full validation</h3>'+table(points,['observation','method','branch','x','Phi','full_validation_status','fallback_used'])+'<h3>离线CI完整性</h3>'+table(contrasts,['contrast_family','arm_a','arm_b','loss','status','eligible_parent_count','parent_count','gain_a_minus_b','paired_percentile_CI95'])+'<p><a href="tables/saved_results/PROBABILITY_DISPLAY_SNAPSHOT.json">数据及状态快照</a> · <a href="tables/saved_results/PROBABILITY_FIGURE_SOURCES.json">逐图来源</a>。缺失NOT_RUN不表示全阶段完成。</p></section>'
    return dict(data=snapshot,html=section)
