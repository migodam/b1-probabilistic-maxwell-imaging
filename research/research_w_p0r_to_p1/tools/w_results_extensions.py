"""Saved-array C/B2 display extension. No physics, Book or queue imports."""
from pathlib import Path
import csv
import json
import math
import numpy as np
import matplotlib.pyplot as plt


def build_extensions(root, sources, errors, finish, table, read_json, atomic, write_csv):
    phase=root/'results/C/phase_diagram'
    def load(path):return read_json(path,root,sources,errors)
    def fail(path,message):errors.append(dict(path=str(path.relative_to(root)),kind='SavedEvidenceError',message=message))
    def rows(stem):
        path=stem.with_suffix('.json')
        if path.exists():
            obj=load(path)
            if isinstance(obj,dict):obj=obj.get('rows',[obj])
            if isinstance(obj,list) and all(isinstance(r,dict) for r in obj):return obj
            fail(path,'Expected list of row dictionaries or dict.rows');return []
        path=stem.with_suffix('.csv')
        if path.exists():
            try:
                def value(v):
                    if v=='':return None
                    if v in ('True','False'):return v=='True'
                    try:return json.loads(v)
                    except(ValueError,TypeError):return v
                with path.open(newline='') as f:result=[{k:value(v) for k,v in r.items()} for r in csv.DictReader(f)]
                sources.append(dict(path=str(path.relative_to(root)),role='saved CSV fallback'))
                return result
            except(OSError,ValueError,TypeError) as exc:fail(path,str(exc))
        return []
    def state(folder):
        if (folder/'RESULTS.json').exists():return 'SAVED_RESULTS'
        if (folder/'FAILURE.json').exists():return 'STOPPED_PARTIAL'
        return 'RUNNING' if folder.exists() else 'NOT_RUN'
    def blank(name,title,status,caption):
        fig,ax=plt.subplots(figsize=(8,3.2));ax.text(.5,.5,status,transform=ax.transAxes,ha='center',va='center',fontsize=17)
        ax.set_axis_off();finish(fig,name,title,caption)
    def saved_npz(path,keys,allowed):
        path=Path(path).resolve()
        if not path.is_relative_to(allowed.resolve()):raise ValueError('saved response path lies outside its named result family')
        with np.load(path,allow_pickle=False) as z:out=[np.asarray(z[k]) for k in keys]
        if not all(np.isfinite(a).all() for a in out):raise ValueError('nonfinite saved array')
        sources.append(dict(path=str(path.relative_to(root)),role='saved response array replay only'))
        return out
    plan=load(phase/'PLAN.json') if (phase/'PLAN.json').exists() else {}
    plan=plan or {};sigma=float(plan.get('sigma',.01))
    subset=plan.get('nested_aperture',dict(tx=[1,2],rx=[3,4,5,6,7,8]))
    phases=rows(phase/'PHASE_TABLE');checks=[]
    for row in phases:
        try:
            a,b=[saved_npz(Path(row[k]),['data'],phase)[0] for k in ('response0','response1')]
            if a.shape!=b.shape or a.ndim!=3:raise ValueError('phase response shape mismatch')
            delta=(b-a)/sigma;full=float(np.linalg.norm(delta));central=float(np.linalg.norm(delta[:,subset['tx'],:][:,:,subset['rx']]))
            if central>full+1e-10*max(1.,full):raise ValueError('nested subset norm exceeds full norm')
            if max(abs(full-row['d_mean_all']),abs(central-row['d_mean_central_nested']))>1e-8*max(1.,full):raise ValueError('array replay disagrees with saved phase distances')
            checks.append(dict(case=row['case'],band=row['band'],loss_scale=row['common_wall_loss_scale'],full_from_array=full,central_from_array=central,replay='MATCH',source0=row['response0'],source1=row['response1']))
        except(OSError,ValueError,KeyError,IndexError) as exc:fail(phase/'PHASE_TABLE.json',str(exc))
    # Only array-replayed rows enter scientific plots; raw records remain in tables.
    key=lambda r:(r['case'],r['band'],r.get('common_wall_loss_scale',r.get('loss_scale')))
    valid={key(r) for r in checks};plot_rows=[r for r in phases if key(r) in valid]
    phase_status=state(phase)
    if plot_rows:
        cases=['cal','joint_contrast3','joint_contrast7'];bands=['narrow','wide']
        fig,ax=plt.subplots(figsize=(9,4.5));x=np.arange(len(cases));width=.34
        for i,band in enumerate(bands):
            vals=[next((r['d_mean_all'] for r in plot_rows if r['case']==case and r['band']==band and r['common_wall_loss_scale']==1.),np.nan) for case in cases]
            bars=ax.bar(x+(i-.5)*width,vals,width,label=band)
            for bar,val in zip(bars,vals):
                if math.isfinite(val):ax.text(bar.get_x()+bar.get_width()/2,val,f'{val:.4g}',ha='center',va='bottom',fontsize=9)
        ax.set_xticks(x,['cal: query air','joint: contrast 3','joint: contrast 7']);ax.set_ylabel('Saved full-aperture mean-response distance / .01');ax.legend()
        ax.set_title('Same fixed pair, different non-nested frequency sets')
        finish(fig,'C_narrow_wide_cal_joint','C：narrow/wide校准与joint','来自PHASE_TABLE并由response0/1实际数组复算。频率列表非嵌套，带宽比较不具有追加观测的单调性；同一材料对未重新优化，Oracle非在线。')
        fig,ax=plt.subplots(figsize=(7,5))
        for band,color in (('narrow','#0369a1'),('wide','#c2410c')):
            rs=[r for r in plot_rows if r['band']==band]
            ax.scatter([r['d_mean_all'] for r in rs],[r['d_mean_central_nested'] for r in rs],label=band,color=color)
        maximum=max(r['d_mean_all'] for r in plot_rows)*1.08;ax.plot([0,maximum],[0,maximum],':',color='gray')
        ax.set_xlabel('Full aperture: 240 complex channels');ax.set_ylabel('Central strict subset: 60 complex channels');ax.legend()
        ax.set_title('Identical response arrays; nested channel selection')
        finish(fig,'C_strict_nested_aperture','C：严格嵌套central/full孔径','从每个相同response数组抽取tx[1,2]与rx[3..8]，60通道严格包含于240通道；无按通道数重新归一化。此处的范数关系来自嵌套子集，不能移用于非嵌套频率。')
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
        for ax,band in zip(axes,bands):
            for case in ('cal','joint_contrast3','fixed_target0_contrast3','joint_contrast7'):
                rs=sorted((r for r in plot_rows if r['band']==band and r['case']==case),key=lambda r:r['common_wall_loss_scale'])
                if rs:ax.plot([r['common_wall_loss_scale'] for r in rs],[r['d_mean_all'] for r in rs],'o-',label=case)
            ax.set_title(band);ax.set_xlabel('Common wall conductivity scale');ax.legend(fontsize=8)
        axes[0].set_ylabel('Full-aperture fixed-pair data distance')
        finish(fig,'C_wall_loss_scale','C：壁体损耗敏感性','同一已保存材料对，仅共同缩放墙体conductivity；目标conductivity仍为.02。曲线是固定对的响应，不表示每个损耗下重新优化、最优下界或在线恢复。缺失设置不插补。')
    else:
        for name,title in [('C_narrow_wide_cal_joint','C：narrow/wide校准与joint'),('C_strict_nested_aperture','C：嵌套孔径'),('C_wall_loss_scale','C：损耗')]:blank(name,title,phase_status,'暂无可由实际数组复算的保存相图点。')

    born_endpoints=rows(phase/'BORN_ENDPOINTS');born_joint=rows(phase/'BORN_JOINT_COMPARISON');born_checks=[]
    for r in born_joint:
        path=phase/'born'/f"JOINT_contrast{int(r['contrast'])}.npz"
        try:
            b0,b1,t0,t1=saved_npz(path,['born0','born1','true0','true1'],phase)
            d1=float(np.linalg.norm(b1-b0)/sigma);d2=float(np.linalg.norm(t1-t0)/sigma);error=float(np.linalg.norm((b1-b0)-(t1-t0))/sigma)
            if max(abs(d1-r['born_joint_d_mean']),abs(d2-r['full_joint_d_mean']),abs(error-r['joint_difference_Born_error_white']))>1e-8*max(1.,d1,d2):raise ValueError('Born/full saved arrays disagree with JSON')
            born_checks.append(dict(contrast=r['contrast'],born=d1,full=d2,error=error,replay='MATCH',source=str(path.relative_to(root))))
        except(OSError,ValueError,KeyError) as exc:fail(path,str(exc))
    if born_checks:
        fig,axes=plt.subplots(1,2,figsize=(10,4.5));x=np.arange(len(born_checks));width=.35
        axes[0].bar(x-width/2,[r['born'] for r in born_checks],width,label='distorted-wave target Born')
        axes[0].bar(x+width/2,[r['full'] for r in born_checks],width,label='full target response');axes[0].legend(fontsize=8)
        axes[1].bar(x,[r['error'] for r in born_checks]);axes[1].set_title('Norm of Born/full pair-difference error')
        for ax in axes:ax.set_xticks(x,[f"contrast {r['contrast']}" for r in born_checks]);ax.set_ylabel('Whitened distance')
        finish(fig,'C_Born_vs_full','C：Born/full实际对照','同一fixed pair和同一平移目标，读取born/JOINT_contrast*.npz四组实际保存预测复算；Born为墙体背景上的目标一阶截断，不是把墙体全部设空气。无新求解。')
    else:blank('C_Born_vs_full','C：Born/full实际对照','RUNNING' if phase.exists() and phase_status!='SAVED_RESULTS' else 'NOT_RUN','未生成完整Born joint数组与对照JSON时不造数；已有endpoint结果仍在下表保留。')

    wide=root/'results/C/wide_reoptimized';wide_record=load(wide/'RESULTS.json') if (wide/'RESULTS.json').exists() else None
    if isinstance(wide_record,dict):
        path=wide/'FINAL.npz'
        try:
            initial,final,reference=saved_npz(path,['initial_data','data','reference'],wide)
            raw=float(np.linalg.norm(initial-reference)/sigma);optimized=float(np.linalg.norm(final-reference)/sigma)
            if max(abs(raw-wide_record['fixed_pair_wide_distance']),abs(optimized-wide_record['reoptimized_pair_wide_distance']))>1e-8*max(1.,raw):raise ValueError('wide FINAL arrays disagree with result')
            fig,ax=plt.subplots(figsize=(7,4));ax.bar(['fixed narrow-derived pair','wide reoptimized retained pair'],[raw,optimized]);ax.set_ylabel('Full-wide data distance');ax.set_title('Oracle existence; target shift remains 0.1')
            finish(fig,'C_wide_reoptimized','C：wide重新优化实际结果','从wide_reoptimized/FINAL.npz实际initial_data/data/reference复算。更新后的墙体不替换早先narrow见证，也不产生不存在性或全局最优结论。')
        except(OSError,ValueError,KeyError) as exc:fail(path,str(exc))
    else:blank('C_wide_reoptimized','C：wide重新优化实际结果',state(wide),'仅在实际RESULTS与FINAL数组存在时展示完成结果；缺失不由相图推断。')

    families=[('native_cal','native_cal_*'),('native_joint','native_joint_direct_*'),('native_cavity','native_cavity_direct_*'),('task_preserving','task_preserving_*')]
    branch_rows=[];task_changes=[]
    for family,pattern in families:
        folders=sorted((root/'results/B2').glob(pattern))
        if not folders:branch_rows.append(dict(family=family,status='NOT_RUN',source=None));continue
        for folder in folders:
            branch_rows.append(dict(family=family,folder=folder.name,status=state(folder),source=str((folder/'RESULTS.json').relative_to(root))))
            records=rows(folder/'RESULTS')
            for r in records:
                retained=r.get('retained_result');r=dict(r,**retained) if isinstance(retained,dict) else r
                task_changes.append(dict(family=family,folder=folder.name,method=r.get('method'),status=r.get('status'),raw_d_cal=r.get('raw_d_cal'),raw_d_T=r.get('raw_d_T'),d_cal=r.get('d_cal'),d_T=r.get('d_T'),edited_observation_distance=r.get('final_distance'),edit_kind=r.get('kind',family),stopping=r.get('stopping',r.get('stop')),identity=r.get('identity'),source=str((folder/'RESULTS.json').relative_to(root))))
            if family=='task_preserving':
                for name in ('PROPOSAL','FINAL','RETAINED_ACCEPTED'):
                    if (folder/f'{name}.json').exists():
                        r=load(folder/f'{name}.json')
                        if isinstance(r,dict):task_changes.append(dict(family=family,folder=folder.name,method=name,status='SAVED_PAIRED_POINT',d_cal=r.get('d_cal'),d_T=r.get('d_T'),source=str((folder/f'{name}.json').relative_to(root)),identity=r.get('identity')))
    for r in task_changes:
        raw,final=r.get('raw_d_T'),r.get('d_T')
        r['task_response_change']='NOT_RECORDED'
        if isinstance(raw,(int,float)) and isinstance(final,(int,float)):
            r['task_response_change']='raw>=3 to final<3' if raw>=3>final else 'retained >=3' if min(raw,final)>=3 else 'raw<3 to final>=3' if raw<3<=final else 'both below3'
    task_array_checks=[]
    for r in task_changes:
        if r['family']!='native_cal' or not r.get('method'):continue
        folder=root/'results/B2'/r['folder'];path=folder/f"{r['method']}_FINAL.npz"
        if not path.exists():continue
        try:
            task,base,proposal=saved_npz(path,['task','task_base','proposal_x'],folder)
            final=float(np.linalg.norm(task-base)/sigma);raw=None;raw_source=None
            for endpoint in sorted((folder/'task').glob('endpoint_*/RESPONSE.npz')):
                with np.load(endpoint,allow_pickle=False) as z:
                    if 'x' in z and np.array_equal(z['x'],proposal):
                        raw=float(np.linalg.norm(z['data']-base)/sigma);raw_source=str(endpoint.relative_to(root));break
            if abs(final-r['d_T'])>1e-8*max(1.,final):raise ValueError('final task response disagrees with JSON')
            if raw is not None and abs(raw-r['raw_d_T'])>1e-8*max(1.,raw):raise ValueError('raw task response disagrees with JSON')
            task_array_checks.append(dict(family=r['family'],method=r['method'],raw_d_T=raw,final_d_T=final,raw_source=raw_source,final_source=str(path.relative_to(root)),final_replay='MATCH',raw_replay='MATCH' if raw is not None else 'NO_MATCHING_SAVED_PROPOSAL_RESPONSE'))
            r['task_array_replay']='MATCH' if raw is not None else 'FINAL_ONLY'
        except(OSError,ValueError,KeyError) as exc:fail(path,str(exc))
    for folder in sorted((root/'results/B2').glob('task_preserving_*')):
        proposal_path=folder/'PROPOSAL.npz';final_path=folder/'FINAL.npz'
        if not final_path.exists():final_path=folder/'RETAINED_ACCEPTED.npz'
        if not proposal_path.exists() or not final_path.exists():continue
        try:
            p,p0=saved_npz(proposal_path,['task','task_base'],folder);f,f0=saved_npz(final_path,['task','task_base'],folder)
            task_array_checks.append(dict(family='task_preserving',method='PROPOSAL to '+final_path.stem,raw_d_T=float(np.linalg.norm(p-p0)/sigma),final_d_T=float(np.linalg.norm(f-f0)/sigma),raw_source=str(proposal_path.relative_to(root)),final_source=str(final_path.relative_to(root)),raw_replay='SAVED_ARRAY_NORM',final_replay='SAVED_ARRAY_NORM'))
        except(OSError,ValueError,KeyError) as exc:fail(final_path,str(exc))
    valid_changes=[dict(family=c['family'],method=c['method'],raw_d_T=c['raw_d_T'],d_T=c['final_d_T']) for c in task_array_checks if c['raw_d_T'] is not None]
    if valid_changes:
        fig,ax=plt.subplots(figsize=(10,4.5));x=np.arange(len(valid_changes));width=.35
        ax.bar(x-width/2,[r['raw_d_T'] for r in valid_changes],width,label='raw task-response difference')
        ax.bar(x+width/2,[r['d_T'] for r in valid_changes],width,label='retained task-response difference')
        ax.set_xticks(x,[f"{r['family']} / {r.get('method')}" for r in valid_changes],rotation=15,ha='right',fontsize=8);ax.set_ylabel('Saved d_T');ax.legend(fontsize=8)
        finish(fig,'B2_task_response_retention','B2：任务响应差是否在corrector后保留','只对保存了raw_d_T与最终d_T的实际结果配对。响应差降至3以下不代表物理目标被删除；直接joint距离不是同目标d_T。Oracle、固定支持与原最小特征语义未定义的限制持续保留。')
    else:blank('B2_task_response_retention','B2：任务响应差保留','NOT_RUN','暂无同时保存raw/final任务响应的完成记录。')
    outputs={'C_phase_table':phases,'C_phase_array_replay':checks,'C_born_endpoints':born_endpoints,'C_born_joint':born_joint,'C_born_array_replay':born_checks,'B2_branch_status':branch_rows,'B2_task_response_changes':task_changes,'B2_task_array_replay':task_array_checks}
    for name,data in outputs.items():write_csv(root/'tables/saved_results'/f'{name}.csv',data)
    snapshot=dict(status='WORK_IN_PROGRESS',new_physics=0,phase_status=phase_status,phase_rows=len(phases),phase_array_checks=len(checks),frequency_lists_nested=plan.get('frequency_lists_nested',False),nested_aperture=subset,born_joint_status='SAVED_RECORDS' if born_checks else 'RUNNING' if phase_status=='RUNNING' else 'NOT_RUN',wide_reoptimized_status=state(wide),wide_reoptimized_result=wide_record,B2_branches=branch_rows,minimum_feature='UNRESOLVED_ORIGINAL_DEFINITION',conditional_fixed_support_geometry=True,oracle_not_online=True,**outputs)
    atomic(root/'tables/saved_results/C_B2_DISPLAY_SNAPSHOT.json',json.dumps(snapshot,ensure_ascii=False,indent=2,allow_nan=False))
    figure_sources=[dict(figure='figures/saved_results/'+name+'.png',json=source,array_replay=array_source,status=status) for name,source,array_source,status in [
        ('C_narrow_wide_cal_joint','results/C/phase_diagram/PHASE_TABLE.json','response0/response1 data',phase_status),
        ('C_strict_nested_aperture','results/C/phase_diagram/PHASE_TABLE.json','same response arrays, declared TX/RX subset',phase_status),
        ('C_wall_loss_scale','results/C/phase_diagram/PHASE_TABLE.json','same fixed-pair responses at saved conductivity scales',phase_status),
        ('C_Born_vs_full','results/C/phase_diagram/BORN_JOINT_COMPARISON.json','born/JOINT_contrast*.npz',snapshot['born_joint_status']),
        ('C_wide_reoptimized','results/C/wide_reoptimized/RESULTS.json','wide_reoptimized/FINAL.npz',snapshot['wide_reoptimized_status']),
        ('B2_task_response_retention','tables/saved_results/B2_task_response_changes.csv','B2_task_array_replay.csv raw/final saved response sources','SAVED_RECORDS' if valid_changes else 'NOT_RUN')]]
    atomic(root/'tables/saved_results/C_B2_FIGURE_SOURCES.json',json.dumps(figure_sources,ensure_ascii=False,indent=2))
    section='<section><h2>C固定材料对与B2扩展 · WORK_IN_PROGRESS</h2><p>Oracle existence、conditional fixed support；minimum_feature原始标量语义UNRESOLVED。narrow/wide频率非嵌套，central/full aperture才是同一数据的严格子集；不作最终科学判定。</p>'+table([dict(phase_status=phase_status,phase_saved_rows=len(phases),array_replayed_rows=len(checks),Born=snapshot['born_joint_status'],wide_reoptimized=snapshot['wide_reoptimized_status'])])+table(phases,['case','band','common_wall_loss_scale','target_shift','d_mean_all','d_mean_central_nested'])+'<h3>B2实际分支与任务响应</h3>'+table(branch_rows)+table(task_changes,['family','folder','method','status','raw_d_cal','raw_d_T','d_cal','d_T','edited_observation_distance','edit_kind','task_response_change','task_array_replay','identity'])+'<h3>Born已保存端点与joint对照</h3>'+table(born_endpoints,['world','target_position','contrast','d_mean_all','d_mean_central_nested'])+table(born_joint)+'<p><a href="tables/saved_results/C_B2_DISPLAY_SNAPSHOT.json">完整数据/图源快照</a>。无新物理、无Book或任务队列执行、无概率推断。</p></section>'
    return dict(data=snapshot,html=section)
