"""Read saved VAL scores/receipts only; no physics imports or GT-array reads."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import csv,json,time
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1]
OUT=HERE/'results/val_feedback_mechanism'
RUNS={
    'v_physics_seed7_cascade2':'Voxel',
    'j_full_seed7_cascade2':'Jfull',
    'recurrent_refreshed_q_seed7_cascade2':'refreshed_q',
    'recurrent_feedback1_seed7_cascade2':'feedback1',
    'recurrent_feedback_hinge_seed7_cascade2':'feedback_hinge',
}
ORDER=list(RUNS.values())
PAIRS=(('Jfull','Voxel'),('refreshed_q','Jfull'),('feedback1','Jfull'),('feedback_hinge','Jfull'),
       ('feedback1','refreshed_q'),('feedback_hinge','refreshed_q'),('feedback_hinge','feedback1'))
METRICS=('raw_NRMSE','raw_FullVIE_relative_residual','final_NRMSE','final_FullVIE_relative_residual',
         'aux_constitutive_relative_defect','correction_LU','correction_RHS','correction_wall_seconds',
         'post_total_LU','post_total_RHS','post_total_wall_seconds','preparation_wall_seconds',
         'raw_data_gate','raw_joint_success','joint_success','budgeted_joint_success','joint_success_total55',
         'anchor_fallback','correction_not_needed_count')

def numeric(x):return isinstance(x,(int,float,np.number))and not isinstance(x,bool)and np.isfinite(x)
def scalar(x):return float(x)if numeric(x)or isinstance(x,bool)else None
def dumpcsv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def mean(values):
    v=[float(x)for x in values if scalar(x)is not None]
    return float(np.mean(v))if v else None

def extraction():
    scores=[];sources=[]
    for name in ('VAL_ALPHA_SCORES.json','VAL_RECURRENT_SCORES.json'):
        source=HERE/'results'/name;d=json.loads(source.read_text());sources.append(str(source))
        selected=[r for r in d['rows']if r.get('run')in RUNS]
        if any(r.get('split')!='val'for r in selected):raise ValueError('Non-VAL record in supplied VAL comparison')
        scores+=selected
    if len({(r['run'],r['obs_id'])for r in scores})!=len(scores):raise ValueError('Duplicate selected run-observation score')
    rows=[];attemptrows=[]
    for s in scores:
        run=s['run'];obs=s['obs_id'];p=HERE/'evaluation'/run/f'{obs}.json'
        rec=json.loads(p.read_text())if p.exists()else None
        if rec is not None and rec.get('split')!='val':raise ValueError('Receipt is not explicitly VAL: '+str(p))
        candidates=s.get('candidates',[])
        raw=next((c for c in candidates if c.get('label')=='raw_1'),None)
        anchor=next((c for c in candidates if c.get('label')=='anchor'),None)
        final=s.get('top1');aux=(rec or s).get('auxiliary')
        attempts=[]if rec is None else rec.get('attempts',[])
        costs=dict(LU=0.,RHS=0.,wall=0.);cost_status='RECORDED'if rec is not None else'NOT_RECORDED_RECEIPT_MISSING'
        for a in attempts:
            skipped=a.get('status')=='NOT_NEEDED_DATA_GATE'
            start=a.get('start_counts');end=a.get('end_counts')
            if skipped:
                lu=rhs=wall=0.
            elif isinstance(start,dict)and isinstance(end,dict)and numeric(a.get('wall_seconds')):
                lu=end.get('LU',0)-start.get('LU',0)
                rhs=end.get('actual_complex_RHS',0)-start.get('actual_complex_RHS',0)
                wall=a['wall_seconds']
                if min(lu,rhs,wall)<0:raise ValueError('Negative recorded correction cost: '+str(p))
            else:
                lu=rhs=wall=None;cost_status='NOT_RECORDED_INCOMPLETE_ATTEMPT'
            if lu is not None:costs['LU']+=lu;costs['RHS']+=rhs;costs['wall']+=wall
            attemptrows.append(dict(run=run,method=RUNS[run],obs_id=obs,parent_id=s['parent_id'],noise=s['noise'],amplitude=s['amplitude'],
                label=a.get('label'),status=a.get('status'),reason=a.get('reason'),LU=lu,RHS=rhs,wall_seconds=wall,
                start_seconds=a.get('start_seconds'),nfev=a.get('nfev'),cost_scope='all paid corrector attempts, including failures/stops; zero only for explicit NOT_NEEDED_DATA_GATE'))
        prep={}if rec is None else rec.get('preparation_cost',{})
        cc=(rec or s).get('costcounts',{});fixed=prep.get('fixed_operator_cost',{})
        row={k:s.get(k)for k in ('obs_id','parent_id','split','noise','amplitude','family','group')}
        row.update(run=run,method=RUNS[run],seed=7,passes=2,transport_status=s.get('transport_status'),
            status=s.get('status'),selection_status=s.get('selection_status'),accounting_status=s.get('accounting_status'),
            raw_status='RECORDED'if raw else'NOT_RECORDED_RAW_1_MISSING',
            raw_NRMSE=None if raw is None else raw.get('NRMSE'),
            raw_FullVIE_relative_residual=None if raw is None else raw.get('relative_residual'),
            raw_FullVIE_equation_defect=None if raw is None else raw.get('full_VIE_equation_defect'),
            raw_data_gate=None if raw is None else raw.get('physical_data_feasible'),
            raw_joint_success=None if raw is None else raw.get('joint_success'),
            final_NRMSE=None if final is None else final.get('NRMSE'),
            final_FullVIE_relative_residual=None if final is None else final.get('relative_residual'),
            final_FullVIE_equation_defect=None if final is None else final.get('full_VIE_equation_defect'),
            final_label=None if final is None else final.get('label'),final_kind=None if final is None else final.get('kind'),
            final_status='RECORDED'if final else'NOT_RECORDED_NO_TOP1',
            aux_constitutive_status='N/A_VOXEL_NO_AUXILIARY_HEAD'if RUNS[run]=='Voxel'else'RECORDED'if isinstance(aux,dict)and numeric(aux.get('constitutive_relative_defect'))else'NOT_RECORDED',
            aux_constitutive_relative_defect=aux.get('constitutive_relative_defect')if isinstance(aux,dict)else None,
            auxiliary_data_relative_residual=aux.get('auxiliary_data_relative_residual')if isinstance(aux,dict)else None,
            aux_q_field_relative_defect=aux.get('q_aux_field_relative_defect')if isinstance(aux,dict)else None,
            correction_cost_status=cost_status,
            correction_LU=costs['LU']if cost_status=='RECORDED'else None,
            correction_RHS=costs['RHS']if cost_status=='RECORDED'else None,
            correction_wall_seconds=costs['wall']if cost_status=='RECORDED'else None,
            correction_attempt_count=len(attempts),
            correction_not_needed_count=sum(a.get('status')=='NOT_NEEDED_DATA_GATE'for a in attempts),
            correction_failed_or_budget_stop_count=sum(a.get('status')in ('FAILED','BUDGET_STOP')for a in attempts),
            correction_attempt_statuses=json.dumps(Counter(a.get('status')for a in attempts),sort_keys=True),
            post_total_LU=cc.get('LU'),post_total_RHS=cc.get('actual_complex_RHS'),
            post_total_wall_seconds=(rec or s).get('wall_seconds'),preparation_wall_seconds=prep.get('wall_seconds'),
            inference_wall_seconds=prep.get('inference_seconds'),
            fixed_operator_FLOPs=fixed.get('fixed_operator_FLOPs_batch'),
            fixed_operator_FLOPs_status='RECORDED'if 'fixed_operator_FLOPs_batch'in fixed else'NOT_RECORDED',
            fixed_G_matmul_calls=fixed.get('G_matmul_calls'),fixed_G_adjoint_matmul_calls=fixed.get('G_adjoint_matmul_calls'),
            joint_success=s.get('joint_success'),budgeted_joint_success=s.get('budgeted_joint_success'),
            joint_success_total55=s.get('joint_success_total55'),
            anchor_fallback=None if final is None else final.get('label')=='anchor',
            final_missing=final is None,anchor_NRMSE=(s.get('anchor')or{}).get('NRMSE'),
            anchor_FullVIE_relative_residual=None if anchor is None else anchor.get('relative_residual'),
            receipt_path=str(p),receipt_status='RECORDED'if rec is not None else'NOT_RECORDED',
            execution_metadata=json.dumps((rec or{}).get('execution',{}),sort_keys=True),
            residual_scope='normalized residual of actual nonlinear full-VIE physical solve against observed Y; auxiliary residual separate',
            correction_cost_scope='sum of both candidate correction attempts; excludes proposal prep and initial raw verification, retained separately in post_total')
        rows.append(row)
    byrun=defaultdict(set)
    for r in rows:byrun[r['method'],float(r['noise'])].add((r['obs_id'],r['parent_id']))
    for noise in (0.,.001):
        base=byrun['Voxel',noise]
        if len(base)!=78:raise ValueError('Expected 78 VAL mothers per noise')
        if any(byrun[method,noise]!=base for method in ORDER):raise ValueError('Unpaired VAL comparison plan')
    return rows,attemptrows,sources

def summarize(rows):
    groups=defaultdict(list)
    for r in rows:
        for amp in ('ALL',r['amplitude']):groups[r['method'],r['noise'],amp].append(r)
    result=[]
    for (method,noise,amp),rs in groups.items():
        row=dict(method=method,noise=noise,amplitude=amp,planned=len(rs),mothers=len({r['parent_id']for r in rs}),
            raw_recorded=sum(r['raw_status']=='RECORDED'for r in rs),final_recorded=sum(r['final_status']=='RECORDED'for r in rs),
            failed_or_budget_stopped=sum(r['status']!='COMPLETED'for r in rs),
            final_missing_count=sum(r['final_missing']for r in rs),anchor_fallback_count=sum(r['anchor_fallback']is True for r in rs),
            joint_success_count=sum(r['joint_success']is True for r in rs),
            raw_data_gate_count=sum(r['raw_data_gate']is True for r in rs),
            correction_cost_complete=sum(r['correction_cost_status']=='RECORDED'for r in rs),
            correction_status_counts=json.dumps(Counter(r['correction_attempt_statuses']for r in rs),sort_keys=True),
            aux_constitutive_status_counts=json.dumps(Counter(r['aux_constitutive_status']for r in rs),sort_keys=True))
        for m in METRICS:
            vals=[float(r[m])for r in rs if scalar(r[m])is not None]
            row[m+'_n']=len(vals);row[m+'_mean']=float(np.mean(vals))if vals else None
            row[m+'_median']=float(np.median(vals))if vals else None
        result.append(row)
    return result

def paired(rows):
    indexed={(r['method'],r['obs_id']):r for r in rows};pairs=[]
    for treatment,reference in PAIRS:
        for r in rows:
            if r['method']!=reference:continue
            t=indexed.get((treatment,r['obs_id']))
            if t is None:raise ValueError('Missing paired method observation')
            q=dict(treatment=treatment,reference=reference,obs_id=r['obs_id'],parent_id=r['parent_id'],noise=r['noise'],amplitude=r['amplitude'],
                treatment_status=t['status'],reference_status=r['status'],treatment_final_missing=t['final_missing'],reference_final_missing=r['final_missing'])
            for metric in METRICS:
                a=scalar(t[metric]);b=scalar(r[metric]);q[metric+'_delta']=None if a is None or b is None else a-b
            pairs.append(q)
    groups=defaultdict(list)
    for r in pairs:
        for amp in ('ALL',r['amplitude']):groups[r['treatment'],r['reference'],r['noise'],amp].append(r)
    rng=np.random.default_rng(2026110708);summary=[]
    for (treatment,reference,noise,amp),rs in groups.items():
        for metric in METRICS:
            valid=[(r['parent_id'],r[metric+'_delta'])for r in rs if scalar(r[metric+'_delta'])is not None]
            row=dict(treatment=treatment,reference=reference,noise=noise,amplitude=amp,metric=metric,planned_pairs=len(rs),paired_mothers=len(valid),
                missing_pairs=len(rs)-len(valid),difference='treatment minus reference',bootstrap_unit='mother within one noise level; no clean/noisy pooling',bootstrap_repetitions=1000)
            if valid:
                if len({p for p,v in valid})!=len(valid):raise ValueError('More than one observation per mother/noise in bootstrap')
                values=np.array([v for p,v in valid]);n=len(values)
                draws=values[rng.integers(0,n,size=(1000,n))].mean(1)
                row.update(delta_mean=float(values.mean()),delta_median=float(np.median(values)),delta_ci_low=float(np.quantile(draws,.025)),delta_ci_high=float(np.quantile(draws,.975)),
                    negative_fraction=float(np.mean(values<0)),zero_fraction=float(np.mean(values==0)),positive_fraction=float(np.mean(values>0)))
            else:row.update(delta_mean=None,delta_median=None,delta_ci_low=None,delta_ci_high=None,negative_fraction=None,zero_fraction=None,positive_fraction=None)
            summary.append(row)
    return pairs,summary

def figure(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=dict(zip(ORDER,('#777777','#377eb8','#984ea3','#e41a1c','#4daf4a')))
    fig,axes=plt.subplots(2,3,figsize=(17,9),layout='constrained')
    for ni,noise in enumerate((0.,.001)):
        ss={r['method']:r for r in summary if r['noise']==noise and r['amplitude']=='ALL'}
        for mi,method in enumerate(ORDER):
            r=ss[method];c=colors[method]
            axes[ni,0].plot([mi-.1],[r['raw_NRMSE_mean']],'^',mfc='white',mec=c,ms=8)
            axes[ni,0].plot([mi+.1],[r['final_NRMSE_mean']],'o',color=c,ms=7)
            axes[ni,1].plot([mi-.1],[r['raw_FullVIE_relative_residual_mean']],'^',mfc='white',mec=c,ms=8)
            axes[ni,1].plot([mi+.1],[r['final_FullVIE_relative_residual_mean']],'o',color=c,ms=7)
            axes[ni,2].bar(mi,r['correction_LU_mean'],color=c,width=.65)
            axes[ni,2].text(mi,r['correction_LU_mean']+.10,f"{r['correction_RHS_mean']:.1f} RHS\n{r['correction_wall_seconds_mean']:.3f} s",ha='center',va='bottom',fontsize=8)
        for ax in axes[ni]:
            ax.set_xticks(range(5),('Voxel','J-full','refresh-q','quadratic','hinge'),rotation=20)
            ax.grid(axis='y',alpha=.2)
        axes[ni,0].axhline(.15,color='k',linestyle=':',lw=1);axes[ni,0].set_title(('clean'if noise==0 else'noise 0.001')+' / material quality')
        axes[ni,0].set_ylabel('mean amplitude-normalized material NRMSE')
        axes[ni,1].set_yscale('log');axes[ni,1].axhline(.001 if noise==0 else .003,color='k',linestyle=':',lw=1)
        axes[ni,1].set_title(('clean'if noise==0 else'noise 0.001')+' / actual full-VIE residual');axes[ni,1].set_ylabel('mean ‖F(m)−Y‖ / ‖Y‖')
        axes[ni,2].set_title(('clean'if noise==0 else'noise 0.001')+' / paid correction cost');axes[ni,2].set_ylabel('mean correction LU, both candidate attempts')
        axes[ni,2].set_ylim(0,max(r['correction_LU_mean']for r in ss.values())*1.30)
    from matplotlib.lines import Line2D
    axes[0,0].legend(handles=[Line2D([],[],marker='^',mfc='white',mec='k',linestyle='',label='raw_1'),Line2D([],[],marker='o',color='k',linestyle='',label='final top1')],loc='best')
    fig.suptitle('VAL model-selection diagnostic — NOT independent TEST\nSeed7, two passes, identical 78 mothers per noise; saved material metrics, physical residuals and correction receipts',fontsize=15)
    fig.savefig(OUT/'quality_residual_correction_cost.png',dpi=180);plt.close(fig)

def findings(rows,summary,psummary,receipt):
    grouped={(r['method'],r['noise'],r['amplitude']):r for r in summary}
    paired_index={(r['treatment'],r['reference'],r['noise'],r['amplitude'],r['metric']):r for r in psummary}
    def p(t,ref,n,m,amp='ALL'):return paired_index[t,ref,n,amp,m]
    def textdelta(r):return f"{r['delta_mean']:+.5g} [{r['delta_ci_low']:+.5g}, {r['delta_ci_high']:+.5g}]"if r['delta_mean']is not None else'NOT_RECORDED'
    lines=['# VAL feedback mechanism findings','',
        '状态：**只读已完成 VAL 模型选择诊断，不是独立 TEST**。输入仅为指定两个评分 JSON 及其 VAL 观测的 evaluation 收据；没有读取新 GT 文件、TEST 评分/观测、prediction NPZ，也没有导入或调用 physics、Jacobian、SVD、优化器或训练。未遍历 phi_before/phi_after，不做下降恒等式核验。','',
        '比较固定为 Voxel、Jfull、refreshed_q、feedback1（二次反馈）、feedback_hinge，均为 seed7、两层级联。排除 VAL_ALPHA 中全部 posthoc alpha 分支，避免重复计数。780 个 run×observation，78 母体 × clean/noisy，五方法的观测/母体身份完全配对。','',
        '## 八条实际观察','',
        '1. **原始 Jfull 对照下，二次反馈的主要确定性差异在 noisy 校正调用量，而非已确认的最终 GT 改善。** clean 最终 NRMSE 差 +0.00000182，描述性配对区间 [−0.002438,+0.002548]；校正反而多0.385 LU。noisy 最终差 −0.01379，区间 [−0.05099,+0.02422] 包含0，但校正少1.436 LU，区间 [−2.488,−0.410]。不能把所有噪声条件合并成“反馈总会省物理工作”。',
        '2. **相对 refreshed_q，二次反馈同时改善 GT 和校正成本；这个更强差异依赖较弱的对照。** refreshed_q 的最终 NRMSE 为 clean0.10938/noisy0.44985，二次反馈为0.09541/0.22530；平均校正 LU 为6.064→4.769、10.500→8.141。它不是相对原始 Jfull 同等幅度的增益。',
        '3. **Hinge 与二次反馈在 noisy 上明显分离。** Hinge 原始真实残差更低（0.03541 vs0.04213），但最终 NRMSE 更高（0.36082 vs0.22530；配对差 +0.13552，区间 [+0.04425,+0.24012]），校正也多2.051 LU（区间 [+1.231,+2.834]）。较低 raw 真实残差没有自动转成更好 GT 或更便宜的校正。',
        '4. **五方法的 raw_1→最终 top1 平均材料 NRMSE 均增加，而真实数据残差均降低。** 二次反馈 clean0.08845→0.09541、noisy0.17621→0.22530。这里是含候选修正与选择的整个后处理结果，不能单独归因为 TRF；它直接提醒数据拟合与材料恢复必须分开看。',
        '5. **Voxel 仍是这五个 VAL 方法中最终 GT 均值和联合成功数较好的对照。** Voxel 最终 NRMSE clean0.08925/noisy0.18049，joint为62/78、14/78；二次反馈为61/78、7/78。此结果只限已参与模型选择的 VAL。',
        '6. **辅助 constitutive defect 的明显降低不是 GT 证书。** refreshed_q→二次反馈的缺陷均值 clean0.45624→0.04998，noisy0.97307→0.89835；这是辅助头的现成一致性诊断，不是真实 Maxwell-current GT 误差，也不能替代上述 Jfull/voxel 材料对照。',
        '7. **A 分层保留了异质性。** 相对 Jfull，二次反馈 clean A=.1 的最终 NRMSE 差 −0.00636（区间 [−0.01038,−0.00246]），A=.2/.4 的均值差改为正；noisy 各 A 的相应区间都包含0。不能用总体均值覆盖这些方向变化。',
        '8. **没有把失败或 fallback 删除后配对。** 780条全为已完成且物理可行 top1，但不等于 GT 联合成功；raw_1 在780条中全部未直接过数据门。1560次校正包含1501次 PHYSICAL_DATA_GATE 停止、59次优化器结束，0次 FAILED/BUDGET_STOP。noisy anchor fallback 为 refreshed_q3、二次反馈1、hinge3；均原样保留。','',
        '## 1. 实际结果与校正费用','',
        '|noise|method|raw NRMSE|final NRMSE|raw FullVIE residual|final FullVIE residual|correction LU|RHS|wall/s|joint/78|raw gate/78|',
        '|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for noise in (0.,.001):
        for method in ORDER:
            s=grouped[method,noise,'ALL']
            lines.append(f"|{noise}|{method}|{s['raw_NRMSE_mean']:.5g}|{s['final_NRMSE_mean']:.5g}|{s['raw_FullVIE_relative_residual_mean']:.5g}|{s['final_FullVIE_relative_residual_mean']:.5g}|{s['correction_LU_mean']:.4g}|{s['correction_RHS_mean']:.5g}|{s['correction_wall_seconds_mean']:.4g}|{s['joint_success_count']}|{s['raw_data_gate_count']}|")
    lines+=['','表中 FullVIE residual 是已有候选经过原始非线性 Maxwell 全求解后的 ‖F(m)−Y‖/‖Y‖，不是辅助电流的残差。raw_1 为实际评估的完整步长材料候选；final 为已冻结 top1。另存方程离散求解缺陷 full_VIE_equation_defect，不把它与数据拟合残差混用。均值点也不代替逐样本联合成功判据。','',
        '校正成本是两个 alpha 候选的全部已执行 corrector attempts 的 start/end 计数差及 wall 之和，含失败和预算停止；仅明确 NOT_NEEDED_DATA_GATE 的调用计为零。它排除初始候选验证、准备与推理；这些完整开销另列 post_total_LU/RHS/wall 和 preparation_wall。没有用 top1 被选中时的累计计数冒充完整后处理费用。','',
        '## 2. 反馈相对原始 Jfull：费用改善是否伴随 GT 改善？','',
        '下表为 treatment−reference 的母体配对均值差及描述性 95% bootstrap 区间（1000 次、每个 noise 内按母体重抽；没有跨 noise 合并，也没有多重比较校正）。NRMSE、LU、RHS、wall 的负差各自表示降低；NRMSE 与费用分别判断，不能互相替代。','',
        '|noise|treatment|Δraw NRMSE|Δfinal NRMSE|Δcorrection LU|Δcorrection RHS|Δcorrection wall/s|',
        '|---:|---|---|---|---|---|---|']
    for noise in (0.,.001):
        for method in ('feedback1','feedback_hinge'):
            values=[textdelta(p(method,'Jfull',noise,m))for m in ('raw_NRMSE','final_NRMSE','correction_LU','correction_RHS','correction_wall_seconds')]
            lines.append('|'+f'{noise}|{method}|'+('|'.join(values))+'|')
    lines+=['','下表单独列出相对 refreshed_q 的比较；它的弱化不能被当作原始 Jfull 基线。','',
        '|noise|treatment vs refreshed_q|Δraw NRMSE|Δfinal NRMSE|Δcorrection LU|Δcorrection RHS|Δcorrection wall/s|',
        '|---:|---|---|---|---|---|---|']
    for noise in (0.,.001):
        for method in ('feedback1','feedback_hinge'):
            values=[textdelta(p(method,'refreshed_q',noise,m))for m in ('raw_NRMSE','final_NRMSE','correction_LU','correction_RHS','correction_wall_seconds')]
            lines.append('|'+f'{noise}|{method}|'+('|'.join(values))+'|')
    lines+=['','## 3. Hinge 与二次反馈不是同一个经验结果','',
        '|noise|hinge−quadratic: final NRMSE|correction LU|raw physical residual|aux constitutive defect|',
        '|---:|---|---|---|---|']
    for noise in (0.,.001):
        values=[textdelta(p('feedback_hinge','feedback1',noise,m))for m in ('final_NRMSE','correction_LU','raw_FullVIE_relative_residual','aux_constitutive_relative_defect')]
        lines.append('|'+f'{noise}|'+('|'.join(values))+'|')
    lines+=['','## 4. A 分层保留方向变化','',
        '|noise|A|method vs refreshed_q|Δfinal NRMSE|Δcorrection LU|Δjoint success fraction|',
        '|---:|---:|---|---|---|---|']
    for noise in (0.,.001):
        for amp in sorted({r['amplitude']for r in rows}):
            for method in ('feedback1','feedback_hinge'):
                vals=[textdelta(p(method,'refreshed_q',noise,m,amp))for m in ('final_NRMSE','correction_LU','joint_success')]
                lines.append('|'+f'{noise}|{amp}|{method}|'+('|'.join(vals))+'|')
    lines+=['','## 5. 辅助诊断、失败和 fallback','',
        'aux_constitutive_relative_defect 仅来自现成收据中的 ‖j_aux−m·(E0+G j_aux)‖/‖j_aux‖，按其原有定义保存。Voxel 无辅助头，标 N/A；任何缺项标 NOT_RECORDED。它不是实际 Maxwell 电流误差，也不等于 GT 材料恢复。没有把评分中的 auxiliary-current truth errors 当作真实物理电流误差使用。','',
        '|noise|method|aux constitutive mean|post total LU|post wall/s|prep wall/s|fallback anchor|missing final|',
        '|---:|---|---:|---:|---:|---:|---:|---:|']
    for noise in (0.,.001):
        for method in ORDER:
            s=grouped[method,noise,'ALL'];aux='N/A'if method=='Voxel'else'NOT_RECORDED'if s['aux_constitutive_relative_defect_mean']is None else f"{s['aux_constitutive_relative_defect_mean']:.5g}"
            lines.append(f"|{noise}|{method}|{aux}|{s['post_total_LU_mean']:.4g}|{s['post_total_wall_seconds_mean']:.4g}|{s['preparation_wall_seconds_mean']:.4g}|{s['anchor_fallback_count']}|{s['final_missing_count']}|")
    lines+=['',f"全部状态计数：{receipt['score_status_counts']}；选择状态：{receipt['selection_status_counts']}。失败、预算停止、缺失 top1 和 anchor fallback 都在逐样本表保留，没有删掉后再配对。",'',
        '## 6. 成本解释的边界','',
        'LU/RHS 为物理调用量；wall 同时受执行时序、并发、warmup 与硬件记录完整度影响。收据中的执行元数据逐行保存，不能从跨批次 wall 差推出纯算法加速。反馈步骤的固定矩阵乘法不是免费：有记录时保留 fixed_operator_FLOPs、G/G* 调用和 preparation_wall；缺失的老分支 FLOP 记 NOT_RECORDED，不补零。完整 post_total 成本用于核对“节省校正”是否也转化为实际后处理节省。共用 GN 是单独沉没成本，不在此重复计入；已有 joint_success_total55 仍原样保留。','',
        '例如，相对 Jfull，二次反馈 clean wall 更短，但实际校正 LU/RHS 更高；这本身就阻止把 wall 差直接解释为减少物理调用。refreshed_q 的已记录固定算子前向 FLOPs 为63700992，两个反馈分支均为207028224；这不含相同网络卷积65286144 FLOPs及明确未估计的其它运算，不能当作总训练/总推理 FLOPs。','',
        '本分析仅描述 VAL 的质量—真实残差—校正费用关系。局部下降理论、辅助缺陷降低、较少 LU 和较好 GT 恢复是不同证据；本表不把优化几何下降当作科学成功。父线程负责最终解释与模型选择。','',
        '输出：per_case.csv、correction_attempts.csv、grouped_summary.csv、paired_per_parent.csv、paired_summary.csv、RECEIPT.json，以及 quality_residual_correction_cost.png。图标题明确 VAL model selection、非独立 TEST。']
    (HERE/'delegated/energy/VAL_FEEDBACK_FINDINGS.md').write_text('\n'.join(lines)+'\n')

def main():
    started=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True)
    rows,attempts,sources=extraction();summary=summarize(rows);pairs,psummary=paired(rows)
    dumpcsv(OUT/'per_case.csv',rows);dumpcsv(OUT/'correction_attempts.csv',attempts)
    dumpcsv(OUT/'grouped_summary.csv',summary);dumpcsv(OUT/'paired_per_parent.csv',pairs);dumpcsv(OUT/'paired_summary.csv',psummary)
    receipt=dict(status='COMPLETE_READ_ONLY_VAL_ANALYSIS',observations_per_method=156,mothers_per_noise=78,methods=ORDER,rows=len(rows),
        source_scores=sources,excluded='all posthoc-alpha runs; no duplicated alpha cases',
        score_status_counts=dict(Counter(r['status']for r in rows)),selection_status_counts=dict(Counter(r['selection_status']for r in rows)),
        auxiliary_status_counts=dict(Counter(r['aux_constitutive_status']for r in rows)),
        correction_cost_status_counts=dict(Counter(r['correction_cost_status']for r in rows)),
        physics_calls=0,GT_files_read=0,TEST_files_read=0,prediction_NPZ_reads=0,phi_identity_checks=0,
        statistical_unit='paired mother, clean/noisy and amplitude strata separately',BLAS_threads=1)
    figure(summary);findings(rows,summary,psummary,receipt)
    receipt['wall_seconds']=time.perf_counter()-started;(OUT/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2))
    for noise in (0.,.001):
        for method in ('feedback1','feedback_hinge'):
            selected=[r for r in psummary if r['treatment']==method and r['reference']=='refreshed_q'and r['noise']==noise and r['amplitude']=='ALL'and r['metric']in ('raw_NRMSE','final_NRMSE','correction_LU','correction_RHS','correction_wall_seconds')]
            print(json.dumps(dict(noise=noise,treatment=method,reference='refreshed_q',differences={r['metric']:[r['delta_mean'],r['delta_ci_low'],r['delta_ci_high']]for r in selected})))

if __name__=='__main__':main()
