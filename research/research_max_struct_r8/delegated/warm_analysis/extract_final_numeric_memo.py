"""Extract existing score/table rows only. No solver, score runner, GT, plots or CI recomputation."""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,math,time
STARTCPU=time.process_time();STARTWALL=time.perf_counter()
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'delegated/warm_analysis'
COHORTS=('IN_DISTRIBUTION','HELDOUT_FAMILY','CROSS_GRID');SHORT=dict(zip(COHORTS,('ID24','OOD12','cross6')))
ARMS=('physical_lbfgs','cie_dual','global_dual','a_plus','a_plus_sourcewise')
CORE=[a+' / '+i for a in ARMS for i in ('uniform','nn')]
ANCHOR=[a+' / nn_affine_anchor' for a in ('cie_dual','a_plus','a_plus_sourcewise')]
def load(p):return json.loads((ROOT/p).read_text())
def csvrows(p):
    out=[]
    with (ROOT/p).open(newline='') as f:
        reader=csv.DictReader(f)
        for r in reader:
            r=dict(r);r['_evidence']=dict(artifact=p,line=reader.line_num);out.append(r)
    return out

def num(x):return float(x) if x not in (None,'') else None
def slim(r):return dict(r)
def fmt(x,n=4):return 'NA' if x is None else f'{float(x):.{n}f}'
def ci(r):return f"{num(r['mean_paired_difference']):+.4f} [{num(r['percentile95_low']):+.4f},{num(r['percentile95_high']):+.4f}]"
def evidence(r):return f"{r['_evidence']['artifact']}:{r['_evidence']['line']}"
original=load('RESULTS/warm_scoring/SCORES.json')['rows'];anchor=load('delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json')['rows']
assert all(r.get('material_legal') for r in original+anchor), 'Failure-category extraction assumes inherited legal-box checks passed'
primary=[r for r in original if r['primary']];extra=[r for r in original if not r['primary']]
assert len(primary)==630 and len(extra)==120 and len(anchor)==126
pco=csvrows('delegated/warm_analysis/PRIMARY_630_cohort_summary.csv');aco=csvrows('delegated/warm_anchor_analysis/CONFIRMATION_cohort_material.csv')
core=[r for r in pco if r['method'] in CORE]+[r for r in aco if r['method'] in ANCHOR]
verification=[]
for r in core:
    scores=[s for s in primary+anchor if (s['cohort'],s['method'])==(r['cohort'],r['method'])]
    n=len(scores);mean=sum(s['material_nrmse'] for s in scores)/n
    assert n==int(r['objects']) and math.isclose(mean,num(r['material_nrmse_mean']),abs_tol=1e-12)
    assert sum(s['material_image_success'] for s in scores)==int(r['material_successes'])
verification.append('39 core method/cohort summary rows independently matched existing SCORES count, mean NRMSE and inherited success count')
pcost=csvrows('delegated/warm_analysis/PRIMARY_630_absolute_cost_by_hardware.csv');acost=csvrows('delegated/warm_anchor_analysis/CONFIRMATION_hardware_costs.csv')
costs=[]
for row in pcost:
    if row['method'] not in CORE or not row['method'].endswith('/ nn'):continue
    cohort,host,method=row['cohort'],row['hardware'],row['method'];ref=method.replace('/ nn','/ uniform')
    u=next(r for r in pcost if (r['cohort'],r['hardware'],r['method'])==(cohort,host,ref))
    nn=[r for r in primary if (r['cohort'],r['hardware'],r['method'])==(cohort,host,method)]
    uniform=[r for r in primary if (r['cohort'],r['hardware'],r['method'])==(cohort,host,ref)]
    assert {r['case_id'] for r in nn}=={r['case_id'] for r in uniform}
    out=dict(cohort=cohort,hardware=host,method=method,reference=ref,paired_objects=len(nn),evidence=[row['_evidence'],u['_evidence']])
    for key in ('solver_cpu_seconds','warm_deployment_cpu_seconds'):
        a=num(row[key+'_mean']);b=num(u[key+'_mean'])
        assert math.isclose(a,sum(s[key] for s in nn)/len(nn),abs_tol=1e-10)
        out[key]=dict(nn_mean=a,uniform_mean=b,ratio_of_means_nn_over_uniform=a/b,mean_paired_extra=a-b)
    costs.append(out)
verification.append('All25 core NN/uniform cohort/hardware CPU mean pairs use identical case IDs; ratios are ratios of means, never cross-hardware averages')
finalci=[r for r in csvrows('delegated/warm_analysis/PRIMARY_NN_minus_uniform_object_paired.csv') if r['method'] in CORE]
budget=[r for r in csvrows('delegated/warm_analysis/PRIMARY_object_matched_budget_paired_bootstrap.csv') if r['axis'] in ('CPU_seconds','LU','complex_RHS') and r['metric'] in ('material_nrmse','original_Phi') and r['method'] in CORE and r['reference']==r['method'].replace('/ nn','/ uniform')]
anchorbudget=[r for r in csvrows('delegated/warm_anchor_analysis/CONFIRMATION_budget_paired_effects.csv') if r['axis'] in ('CPU_seconds','LU','complex_RHS') and r['metric'] in ('material_nrmse','original_Phi')]
training=load('delegated/warm_analysis/TRAINING_COST_ACCOUNTING.json');amort=csvrows('delegated/warm_analysis/PRIMARY_training_amortization_mean_paired_CPU.csv')
assert math.isclose(training['original_physics_network']['cpu_seconds'],402.558889,abs_tol=1e-9)
portco=csvrows('delegated/warm_analysis/PAID_PORTFOLIOS_cohort_summary.csv');portcost=csvrows('delegated/warm_analysis/PAID_PORTFOLIOS_absolute_cost_by_hardware.csv')
portfolios=[]
for r in portco:
    if not r['method'].endswith(' / nn / retained_physics'):continue
    c=r['cohort'];u=next(x for x in portco if (x['cohort'],x['method'])==(c,r['method'].replace('/ nn /','/ uniform /')))
    d=next(x for x in portco if (x['cohort'],x['method'])==(c,'physical_lbfgs / nn'))
    portfolios.append(dict(cohort=c,method=r['method'],objects=int(r['objects']),mean_nrmse=num(r['material_nrmse_mean']),successes=int(r['material_successes']),
        minus_equal_candidate_uniform_nrmse=num(r['material_nrmse_mean'])-num(u['material_nrmse_mean']),
        minus_single_NNdirect_nrmse=num(r['material_nrmse_mean'])-num(d['material_nrmse_mean']),
        single_NNdirect_successes=int(d['material_successes']),uniform_portfolio_successes=int(u['material_successes']),evidence=[r['_evidence'],u['_evidence'],d['_evidence']]))
portpaid=[]
for r in portcost:
    if not r['method'].endswith(' / nn / retained_physics'):continue
    c,h=r['cohort'],r['hardware'];u=next(x for x in portcost if (x['cohort'],x['hardware'],x['method'])==(c,h,r['method'].replace('/ nn /','/ uniform /')))
    d=next(x for x in portcost if (x['cohort'],x['hardware'],x['method'])==(c,h,'physical_lbfgs / nn'))
    a,b,dd=[num(x['warm_deployment_cpu_seconds_mean']) for x in (r,u,d)]
    portpaid.append(dict(cohort=c,hardware=h,method=r['method'],objects=int(r['objects']),warm_CPU_mean=a,
        uniform_two_candidate_warm_CPU_mean=b,NNdirect_warm_CPU_mean=dd,extra_vs_NNdirect=a-dd,ratio_vs_NNdirect=a/dd,
        ratio_vs_uniform_two_candidate=a/b,evidence=[r['_evidence'],u['_evidence'],d['_evidence']]))
portci=csvrows('delegated/warm_analysis/PORTFOLIO_paired_bootstrap_vs_equal_two_candidate_uniform.csv')
fixed=csvrows('delegated/warm_analysis/FIXED8_120_EXTRA_WITH_MATCHED_CONTROLS_cohort_summary.csv')
fixedmethods=[a+' / '+i for a in ('cie_single','cie_dual','a_plus') for i in ('uniform','uniform1','born','dpcs','random','nn','nn_joint','nn_born_features')]
fixed=[r for r in fixed if r['method'] in fixedmethods]
survival=csvrows('delegated/warm_mechanism_analysis/COHORT_SURVIVAL_COUNTS.csv')
survival=[r for r in survival if r['arm'] in ('cie_dual','global_dual','a_plus','a_plus_sourcewise') and r['initialization'] in ('nn','nn_joint','nn_affine_anchor')]
conditions={r['case_id']:r for r in load('data_public/warm_test/test_manifest.json')['records']}
failure_cases=[];failure_counts=[]
for c in COHORTS:
    rr=[r for r in primary if r['cohort']==c and r['method']=='physical_lbfgs / nn'];categories=Counter()
    for r in rr:
        if r['material_image_success']:continue
        category=('both_NRMSE_and_IoU' if not r['material_nrmse_pass'] and not r['iou_pass'] else 'NRMSE_only' if not r['material_nrmse_pass'] else 'IoU_only')
        categories[category]+=1
        failure_cases.append(dict(case_id=r['case_id'],cohort=c,hardware=r['hardware'],material_nrmse=r['material_nrmse'],otsu_iou=r['otsu_iou'],failure_components=category,
            public_condition=conditions[r['case_id']]['condition'],declared_amplitude=r['amplitude'],corrector_stop=r['corrector_stop'],projected_gradient_inf=r.get('projected_gradient_inf'),
            evidence=dict(artifact='RESULTS/warm_scoring/SCORES.json',row_index=original.index(r),task_id=r['task_id']),condition_evidence=dict(artifact='data_public/warm_test/test_manifest.json',case_id=r['case_id'])))
    failure_counts.append(dict(cohort=c,objects=len(rr),material_success=sum(r['material_image_success'] for r in rr),failure_categories=dict(categories),case_ids=[x['case_id'] for x in failure_cases if x['cohort']==c]))
# Public failure conditions for all core methods, never infer hidden family names from geometry.
all42=[]
for case in conditions:
    direct=next(r for r in primary if r['case_id']==case and r['method']=='physical_lbfgs / nn')
    all42.append(dict(case_id=case,cohort=direct['cohort'],condition=conditions[case]['condition'],declared_amplitude=direct['amplitude'],
        methods=[dict(method=r['method'],status=r['status'],nrmse=r['material_nrmse'],iou=r['otsu_iou'],material_success=r['material_image_success']) for r in primary+anchor if r['case_id']==case and r['method'] in CORE+ANCHOR]))
# Label discordance uses only the inherited existing success field; no new threshold.
label_discordance=[]
for cohort in COHORTS:
    direct={r['case_id']:r for r in primary if r['cohort']==cohort and r['method']=='physical_lbfgs / nn'}
    for method in [a+' / nn' for a in ('cie_dual','global_dual','a_plus','a_plus_sourcewise')]+['a_plus / nn_joint']+ANCHOR:
        rr=[r for r in primary+anchor if r['cohort']==cohort and r['method']==method]
        losses=[r['case_id'] for r in rr if direct[r['case_id']]['material_image_success'] and not r['material_image_success']]
        gains=[r['case_id'] for r in rr if not direct[r['case_id']]['material_image_success'] and r['material_image_success']]
        label_discordance.append(dict(cohort=cohort,method=method,objects=len(rr),direct_success_aux_failure_case_ids=losses,
            aux_success_direct_failure_case_ids=gains,direct_success_aux_failure_count=len(losses),aux_success_direct_failure_count=len(gains),
            evidence=[dict(artifact='RESULTS/warm_scoring/SCORES.json' if r in primary else 'delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json',task_id=r['task_id']) for r in rr],
            selection='Inherited saved material_image_success labels only; no newly chosen error tolerance'))
manifest=load('configs/warm_all_manifest.json')['tasks']+load('configs/warm_anchor_confirmation_manifest.json')['tasks'];taskmap={t['task_id']:t for t in manifest}
case27=[]
for row in primary+anchor:
    if row['case_id']!='r8_test_0027' or row['method'] not in CORE+ANCHOR+['a_plus / nn_joint'] or row['initialization']=='uniform':continue
    d=Path(taskmap[row['task_id']]['output_dir']);d=d if d.is_absolute() else ROOT.parents[1]/d
    rec=json.loads((d/'RECEIPT.json').read_text());selected=rec.get('auxiliary_selected_stage')
    pre=next((c.get('objective') for c in rec.get('auxiliary_candidate_scores',[]) if c.get('stage')==selected),None)
    corr=rec.get('corrector') or {}
    case27.append(dict(case_id=row['case_id'],cohort=row['cohort'],method=row['method'],selected_stage=selected,
        pre_corrector_Phi=pre,final_Phi=row['original_Phi'],material_nrmse=row['material_nrmse'],otsu_iou=row['otsu_iou'],
        inherited_material_success=row['material_image_success'],corrector_stop=corr.get('stop'),
        corrector_accepted_iterations=corr.get('iterations'),corrector_evaluations=corr.get('evaluations'),
        corrector_receipt_success=corr.get('success'),solver_CPU_seconds=row['solver_cpu_seconds'],
        evidence=[dict(artifact=str((d/'RECEIPT.json').relative_to(ROOT)),keys=['auxiliary_selected_stage','auxiliary_candidate_scores','corrector']),
                  dict(artifact='RESULTS/warm_scoring/SCORES.json' if row in primary else 'delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json',task_id=row['task_id'])]))
limitations=[
 'Final endpoints have unequal actual costs; NNdirect/uniform are the same algorithm with6000-step ceiling, auxiliary flows have3000 auxiliary plus3000 corrector ceilings.',
 'Matched-budget endpoints use min(finalA,finalB) independently per object and cost axis, sparse paid Phi incumbents, unused allowance and checkpoint gaps; they are at-most allowance, not exact equal expenditure or a GT iterate oracle.',
 'Hardware partitions are different fixed cases: ID Mac10/XINAN14, OOD Mac6/XINAN6, cross-grid Mac6. Absolute CPU is never pooled across hardware or interpreted as monetary cost.',
 'Global A+ anchor vs nn_joint is the available truncated physical-current control. Dual/sourcewise have no analogous nn_joint control; anchor changes their initialization and affine family together. No42-object uniform-anchor arm exists.',
 'Fixed8 has ID4/OOD2/cross2 and120 additional rows; its extra initializers/Born-only network cannot be assigned a42-object denominator. Training costs of Born/domain repair are not the original physics-network numerator.',
 'All876 solver predictions were sealed; a solver can terminate successfully while failing material metrics. NRMSE/IoU and noise discrepancy are separate, and smooth support topology is diagnostic.',
 'Cross-grid36→24 remains within the same discrete scalar-TM physical family; it is not independent continuum-PDE validation.',
 'Strict signed pre/post flags include possible floating-point micro-differences; not every counted case is a substantive material/basin failure. Direct ceiling6000 vs auxiliary corrector3000 also prevents interpreting these counts as one common fixed correction map.',
 'Stage−1 means the initial candidate entered correction, not that the final image was unchanged. Lower prePhi followed by worse endpoints is descriptive; it does not identify a basin or establish causality.',
 'Portfolios charge both complete solvers and share measured NN prefixes once. They were selected by public Phi, not material error; no partial-branch budget schedule was frozen, so no invented equal-budget portfolio trajectory is claimed.',
 'Existing bootstrap intervals are descriptive object-paired percentile95% intervals; no significance claim, threshold search or new resampling was performed here.'
]
result=dict(scope=dict(original_primary=630,primary_objects=42,primary_methods=15,extra_tasks=120,extra_objects=8,extra_cohorts=dict(Counter(conditions[c]['cohort'] for c in {r['case_id'] for r in extra})),anchor_tasks=126,solver_statuses=dict(Counter(r['status'] for r in original+anchor))),
    verification=verification,core_material_summary=core,same_algorithm_NN_vs_uniform_final_CIs=finalci,same_hardware_CPU_means_ratios=costs,
    selected_existing_matched_budget_CIs=budget,anchor_existing_matched_budget_CIs=anchorbudget,training_original_physics=training['original_physics_network'],
    all_existing_training_amortization=amort,portfolio_material_differences=portfolios,portfolio_hardware_costs=portpaid,portfolio_existing_CIs=portci,
    fixed8_initialization_and_Born_feature_rows=fixed,stage_survival_and_pre_post_counts=survival,NNdirect_failure_category_counts=failure_counts,
    NNdirect_failed_case_conditions=failure_cases,all42_core_outcomes_and_public_conditions=all42,
    inherited_material_label_discordance=label_discordance,predeclared_case0027_stop_and_endpoint=case27,limitations=limitations,
    extraction_cost=dict(cpu_seconds=time.process_time()-STARTCPU,wall_seconds=time.perf_counter()-STARTWALL,LU=0,complex_RHS=0,new_physics=0,new_training=0,new_GT_reads=0,bootstrap_recomputed=False,plots_overwritten=False))
(OUT/'FINAL_NUMERIC_MEMO.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

md=['# 冻结结果数字备忘录','', '只摘录既有评分/汇总/置信区间；未重跑评分、实验、训练或bootstrap，未读新GT，也未覆盖图。精确原表行号、所有硬件层及逐对象公开条件见同名JSON。以下均不作科学胜负判断。','',
 '原主比较630=42×15；额外120仅固定8对象（ID4/OOD2/cross2）；anchor126=42×3独立扩展。876个回执均为SEALED_PUBLIC_PREDICTION，不能据此称材料恢复成功。', '',
 '## 最终材料：均值NRMSE；继承复合成功数/分母','', '|方法|ID24|OOD12|cross6|','|---|---:|---:|---:|']
for method in CORE+ANCHOR:
    rr=[next(r for r in core if (r['cohort'],r['method'])==(c,method)) for c in COHORTS]
    md.append('|'+method+'|'+'|'.join(f"{num(r['material_nrmse_mean']):.5f}; {r['material_successes']}/{r['objects']}" for r in rr)+'|')
md+=['','NRMSE/IoU定义及成功门槛完全继承；表不合并三个总体。源：PRIMARY_630_cohort_summary.csv和CONFIRMATION_cohort_material.csv；JSON逐单元保留原CSV物理行号。','',
 '## NNdirect vs uniform：相同原物理算法、同硬件实际CPU','', '|cohort/硬件(n)|solver CPU均值 uniform→NN|warm CPU均值 uniform→NN|warm NN/uniform|','|---|---:|---:|---:|']
for r in costs:
    if r['method']!='physical_lbfgs / nn':continue
    a,b=r['solver_cpu_seconds'],r['warm_deployment_cpu_seconds']
    md.append(f"|{SHORT[r['cohort']]}/{r['hardware']}({r['paired_objects']})|{a['uniform_mean']:.3f}→{a['nn_mean']:.3f}|{b['uniform_mean']:.3f}→{b['nn_mean']:.3f}|{b['ratio_of_means_nn_over_uniform']:.3f}|")
md+=['','CPU比为配对同对象的**均值之比**。原feature/inference字段为Mac实际计时；deployment_prefix为独立同机计时；actual_workflow可包含Mac预处理+XINAN求解。所有25个核心NN/uniform硬件层均值/比值见JSON，未混硬件平均。', '',
 '## NNdirect相同at-most预算端点：NRMSE差 NN−uniform [95%区间]','', '|cohort/硬件(n)|CPU allowance|LU allowance|RHS allowance|','|---|---|---|---|']
for cost in costs:
    if cost['method']!='physical_lbfgs / nn':continue
    cc,hh=cost['cohort'],cost['hardware'];rr=[next(r for r in budget if (r['cohort'],r['hardware'],r['method'],r['axis'],r['metric'])==(cc,hh,'physical_lbfgs / nn',axis,'material_nrmse')) for axis in ('CPU_seconds','LU','complex_RHS')]
    md.append(f"|{SHORT[cc]}/{hh}({cost['paired_objects']})|"+'|'.join(ci(r) for r in rr)+'|')
md+=['','预算是每对象min(两方法最终实付成本)，按该预算内原Phi incumbent取已保存状态；停止方法不续算。不是精确等支出；unused和checkpoint间距来自现有端点表。所有原NN/均匀及anchor三方法对原NN/NNdirect/global nn_joint的CPU/LU/RHS区间已逐行摘录JSON。例：ID/XINAN global-anchor−globalNN在CPU预算NRMSE差+0.11044 [+.00077,+.23583]；同组anchor−nn_joint为−.00828 [−.07402,+.06107]。不把区间包装成显著性或普遍机制。', '',
 '## 原physics网络训练摊销','', '固定分子402.558889 Mac CPU秒=生成1536训练11.396809 + 192验证3.472078 + 全120epoch训练387.690002；生成总体已含feature/setup，未重复加。Born/domain费用不计入此部署分子。NNdirect五个cohort/硬件层平均CPU节省全≤0，均无有限CPU摊销优势。其他支路的全部25条均值摊销记录在JSON；Mac为CPU记账，XINAN仅给(cMac/cXINAN)乘系数，未指定价比不能称经济盈亏平衡。','',
 '## 两完整支路portfolio：最终均值NRMSE；成功数','', '|NN physical+辅助支路|ID24|OOD12|cross6|','|---|---:|---:|---:|']
for arm in ('cie_single','cie_dual','global_dual','a_plus','a_plus_sourcewise'):
    rr=[next(r for r in portfolios if (r['cohort'],r['method'])==(c,'portfolio / '+arm+' / nn / retained_physics')) for c in COHORTS]
    md.append('|'+arm+'|'+'|'.join(f"{r['mean_nrmse']:.5f}; {r['successes']}/{r['objects']}" for r in rr)+'|')
md+=['','例如global A+ portfolio相对单NNdirect均值差ID −.01893、OOD +.18253、cross +.14450；相对等候选数uniform portfolio均值差ID −.17071、OOD −.29260、cross −.09927。不是GT最优候选，而是已封存原Phi选择。全部5种组合的材料差和旧CI在JSON。','', '|global A+ portfolio硬件层|NNdirect warm CPU|NN两支路warm CPU|uniform两支路warm CPU|NN两支路/NNdirect|','|---|---:|---:|---:|---:|']
for r in portpaid:
    if r['method']=='portfolio / a_plus / nn / retained_physics':md.append(f"|{SHORT[r['cohort']]}/{r['hardware']}|{r['NNdirect_warm_CPU_mean']:.3f}|{r['warm_CPU_mean']:.3f}|{r['uniform_two_candidate_warm_CPU_mean']:.3f}|{r['ratio_vs_NNdirect']:.3f}|")
md+=['','两条完整solver成本均支付，共享NN prefix只计一次。没有冻结部分支路执行调度，因此没有等预算portfolio轨迹。', '',
 '## 固定8对象初始化/Born-only：仅global A+示例，全部3支路×初始化在JSON','', '|初值/网络|ID4 均值;成功|OOD2 均值;成功|cross2 均值;成功|','|---|---:|---:|---:|']
for init in ('uniform','uniform1','born','dpcs','random','nn','nn_joint','nn_born_features'):
    rr=[next(r for r in fixed if (r['cohort'],r['method'])==(c,'a_plus / '+init)) for c in COHORTS]
    md.append('|'+init+'|'+'|'.join(f"{num(r['material_nrmse_mean']):.5f}; {r['material_successes']}/{r['objects']}" for r in rr)+'|')
md+=['','Born-only网络是另一个训练控制，不能据ID4的结果宣称42对象优势。解析Born和Born-only NN是不同初值。', '',
 '## 初值进入corrector及“更低prePhi、但finalPhi与材料皆劣于NNdirect”','', '|方法|ID stage−1 /双劣|OOD stage−1 /双劣|cross stage−1 /双劣|','|---|---:|---:|---:|']
for arm,init in [('cie_dual','nn'),('global_dual','nn'),('a_plus','nn'),('a_plus_sourcewise','nn'),('a_plus','nn_joint'),('cie_dual','nn_affine_anchor'),('a_plus','nn_affine_anchor'),('a_plus_sourcewise','nn_affine_anchor')]:
    rr=[next(r for r in survival if (r['cohort'],r['arm'],r['initialization'])==(c,arm,init)) for c in COHORTS]
    md.append('|'+arm+'/'+init+'|'+'|'.join(f"{r['stage_minus1_selected_count']}/{r['all_objects']} ; {r['lower_pre_Phi_but_both_final_Phi_and_material_worse_than_NNdirect_count']}/{r['all_objects']}" for r in rr)+'|')
md+=['','stage−1只指初始候选进入校正，非最终图像保持不变。“双劣”为严格符号比较，会包含浮点微差，不能称每一计数都是实质basin失败。NNdirect上限6000 vs辅助corrector3000也使计数不是统一固定C映射的证明；不另设差值阈值。','',
 '### 更具体的继承成功标签分歧：direct成功/辅助失败；反向增益','', '|辅助方法|ID24 丢失/新增|OOD12 丢失/新增|cross6 丢失/新增|','|---|---:|---:|---:|']
for method in [a+' / nn' for a in ('cie_dual','global_dual','a_plus','a_plus_sourcewise')]+['a_plus / nn_joint']+ANCHOR:
    rr=[next(r for r in label_discordance if (r['cohort'],r['method'])==(c,method)) for c in COHORTS]
    md.append('|'+method+'|'+'|'.join(str(r['direct_success_aux_failure_count'])+'/'+str(r['aux_success_direct_failure_count']) for r in rr)+'|')
md+=['','每个丢失/新增对象ID集合已列入JSON；完全继承原材料成功标签，无新门槛。','',
 '### 预定对象0027：停止原因与实际accepted次数','', '|方法|prePhi→finalPhi|NRMSE / IoU|accepted / evaluations|corrector stop|','|---|---:|---:|---:|---|']
for r in case27:
    md.append(f"|{r['method']}|{r['pre_corrector_Phi']:.9g}→{r['final_Phi']:.9g}|{r['material_nrmse']:.6f} / {r['otsu_iou']:.6f}|{r['corrector_accepted_iterations']} / {r['corrector_evaluations']}|{r['corrector_stop']}|")
md+=['','0027的NNdirect在1814 accepted步按PGTOL正常停止（低于3000）；原四个NN辅助方法均达到3000 accepted上限且材料标签失败。它比微差符号计数更具体，但辅助若继续迭代会怎样仍未运行；不能据此断言不可逃离盆地。','',
 '## NNdirect未达继承材料复合门槛：公开条件及组件','', '|cohort|未达门槛|仅NRMSE|仅IoU|两者|','|---|---:|---:|---:|---:|']
for r in failure_counts:
    d=r['failure_categories'];md.append(f"|{SHORT[r['cohort']]}|{r['objects']-r['material_success']}|{d.get('NRMSE_only',0)}|{d.get('IoU_only',0)}|{d.get('both_NRMSE_and_IoU',0)}|")
md+=['','全部25个NNdirect未达门槛对象的case_id、contrast、源数、孔径、噪声、生成网格、NRMSE/IoU、corrector stop和PG norm在JSON逐行列出；另含全42×核心方法结果与公开条件。不从图像推断未公开的具体家族。', '', '## 缺口与解释边界','']
md += ['- '+x for x in limitations]
md+=['','证据行：每个JSON摘录原CSV `_evidence.artifact/line`；derived均值差/比值列出两个或三个原表行。原score逐case证据以rows索引和唯一task_id定位。所有CI均来自已有分析输出，未重新采样。']
(OUT/'FINAL_NUMERIC_MEMO.md').write_text('\n'.join(md)+'\n')
print(json.dumps(dict(status='EXTRACTED_VERIFIED_EXISTING_TABLES',core_rows=len(core),CPU_pairs=len(costs),fixed8_cohort_counts=result['scope']['extra_cohorts'],NNdirect_failure_cases=len(failure_cases),outputs=['FINAL_NUMERIC_MEMO.md','FINAL_NUMERIC_MEMO.json'],new_physics=0,new_GT_reads=0)))
