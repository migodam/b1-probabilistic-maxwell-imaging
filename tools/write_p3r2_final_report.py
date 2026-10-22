"""Write the research interpretation from closed, independently scored outputs.

No solver, training, estimator reselection, or raw truth is used here. This
script edits final prose/metadata only; original experiments remain immutable.
"""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr2-risk-posterior-20261023'


def read(path):return json.loads((OUT/path).read_text())
def put(name,text):(OUT/name).write_text(text.strip()+'\n')
def fmt(x):return f'{x:.9f}'
def pct(x):return f'{100*x:.2f}%'
def interval(x):return f'[{pct(x[0])}, {pct(x[1])}]'
def table(head,rows):return '| '+' | '.join(head)+' |\n|'+'|'.join('---' for _ in head)+'|\n'+'\n'.join('| '+' | '.join(str(x) for x in r)+' |' for r in rows)


def main():
    for path in ('FRESH_CHECKPOINT.json','POST_SCORING_CHECKPOINT.json','MECHANISM_CHECKPOINT.json'):
        if read(path).get('status')!='COMPLETE':raise RuntimeError('All independent experiments must close before final interpretation')
    first=read('L/evaluation/RESULTS.json');traj=read('L/trajectory_evaluation/RESULTS.json')
    confirm=read('L/mechanism_confirmation/evaluation/RESULT.json');oracle=read('L/offline_geometry/RESULTS.json')
    saturation=read('L/saturation_risk/RESULTS.json');phase=read('L/material_risk_mechanism/RESULTS.json')
    cdf=read('B/jobs/B_field2_direct_CDF_DEV030_r0/RESULT.json');bracket=read('C/jobs/C_projector_brackets/RESULTS.json')
    mala=read('B/jobs/B_metric_pair_DEV030/RESULT.json');ledger=read('RESOURCE_LEDGER.json')
    mse=first['summaries']['fullgrid_MSE'];mae=first['summaries']['fullgrid_MAE'];primary=first['primary_comparison']
    cp=traj['co_primary_vs_strong_probit'];rows=oracle['rows']
    rate=primary['relative_risk_reduction'];ci=primary['relative_reduction_CI_97_5']
    first_table=table(['材料估计','全网格 complex MSE','complex MAE'],[
        (label,fmt(mse[key]['mean']),fmt(mae[key]['mean'])) for key,label in (
            ('exactprior_material_mean','原先验材料均值'),('native','native TRF 点（96-state cap）'),
            ('probit','强 Probit 点（再64-state cap）'),('matched_point176','普通额外16-state对照'),
            ('direct','measurement-only 监督 MLP'),('fused','冻结的局部几何 CNN＋材料矩形投影'))])
    family_labels={'full_geometry':'完整灵敏度＋回投影','local_anchor':'仅基础图＋属性','capacity':'同参数、附加通道置零',
        'permuted':'附加图固定空间置乱','sensitivity_only':'仅训练谱灵敏度图','backprojection_only':'仅残差回投影图',
        'decoder_random_sensitivity':'仅随机谱 decoder 灵敏度图'}
    confirm_table=table(['第二批固定家族（三 seed 风险均值）','全网格 MSE'],
        [(family_labels[k],fmt(v)) for k,v in confirm['family_mean_risk'].items()])
    contrast_table=table(['第二批事前冻结的主比较','相对风险降低','98.333% 配对区间','实际裁决'],[
        (family_labels[r['candidate']]+' / '+family_labels[r['baseline']],pct(r['relative_risk_reduction']),
         interval(r['relative_reduction_CI']),r['descriptive_verdict']) for r in confirm['primary_comparisons']])
    interpretations=[]
    for r in confirm['primary_comparisons']:
        name=family_labels[r['candidate']]+' 相对 '+family_labels[r['baseline']]
        if r['difference_CI'][1]<0:sentence=name+'在这批确认对象上风险较低。'
        elif r['difference_CI'][0]>0:sentence=name+'没有增益，后者在这批确认对象上风险更低。'
        else:sentence=name+'的调整后区间跨零，独立增量没有建立。'
        interpretations.append(sentence)
    mechanism_text=' '.join(interpretations)
    gtweak={str(r):float(np.mean([v['prior_aware_weak_energy_fraction'][str(r)] for v in rows])) for r in (.02,.04,.08)}
    median=lambda key:float(np.median([v[key] for v in rows]))
    forward=read('L/jobs/locked_decision_forward/RESULTS.json')['rows']
    fwd=[]
    for name in ('classical_postprocess','fused','dynamic_curvature__arc_length_trapezoid','projected_direct'):
        selected=[r['Phi'] for r in forward if r['method']==name]
        fwd.append((name,f'{np.median(selected):.3f}',f'{sum(x<=282.8737921622988 for x in selected)}/256'))
    forward_table=table(['风险动作（真实 Full forward 检查）','Φ 中位数','noise-ball 内数量'],fwd)
    stage_table=table(['方向','新增 complex RHS','state attempts','已结算受控 wall 秒'],
        [(k,v['RHS'],v['states'],f"{v['controlled_wall_seconds']:.2f}") for k,v in ledger['stages'].items()])
    total_rhs=sum(v['RHS'] for v in ledger['stages'].values())
    first_family=[]
    manifest=read('METHODS_FROZEN_FOR_LOCKED.json')
    for family in ('local_anchor_cnn','local_geometry_cnn','local_capacity_cnn','local_permuted_geometry_cnn',
                   'weak_decoder','weak_strong_decoder','random_decoder','all_decoder_matched','u2_decoder'):
        vals=[mse['ablation__'+e['name']]['mean'] for e in manifest['ablation_checkpoints'] if e['name'].startswith(family+'_n512_')]
        if len(vals)!=3:raise ValueError('All fixed training seeds required')
        first_family.append((family,fmt(float(np.mean(vals))),' / '.join(fmt(x) for x in vals)))
    family_table=table(['首批次级方法','三训练 seed 平均 MSE','三个固定 seed 分列 MSE'],first_family)
    checks=read('AUDIT/locked_statistics/RESULTS.json')
    summary=dict(status='RESEARCH_COMPUTATION_COMPLETE',target_id='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',
        model_scope='2D scalar TM/VIE24, original nonlinear66-dimensional P1, same prior/noise/acquisition',
        first_fresh_parents=256,independent_confirmation_parents=256,train_parents=32768,validation_parents=2048,
        physics_context_train_subset=512,physics_context_validation_subset=128,
        primary_first_cohort=primary,second_cohort_confirmations=confirm['primary_comparisons'],
        first_cohort_material_MSE={k:v['mean'] for k,v in mse.items() if not k.startswith('ablation__')},
        independent_statistical_audit_status=checks.get('status'),
        statistical_exception='Strict secondary MAE ROI scaling alarm explained by storedfloat32 background; original alarm and raw scores preserved',
        FULL_posterior_qualified=False,posterior_probabilities=None,
        trajectory_arc_primary=cp,trajectory_scope='candidate-only endpoints and supervised risk actions, never posterior weights',
        decision=dict(supervised_material_risk='GO_WITHIN_FROZEN_P1',FULL_posterior='UNQUALIFIED',
            trajectory_risk_integration='NO_GO_CURRENT_POLICY',A17_dense_TM24='NO_INDEPENDENT_INCREMENT_ESTABLISHED',
            strong_TAP_claim='NOT_ESTABLISHED; mechanism specificity and external generalization still needed'),
        resource_snapshot_sequence=ledger['sequence'],snapshot_new_RHS=total_rhs,
        resource_scope='Snapshot before final report/packaging close; final ledger and separate publication receipt authoritative')
    (OUT/'RESULTS_AND_METRICS.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    put('START_HERE.md',f'''# P3L-R2：真实材料风险改善，完整后验仍未解决

本轮完成 32,768 个独立训练材料、2,048 个 validation 材料，以及两批分别冻结的 256 个全新测试材料。沿用原 P1 与真实 24×24 TM Maxwell，不训练新的三维模型，不把旧128测试对象冒称新盲测。

**首批主要结果：** 材料 MSE 从 **{fmt(mse['probit']['mean'])}** 降为 **{fmt(mse['fused']['mean'])}**，降低 **{pct(rate)}**；两项主比较调整后的 **97.5% 区间 {interval(ci)}**。217/256对象改善，39个变差，全部保留。普通额外16-state优化没有解释该收益。

这是一项**监督材料风险**结果，不是已认证 posterior mean。学习图仅138/256通过原noise-ball，不能作为256个新物理世界；66维 FULL posterior 仍 UNQUALIFIED。弧长轨迹平均未改善材料风险。

[完整中文报告](report.html) · [固定对象图库](gallery.html) · [原始 Markdown 报告](RESEARCH_REPORT_ZH.md) · [裁决](GO_NO_GO.md)

![首批独立材料风险](figures/main_v2/PRIMARY_RISK.png)

第二批独立机制确认采用21模型、7家族、3个训练seed，另有3项预注册主比较。{mechanism_text}

[监督风险](SUPERVISED_RISK_RESULTS.md) · [网络与消融](LEARNED_MODEL_COMPARISON.md) · [非Gaussian推断](NON_GAUSSIAN_POSTERIOR_REPORT.md) · [几何轨迹](GEOMETRY_TRAJECTORY_RESULTS.md) · [两批泛化合同](LOCKED_GENERALIZATION.md) · [数学定义](THEORY_AND_METHODS.md) · [独立审查](AUDIT/locked_statistics/FINAL_REVIEW.md) · [失败](FAILURE_ANALYSIS.md)

[GPT Pro 独立研究交接](GPT_PRO_HANDOFF/START_HERE.md) · [可复现入口](REPRODUCE.md) · [环境](ENVIRONMENT.json) · [费用](RESOURCE_SUMMARY.md) · [完整证据 Release]({RELEASE})

源码、表图与说明在原B1仓库；较大端点、链、检查点和完整费用账本在证据Release。本轮科学资格与公开可访问状态分别记录。''')
    put('SUPERVISED_RISK_RESULTS.md',f'''# 直接材料风险学习的实际结果

本轮没有等待MCMC teacher。独立原先验材料与原Maxwell noisy observation形成联合模拟数据，实际复数材料MSE直接监督风险估计器。3种数据规模4096/16384/32768、多个网络容量和三个训练seed均实际运行；模型选择只使用既有validation，不能按测试对象换checkpoint。

measurement-only MLP随独立parent增加而改善，但首批新对象上仍弱于物理点估计：MSE为{fmt(mse['direct']['mean'])}。大网络不自动有优势；仅512train物理上下文的小型共享局部CNN才取得主要风险改善。训练上下文是native48+probit48，测试强上下文native96+probit64；这个budget/covariate shift明确保留。

{first_table}

主模型5842参数，由验证集冻结seed33及材料矩形投影；alpha=1，实质就是选定CNN风险动作，不是测试后再融合。MSE下降{pct(rate)}、97.5%区间{interval(ci)}。全网格MAE也从{fmt(mae['probit']['mean'])}降到{fmt(mae['fused']['mean'])}。训练和研究费用不为免费；线上各方法共享同一实付物理重建，网络本身不新增Maxwell求解。

217人改善、39个对象变差。总体风险成功不保证每个人正确，也不保证任意新prior/噪声/采集几何泛化。第一个256-parent主结果不因第二批某个控制更好而重选。

## 风险收益不是后验资格

平方损失下 E[χ|y] 的最优性是数学事实；有限网络在新对象上风险较低，只支持一个较好的条件材料估计器，不能证明它已经等于该条件均值，更不能赋予候选概率。矩形投影包含所有原P1真值，逐材料平方误差不增加；它并不保证输出在单个decoder image中。

{forward_table}

上述是风险图本身的Full prediction，不是posterior predictive分布。强物理点全部位于noise-ball，学习风险图只有138/256位于其中。因此需要物理候选时应继续保留基础点；不能把学习图标为新增合法材料解。

## 后验评分后的机制归因（探索性）

在GT材料相分层，误差下降来自宿主、过渡和纯夹杂三个区域，均按parent重采样，未把像素当独立对象。按基础点的decoder相分层，{pct(saturation['point_decoder_phase']['transition']['fullgrid_risk_gain_contribution']['mean']/saturation['mean_risk_gain']['mean'])}的总MSE收益来自基础点仍处于过渡带的单元。它更符合修正基础成像的灰化/风险偏差，不是已经证明的饱和nullspace探索。

学习修正的{pct(saturation['correction_normal_energy_fraction'])}能量在基础点先验线性切空间之外，但该部分只贡献{pct(saturation['outside_tangent_risk_gain']['mean']/saturation['mean_risk_gain']['mean'])}的风险改善；大部分收益在切空间投影分量中。不能用“所有收益都来自曲率”解释本实验。此归因未改变任何在线方法。

![实际风险归因](L/saturation_risk/SATURATION_RISK.png)''')
    put('LEARNED_MODEL_COMPARISON.md',f'''# 学习模型、容量与几何控制

先前P3-L的15次训练全部update0最佳；本轮实际有多个小型模型选择正更新。它说明此前特定训练失败可被改变，不能推导所有NN无效；同样不能把本轮小CNN成功泛化成任意架构都有效。数据规模、模型容量和共享局部归纳偏好均有实际作用。

{family_table}

首批几何CNN相对仅基础图CNN、等参数零通道及固定空间置乱，三seed平均风险改善约3.45%/3.49%/3.47%。这些是次级比较，不是三个额外测试集。零通道控制只匹配名义参数量；固定位置置乱更直接检验空间对应。随机decoder比较的增量很小，调整区间跨零，不得宣称所有表示上几何都优越。

## 第二批独立定向确认

发现固定frame灵敏度图可因式分解后，加入两个通道去除控制和随机谱decoder控制。全部参数、旧512/128数据、100epoch及训练规则保持；两个控制六run和随机控制三run均只读TRAIN/validation。新21检查点在新256观测前冻结。

第二批为raw unprojected、无blend的公平家族比较，不能直接冒称第一批选定的projection动作。每个parent先平均三个固定训练seed的损失，不是先把三张图ensemble，不是n=768。

{confirm_table}

{contrast_table}

{mechanism_text}

不论显著性如何，随机控制只检验一个预定frame；它不能证明对全部随机子空间或所有Maxwell几何的定理。不同通道带来的有限样本学习价值，也不等于取得了新增测量信息。86个不同拟合/基线检查点中，首批实际评价77个、第二批评价其21个；完整范围见manifest，不按第二批结果换主模型。

## 归一化坐标密度

第一批固定10维informative density NLL为{first['density']['density']['mean_nll']:.6f}，自己的原先验基线为{first['density']['density']['unit_prior_mean_nll']:.6f}；risk-weak10维为{first['density']['weak_density']['mean_nll']:.6f}，自己的先验为{first['density']['weak_density']['unit_prior_mean_nll']:.6f}。两坐标不同，不能直接按原始NLL排胜负。

90%逐分量宏覆盖分别为{pct(first['density']['density']['mean_component_coverage90'])}和{pct(first['density']['weak_density']['mean_component_coverage90'])}，风险弱坐标仍显示小幅欠覆盖。它们不是10维联合覆盖率，更不是66维FULL或材料区间校准。把学到10维与独立原先验56维相乘，再合法Full-MH修正，仍没有解决高维耦合。''')
    put('NON_GAUSSIAN_POSTERIOR_REPORT.md',f'''# 非 Gaussian 推断：新机制、实际失败与条件参考

本轮并非重复上一轮单一GN Gaussian。实现完整支撑多中心三角Student proposal、逐照明源渐进同化SMC、真正state-dependent Full-J metric MALA、source-adjoint frozen-metric MALA、归一化learned proposal加exact MH，以及两个不同的Full-Maxwell二维条件目标。所有被拒提案、重复状态、祖先和失败费用保存。

## 完整66维 P1：全部仍 UNQUALIFIED

|方法|实际结果|新增RHS|含义|
|---|---|---:|---|
|多中心非线性Student独立MH|4096 production提案全部拒绝|32800|合法proposal与完整支撑不保证mixing|
|逐源SMC，DEV000|到达8源；ESS64但仅1个原先验祖先|199168|不能用最终weightESS掩盖genealogical collapse|
|逐源SMC，DEV030|到达8源；ESS57.19但仅1个祖先|123392|不是合格FULL reference|
|真实动态Full-J metric MALA|最小有限material ESS4.14，最大Rhat3.00|102400|刷新几何并未解决混合|
|冻结metric＋真实source adjoint梯度|最小有限ESS5.61，最大Rhat1.86|102400|高接受率也不意味着准确|
|learned10＋prior56，两个frame的完整校正|大量/全部拒绝，常数函数未评估|147456|低维NLL改善不能补齐56维conditional|

动态/冻结MALA各有4条链，production RHS均81920；能做到的transition数量不同，必须按实际RHS/MCSE比较。动态材料均值RMS MCSE约0.00512、冻结约0.01037，但两者都缺混合资格，不能据此宣布后验准确度胜出。常数链的MCSE0不是误差为零。

被拒Student提案的中位真实ΔΦ约40203，而latent线性预测约10.98。decoder饱和区跨越、材料secant和Maxwell余项共同导致失真，余项有显著负交叉项，不可将范数分数当作可加因果百分比。逐源SMC改善部分能量并没有获得可靠祖先覆盖；完整后验失败并未阻断监督和轨迹实验。

## 原属性方框 SLICE

固定64维材料场，仅在原两个uniform属性完整[0,1]²上积分。8/12/20阶分区积分与4条MH独立对照。两个DEV都达到该条件积分的预设数值细化停止标准；material mean链/积分RMS差约1.24e-5和5.70e-5，属性Rhat约1.0015–1.0040。整体posterior_qualification仍false，因为这不是积分掉其余64维的FULL marginal。

早期全局GH/GL欠分辨结果保留，没有覆盖。不能把近边界uniform属性当无界Gaussian或随手乘半Gaussian常数。

## 真正非线性 Gaussian field2 SLICE

在DEV030 paid anchor固定其余62个正交Gaussian坐标和两属性，absolute a∼N(0,I2)，使用真实非线性decoder与Full likelihood。它可以检验材料场Jensen effect，不是完整posterior。

GH9/17/33/65/97出现CDF节点混叠；后续对**同一目标**用t=NormalCDF(a)转换整单位方框，原Gaussian prior与Jacobian精确抵消，再用64块GL8/12/20直接积分CDF。未截断掉尾部，也不加Gram determinant。

最终logZ={cdf['final_logZ']:.9f}，材料Jensen gap RMS={cdf['final_decoder_Jensen_gap_material_RMS']:.9f}。12→20的材料均值RMS差为{cdf['rows'][-1]['material_mean_successive_RMS']:.3g}，切点CDF最大差约{max(cdf['rows'][-1]['direct_marginal_CDF_successive_max_difference']):.3g}；但logZ差{cdf['rows'][-1]['absolute_logZ_successive_difference']:.3g}仍高于冻结的1e-6，未达到连续两次全标准，**UNQUALIFIED标记保留**，没有为取得PASS放宽门槛。

GH97节点CDF与直接积分的最大差约0.276/0.279；不能将此误说成MCMC全部漏模式。已有pCN pooled ECDF与新积分差约0.0153/0.0164，material mean差约0.000481，位于其有限MC误差尺度。条件Jensen gap不等于总体成像收益，也不能解释CNN全部风险改善。

![同一条件目标的直接CDF](figures/main_v2/FIELD_SLICE_DIRECT_CDF.png)

## 明确未完成

没有合格66维FULL posterior、全材料可信区间、全局漏质量界或后验校准；不输出候选posterior probability。没有把二维SLICE称为FULL marginal。下一步应研究保留真实强弱耦合和非Gaussiandecoder边界的conditional transport，而非继续加同一个GN proposal的链长。''')
    diag=traj['actual_candidate_diagnostics']['dynamic_curvature']
    bracket_table=table(['已暴露anchor','h','normal leakage','材料回路/bracket相对误差'],[
        (p['parent'],r['h'],f"{r['normal_leakage_fraction']:.4f}",f"{r['loop_material_bracket_relative_error']:.4f}") for p in bracket['parents'] for r in p['rows']])
    put('GEOMETRY_TRAJECTORY_RESULTS.md',f'''# 局部几何与真实有限轨迹

所有轨迹均从只依赖观测的paid physics anchor出发，GT不参与提案、方向、修正或停止。Gaussian64方向使用材料RMS归一化的prior-aware广义谱；两个bounded属性固定。弱维数4，任务步长0.02/0.04，signed progress至少半步，decoder先验保持合法。

实际比较frozen、每步dynamic、curvature threshold、普通随机方向、contextual bandit和小型model-based horizon3 policy。bandit/MBRL用明确TRAIN轨迹风险监督；没有把访问频率叫posterior。PPO/SAC未扩展：已有行为数据只覆盖8个TRAIN parent，简单策略没有独立风险收益，不能以更大policy掩盖这个限制；这与FULL reference未合格无关。

## 首批256个对象的冻结主轨迹

256条曲率轨迹全部完成，共{diag['selected_endpoints_observed']}个被选端点；{diag['selected_prior_support_legal_count']}个prior合法、{diag['selected_radius_retained_count']}个保持实际步长、{diag['selected_noise_ball_count']}个noise-ball内、{diag['selected_noise_ball_and_delta_Phi_le_2_count']}个同时noise-ball内且相对固定基础点ΔΦ≤2。完整predictor/corrector共{diag['measured_predictor_corrector_endpoints']}次实付端点，incremental{int(diag['incremental_step_RHS'])}RHS。

因此存在任务尺度的有限噪声、prior-supported材料移动。但ΔΦ≤2是启发式likelihood竞争诊断，不是95%credible region；有限Gaussian prior penalty不等于某个区域有重要posterior mass。相邻端点相容也不是整条连续连接路径的证明。

弧长风险动作MSE={cp['candidate_mean_MSE']:.9f}，相对强点的风险降低为{pct(cp['relative_risk_reduction'])}，两主比较调整后的97.5%区间{interval(cp['relative_reduction_CI_97_5'])}，没有达到风险改进。额外16-state普通优化对照也不能将其变为成功。开发集约9.95%的弧长平均改善未复现在新对象上。

监督凸组合次级风险约0.01009647，区间跨零；不得事后替换失败的弧长主比较。前32人的五臂实际对照和全部失败保留；随机方向只7/256被选端点同时通过两个数据判据，说明几何对**数据可行移动**有作用，但不能据此宣布风险或posterior胜利。

## GT几何只能离线解释

新256个对象的GT−point材料RMS中位数{median('material_RMS'):.5f}；完整prediction距离中位数{median('full_whitened_prediction_distance'):.3f}。latent线性响应中位数{median('latent_linear_response_norm'):.2f}，材料secant线性响应{median('material_linear_response_norm'):.3f}，Maxwell余项{median('Maxwell_response_remainder_norm'):.3f}。对这些GT secant，decoder非线性是更大的线性化障碍，不能把它全叫Maxwell feedback。

任务尺度0.02/0.04/0.08的prior-aware弱能量平均分别为{pct(gtweak['0.02'])}/{pct(gtweak['0.04'])}/{pct(gtweak['0.08'])}。这不是此前单一obs_A的79.46%的普遍复现；阈值/先验度量不同也须分开。有限GT差量投影在弱空间不意味着可沿在线路径找到GT。

## 真实P1动态弱分布的Frobenius诊断

两个paid DEV anchor，固定rank4 projector定义材料RMS归一化场；两步长中央差分并做两种反向Euler回路，总48states/1920RHS，非toy代数。

{bracket_table}

DEV000两尺度材料bracket相对变化0.00236，DEV030约0.0236；probe材料rank保持64，弱/其余谱gap约0.28–0.30。回路/bracket误差随h减半而减小。这是该具体弱分布局部不对合的数值证据，提示不能默认把dynamic弱方向当单张坐标面。

它不证明所有P1区域不可积，不否定finite-noise不确定性，不证明不存在更高维路径或精确fiber。零宽轨迹的弧长law也不是原66维Lebesgue/Gaussian posterior law。Trajectory tubes、重叠质量与FULL posterior积分仍NOT_RUN，不伪造权重。

![轨迹风险与成本](figures/final/08_trajectory_risk.png)''')
    put('LOCKED_GENERALIZATION.md',f'''# 两批独立新对象与冻结合同

## 第一批：主要成像比较

原旧P3-L128个locked parent已暴露，仅作历史证据。本轮第一批R2_locked_00000–00255在第一份新观测产生前冻结，source revision `{first['source_revision']}`，77检查点/主融合/轨迹/失败fallback/额外16-state普通优化/统计规则均写入METHODS_FROZEN_FOR_LOCKED.json。256个在线点、全部神经输出和轨迹闭合后，独立evaluator才打开GT。

两项主比较以parent配对bootstrap2000次、seed2026102292，给出97.5%双侧区间。主风险降低{pct(rate)}，区间{interval(ci)}；弧长动作没有改善。无parent删除，无point/context/model fallback。三训练seed不能被当三倍独立parent。

## 第二批：机制的独立确认

根据第一批无法拆分通道的明确缺口，以及P1灵敏度恒等式，另建两份protocol设计记录和MECHANISM_CONFIRMATION_FROZEN.json。21个checkpoint及source revision `{confirm['source_revision']}` 在第二批第一份观测前登记。新IDs00256–00511使用原master seed规则下新的独立parent/noise，不是旧数据的新命名。

只使用旧512train/128validation训练9个新控制，未使用第一批GT作新标签，也未按第二批选winner。三项主比较用固定seed2026102393、2000次parent bootstrap、98.3333%区间。每family是三个固定训练模型的平均**损失**，不是ensemble图。

{contrast_table}

第二批估计均为raw/no-blend，不改变第一批投影融合主结果；两批不混作一个事前512人主检验。确认对象完成全部在线估计与物理子作业后才读取GT，全部256个对象及fallback进入分母。成本仍在原L3m/总20m预算中，不追加授权。

## 外推范围

仅验证同一合成P1 prior、原噪声和acquisition下的新parent。没有新噪声/新频率/新prior泛化，没有真实测量或三维vector-Maxwell外推。独立测试成功不消除模型误差或模型外材料质量。''')
    put('FAILURE_ANALYSIS.md',f'''# 保留的失败与局限

- **FULL probability失败仍在。** 原P3-L11组参考不合格、SMC祖先退化与Gaussian平均严重偏差保留。本轮新的Student、sourceSMC、动态MALA及corrected learned proposal也没有FULL资格；没有把它们升级为teacher或概率权重。
- **高接受率/高ESS误导。** sourceSMC终点ESS近64但只有一个祖先；constant chains的MCSE0是停滞而非准确。所有undefined函数保持未评估。
- **条件积分数值限制。** field2直接CDF显著减少GH节点混叠，但冻结的严格logZ细化标准未通过，没有改容差取得PASS。二维SLICE不等于66维边缘化。
- **轨迹开发收益未泛化。** 第一批主弧长估计风险反而增加1.96%，policy终点更差；次级凸组合也无可靠增量。它们产生可行材料移动，不等于材料风险或后验改进。
- **监督模型也会失败。** 大型物理CNN、普通残差MLP及measurement-only网络并未自动胜过物理基线；39/256个对象的主学习图更差，不删失败。所有候选网络只按既有validation选checkpoint。
- **材料风险图不都是新解。** 主风险图仅138/256处于原noise-ball。它是含原先验材料的凸矩形中的风险动作，不保证单个P1decoder、posterior sample或数据相容世界。
- **控制执行错误保留。** 一次对照误传200epochs在参数检查处退出，0新RHS、0训练；新run按原100epochs补齐，无删除重置。一次旧轨迹风险报告metadata导出失败亦保留，恢复只导出已冻结权重、未重训。
- **严格MAE检查警报保留。** 次级38个模型因保存complex64背景，在320个ROI外单元有常数约1.74e-9误差，导致全网格MAE比ROI×256/576高约9.65e-10；原FAIL、容差与评分未改。主MSE/主模型MAE不受影响，独立审查已解释，不能写“全部检查PASS”。
- **归因边界。** 几何通道的增量不等于全部25.27%收益。固定训练谱只通过空间包络进入该局部CNN；随机frame对照与第二批结论必须保留。有限训练seed和单prior的正面结果尚不足以作强IEEE TAP新颖性宣言。
- **未执行不伪造。** PPO/SAC、大型diffusion、三维/新Mie、完整trajectory tube概率、合格FULL材料credible intervals与A17新算力胜利均NOT_RUN或NOT_ESTABLISHED。资源充足不是继续消耗失败机制的理由。

失败对象、逐parent风险、全部拒绝、条件链、ancestor history、实际训练history及未成功作业见独立证据Release。未用GT进行任何在线修正或失败替换。''')
    put('GO_NO_GO.md',f'''# 本轮研究裁决

**GO：继续小型 physics-context 材料风险估计；FULL posterior 和轨迹风险不随之自动通过。**

{table(['路线','裁决','直接证据'],[
('监督材料风险','GO_WITHIN_FROZEN_P1',f'256fresh MSE降低{pct(rate)}，调整97.5%CI{interval(ci)}；强Probit及普通额外16-state对照保留'),
('特定几何通道机制','CONDITIONAL_GO / 以第二批分项结果为限','21模型、第二批256个对象、3项调整主比较；不把包络收益说成在线feedback'),
('完整66维posterior','UNQUALIFIED','多种独立机制仍混合/祖先/conditional耦合不足'),
('二维条件SLICE','数值机制证据','属性分区较好；field2CDF部分收敛，严格全标准未通过'),
('弧长/当前policy风险积分','NO_GO_CURRENT_POLICY','冻结fresh主比较无收益，不能换次级估计救结论'),
('A17 dense TM24 加速','NO_INCREMENT_ESTABLISHED','复用已封账负证据，没有新优势假定'),
('强IEEE TAP算法论文','NOT_ESTABLISHED','有风险结果和真实几何现象，仍需Maxwell特化增量、跨acquisition/模型及更强prior覆盖证据')])}

{mechanism_text}

下一阶段应以“低容量材料风险修正及其空间物理上下文”为主线，保留最强普通点/局部CNN/decoder-envelope基线。先明确训练谱到底提供何种独立信息，再决定是否需要新的物理输运模型。停止用同一失败Gaussian强方向补齐、只加链长、只增大policy/网络容量寻找正例。

真正Bayesian posterior研究可继续，但需要能保留剩余维conditional耦合的新机制，而非重复10维密度×prior56。未来研究不是本轮自动追加的新预算。

本裁决不把单个GT误差、低残差材料数量或数学恒等式当作概率证据。全部原P0/P1/P2/P3-L裁决和预算保持原样。''')
    put('RESEARCH_REPORT_ZH.md',f'''# P3L-R2：从失败的 Gaussian 平均到可泛化的材料风险估计

## 本轮回答了什么

这轮获得了明确的材料成像风险改善，也保留了完整后验和轨迹方法的失败。首批256个全新独立对象上，验证集冻结的5842参数局部物理上下文CNN，将强Probit点的全网格复材料MSE从 **{fmt(mse['probit']['mean'])}** 降到 **{fmt(mse['fused']['mean'])}**，下降 **{pct(rate)}**；两项预注册主比较校正后的97.5%区间为 **{interval(ci)}**。全部对象、39个变差对象和所有失败均进入分母。

这项正面结果属于 joint-simulation 监督材料风险；没有合格66维FULL posterior，也未证明学习图等于真实posterior mean。主学习图仅138/256通过原noise-ball，不能把它们称为256个新物理材料世界。冻结弧长轨迹平均风险增加约1.96%，没有获得预期增量。

![全新测试的材料风险](figures/main_v2/PRIMARY_RISK.png)

## 从上一轮失败中改变了什么

原P3-L完成800个独立parent、128个当时全新的测试对象和约837万RHS，但11组FULL参考全部未合格，Gaussian平均MSE0.061502、Probit点0.011031；15次神经训练全部选择update0。这些是封账的具体算法失败，不能推出真实Bayes均值不如MAP或材料概率不可学习。

本轮把三个目标独立推进：A用真实材料MSE监督条件估计；B使用新的非Gaussian inference并保留概率资格；C实际执行非线性有限轨迹和风险策略。没有让B的失败关闭A/C，也没有将A的成功转写为B的概率正确性。

## 固定的物理与概率合同

沿用P1原64维相关Gaussian场＋2维uniform属性、C²饱和decoder、24×24二维标量TM/VIE、256-cell ROI、8 sources×32 complex receivers、原背景与acquisition。σ固定0.00119211798514519，实虚噪声各σ²/2；Φ=||F−y||²/σ²。prior与noise没有为制造成功而放宽。所有Full是这个离散二维模型，不是3D vector或真实测量认证。

原native目标为Φ+||z||²/2；probit将u=NormalCDF(η)，原prior恰推送为Gaussian66，但该坐标密度的mode与native MAP不同。我们称它“Probit坐标点”，不混淆参数化依赖的MAP。GT仅用于生成、TRAIN监督与闭合后的独立离线评分，未用于在线init、提案、排序、stopping或fallback。

新32,768train、2,048validation；只有512/128子集建立physics contexts，避免为所有训练材料支付96-step反演。原来暴露过的对象仅作开发诊断。第一批新256个对象比较主风险和轨迹；随后另以固定新256个对象确认通道机制。两批分别冻结、分别评分，不混作一次事前512个对象检验。

## 材料成像的主要实际结果

{first_table}

普通额外16-state物理优化没有降低材料风险，控制了简单多迭代解释。measurement-only网络有数据规模增益，但仍弱于强点；大物理CNN未必优于小局部卷积。低容量、空间共享和合适物理上下文，而非模型名字或更多参数，解释了最直接的工程差异。

217/256个对象的风险改善，39个对象变差。主模型在原材料矩形上的投影具有逐真值平方误差不增性质，但矩形大于decoder image。网络和轨迹材料平均都作为**风险动作**评价，不能当posterior samples、事件概率或合法新解。

## 几何究竟贡献多少：两批对照

首批三固定训练seed中，几何局部CNN的MSE均值0.00766939，仅基础图CNN为0.00794370，额外收益约3.45%；同参数零通道和固定空间置乱也不能解释这点增量。但大部分主改善已经出现在小型局部anchor CNN中，所以不能说全部25.27%来自弱空间。

进一步代数检查发现：固定frame V下的sensitivity map恰为 `|δ| h'(Az/τ)/τ · ||A_i V||`。当前材料只通过decoder导数进入；Maxwell只通过冻结的TRAIN平均谱形成固定空间包络。这一事实要求拆开sensitivity/backprojection，并用独立J的随机包络作强对照。

第二批21个模型、7家族、三seed，全部checkpoint和3项比较在新观测前冻结；训练仍只读旧512/128。结果如下，未按新测试换主模型：

{confirm_table}

{contrast_table}

{mechanism_text}

![独立机制确认](figures/confirmation/MECHANISM_CONFIRMATION.png)

这个确认所能支持的，是具体输入通道与固定框架的有限样本风险作用；它不自动证明动态Maxwell几何、A17反馈或普遍采样定理。原P1训练先验之外的泛化仍缺证据。

![首批几何与合法decoder控制](figures/main_v2/GEOMETRY_CONTROLS.png)

## 风险改善发生在哪里

独立离线归因显示，基础点的过渡相单元贡献约{pct(saturation['point_decoder_phase']['transition']['fullgrid_risk_gain_contribution']['mean']/saturation['mean_risk_gain']['mean'])}的总风险改善；多数真实纯相在基础图中仍灰化。学习动作更接近纠正这类材料风险偏差，而不是发现一个已经认证的nullspace世界。

在材料RMS正交切空间分解下，约{pct(saturation['correction_normal_energy_fraction'])}的修正能量在局部prior tangent之外，这部分贡献约{pct(saturation['outside_tangent_risk_gain']['mean']/saturation['mean_risk_gain']['mean'])}的风险收益。大部分收益位于切空间投影；不能把全部改善解释为decoder外在曲率。这是评分后的探索性解释，不用于改变模型。

## 概率路线：独立推进，但不降低资格

新的三角Student完整MH4096次production提案全拒；逐源SMC确实到达全部8源，但终点ESS64/57仍只有1个先验祖先。动态Full-J metric MALA和冻结source-adjoint MALA在同102400RHS下都不合格，最大material Rhat约3.00/1.86。两个learned10＋prior56 exact-correction控制仍不能补齐高维conditional；MCSE0停滞函数保留为未评估。

低维参考也产生了具体新证据。原属性box的两个SLICE分区积分达到其细化停止标准，并与四链核对。真正field2非Gaussian SLICE从GH节点CDF改为同目标的prior-CDF整方框分区积分，切点CDF更稳定；但严格logZ差标准未通过，原UNQUALIFIED保留。条件Jensen gap约0.000749材料RMS，不能转写为FULL后验或解释全部成像收益。

固定10维联合模拟密度在首批新对象上NLL明显优于自身无条件prior；90%逐分量宏覆盖约89.02%/88.48%，后者仍有小幅欠覆盖。它们是投影坐标概率诊断，不是66维或材料可信区间校准。

## 几何与轨迹：可行移动不等于风险积分

动态曲率256条轨迹全部完成，2048个合法选中端点保留0.02/0.04实际RMS步长；1934个同时处于noise-ball且相对基础点ΔΦ≤2。确有任务尺度、prior-supported有限噪声材料运动，但这不是重要posterior mass或整条连续等数据路径证明。

主弧长动作风险却更差，调整区间不支持改善；次级凸组合和简单策略也未建立稳健收益。训练risk policy不需要MCMC teacher，我们已实际训练bandit及model-based规划并在线运行，失败没有归咎于概率reference门槛。

两个真实P1 anchor的rank4弱projector bracket诊断进一步显示normal leakage约0.93/0.82，两个差分尺度和反向回路一致到各自数值误差。它是当前动态弱分布非对合的局部数值证据，不能默认把它当单张流形；也不否定更高维finite-noise区域。48个states和1920RHS全计费。

GT−point离线弱能量随尺度显著变化：0.02/0.04/0.08平均{pct(gtweak['0.02'])}/{pct(gtweak['0.04'])}/{pct(gtweak['0.08'])}。原单对象79.46%不是全体事实，更不是在线找到GT的成功证据。

## 完整费用、内存和复现

{stage_table}

以上是写报告时的资源快照，共{total_rhs:,}RHS；最终打包/分析结算及完整追加式ledger另见[资源报告](RESOURCE_SUMMARY.md)。总硬上限20,000,000RHS、24h受控CPU、24GPU-hours未扩张，GPU实际0，不使用付费外部服务。计入truth、forward、receiver adjoint、source gradient、失败、拒绝、line-search和训练；历史已付background/setup与cache复用有单独记录。

只运行一个物理worker，transfer cache最多4个，原始endpoint立即持久化。可重复材料网络重放已在一份新对象上实际得到0差、0RHS；它只证明保存上下文下的重放，不把完整物理生成或部署setup说成免费。原A30外部源码版本、环境、依赖和大数组恢复方式见[复现说明](REPRODUCE.md)。

## 最终判断与下一步

本轮是 **监督材料风险 GO，FULL posterior UNQUALIFIED，当前轨迹风险积分 NO-GO**。A17在dense TM24没有已建立独立成本优势，不强行整合；已有负结果继续有效。最值得继续的是低容量材料风险模型及其可解释的decoder/空间物理上下文，而非重复失败Gaussian补齐或盲目扩大网络/policy。

目前有真实风险结果、定向控制和P1几何现象，但强IEEE TAP贡献尚未成立。需要证明哪些Maxwell特化设计超越普通重建后处理、decoder几何和已有amortized方法，并验证独立采集/噪声/材料支持。不能把CNN、LIS、Frobenius或Schur恒等式本身当原创贡献。

三个具体新理论问题及必要原始证据已整理为[GPT Pro Handoff](GPT_PRO_HANDOFF/START_HERE.md)。它们围绕：材料风险方向与信息方向的差别；饱和decoder下非Gaussian强弱conditional耦合；不对合弱分布如何支持合法、有效且值得计算的探索/积分。下一步须据真实证据重新设计，不自动续烧本轮预算。

[完整固定图库](gallery.html) · [学习控制](LEARNED_MODEL_COMPARISON.md) · [概率报告](NON_GAUSSIAN_POSTERIOR_REPORT.md) · [几何报告](GEOMETRY_TRAJECTORY_RESULTS.md) · [冻结泛化](LOCKED_GENERALIZATION.md) · [失败与限制](FAILURE_ANALYSIS.md) · [数学及文献范围](THEORY_AND_METHODS.md)''')
    handoff=OUT/'GPT_PRO_HANDOFF';handoff.mkdir(exist_ok=True)
    (handoff/'START_HERE.md').write_text('''# P3L-R2：给GPT Pro的实际研究交接

请先读[完整中文报告](../RESEARCH_REPORT_ZH.md)、[数学定义](../THEORY_AND_METHODS.md)、[首批独立审查](../AUDIT/locked_statistics/FINAL_REVIEW.md)、[第二批控制](../LEARNED_MODEL_COMPARISON.md)，再读本目录的三个问题。不要把本包当自动执行授权。

本轮给出了监督材料风险的正面结果，同时保留FULL probability与轨迹平均的失败。需要研究具体机制，而不是再做一轮泛化综述或扩大相同失败MCMC。

[真实新现象](ACTUAL_FINDINGS.md) · [三个最值得推导的问题](THREE_THEORY_QUESTIONS.md) · [证据索引](EVIDENCE_INDEX.json) · [复现说明](../REPRODUCE.md)
''')
    (handoff/'ACTUAL_FINDINGS.md').write_text(f'''# 可用于理论推理的实际现象

1. 首批256fresh，固定局部CNN风险降低{pct(rate)}，调整CI{interval(ci)}；不是后验reference结论。第二批256fresh的21模型拆分控制结果单独保留，未换第一批winner。
2. P1灵敏度图精确因式分解为decoder导数×固定空间包络。{mechanism_text} 这要求区分Maxwell frame、decoder几何与卷积空间归纳偏好。
3. 约{pct(saturation['point_decoder_phase']['transition']['fullgrid_risk_gain_contribution']['mean']/saturation['mean_risk_gain']['mean'])}风险收益位于基础点过渡相；修正13.03%能量在prior线性切空间外，但仅约5.83%收益来自该部分。不能简单说NN只在nullspace填充或曲率解释全部风险。
4. 66维Student、sourceSMC、动态MALA、learned10×prior56都没得到合格FULL reference；projection NLL有改进但全维校正仍失败。新Proposal的真实ΔΦ与latent线性预测差数千倍，decoder/Maxwell余项存在抵消。
5. 相同field2 SLICE的直接CDF积分与GH节点CDF差约0.28，已有pCN与直接CDF差约0.016；不能将坏quadrature误差全归于sampler。Jensen gap约7.49e-4，严格logZ refinement仍未过。
6. rank4动态弱分布在两个真实P1 anchor上表现为非对合，normal leakage约0.93/0.82，差分和反向回路有数值对照；不是toy反例，也不是全P1定理。
7. 1934/2048个所选轨迹端点在固定任务尺度和数据判据下合格，但冻结弧长平均风险未改善。数据可行移动、材料风险动作与正确后验积分是不同对象。

所有数学身份、真实Maxwell结果、探索性离线归因和未执行猜想分开；本包不提供FULL posterior probabilities，不提供原P1之外的mass。
''')
    (handoff/'THREE_THEORY_QUESTIONS.md').write_text('''# 三个具体的新理论问题

## 1. 材料风险方向与Maxwell信息方向如何统一，而不是事后冠名

已知真实risk改善、等容量/空间置乱控制，以及第二批的Maxwell训练谱与随机decoder-envelope控制。请推导固定观测函数下，何种可计算的局部风险算子能够优于仅按J最小奇异值分割；明确哪些信息只是decoder导数/位置编码，哪些需要真实Maxwell谱。应给出可证伪的必要对照，而不是声称LIS首次用于逆散射。还需解释为何局部材料输出优于受限latent MLP、而大部分真实风险收益又位于prior线性切空间。

## 2. 饱和decoder与强弱耦合能否形成可归一化的非Gaussian conditional transport

10维proper NLL改善不能补齐余下56维；多中心/full-support/高weightESS都不足以保证探索。请利用实际被拒提案的decoder secant/Maxwell余项、sourceSMC祖先丢失、field2连续CDF与材料矩，构造有明确密度/Jacobian及合法Full校正的提案。需要说明它解决什么条件耦合，如何避免强方向Gaussian体积失配，并预言一个比已有控制更高的material ESS/RHS或更低的MCSE，不能只加flow层数。

## 3. 非对合弱分布下，什么探索或概率积分值得支付物理费用

真实P1的rank4弱场bracket有显著normal leakage；有限噪声路径仍能走，但弧长平均无风险收益。请研究以风险为目标的横向补偿、tube/多chart测度或随机核：什么时候必须加入弱分布以外的方向，什么时候局部采样的统计收益小于一次监督风险动作？提出含decoder余项、Maxwell余项、谱gap/rotation及统计误差的可检验成本判据。任何不变分布或probability claim均须给出参考测度、反向概率/重叠权重；零宽trajectory occupation不是后验。

以上都需要新的实质数学或明确反例，不要求维护A17、Gaussian、单流形或某个NN架构。未来实验与训练必须有新阶段授权，不能自动复用已封账余额。
''')
    evidence=dict(scope='Actual R2 evidence, no new budget authority',sources={
        'first_frozen':'../METHODS_FROZEN_FOR_LOCKED.json','confirmation_frozen':'../MECHANISM_CONFIRMATION_FROZEN.json',
        'first_risk':'../L/evaluation/RESULTS.json','confirmation':'../L/mechanism_confirmation/evaluation/RESULT.json',
        'saturation_risk':'../L/saturation_risk/RESULTS.json','oracle_geometry':'../L/offline_geometry/RESULTS.json',
        'brackets':'../C/jobs/C_projector_brackets/RESULTS.json','same_field_slice':'../B/jobs/B_field2_direct_CDF_DEV030_r0/RESULT.json',
        'full_MALA':'../B/jobs/B_metric_pair_DEV030/RESULT.json','new_ledger':'../RESOURCE_LEDGER.jsonl'},
        probability_scope='SUPERVISED_RISK or explicitly conditionalSLICE; FULL remains UNQUALIFIED',
        external_dependency='A30 frozen backend commit717d3d26a290248bd2b9aeb55f3d03694b407b84',
        raw_evidence_release=RELEASE)
    (handoff/'EVIDENCE_INDEX.json').write_text(json.dumps(evidence,indent=2,ensure_ascii=False)+'\n')
    print('Final research documents written from closed outputs; no model or raw score changed')


if __name__=='__main__':main()
