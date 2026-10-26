"""Render the completed independent confirmation from immutable score receipts."""
from pathlib import Path
import json,time
P=Path(__file__).resolve().parents[1]
def load(name):return json.loads((P/name).read_text())
def write(name,text):
    f=P/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text(text.strip()+'\n')
NAMES={'physics_probit':'Physics Probit','strong_cnn':'R4 Strong CNN','k32':'原 Gaussian K32','deterministic_matched':'确定性 Q64 + matched head','distilled_w128':'214k 蒸馏 W128'}
def table(means):
    return '|冻结动作|Full-grid complex MSE|\n|---|---:|\n'+'\n'.join(f'|{NAMES[m]}|{means[m]:.9f}|' for m in NAMES)
def main():
    st=load('RESULTS_AND_STATISTICS.json');sc=load('raw/offline_scores/SCORING_SUMMARY.json')
    cost=load('ONLINE_AND_END_TO_END_COST.json');ledger=load('RESOURCE_LEDGER_SNAPSHOT.json');var=load('ENSEMBLE_VARIANCE_DECOMPOSITION.json')
    raw=st['raw_primary'];m=raw['method_means'];ret=raw['retention_ratio'];lo,hi=ret['bounds']['ci_two_sided_98_75'];c=raw['comparisons']
    clip=st['clipped_secondary'];rn=ret['estimate']*100
    gain=lambda a,b:100*(m[a]-m[b])/m[a]
    detgain=gain('k32','deterministic_matched');studentgain=gain('strong_cnn','distilled_w128')
    root_summary=f'''# R5-B 独立确认：真实风险收益成立，80%保留率尚未确认

**512 个全新 Maxwell 对象、冻结模型、原 Native96 / Probit64、完整预测封存后评分，已全部完成。214272 参数蒸馏器相对 Strong CNN 的 MSE 降低 {studentgain:.2f}%，有配对统计支持；但增量保留率仅为 {rn:.2f}%，调整后的区间为 [{100*lo:.2f}%, {100*hi:.2f}%]，不能统计确认至少80%。**

{table(m)}

确定性 Q64 + matched head 比 K32 的 MSE 进一步低 {detgain:.2f}%，四比较校正后区间排除零。两批均同向。蒸馏器的 batch64 provider 费用降低97.97%，完整学习链降低46.32%；包含平均1.7303秒物理锚点时，推导的全系统节省仅0.0117%，不能称整个 Maxwell 系统快约49倍。

[完整中文报告](RESEARCH_REPORT_ZH.md) · [A批](NEW_TEST_A_RESULTS.md) · [B批](NEW_TEST_B_RESULTS.md) · [512配对统计](POOLED_512_STATISTICS.md) · [风险比较](MATERIAL_RISK_COMPARISON.md) · [完整成本](ONLINE_AND_END_TO_END_COST.md) · [裁决](GO_NO_GO.md) · [Pro交接](GPT_PRO_HANDOFF/START_HERE.md)

![固定独立对象，非最佳样本](figures/confirmation/FIXED_PARENT_MATERIAL_A.png)

实际物理费用 **1,999,504 / 3,300,000 RHS**，56,734 次状态/LU，未新增训练。原 R4/R5、旧验证集结果和原 CPU/MPS FAIL 不改写。新图是监督材料风险动作，不是 posterior 样本；本轮没有 FULL posterior、校准概率或新 Maxwell 专属机制的成功声明。
'''
    write('START_HERE.md',root_summary)
    for g in ('A','B'):
        x=st['batches_descriptive'][g]['raw'];r=x['retention_ratio'];a,b=r['bounds']['ci_two_sided_98_75'];con=x['comparisons']['retention_80pct']
        write(f'NEW_TEST_{g}_RESULTS.md',f'''# TEST-{g}：256个全新独立对象

本批原先只冻结seed、未生成GT或观测；本次实际运行原P1 Full Maxwell，并完成完整Native96+Probit64。使用与另一批完全相同的冻结模型和配置，所有预测封存后才评分。没有删除异常对象、修改checkpoint或重训。

{table(x['method_means'])}

增量保留率为 **{r['estimate']*100:.4f}%**；98.75%配对bootstrap区间为 [{a*100:.2f}%, {b*100:.2f}%]。80%线性contrast为 {con['estimate']:.8g}，单侧98.75%下界为 {con['lower_one_sided_98_75']:.8g}，不能确认达到80%。

这是预冻结的独立复验描述。pooled512的四比较family是主要推断范围；不能从两批中挑显著者作为额外主结论。确定性方案相对K32、蒸馏器相对Strong的方向与另一批一致。蒸馏器相对K32点风险略高，不将“差异未显著”解释为等价证明。

逐parent、逐seed的MSE/MAE、ROI风险、被动裁剪、图像ensemble与fallback在证据包 `raw/offline_scores/PARENT_METRICS.csv/.npz`；正式预测在 `raw/predictions/{g}/`。实际图采用事先固定的parent0/1、seed2026102601，非GT选优。

![固定新对象](figures/confirmation/FIXED_PARENT_MATERIAL_{g}.png)
''')
    comps='|主比较（正值有利于候选）|配对MSE差|调整后98.75%双侧区间|\n|---|---:|---|\n'
    labels={'deterministic_vs_k32':'K32 − deterministic','distilled_vs_k32':'K32 − distilled','distilled_vs_strong':'Strong − distilled','retention_80pct':'(Strong−distilled) − .8(Strong−K32)'}
    for k,v in c.items():
        a,b=v['ci_two_sided_98_75'];comps+=f'|{labels[k]}|{v["estimate"]:.8g}|[{a:.8g}, {b:.8g}]|\n'
    write('POOLED_512_STATISTICS.md',f'''# 512-parent 预冻结配对统计

统计单位是512个独立parent，不是1536个seed/parent组合。主动作先对三个冻结seed的对象级loss取平均；ensemble图像平均的风险另列。两批各256、pooled按批分层重采样，20,000次，共用分子分母和全部方法的parent抽样。四比较Bonferroni family-wise alpha=.05，双侧区间98.75%，retention单侧下界98.75%。脚本、RNG和loss数组完整提供。

{table(m)}

{comps}

K32新增改善分母为 {ret['denominator_estimate']:.9g}，调整区间排除零，20,000次pooled抽样中非正分母为0，故比率可解释。保留率 **{rn:.5f}%**，双侧98.75%区间 **[{lo*100:.3f}%, {hi*100:.3f}%]**，单侧下界 {ret['bounds']['lower_one_sided_98_75']*100:.3f}%。正式门以线性contrast为准：下界 {c['retention_80pct']['lower_one_sided_98_75']:.9g}<0，**HOLD**。

点估计略高80%，不等于统计证明≥80%。A为78.535%，B为82.106%，80%contrast在两批符号相反。总体比率按合并风险差计算，不是简单平均两个批次的百分比。蒸馏相对Strong的调整区间为正；det相对K32也为正；蒸馏相对K32区间跨零，既没有证明更优，也没有证明与K32等价。

被动裁剪为预声明次要敏感性：retention点估计 {clip['retention_ratio']['estimate']*100:.3f}%，低于80%，同样不支持该门。不能改用ensemble或只报raw点估计救结论。所有区间为有限parent的bootstrap推断，不是精确有限样本概率保证，也不包含无限训练seed的不确定性。

![主比较与保留率](figures/confirmation/RETENTION_AND_PAIRED_CONTRASTS.png)
''')
    detail='|冻结方法|Raw MSE|Complex MAE|ROI MSE|三seed图像ensemble MSE（次要）|\n|---|---:|---:|---:|---:|\n'
    for method in NAMES:
        q=sc['raw_seed_mean'][method];detail+=f'|{NAMES[method]}|{q["mse_full"]:.9f}|{q["mae_full"]:.9f}|{q["mse_roi"]:.9f}|{sc["raw_ensemble"][method]["mse_full"]:.9f}|\n'
    frozen_det=sc['raw_seed_mean']['det_frozen']['mse_full'];direct=m['k32']-frozen_det;total=m['k32']-m['deterministic_matched']
    write('MATERIAL_RISK_COMPARISON.md',f'''# 材料风险：泛化收益与门槛分别判断

{detail}

原K32相对Strong CNN的MSE改善 {gain('strong_cnn','k32'):.3f}%；确定性匹配方案 {gain('strong_cnn','deterministic_matched'):.3f}%；蒸馏器 {studentgain:.3f}%。相对同原R4配置的Probit点，deterministic与distilled分别降低 {gain('physics_probit','deterministic_matched'):.2f}%、{gain('physics_probit','distilled_w128'):.2f}%。该点中338/512报告成功，其余按原规则保留有限最佳点，不能称全部优化充分。这些是新对象上的监督材料风险，不是后验正确性的资格。

确定性 matched 比K32低 {detgain:.3f}%，配对校正区间支持正改善。额外冻结头控制 det_frozen 的MSE为 {frozen_det:.9f}，保持原head仅换输入已解释观察总改善的 {100*direct/total:.2f}%；matched-head额外差为 {frozen_det-m['deterministic_matched']:.9g}。该分解为描述，不能给不同训练seed/head制度的差异赋予单因素因果结论。K32方差为ddof0、确定性为population variance，故输入变化包含MC噪声及31/32方差尺度差，不是仅一项“积分精确性”。

三seedensemble的蒸馏MSE约与K32ensemble接近，但主门固定为seed平均loss。它们不是可互换的口径；不能用ensemble二次宣布80%主门通过。预测不是prior合法单一世界，passive-clipped也不等于原P1 decoder支持。未对每份学习图额外求Maxwell，不给它们赋noise-ball或posterior标签。

![风险与独立批次](figures/confirmation/RISK_AND_PAIRED_BATCH_EFFECTS.png)
![Raw、clipped和ensemble分开](figures/confirmation/RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE.png)
''')
    phase='|Anchor相区|ROI占比|Strong每单元MSE|K32|Det matched|Distilled|\n|---|---:|---:|---:|---:|---:|\n'
    for key in ['host_plateau','transition','inclusion_plateau']:
        vals=[sc['raw_seed_mean'][a]['anchor_phase'][key]['cell_weighted_mse'] for a in ['strong_cnn','k32','deterministic_matched','distilled_w128']]
        phase+=f'|{key}|{sc["phase_roi_proportions"][key]*100:.3f}%|'+'|'.join(f'{v:.8f}' for v in vals)+'|\n'
    contribution=lambda a,b:100*(sc['raw_seed_mean'][a]['anchor_phase']['transition']['full_grid_mse_contribution']-sc['raw_seed_mean'][b]['anchor_phase']['transition']['full_grid_mse_contribution'])/(m[a]-m[b])
    write('PHASE_RISK_ANALYSIS.md',f'''# Anchor相区误差机制

相区完全由原physics anchor场值/halfwidth在ROI内按≤−1、(−1,1)、≥1划分，没有用GT材料给在线算法分类。

{phase}

过渡相占ROI **{sc['phase_roi_proportions']['transition']*100:.2f}%**，不是薄薄一圈边界。Strong基准几乎全部风险已位于过渡相；因此大部分改善集中于它，不足以证明特殊Maxwell传播机制。Strong→distilled增量的 **{contribution('strong_cnn','distilled_w128'):.2f}%**、K32→deterministic增量的 **{contribution('k32','deterministic_matched'):.2f}%** 来自该区。百分比可受其他区负改善抵消，不能当像素成功率。

每单元风险和总贡献不同：前者除该相区实际单元数，后者除576并按parent平均。空相区有显式valid mask，不按零误差混入。纯inclusion仅约2.39%，其中属性估计与anchor相区误分类耦合，未把该组均值当真实posterior材料属性。

![相区占比、每单元风险与全图贡献](figures/confirmation/ANCHOR_PHASE_RISK_AND_GAIN.png)
''')
    write('DEVICE_NUMERICAL_REPRODUCIBILITY.md','''# 数值复现：canonical CPU与原设备失败

本轮全部正式预测使用预冻结单线程CPU、batch64、原float32网络和complex128物理。旧VAL的16个动作/seed组合，两个全新CPU实例逐数组重复一致；deterministic和distilled主要均值精确重现，Strong/K32与旧封存值差≤9.32e−10。优化K32与未优化CPU入口预测完全一致，原随机合同不变。

原R5 CPU/MPS逐像素3e−6容差曾FAIL，最大3.94517929e−6。该失败没有放宽、覆盖或改称PASS。独立CPU复现解决的是本轮canonical确定性，不替另一个设备的原失败翻案。本次受控运行环境MPS/CUDA不可用，新的MPS推断计时为NOT_AVAILABLE，不代表设备物理上不存在GPU。

新两批预测先以原字节和BLAKE2b封存，再评分。评分与最终QA从已存数组重算MSE/MAE一致，不用重新生成的容差图替换已评分图。详见CHECKPOINT_REPRODUCTION、32份prediction seal及FINAL_INDEPENDENT_QA。
''')
    gap=var['raw']['distilled_minus_k32'];vtable='|冻结方法|Seed平均risk|Ensemble risk|三checkpoint分歧项|\n|---|---:|---:|---:|\n'
    for method in NAMES:
        q=var['raw']['methods'][method];vtable+=f'|{NAMES[method]}|{q["seed_mean_risk"]:.9f}|{q["ensemble_risk"]:.9f}|{q["three_checkpoint_dispersion"]:.9g}|\n'
    write('ENSEMBLE_RISK_DECOMPOSITION.md',f'''# 封存后的解释：风险缺口与三个模型的分歧

这是预声明ensemble输出上的后验解释，不是新主检验，也没有重新训练、选择seed或改变图像。对复材料平方损失，有限集合恒等式为：seed平均risk = ensemble risk + seed对ensemble的均方分歧。仅是标准Hilbert空间平方范数恒等式，不是新Maxwell定理。

{vtable}

蒸馏−K32主risk缺口为 {gap['seed_mean_gap']:.9g}，分解为 ensemble缺口 {gap['ensemble_gap']:.9g} 与分歧缺口 {gap['dispersion_gap']:.9g}。在此固定三模型与512对象上，约 **{gap['fraction_gap_due_to_dispersion']*100:.2f}%** 缺口来自分歧项差。恒等式直接从材料数组核验，最大舍入误差 {var['raw']['max_identity_error']:.3g}。

这支持下一步检查student/provider+head训练稳定性和误差相关性，不能断言214k表示容量已充分，也不是无限训练seed总体方差的无偏估计。三seed ensemble支付三个动作的费用；K32的Gaussian矩也有共享权利，当前单seed计时不能直接乘三当最优ensemble部署成本。不能用该解释将未通过的80%非劣性主门改成PASS。
''')
    write('FAILURE_ANALYSIS.md',f'''# 失败、未确认和范围限制

1. **80%保留率统计门未通过。** pooled点估计{rn:.2f}%仅稍过门；单侧contrast下界为负；A批低于80%，clipped总体为{clip['retention_ratio']['estimate']*100:.2f}%。保留原VAL92.96%的开发事实，但不能再将其当独立泛化估计。新结果并不证明真实保留率低于80%，只表示不足以确认。
2. **蒸馏与K32等价性未建立。** 点MSE略高、差异区间跨零；缺乏显著差异不是非劣性或等价证书。相对Strong的改善有支持，故也不应把整个蒸馏模型称为无效。
3. **全系统加速不明显。** provider节省约98%是真，但平均1.73秒锚点主导，batch64全系统推导节省仅0.0117%。离线teacher与学生训练不能忽略，CPU/GPU秒不能直接合价。
4. **确定性更准确不等于更快。** canonical Q64完整学习链约为K32的2.46倍；网络部分较慢，端到端额外约0.037%。原Q16图不能冒用为Q64成本。
5. **数值/优化限制保留。** 原CPU/MPS门仍FAIL；338/512 Probit报告scipy success，174保留有限最佳点，未删除或看GT重启。无全局优化、posterior或所有anchor驻点资格主张。
6. **具体机制未完全归因。** 共同P1统计合同与历史一致，但有限样本波动、开发选择乐观和student/head分歧仍可能影响VAL到TEST差异。本轮没有证据证明分布改变或唯一归因teacher失配。旧VAL已选择过；禁止回头改模型再称独立。
7. **Maxwell专属增量未建立。** 过渡相占比大且风险原本集中，不能把区域改善当作传播信息不可替代性；上一轮RRR控制负结果继续保留。

实际数值求解接口失败0，未隐藏failed proposal或fallback成本。事后最差图是解释性GT选例，明确标注OFFLINE_ERROR_SELECTED，不能充当在线选图或新方法。微小fontcache警告后使用临时缓存正常渲染，未影响统计。

![事后误差选例，仅解释](figures/confirmation/OFFLINE_ERROR_SELECTED_DIAGNOSTICS.png)
''')
    decisions=dict(GO_1=dict(risk_gain='GO',retention_at_least_80pct='HOLD_NOT_STATISTICALLY_CONFIRMED',point=ret['estimate'],adjusted_ratio_CI=ret['bounds']['ci_two_sided_98_75'],contrast_lower=c['retention_80pct']['lower_one_sided_98_75']),
      GO_2='GO_DETERMINISTIC_MSE_GAIN_WITH_HIGHER_INFERENCE_COST',GO_3='GO_MODULE_COMPRESSION_NO_MEANINGFUL_END_TO_END_MAXWELL_SPEEDUP',
      GO_4='STRONG_IEEE_TAP_NOT_YET_QUALIFIED',GO_5='NO_ARCHITECTURE_SEARCH; versioned acquisition/noise and training-stability investigation next',
      posterior_claim=False,new_training=0,physics_RHS=ledger['RHS'],locked_parents=512,primary_family_unchanged=True)
    write('GATE_SUMMARY.json',json.dumps(decisions,indent=2))
    write('GO_NO_GO.md',f'''# 科研裁决

|问题|裁决|真实依据|
|---|---|---|
|GO-1：统计蒸馏是否有效|**风险泛化GO；80%保留率HOLD**|相对Strong MSE降低{studentgain:.2f}%，两批同向且校正区间为正；保留率{rn:.2f}%但区间[{lo*100:.2f}%,{hi*100:.2f}%]，非劣性下界<0|
|GO-2：确定性材料矩|**GO，准确度优先选项**|相对K32降低{detgain:.2f}%MSE，独立两批同向；Q64推断更慢，共同物理成本主导时额外全系统费用很小|
|GO-3：真正计算加速|**模块GO；整套Maxwell明显加速NO-GO**|batch64 provider降97.97%、完整学习链降46.32%；anchor-inclusive推导仅降0.0117%，另有历史与学生离线成本|
|GO-4：强IEEE TAP|**尚不GO**|获得独立可复核风险确认；尚无Maxwell信息不可替代性、跨采集/噪声、正确概率或新的物理计算机制证据|
|GO-5：下一步|**停止本轮调参；有条件转向稳健性与误差相关性**|固定现有模型，在另起版本/授权的noise/layout/contrast上检验；研究有限三seed分歧而非未经验证直接加大网络|

核心回答：214k蒸馏器能以显著更低的附加在线模块费用，在真正未知对象上复制一部分有统计支持的成像收益；但**尚不能确认稳定保留至少80%原Fusion增量，也不支持整个物理系统大幅加速**。该混合结论不能压成全PASS或全失败。
''')
    report=f'''# Research P3L-R5-B-Confirm：冻结独立泛化完整报告

{root_summary.split('[完整中文报告]')[0].split('\n',1)[1]}

## 1. 本轮真正解决了什么

原R5压缩结果来自被用于容量选择的VAL256。此次不训练、不搜索架构，不将已有测试改名；从原先冻结但尚未生成数据的两批清单实际生成512份原P1 Maxwell观测。独立成本探针8个只计费用，不评分。模型、prior、噪声、source/receiver、目标函数、Native96/Probit64和四项统计比较先冻结。模型推断与GT评价隔离，两批32份模型/seed文件全部封存后才评分。

数学和物理对象没有换成更容易的模型。原24×24被动TM/VIE、Gaussian64×Uniform2和C²饱和场原样，source-major复数据，sigma=0.00119211798514519。直接复用原完整优化器，不按噪声球提前停，也不按GT决定初始化或重启。338个Probit报告成功，174个保留有限最佳端点；所有对象进入总体。矩阵残差最大1.64e−15，不等同于全局逆问题已解决。

## 2. 主结果和风险定义

{table(m)}

平方损失采用复材料实Hilbert范数，除以全网格576。每parent先对三个冻结seed的loss平均，统计单位为512个独立parent。不是Best-of-K，也不是用GT挑一个seed。ROI、MAE、raw/clipped、ensemble与相区分层另表保存。[完整风险表](MATERIAL_RISK_COMPARISON.md)

两批retention分别78.535%与82.106%；pooled为{rn:.3f}%。K32增量分母稳定为正，但20,000次配对bootstrap四比较调整后的保留率区间为[{lo*100:.2f}%,{hi*100:.2f}%]。正式80%contrast下界为−5.43036e−5，因此HOLD。不能引用原VAL92.96%，也不能拿ensemble口径把本门改成PASS。相对Strong的蒸馏风险改善和相对K32的确定性改善均有调整区间支持。[统计细节](POOLED_512_STATISTICS.md)

![独立风险和配对效应](figures/confirmation/RISK_AND_PAIRED_BATCH_EFFECTS.png)

## 3. 输入统计与误差机制

保持旧head只换确定性矩，MSE为{frozen_det:.9f}，已解释此处确定性总改进的{100*direct/total:.1f}%；新matched head额外变化较小。这更支持“材料统计表示对局部估计有用”，不是新posterior正确性。K32的有限MC、population variance尺度和非线性head都参与，不以理论积分更准确直接推出MSE必优。

过渡相占anchor ROI的78.28%，Strong→distill增量约{contribution('strong_cnn','distilled_w128'):.2f}%来自它；纯相占比、每单元风险和总贡献一起报告。[相区分析](PHASE_RISK_ANALYSIS.md) 该现象不等于Maxwell专属传播机制，旧物理特征控制的负结果未被覆盖。

封存后另发现一个可精确核对的解释：蒸馏与K32的seed平均risk差为5.37095e−5，其中ensemble差仅2.76567e−6，余下5.09439e−5来自固定三checkpoint分歧，约占94.85%。这是标准有限集合平方损失恒等式，不是新定理；它提示训练稳定性/误差相关性值得研究，不证明容量足够，也不能用次要ensemble结果救主非劣性门。[分解](ENSEMBLE_RISK_DECOMPOSITION.md)

## 4. 成本：模块和系统必须分开

真实物理账本合计1,999,504 RHS：新测试truth4096、主512强锚点1,960,400、独立probe35,008；总预算3.3M，未使用部分保留。56,734次新状态/LU，物理接口失败0；反演收敛标志另列，未隐藏cap。没有新训练或云计算。

CPU1、batch64、每parent的p50：优化K32 provider207.606µs、全predict438.259µs；蒸馏4.209µs、235.236µs；Q64为847.091µs、1079.990µs。均有30次原始重复和p10/p90，不用重复encoder的低效API。Native+Probit平均1.7303秒，所以warm全系统推导节省仅0.0117%；这不是联合实测p50，更不是约49倍Maxwell系统加速。

教师数据、278,528 RHS的原joint生成、4,928,104 RHS的昂贵context、原Gaussian训练、29.15秒teacher缓存与约82.76秒三seed学生/head MPS操作等历史成本单列；不把GPU秒与CPU秒当等价货币。只回收teacher CPU的一项batch64回收点约143,596次部署，不是总经济break-even。[完整成本](ONLINE_AND_END_TO_END_COST.md)

## 5. 正确性、独立性和失败

37,240个历史ID和74,100历史RNG状态已做元数据审计，A/B和probe不重叠；清单精确沿用，不挑新样本。历史源码与ad71a86e/85b7固定版本核对一致。checkpoint和预测内容采用BLAKE2b身份核对，旧R4/R5未改写。原CPU/MPS3e−6 FAIL仍保留；本轮CPU复现和数值一致性独立通过。新物理端点立即无损持久化，没有首批状态数组丢失。

评分后的独立QA复算全部512×3×6动作的MSE/MAE，与已存分数一致；流式重放317,185账本行与实际solver计数一致。未认证FULL posterior、学习图Maxwell相容性、跨采集鲁棒性或Maxwell不可替代增量。所有这些缺口保留，不因论文目标改门。[失败与范围](FAILURE_ANALYSIS.md)

## 6. 最终判断

风险泛化成立，但≥80%保留率不够确定。确定性矩具有可复验MSE收益，Q64略增附加费用；蒸馏器具有很强的模块压缩，但现有流程端到端由物理锚点主导。现阶段证据更支持“全局功能统计+局部材料估计的风险—成本权衡”，尚不足以支持强IEEE TAP独立Maxwell新算法。

后续先保留当前冻结候选与negative evidence，在另起noise/layout/contrast协议中检验物理依赖和误差稳定性；不在本轮重新挑checkpoint或加大网络。给Pro的问题是：哪些统计必须保留、种子分歧与期望风险怎样权衡、怎样实际降低共同物理锚点费用。[科研裁决](GO_NO_GO.md) · [Pro交接](GPT_PRO_HANDOFF/START_HERE.md)
'''
    write('RESEARCH_REPORT_ZH.md',report)
    write('GPT_PRO_HANDOFF/START_HERE.md',f'''# R5-B-Confirm：给下一轮Pro的真实证据

本轮是完成的冻结独立确认，不是新训练。两批新256对象，原P1、Native96/Probit64、3seed、GT隔离，1,999,504 RHS，0训练。请先读完整报告、POOLED_512_STATISTICS、风险—成本和有限seed分解。

核心事实：Strong/K32/det/distilled MSE分别0.007221538/0.006943309/0.006854734/0.006997018。蒸馏相对Strong改善{studentgain:.2f}%，但retention{rn:.2f}%的调整区间[{lo*100:.2f}%,{hi*100:.2f}%]，80%门HOLD。raw/clipped、A/B、seedmean/ensemble全部分开。

三个最值得继续推理的问题：

1. **功能蒸馏的风险目标。** 给定Gaussian近似矩teacher不是正确posterior，怎样直接约束“共享local head作用后的材料risk”，而不是仅矩空间误差？哪些矩/空间correlation不可省？下一轮必须用新的开发和独立验证，不拿本512继续调参后称盲测。
2. **有限seed分歧的可预测性。** 标准复Hilbert恒等式显示，student与K32当前主risk缺口约94.85%由三个固定模型的分歧差解释；ensemble gap仅2.77e−6。怎样在不借GT选seed、不提高未知对象费用过多的条件下控制这一项？这只是有限三模型现象，不是训练分布定理或本轮80%门通过。
3. **真正的物理计算瓶颈。** provider降约98%但整套原流程仅约0.012%推导节省。能否保持材料风险与原观测信息，减少1.73秒共同Native+Probit锚点，而不是继续压缩几微秒？需要同精度、同信息、同先验成本对照，不能降低anchor后冒称原确认。

已支持：独立材料风险增益、确定性统计替换、模块费用压缩。未支持：≥80%统计资格、FULL posterior、Maxwell专属增量、强TAP论文核心或49倍系统加速。不要把这些现象重新整理成新综述，优先提出可证伪的风险/训练稳定性/物理费用机制。
''')
    # Machine-readable delivery summary; authoritative raw journal remains separate.
    summary=dict(status='COMPLETE',date_display='2026-10-27 Asia/Singapore',new_parents=512,probe_parents=8,
        data_and_model_frozen_before_new_GT_scoring=True,main_scores=m,retention=ret,comparisons=c,
        decisions=decisions,physical_ledger=ledger,new_training_updates=0,new_GPU_hours=0,
        no_new_physics_in_analysis=True,source_identity='configs/FREEZE.json',
        no_unqualified_posterior_claim=True,UTC=time.time())
    write('RESULTS_SUMMARY.json',json.dumps(summary,indent=2))
    write('RESOURCE_LEDGER.json',json.dumps(dict(status='CLOSED',hard_RHS_cap=3300000,unused_RHS=3300000-ledger['RHS'],physical=ledger,
          raw_authority='RESOURCE_LEDGER.jsonl (evidence archive)',timing='audit/timing/CONFIRMATION_TIMING.json',
          historical_not_billed_as_new='audit/HISTORICAL_OFFLINE_COST_SCOPE.json',analysis_and_predictions='ONLINE_AND_END_TO_END_COST.json',
          note='Nested physics/job timers must not be added; not all interactive/orchestration work has isolated CPU timing; observed process peaks are not continuous host aggregate peaks',
          new_training_updates=0,new_GPU_hours=0),indent=2))
    print(json.dumps(dict(status='REPORTS_RENDERED_FROM_ACTUAL_RECEIPTS',retention_percent=rn,new_RHS=ledger['RHS'])))
if __name__=='__main__':main()
