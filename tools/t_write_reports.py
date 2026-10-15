"""Write substantive frozen Research T reports from actual saved results."""
from pathlib import Path
import json,csv
import numpy as np
from research_t.target import OUT

WEB='https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_t_20261015/'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-t-20261015'
def write(name,text):(OUT/name).write_text(text.strip()+'\n')

def main():
    m=json.loads((OUT/'METRICS.json').read_text());o=json.loads((OUT/'offline/AGGREGATE_METRICS.json').read_text())
    archive=m['archive'];res=m['resources'];k=m['offline_best_of_K']
    assert all(abs(r['best_MAE']-r['base_MAE'])<1e-14 for r in k)
    rows=list(csv.DictReader((OUT/'PER_EDIT_METRICS.csv').open()))
    case='p250065_ring_a1_snr20';ex=[r for r in rows if int(r['version'])==1 and r['case']==case and int(r['edit_id'])==0][0]
    table='|轮次|基础规则／主要变化|raw通过|补偿后通过|独立parent有成功|\n|---|---|---:|---:|---:|\n'
    for v,label in [(1,'最低残差base＋普通全材料补偿'),(2,'同提案＋非core正材料相保护'),(3,'最小N3改动的合格base；固定两个编辑')]:
        a=m['iterations'][v-1];table+=f'|{v}|{label}|0/{a["proposals"]}|{a["invalid_to_valid"]}/{a["proposals"]}|4/4|\n'
    write('START_HERE.md',f'''# Research T：Maxwell物理、基础反演与补偿成核

**已实际执行并封口。8/8条件有合格默认基础材料；40个固定结构提案产生32个合格补偿端点记录。数据相容成功，材料精度没有改善。**

{table}

跨轮按相同观测、固定材料距离与相位门槛构建archive，48个保存材料端点（16基础变体＋32补偿结果）全部保留。每份ring观测9／10／12／9个；每份weak观测2个。数量不是独立连通解分支、后验模式或完整覆盖认证。

- [中文研究报告](RESEARCH_REPORT_ZH.md) · [成像与散射图库](offline/index.html) · [裁决](GATE_DECISION.md)
- [G物理交接](T_PHYSICS_HANDOFF/README.md) · [默认base与48份真实材料索引](T_PHYSICS_HANDOFF/CANDIDATE_INDEX.json)
- [逐端点指标](PER_EDIT_METRICS.csv) · [离线成像误差](offline/OFFLINE_SCORING_REPORT_ZH.md) · [完整成本](RESOURCE_SUMMARY.md)
- [原始数组与可复现证据Release]({RELEASE})。下载解压至仓库根目录，恢复相同路径；读取结果无需重跑。

**必须同时读的负结果：**最低残差base的平均complex MAE为0.147656，原N3为0.041756；改为最小材料改动选择后为0.050757。补偿候选的MAE约0.49–0.51，Best-of-4包含基础图时改善为0。通过噪声球不能代替成像准确、posterior质量或统计不可区分。

新授权消费192/192 full-state、5832/6144 complex RHS；受控求解jobs共213.213秒。GPU0，新3D0，旧余额0。剩余312RHS不能越过已耗尽的state预算。历史A30/A/B、Joint P2和A17冻结稿保持原样。

范围为历史暴露的四parent×weak/ring，24×24二维scalar TM；不是独立盲验证、连续网格证书或一般三维Maxwell成像。所有概率为null/candidate_only。公开状态与科研结论分别记录。
''')
    write('RESEARCH_REPORT_ZH.md',f'''# Research T 独立研究报告：能补偿出新孔洞，但未改善成像

## 第一屏结论

1. **基础解：**固定8条件均有合法、数据相容的默认基础材料，原N3只有2/8。原32条件零物理审计为2/32。两套基础选择各8个变体已保存。
2. **成核：**三轮共40个有限编辑，raw全部不相容；32个经实际full-TM端点验证后相容，13cell严格空气核及闭合材料guard均保留。第一轮15/16成功，覆盖四个独立parent。
3. **区别：**同观测材料之间满足固定0.03归一化L2和至少6cell相差的archive共48个端点；这是有限栅格材料集合，不是48个拓扑分支证明。
4. **精度：**所有轮次Best-of-up-to-4含基础图时材料MAE改善为0；未建立更准确重建或更快搜索。没有运行multi-start、QD或CEM竞赛，不能对它们作优越性结论。
5. **A17：**独立增量NOT_RUN。普通全材料Jacobian已使用完整resolvent；移除它制造“无feedback”对照会改变物理。没有证据把A17列为当前生成器必要组件。
6. **交付：**G可直接复用8个minimum-change默认base、共同corrector/validator、48份真实材料与预测、完整约束及独立预算接口。
7. **裁决：**继续受限二维补偿机制与G的独立生成研究；停止最低残差即最佳材料、linear score即非线性排除、候选数即posterior质量三类主张。当前不足以支持强IEEE TAP原创算法论文。

## 为什么先审计

收到的ZIP实际包含执行说明、T↔G接口和一个manifest示例三份文件，均原样保存于pro_original；未假定还有未提供的大型Pro实证库。原A30源码、A17理论、Research A反证与其已有物理数组为主要证据。E0没有求解、分解或构造新Green矩阵。

32个历史条件保存的clean response、current、本构与场方程一致，最大相对误差1.97e-15，噪声metadata重放误差0。原N3缓存与pool只有质量打包带来的4.45e-16差；Research A一步更新后的base差0.03–0.07，不能混用残差。历史图的abs(chi-background)是增量支持，不能当空气；p250065_weak的历史三个阈值均7个增量组件、0洞。32真值均无严格χ=0的单元。Sparse64的35.69%是真值目标离线拟合容量，不是测量到材料的在线恢复。

四个真值噪声实现落在95%噪声球外，其中p250065_weak真值残差24.090925，阈值23.785449。这不是现有模型错配证据。所有实验保持原sigma和阈值；也不能把某一次truth不在置信球解释为其物理非法。

## 冻结的物理与资格

原24×24二维scalar TM/VIE，8照明、32接收、576个复材料，单一chi=epsilon_r−1共享全部源；Rechi∈[-.5,4]、Imchi∈[0,.5]，已知gain=1。原complex RMS噪声sigma，逐源Re32/Im32并乘sqrt(2)/sigma，共512实通道，tau=sqrt(chi²_512(.95))。没有Mie台架的实噪声约定移植、自由gain、噪声放宽或真值在线输入。

普通全材料1152实自由度与Haar64从相同N3开始，最大4更新／6新端点、trust0.6、LM damping0.001，两臂每条件均一次新端点就达到0.95tau。资格阶段16states/448RHS。数据相容恢复说明历史失败不能当不存在可补偿材料；它没有证明材料更接近真值。

|基础|合格条件|平均残差|complex MAE|Re MAE|Im MAE|
|---|---:|---:|---:|---:|---:|
|原N3|2/8|29.139042|0.041756|0.039688|0.009997|
|full GN|8/8|20.148396|0.167560|0.141362|0.062463|
|Haar64|8/8|21.162010|0.054811|0.048762|0.017934|
|最低残差选择|8/8|20.136726|0.147656|0.123295|0.057465|
|最小N3位移选择|8/8|21.927005|0.050757|0.045796|0.015730|

最后一行受离线汇总负结果启发后提出，但每case选择只比较合法相容端点对N3的材料L2；没有用truth选图。这是exposed development修订，不是独立确认。

## 三轮实际机制实验

{table}

固定位置(7,7)/(7,16)/(16,7)/(16,16)，空气核半径2cells，共13cells；guard宽1cell、Rechi≥.30。core为严格chi=0，guard与core均为corrector硬约束。guard被提高的材料也是提案的一部分，并非免费校正。所有core在原base的13cells都不是primary近空气，避免把原孔重命名为新孔。

v1采用全材料普通GN，最多3更新／3新端点；v2仅增加非core原Rechi≥.15处保持≥.15；v3只改用minimum-change base并在所有四parent固定前/第三位置执行，其他物理与corrector不变。raw与中间失败全保存；不得把v3换base说成与v1同一raw。

v2使core外新增near-air细胞均值从29.44降至11.06，但成功率从15/16降到10/16，材料MAE没有改善。v3仍7/8成功，说明第一轮并非仅能从最差full-base成核；同时其合法候选MAE均值仍0.488023。32个通过记录跨轮相关，不能当32个独立重复。

## 一个可核验实例

同一p250065_ring观测、gain=1：基础残差20.050，固定edit0 raw残差{float(ex['raw_residual']):.6f}；三步有界补偿后{float(ex['corrected_residual']):.6f}≤23.785449。13cell core=0、闭合guard≥.30，归一化材料L2={float(ex['normalized_material_L2']):.6f}。精确有限secant恒等式在全部40个端点最大相对缺陷4.10e-15。

![真实材料、基础、raw与补偿端点](offline/figures/mechanism_v1_p250065_ring_a1_snr20_edit00.png)

它证明“有限结构编辑＋合法全材料补偿”能恢复冻结的观测相容；它不证明truth有该孔。该例MAE由基础0.1528变成0.5093。v1通过候选与base预测差距仍6.24–11.36白化单位，v3为10.07–14.27；即使双方在噪声球内，也不能称同一分布、统计不可区分或高posterior mass。概率字段保持null。

## 评分、A17与数值边界

约束linear compensation给局部trust box的可行上界U与合法dual下界L；没有uniform非线性误差界。v1有13个U>tau但后续三步full-nonlinear成功，因此用U>tau直接删除会误排有效端点。U排序只作付清所有raw/32-RHS评分成本的离线顺序重放；未建立独立online screen加速或ordinary multistart优势。

实际24网格单侧差分步长1e-4/1e-5的Jacobian相对误差1.07e-4/1.07e-5，呈一阶收敛；独立组装NumPy直接解与SciPy端点一致。两者仍共享离散算子/LAPACK家族，不是独立物理模型或连续误差证书。构成一致性和backward residual约1e-15并不提供输出误差上界。所有候选为FULL_VALIDATED_NUMERICAL，严格bound状态UNKNOWN_NUMERICAL_ERROR。

A17 paired law针对固定材料注入及指定Galerkin增广的局部二次问题。这里没有遗漏反馈的current ROM。新的core/guard改变active face；不能照搬无约束双侧步差定理。等成本paired/Krylov额外消融未运行，原A17局部GN结论不受改写。

## 资源、失败与后续

192states、5832complex RHS、213.213秒受控求解wall；详细split和零物理后处理在RESOURCE_SUMMARY。一次离线CSV失败保留，wall未记录为null，没有假装零费用。无大型NN、GPU、新3D或旧预算借用。新solver消费封口。

下一阶段先把minimum-change／discrepancy regularization与结构约束结合，保护材料质量；不要用更大的生成器弥补可行集合过宽。G可在同一后端独立授权下比较生成策略；T本轮不重复。非同心3D与强TAP理论声明均NOT_RUN/证据不足。有限可行端点不解决唯一性，不证明连通性或完整解族覆盖。
''')
    write('E0_AUDIT.md','''# 零新增物理E0审计

完整逐条件结果见[e0/E0_DECISION](e0/E0_DECISION.md)、[root-cause矩阵](e0/ROOT_CAUSE_MATRIX.md)、[原始审计JSON](e0/E0_AUDIT.json)、[mask语义](e0/MASK_SEMANTICS.md)。32条件均为历史exposed development。

原N3为2/32相容，固定8为2/8；Research A一步GN后base为3/8。真值cached clean完全匹配，噪声重放误差0，物理恒等式缺陷≤1.97e-15，未发现MODEL_MISMATCH证据。原图增量不等于空气，p250065_weak不是真空多孔真值；历史阈值为7增量组件0孔洞。32个真值没有严格χ=0单元。Sparse64离线目标容量没有在线恢复器。

分类：原不相容base记INVALID_BASE；增量暗区误读为空气记TOPOLOGY_MISINTERPRETATION；更完整corrector已找合法base，支持OPTIMIZATION_FAILURE是历史瓶颈之一；旧低秩图表失败不排除REPRESENTATION_FAILURE。后来的相容不同端点属于POSSIBLE_AMBIGUITY的数值见证，不是独立posterior证据。没有把高残差直接归为MODEL_MISMATCH。

E0 source/cache可复原，未阻断本轮独立192state资格实验。新增states/RHS/LU=0，两个成功审计脚本wall合计4.949秒。缺失缓存与Research A不同端点的来源差异均保留，未为历史缺项补跑。
''')
    write('MECHANISM_TEST.md',f'''# 有限结构编辑与补偿

{table}

三轮raw均0通过，所有成功均满足原观测tau、全球被动材料盒、共享gain1、本构/场检查及core/guard。core13cell，物理半径2/24；守护层1cell。材料距离为体积RMS/固定1.2，phase差≥6cell、primary absχ≤.1，另报.075/.125敏感性。这些是离散几何定义，不是连续persistent homology。

v1 raw残差265.99–394.47，经补偿后通过端点为20.24–23.79内；单个失败27.691。v2有6失败；v3有1失败。所有失败保留完整数组和line-search/KKT记录。

[逐编辑成本与物理距离](PER_EDIT_METRICS.csv) · [same-observation配对材料/观测距离](PAIRWISE_VALID_MATERIALS.csv) · [图与离线误差](offline/index.html)。每条成功core的13个格点在原base都不是primary近空气；纠正器不能用填孔或开通气道冒充成功。

成功伴随全域补偿与额外相变化；材料质量没有提高。噪声球相容与低Bayesian evidence可能并存，不报告任何posterior probability。没有验证全部中间路径相容或证明disconnected branch。
''')
    write('THEORY_REFINEMENTS.md',r'''# 理论修订：材料注入、局部评分与有限成核

原离散Maxwell方程为 Lχ c = diag(χ)e_inc，Lχ=I−diag(χ)G；e=e_inc+Gc，y=Sc。c/e均沿用原质量归一化，χ是raw susceptibility。所有照明同一χ。

## 精确有限注入不是任意dark current

对于两个已求得的合法材料χ0、χ1，令Δχ=χ1−χ0，则

`L1(c1−c0)=diag(Δχ)e0`，故 `F1−F0=S L1^{-1} diag(Δχ)e0`。

这来自原方程相减，不是原创恒等式。相较Jacobian `J0d=S L0^{-1}diag(d)e0`，有限编辑使用新材料resolvent L1。40个实际端点的恒等式缺陷最大4.10e-15；这验证本轮injection与场/材料合同，不是continuum认证。没有将任意receiver-dark current当成材料候选。

## 有界线性补偿的正确范围

对当前状态及单个trust box，`min_{l<=d<=h} ||r+Ad||`。任意||v||≤1给合法下界

`L(v)=max(0, v^T r + Σ_j min(l_j (A^T v)_j, h_j (A^T v)_j))`。

任一合法d给上界U。上下界对这个凸线性问题成立。只有额外取得对同一整个box的uniform nonlinear remainder ε，才可在该box用L−ε>tau排除或用某个实际试点加误差界支持接受。多步更新改变box/resolvent，当前L不能排除后续路径；U>tau只是该attained点失败，甚至不是线性最优下界。

实际反例：v1有13个U>tau但三步非线性通过。不得把单步upper、未收敛primal残差、局部dual或profile minimum升级为全材料不可能性或posterior odds。

## 成核结构必须是约束，而非展示掩膜

core固定χ=0但不保护周围材料，不能保证封闭孔；guard保持Reχ≥.30且闭合，使所有core在8连通空气定义下被包围。最低6cell及2cell半径排除单像素阈值收益。仍有额外非core相变化，v2保护使其减少但损失成功率，说明补偿自由度有实质作用。

## 数据相容和成像质量不同

减小||Wr||未必减小材料误差，尤其1152实材料自由度对应512观测并有弱方向。观测数本身不是有效秩定理，但本轮full vs Haar的实证明确暴露过拟合。minimum-change规则只在已付费合格端点集合中选最靠近N3者；没有声称求得全局最小改动解。

## A17的作用边界

A17 paired law是在固定材料注入、残差、先验和指定Galerkin差异下刻画局部二次步差。全材料exact adjoint的S Lχ^{-1}本身包含完整feedback。增设空气core/guard会改变活动约束；只有固定free face后重新限制H/J才有相应局部线性结论，不能沿用无约束公式认证有限拓扑变化。当前没有相同成本paired/Krylov付费比较，独立增量NOT_RUN。未改写A17冻结稿。

## 能声明与不能声明

可以声明：特定exposed离散TM目标下，实际构造了多个材料上显著不同、满足同一数据球和被动盒的端点。不能声明：统计不可区分、后验显著、完整覆盖、连续可行路径、全局不唯一性解决、网格收敛、一般三维效果或算法原创定理。标准GN、有限resolvent恒等式与约束LS本身不构成新颖性。
''')
    write('METHOD_COMPARISON.md','''# 方法比较与未测部分

|方法|实际运行|可支持的结论|
|---|---|---|
|同N3全材料GN／Haar64资格|两臂8case各一次有效更新|均获得8/8数据相容；Haar较少材料偏离|
|相同有限编辑，无补偿|40个实际raw端点|0/40相容|
|普通受约束full GN|三轮40条轨迹|32条相容且保住结构；成像MAE无改善|
|额外正材料相保护|相同16raw与步数|10/16 vs15/16；额外near-air减少，精度无收益|
|minimum-change合格base|同4parent固定8编辑|7/8；baseline质量较好但候选MAE仍大|
|constrained score排序|同paid端点轨迹离线重放|有完整取得score成本；不是新online加速试验|
|A17 paired／同成本Krylov|NOT_RUN|未建立独立机制增量；完整普通Jacobian已有feedback|
|QD/CEM／全搜索multistart|NOT_RUN，本任务属于G分工|不能宣称T胜过其搜索效率|
|新3D／GPU NN|NOT_RUN／不授权|无一般3D外推|

RANKING_COST_REPLAY.csv中U排序先支付每个raw的1state/8RHS与32-adjoint初始评分，再支付被修正端点的剩余成本；不会把已知成功顺序当免费oracle。每个编辑实际约128RHS，其中首次Jacobian可由score与corrector共享，另有两次refresh。score的额外小矩阵优化wall未独立分离；所有pipeline wall包含在总job费用。

OFFLINE_BEST_OF_K.json按基础图再固定edit次序，K=1/2/4/8；不足K标实际数量，不复制补齐。所有case含base的Best-of-up-to-K均无MAE收益。完整MAE、support、ROI、Re/Im、边界和几何结果见offline/MATERIAL_METRICS.csv。truth只作离线评分。
''')
    write('ALGORITHM_ABLATIONS.md','''# 实际消融

1. 同N3、同物理的full与Haar64：数据阈值两臂通过，但full MAE更差；这是材料表示/正则自由度的重要负结果。
2. 40个完全相同编辑的无补偿版本均失败；普通nonlinear补偿有32条通过，不能归因于只做一个免费几何编辑。
3. v1与v2相同raw和三步预算，只增加非core正材料相保护。额外新增near-air均值29.44→11.06，而通过15→10；非core自由度参与补偿。
4. v3更换minimum-change base，其余corrector不变，固定8个提案7个通过。它保留现象但不是独立盲测，也没有提高材料MAE。
5. constrained screen只是付费轨迹顺序重放；A17独立增量、G的adaptive/QD、真正3D均NOT_RUN，不能从合成代数通过推断这些增量。
''')
    write('GATE_DECISION.md','''# Research T 最终科研裁决

|问题|状态|证据|
|---|---|---|
|原物理合同、noise和cache可复原|PASS|E0全32条件，0新增solve；实际24网格资格通过|
|同一target获得合格基础材料|PASS|8/8默认base；最小N3改动规则不读取truth|
|有限空气核＋guard经合法补偿成核|PASS_NUMERICAL_RESTRICTED|15/16、10/16、7/8；四parent均有成功|
|材料显著不同|PASS_DISCRETE_ARCHIVE|48保存端点，same-observation固定门槛；非连续分支证明|
|真实成像MAE／Best-of-4改善|NO_GO|所有轮含base的Best-of-K无改善；候选MAE约.49–.51|
|constrained score加速／优越性|NOT_ESTABLISHED|仅付清成本的顺序重放，无独立online加速试验|
|A17独立增量|NOT_RUN|普通完整Jacobian已有resolvent，无未校准paired比较|
|完整posterior／统计不可区分|NOT_ESTABLISHED|candidate_only，probability=null；噪声球不等于显著posterior|
|强IEEE TAP原创理论／算法贡献|INSUFFICIENT_EVIDENCE|标准机制的受限pilot；缺少独立泛化与同成本新算法优势|
|G物理交付|READY|8默认base、48真实bundle、共同corrector/validator、独立预算|
|继续T新物理|STOP_BUDGET_EXHAUSTED|192/192states；剩312RHS不可借用|
|新3D/GPU|NOT_RUN|未授权|

继续受限补偿机制及材料质量约束；不要继续把minimum residual、候选数量或local linear score包装为概率/成像成功。向G交付稳定后端，G用其独立真实授权开展搜索。所有历史Gate和论文未改。
''')
    write('FAILURE_CASES.md','''# 保留的失败与限制

- 原N3相容失败30/32，其中固定子集6/8；Research A一步base5/8失败。它们是无效base，不是“无其他解”证明。
- full-space GN平均MAE0.167560，高于N3 0.041756；最低残差base选择0.147656。修复相容没有修复成像。
- v1：p250069_ring/edit1最终27.691464>tau；v2六个、v3一个失败；40个raw都失败。逐端点路径/数组和物理费用未删除。
- v2减少额外空气却降低成功；v3换base不消除补偿后材料误差。所有Best-of-K无改善，不继续调参制造正结果。
- 32truth无严格真空单元。生成空气核是另一种数据相容解释，不能描述为恢复真实孔洞。
- E0有4/32 truth-noise球外实现；这不等于model mismatch，原sigma/tau未改。
- 离线最终CSV首次刷新因重复字段TypeError失败；原始异常和wall=null/UNRECORDED在offline/FAILED_REFRESH_RECEIPT.json。后续修复没有新物理。
- 物理求解数值异常0不等于无失败：raw、失败中间端点、最终不相容材料均已计费。未记录的历史setup/单次后处理wall保持未知，未写成0。
- 没有uniform output error bound、continuum mesh验证、完整覆盖、same-cost A17增量、独立生成泛化或3D证据。
''')
    write('NEXT_STAGE_PROPOSAL.md','''# 下一阶段建议：先保护材料质量，再谈新生成器

本阶段full-state授权已耗尽；以下均未执行，不自动申请或启动新物理。

优先假设：在相同噪声球内，用最小材料变化/适当空间正则限制corrector，可以在保留严格core/guard时减少不必要的全域材料误差。应比较discrepancy-constrained correction和同成本普通GN，固定允许支持，不把正则当免费prior或改变数据噪声。若结构无法保留，应如实记录可行集/表示限制，不训练更大的网络。

G使用本交接后端、默认minchange base与真实成本，独立比较其生成机制。T不重做QD/CEM。若新事件只有很差成像且很低合理prior支持，优先转向event-aware追加测量或约束集合，而非宣称发现重要posterior模式。

A17 paired只有在明确低秩修正成本瓶颈并建成固定活动face的等成本Krylov对照后再考虑；本阶段没有证据支持强行加入。非同心3D目前不具备成像收益依据，暂不扩展。新预算由用户另行明确授权，不继承剩312RHS或其他旧stage余额。
''')
    write('RESOURCE_SUMMARY.md',f'''# 完整成本与独立授权

用户Research T独立overlay：192full-state、6144complex RHS、7200秒受控CPU；GPU0、新3D0。不叠加Pro建议，不借A30/A/B/Joint P2余额。

|阶段|full-state|complex RHS|受控job wall|
|---|---:|---:|---:|
|资格Q，两臂8case|16|448|5.217s|
|机制v1，16编辑＋base重分解|68|2112|80.830s|
|机制v2，同16编辑|68|2112|83.836s|
|机制v3，8编辑|36|1088|以JOB_END原值为准|
|实际Jacobian/独立直接解核验|4|72|以JOB_END原值为准|
|**合计**|**{res['new_full_state_attempts']}**|**{res['complex_full_solver_RHS']}**|**{res['controlled_execution_wall']:.6f}s**|

RHS独立核算：192×8=1536 forward；8个Q初始Jacobian＋40编辑×3刷新＋1核验Jacobian=4128adjoint；21次缓存LU核对×8=168，合计5832。没有免费adjoint或失败重试。

第一瓶颈state预算已耗尽，剩312RHS不可单独支付新fullstate。RUN_CLOSED.json及唯一Book在求解前拒绝超额；历史ledger原样只读。新账本COST_LEDGER.jsonl包含所有STATE/RHS/ACTION/GN试点/JOB_BEGIN/END，global snapshots不能逐bundle再求和。

独立零物理E0成功wall4.949秒；离线三个成功后处理累计19.250秒；导出/分析各receipt另列。一项失败后处理wall未记录（null），所以213.213s是完整计费solver jobs wall，不是全部人工实现与文件处理end-to-end wall。cached Green/S矩阵运算、保存背景setup/历史摊销未知额外列为null，不冒写0。

未运行同配置冷/暖重复物理计时；各新进程job是实际测得的流水线wall，已有效支付cache读取和refactor。旧bundle可零物理重放，但I/O快不能称warm correction加速。没有用不同轮算法差异冒充cold/warm对照。
''')
    write('README.md',f'''# Research T 复现入口

[独立中文研究报告](RESEARCH_REPORT_ZH.md) · [图库]({WEB}) · [完整证据]({RELEASE})。

代码在仓库src/research_t/，Python使用现有NumPy/SciPy，CPU BLAS单线程。没有新增包安装、GPU训练或重写Maxwell solver。原A30 commit717d3d2、其frozen A29 commit5d5e31f复用；source、原始输入子集和所有NPZ在Release恢复相同路径。

读结果不需要求解：先看METRICS.json、PER_EDIT_METRICS.csv、offline/和T_PHYSICS_HANDOFF/CANDIDATE_INDEX.json。将证据zip解压至仓库根目录后可运行保存数组检查：

```sh
PYTHONPATH=src python research/research_t_20261015/T_PHYSICS_HANDOFF/validate_handoff.py research/research_t_20261015/T_PHYSICS_HANDOFF/candidates/p250065_ring_a1_snr20__minchange_base
```

它不调用solver、不检查SHA、不建立新Gate；只校验真实保存材料与其预测、白化残差和合同。需要重画可用tools/t_offline_results.py；真值只在该离线脚本中。tools/t_analyze.py重算已付费指标，tools/t_export_handoff.py导出已付费端点；不会新建Maxwell状态。

实际已执行顺序（历史记录，不在封口run上重跑）：

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run qualification
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism2
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.run mechanism3
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python -m research_t.verify_actual
```

各付费阶段代码冻结commit为c752aa6、439773b、1e1b613、a42845f；原协议与收到的Pro包保留。默认Targets通过A30_ROOT指向现有A30或证据中的T_PHYSICS_HANDOFF/vendor_backend/a30_root。任何新物理都需要新的真实用户授权与独立enforcing book，不可删除RUN_CLOSED/重命名job重置预算。

原始cache、truth和offline摘要只作历史/离线证据。G的在线接口加载sanitized背景和默认base，不加载teacher J、truth或oracle step。具体用法见T_PHYSICS_HANDOFF/README.md；它不是另一套Maxwell solver。
''')
    prior=(OUT/'ITERATION_LOG.md').read_text().replace('结果待实际执行后补写。','实际10/16通过，core外新增near-air均值29.44→11.06，但合法端点MAE未改善；累计152states/4672RHS/169.883秒。减少额外相变化要付出可补偿性代价。')
    if '## 3：' not in prior:prior+='''\n## 3：minimum-change合格基础\n\n离线聚合显示最低残差选择MAE0.147656，原N3为0.041756，Haar64为0.054811。仅改为已合法相容材料中最小N3位移选择（不读取truth），全8默认base MAE为0.050757。四parent固定indices0/2、同ordinary corrector共7/8通过；合法候选MAE0.488023，Best-of-K仍无收益。不是独立盲验证。新增36states/1088RHS，累计188/5760。生成代码1e1b613。\n\n## 实际数值资格与封口\n\n另付4states/72RHS，对真实24网格作两步长Jacobian差分与独立组装直接解；全部40端点有限注入恒等式缺陷≤4.10e-15。累计192states/5832RHS，213.213秒求解job，停止所有新物理。A17独立增量未测，G搜索未重复。代码a42845f。\n'''
    write('ITERATION_LOG.md',prior)
    # First screen is a scientific decision page; full image gallery stays linked.
    content=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research T · Maxwell补偿成核</title><style>body{{font-family:system-ui,-apple-system,"PingFang SC",sans-serif;background:#f3f6f8;color:#172b38;margin:0}}main{{max-width:1100px;margin:auto;padding:38px 24px}}section{{background:white;border:1px solid #d7e1e8;border-radius:12px;padding:24px;margin:24px 0}}h1{{font-size:34px}}p,li{{line-height:1.8}}img{{max-width:100%;height:auto}}a{{color:#1467a0}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px}}.n{{font-size:30px;font-weight:700}}table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #ddd;padding:12px;text-align:left}}.warning{{border-left:5px solid #bb7434}}</style><main><nav><a href="../../">B1总入口</a> · <a href="offline/index.html">全部材料图库</a> · <a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_t_20261015">报告与代码</a> · <a href="{RELEASE}">完整证据下载</a></nav><h1>Research T：结构可补偿，成像未改善</h1><p>同一份二维Maxwell观测允许多个显著不同的合法材料端点。本轮修复基础解、检验严格空气核的非线性补偿，并公开保留材料误差增大的负结果。</p><div class="grid"><section><div class="n">8 / 8</div>合格默认基础材料</section><section><div class="n">15 / 16</div>第一轮无效编辑经补偿合格</section><section><div class="n">0%</div>含基础图Best-of-4 MAE改善</section><section><div class="n">192 / 192</div>full-state预算已封口</section></div><section class="warning"><h2>通过物理残差，不等于重建正确</h2><p>原N3平均材料MAE 0.041756；最低残差基础图0.147656；最小改动的合格基础图0.050757。补偿后的候选MAE约0.49–0.51，没有成像改善。概率为null；48个保存端点不是48个已证明独立的拓扑分支。</p></section><section><h2>三轮真实实验</h2><table><tr><th>实验</th><th>无补偿</th><th>补偿后合格</th><th>结论</th></tr><tr><td>普通全材料补偿</td><td>0/16</td><td>15/16</td><td>可保住严格空气核与包围材料</td></tr><tr><td>额外保护原有正材料相</td><td>0/16</td><td>10/16</td><td>额外空气减少，但成功率下降</td></tr><tr><td>最小改动合格基础图</td><td>0/8</td><td>7/8</td><td>现象保留，MAE仍无收益</td></tr></table></section><section><h2>同一测量的一个真实见证</h2><p>p250065_ring：固定有限编辑使残差385.079；有界补偿后22.791，小于原阈值23.785。13个格点的严格空气核与闭合材料guard均保留；材料MAE从基础0.1528变成0.5093。</p><img src="offline/figures/mechanism_v1_p250065_ring_a1_snr20_edit00.png" alt="离线真值、基础、raw编辑与实际补偿后材料，实部和损耗使用相同色标"><p><a href="offline/index.html">浏览8条件基础图、40条编辑、散射拟合与全部失败 →</a></p></section><section><h2>可以交给G什么</h2><p>8个minimum-change默认base、48份材料/观测/预测/内部场交接包，原24×24 TM后端、共同corrector/validator和独立预算接口。每份ring观测保留9／10／12／9端点，每份weak观测2端点。A17独立增量与G的搜索算法比赛未运行。</p><p><a href="T_PHYSICS_HANDOFF/README.md">交接说明</a> · <a href="T_PHYSICS_HANDOFF/CANDIDATE_INDEX.json">机器可读索引</a> · <a href="METRICS.json">指标JSON</a> · <a href="RESEARCH_REPORT_ZH.md">完整中文报告</a> · <a href="GATE_DECISION.md">科学裁决</a></p></section><p>受限、历史exposed二维scalar TM pilot；无一般三维、posterior或完整覆盖认证。192states、5832complex RHS、213.213秒受控求解wall；GPU0。继续材料质量约束与补偿机制，强IEEE TAP算法贡献证据不足。</p></main></html>'''
    write('index.html',content)
if __name__=='__main__':main()
