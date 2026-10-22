"""Build explanatory reports from measured artifacts; no Maxwell calls."""
from pathlib import Path
import json,sys,shutil,html
import numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'tools'))
import summarize

def main():
    summarize.main()
    load=lambda p:json.loads((root/p).read_text())
    cfg=load('w_config.json');g=load('results/GATE_DECISION.json');cost=load('results/COST_SUMMARY.json')
    a=load('results/p0a/results.json');grid=load('results/qualification_source_average/WALL_GRIDS.json')
    dd=load('results/qualification_source_average/DERIVATIVES.json');fb=load('results/four_blocks/RESULTS.json');pair=load('results/REFERENCE_PAIR.json')
    sigma=.01
    four_error=max(r['full_field_equivalence_relative'] for r in fb)
    born_min=min(r['born_induced_relative_error'] for r in fb);born_max=max(r['born_induced_relative_error'] for r in fb)
    beta_min=min(r['core_sigma_min'] for r in fb);beta_max=max(r['core_sigma_min'] for r in fb)
    rows='\n'.join(f"|{r['seed']}|{r['mesh']['wide_cal_mesh_whitened']:.6f}|{r['mesh']['wide_task_mesh_whitened']:.6f}|{r['mesh']['narrow_cal_mesh_whitened']:.6f}|{r['mesh']['narrow_task_mesh_whitened']:.6f}|" for r in grid if r['nx']==192)
    bands='\n'.join(f"|{r['fractional_bandwidth']*100:g}%|{r['wall_distance']:.6f}|{r['joint_shifted_target_distance']:.6f}|" for r in a['bandwidth_sweep'] if 0<r['fractional_bandwidth']<=.05)
    pairrows='\n'.join(f"|{r['band']}|{r['d_cal']:.6f}|{r['d_fixed_target_total']:.6f}|{r['d_fixed_target_induced']:.6f}|0|" for r in pair)
    logs=list((root/'results/jobs').glob('*/OUTPUT.log'))
    assert any('Ran 20 tests' in p.read_text() and '\nOK' in p.read_text() for p in logs)
    report=f'''# Research W P0：未知不规则墙体的传播歧义验证

日期：2026-10-22。研究状态：**P0-A PASS；P0-B UNRESOLVED（二维前向精度未达标）**。

## 1. 先说本轮能够回答什么

Pro 的一维反例被独立正确复现，包括共同参考面、完整 PEC 多次反射级联和带宽负结果。二维确实执行了有限尺寸、非均匀、有损、不规则墙体的多频多源完整多次散射计算，保存了材料与复场数组；没有用 Born 或手绘场图替代全波。

但二维平台尚不足以支撑本题要求的强歧义判决。最细两层192×144与128×96网格之间，完整多频数据的白化差仍为 **{g['grid_whitened_min']:.3f}–{g['grid_whitened_max']:.3f}**，超过原协议的0.1。这个误差已与“两个世界距离≤1”的关键门槛同阶。不能因为线性系统残差小、伴随正确或者图片看起来正常，就去宣布数据近等价的物理反例。

因此本轮停止在资格审核：普通Full-J、随机、multistart FWI和task-aware见证搜索均未启动；calibration与joint的见证不存在“通过”记录。**这不是证明歧义不存在，不是二维宽带唯一性结论，也不是Research W物理假说NO-GO。** 未耗尽总预算；停止原因是冻结的前向资格阶梯未通过，不能用剩余搜索额度绕过Gate1。

## 2. 原资料与数学审计

ZIP实际包含26个文件：52页PDF、完整LaTeX章节、解析脚本、原始JSON/CSV及图。输入解压到只读 `source_package/`，清单只记录路径和大小；没有执行哈希校验或备份流程。父线程核对核心数学，Codex子代理承担来源、实现和图表整理；没有调用DeepSeek。

**纤维判据。** F相同必然H相同，当且仅当存在集合意义的H=h∘F。可微或局部Lipschitz因子化才给出ker DF⊂ker DH的必要性。原文已经写出可微h的限定；本轮补充F(x)=x³、H(x)=x说明，单射且连续可辨并不自动满足零点的导数核包含。邻域恒秩、每点核包含和局部连通纤维片，才允许积分成局部因子化；不能把一个小奇异值当成有限非线性纤维。

**任务歧义模数。** 找到一个合法非线性端点对可以给ωT(δ)下界；未找到不能给其全局上界。小观测导数、大任务导数只能提出初猜，不是物理见证或后验质量。

**四块公式。** 对L=L_W−UVU*，d=RGQ、A=RGU、B=U*GQ、C=U*GU，有F=d+AV(I−CV)⁻¹B。非交换顺序正确，V必须包含本构、频率与离散质量因子。它属于标准resolvent/T-matrix结构，不是新原创代数。相关逆不存在或接近不稳定时，不能套用统一误差界。

**后验界。** 原联合prior不变、likelihood相同且势能误差uniform受控时，保守TV≤tanh ε成立。替换墙代表时若丢掉目标/仪器条件依赖，结论不再自动成立。它是未来P1的条件工具，本轮没有建立后验、事件概率或压缩证书。

详细条件和状态标签见 [THEORY_AUDIT.md](THEORY_AUDIT.md)。没有发现一维相位、完整目标级联或四块次序的核心错误；机器精度下原归档与本次独立实现的数值分别保留。

## 3. 一维回归：反例成立，特定宽带补偿失败

采用(n,d)=(2,1/4)、(3,1/6)，同前表面，固定后参考面D=1；PEC位置L₁=1.5、L₂=1.5−1/12。行波界面矩阵、独立Fresnel求和、完整PEC边界三条计算路径相互核对。

|量|本次独立实现|
|---|---:|
|中心最大墙反射模|{a['tests']['max_center_reflection_abs']:.6g}|
|共同参考面双程相位差|{a['tests']['center_task_phase_difference_rad']:.12f} rad = π/3|
|同位置目标双程响应差模|{a['tests']['center_task_response_difference_abs']:.12f}|
|PEC位置补偿后中心联合差|{a['tests']['joint_center_data_difference_abs']:.6g}|
|最大无损幺正性缺陷|{a['tests']['maximum_lossless_unitarity_defect']:.6g}|
|直接PEC边界与完整级联最大差|{a['tests']['maximum_direct_pec_cascade_difference']:.6g}|

每档带宽9点，噪声每实/虚正交分量标准差0.01。白化均值距离为||f₁−f₂||/0.01，**不除以通道数的平方根**。

|相对全带宽|墙校准距离|指定联合补偿对距离|
|---|---:|---:|
{bands}

0.5%下校准0.88714仍较难分，而同一联合补偿对已达到3.03825。目标存在的完整数据确实能击破这组特定补偿。这里的位移1/12λ₀也小于P0-B另行声明的0.1λ₀任务阈值，一维回归不替代二维joint门槛。不同带宽重新布9点，不是嵌套采集，不能要求这些距离严格单调。

![一维带宽回归](figures/01_1d_bandwidth.png)

图1只来自独立解析数组，不是有限横向墙体证据。

## 4. 二维模型、开放空间和冻结采集

复用A2的二维TM出射Green VIE数学，适配矩形域和单侧天线。时间为e⁻ⁱωt，g=iH₀⁽¹⁾/4，电流满足j=χ(Einc+Gj)。直接复数LU求解完整耦合，没有截断多次散射。χ=εr−1，本轮二维直接χ本构不是三维CM/RR极化率。

域为[-2,2]×[0,3]，天线y=−0.3，4发12收，孔径3λ₀；宽带0.8/0.9/1/1.1/1.2，窄带0.9975/0.99875/1/1.00125/1.0025。两带共9个不同频率。频率之间共用ε∞与σ₀，εr(f)=ε∞+iσ₀/f；源谱和几何固定，不能每频率拟合自由增益。

噪声协方差为0.0001I，单位线源强度固定1；没有依据墙体结果重新缩放噪声。每带240个复观测，480实维。噪声样本的总范数通常远大于1，并不意味着两个拟合世界就相距≤1。

查询T=[−1.2,1.2]×[1.6,2.8]只是计算分块；本轮calibration fixture明确将其实际置为空气。固定12位置、半径0.12、对比3/7的介质目标类。joint任务原定为T内总正介电对比的质心差≥0.1λ₀，并约束质量相对差≤20%。没有执行仪器nuisance压力测试或未知目标的joint补偿搜索。

**开放边界的实现说明：** VIE直接使用无限均匀外部的Sommerfeld出射Green，不存在PML厚度或人工外边界参数。因此没有假称完成PML增厚实验；使用独立圆柱散射解和网格细化审核这一表示。其辐射条件成立不代表内部空间离散足够准确。

![真实墙体材料](figures/02_geometry.png)

图2为固定seed的实际材料。两者有倾斜外形、变厚度、非均匀介电率、合法损耗和不规则空气区。按Re ε>1.1、前景8/背景4连接作离线描述，第一个空气区与外界连通、封闭hole数0，第二个为封闭hole数1；不把两者都宣称为封闭空腔，也不在看到结果后更换seed。该拓扑统计不是在线输入或门槛。

![实际复场](figures/03_first_source_field.png)

图3是f/f₀=1、第一个源的保存离散总场幅度/相位，分别为无目标和同一固定目标状态。它不是重建图，不表示已经知道墙后目标。

## 5. 二维资格为何没有通过

先做公开benchmark，再保留原8000 RHS/2000分解/7200 CPU秒预算。使用32、64、96、128、192的预先登记方格分辨率阶梯；细网格不是把32×48逆表示误当足够准确的正演。

常规离散一致性结果良好：四步长材料差分的最佳相对误差约{min(r['relative_error'] for r in dd['finite_differences']):.3g}，实伴随dot缺陷约{dd['real_adjoint_relative_error']:.3g}，线性后向残差约10⁻¹⁵。20项代数/合同测试通过。这些只能说明对当前离散矩阵求导和求解基本一致。

独立均匀介质圆柱Mie参考揭示了尚未解决的空间误差：原实现从nx32到128的相对散射误差约84.1%、22.8%、9.91%、7.03%；匹配源平均后nx128仍约7.00%。误差随细化下降，但没有达到足以支撑本题强均值距离门槛的精度。

对两个真实不规则墙体，修订源接口后192×144相对128×96的**完整数据向量**白化差如下：

|seed|宽带cal|宽带固定目标total|窄带cal|窄带固定目标total|
|---|---:|---:|---:|---:|
{rows}

所有值都超过0.1，差约为目标容差的9.5–12.3倍。没有连续PDE严格误差上界；相邻网格差本来就只是数值一致性证据，如今连这一要求都没有通过。不能以某个粗网格d=0.999宣布GO，也不能把网格相近直接称连续Maxwell认证。

![独立参考与网格收敛](figures/04_convergence.png)

源接口和内存修复均保留原记录：接收单元积分与入射中心值不匹配导致交换源/接收的相对差3.49×10⁻⁷；匹配入射单元平均后降至7.20×10⁻¹⁶。但这不改善主要网格误差，故不能把互易性修复当成资格通过。原wall job峰值RSS2.05GB；分块、原位LU后降至0.572GB，同材料无损工程改动的数据对比完全一致。未使用jitter、伪逆或降低精度换取内存。

现有迹象指向空间系数/场近似精度，而非LU未收敛或伴随符号错误；尚未独立量化几何面积平均、常数单元场近似各占多少，不能给出未经测量的误差百分比分解。

## 6. 两个固定参考世界的数据：不是成功见证

以下全部由已保存的修订版最细网格数组重算，没有新增求解。两个seed并非优化得到的数据近等价端点。同一固定目标位置在两世界相同，任务质心差为0。

|带宽|两墙cal距离|同目标完整total距离|目标诱导响应差|目标质心差|
|---|---:|---:|---:|---:|
{pairrows}

目标诱导差是[(F₁(X)−Fcal₁)−(F₂(X)−Fcal₂)]的实验机制诊断；calibration场本来就是额外公开fixture，不能拿这一相减充当joint在线输入。原定义dT使用完整目标响应，该表同时给total与induced，避免混用。虽然诱导差达到若干噪声单位，cal距离约36–38，完全不满足近等价条件；更没有非零目标任务差。且网格仍未qualified，不能拿表中数值作为认证物理下界。

![固定参考对的带宽比较](figures/05_reference_pair_bandwidth.png)

## 7. 四块传播：同一离散系统中的有限验证

在固定两个seed、f/f₀=0.8和1.2、同一目标位置两种固定对比上，额外完成8个实际离散全波对比。四块只覆盖该目标支撑的21个单元端口，**没有计算整个T的完整Green签名**。

用环境LU求任务端口后，四块重建与直接含目标全波数据的最大相对差为{four_error:.3g}。Born目标诱导响应相对误差为{born_min*100:.1f}%–{born_max*100:.1f}%。这验证了非交换公式和内部多次反馈在该离散实例中不能省略；不是未知墙体歧义、ROM优势或连续PDE精度证明。

核心最小奇异值在{beta_min:.3f}–{beta_max:.3f}。没有找到接近奇异的冻结实例，不额外调材料或频率制造“近共振”图。每个环境四块setup需要21个复端口RHS，而直接一个新目标前向为4个源RHS；再加小核因子与4 RHS。不能把21维小核、缓存后的局部计时或恒等式精确当作总成本下降。

![四块离散核验](figures/07_fourblock_diagnostic.png)

## 8. 费用、失败与没有运行的内容

截至本报告生成的计费快照：实际完整离散状态{cost['native_state_evaluations']}次；物理因子分解{cost['factorizations']}次，其中含8个实际目标小核分解；complex RHS为{cost['rhs']}。源前向{cost['by_rhs_role']['forward_sources']} RHS、接收转置{cost['by_rhs_role']['receiver_transpose']} RHS、四块端口{cost['by_rhs_role']['fourblock_port_rhs']} RHS、小核{cost['by_rhs_role']['fourblock_target_core_rhs']} RHS。全部在执行前预约；同材料、不同实现的重复和修复同样计费。

validation使用{cost['by_phase']['validation']['rhs']}/1600 RHS；repair使用{cost['by_phase']['repair']['rhs']}/800 RHS；search和endpoint搜索预算都为0。全局没有达到预算上限。实际包装作业wall合计约{cost['job_wall_seconds']:.1f}秒、process CPU约{cost['job_process_cpu_seconds']:.1f}秒；对包装前的解压、源码/文献审读及环境尝试另保守计费600秒，**不是测量值**。账本收费合计约{cost['charged_cpu_seconds']:.1f}秒。报告、图表与发布尾项仍追加到主账本，最终数值以提交中的账本为准；模型推理不当成本地物理CPU用量。

保留的失败包括系统Python缺NumPy、pytest未安装导致一次测试启动失败、初版源/接收离散不对称、较高内存峰值及两seed的资格失败。使用已有科学环境和标准库unittest后20项测试通过；没有修改冻结的旧项目环境或删除失败记录。

![部分计时与误差](figures/06_partial_cost.png)

图7只显示保存的状态计算部分计时；没有完整推理、ROM加速或time-to-image结论。全部峰值、CPU/wall、LU/RHS及日志见 [物理账本](results/PHYSICS_LEDGER.jsonl) 和 [失败记录](results/FAILURE_LEDGER.jsonl)。公开物理调用口径不把人工矩阵代数单元测试伪造成材料状态；其CPU费用仍计入。

## 9. 七个科学问题逐项回答

1. **一维是否复现？** 是。相位、共同参考面、完整PEC边界/级联、幺正性、白化及0.5%带宽负结果均通过独立核验。
2. **有限不规则二维墙是否有合格calibration见证？** 尚未回答。只完成实际参考场和资格计算；Gate1失败，未搜索，不能报告“找不到所以不存在”。
3. **目标未知的joint歧义是否持续？** 未检验。没有把calibration或同一固定目标的比较冒充未知目标补偿。
4. **什么采集消除歧义？** 一维指定补偿对会被0.5%带宽击破。二维没有合格的歧义对，不能判断其宽带/孔径/目标自照明的消歧作用；没有运行孔径变化。
5. **task-aware搜索同成本更好吗？** 未建立。四种方法均为NOT_RUN，不存在成功率、同成本Pareto或算法优越性结论。
6. **四块/A17比单程或Born揭示什么？** 在同一离散系统里四块恢复完整内部反馈，Born在固定强目标上误差很大；但空间精度不足，量化只属离散诊断。没有执行A17新ROM，Schur同空间等价不构成新算法或加速证明。
7. **P1、测量设计还是停止？** P1-A/P1-B均HOLD，先解决二维前向资格。现在也不能宣布REDIRECT_TO_ACQUISITION_IDENTIFIABILITY，因为尚未有被宽带joint击破的合格二维见证。停止本轮搜索扩展，不停止尚未被检验的物理问题。

## 10. 下一轮可以做的最小工作

优先以这两个已冻结实例做更准确的前向离散核验，比较高阶/界面相容方法或独立FEM/BEM参考与当前VIE；保持同一源、噪声、世界及0.1门槛。按目前近似二阶下降趋势，简单加密可能需要数倍线性分辨率，直接活跃稠密LU的内存与费用会快速增加。这是工程预测，不是严格误差外推，也没有在本轮付费展开。

只有前向资格通过后，才按冻结数据和同成本普通强对照恢复P0-B搜索。若随后发现仅窄带calibration成功、宽带joint失败，再转向采集可辨识性；若joint见证成立，先做可独立积分的低维P1-A，并以同prior/likelihood direct joint方法作参照。P1-B必须另做单点/多解×joint/progressive的2×2消融。本轮两者、3D、大网络和概率权重均没有启动。

## 11. 阅读与复现

- 入口：[START_HERE.md](START_HERE.md)；冻结配置：[w_config.json](w_config.json)。
- 原始指标：[RESULTS.csv](RESULTS.csv)；机器判决：[results/GATE_DECISION.json](results/GATE_DECISION.json)。
- 实施细节：[PHYSICS_VALIDATION.md](PHYSICS_VALIDATION.md)、[PHYSICS_PROTOCOL.md](PHYSICS_PROTOCOL.md)。
- 代码：[src/research_w](src/research_w)；实际命令与CPU凭据：[results/jobs](results/jobs)。
- 全部图只读实际缓存；PNG/SVG及作图CSV位于 [figures](figures)。

最接近文献边界：单侧通过未知墙的应用可参考[Watson 2024](https://doi.org/10.1088/1361-6420/ad5b83)；反射/透射多模歧义已有[Waagaard–Skaar 2007](https://doi.org/10.1137/060672224)；部分边界数据的连续唯一性条件见[Brown–Marletta–Reyes 2016](https://doi.org/10.1016/j.jde.2016.01.002)。本轮仅用ScholarQA核对5条文献身份/可用摘要，未做穷尽新颖性检索；摘要缺口和与A17/B1/C1等既有代数的重叠保留在资料审计中。不宣称二维一般Maxwell唯一性、后验校准、首次反馈公式或IEEE TAP论文级算法优势。
'''
    (root/'RESEARCH_REPORT_ZH.md').write_text(report)
    start='''# Research W P0：执行入口

**结论：P0-A PASS；P0-B UNRESOLVED / UNQUALIFIED_SOLVER。**

一维完整Maxwell反例及其带宽负结果正确复现。二维已经运行两个有限、不规则、非均匀、有损墙体的真实多频全波状态，但最细网格白化差0.951–1.229，大于冻结要求0.1。故停止四种见证搜索；calibration/joint没有合格反例，P1没有授权执行。这不是歧义不存在的证明。

- 完整中文解释：[RESEARCH_REPORT_ZH.md](RESEARCH_REPORT_ZH.md)。
- 状态与原因：[GATE_DECISION.md](GATE_DECISION.md)。
- 最重要图片：[真实几何](figures/02_geometry.png)、[复场](figures/03_first_source_field.png)、[未通过的收敛](figures/04_convergence.png)、[四块有限诊断](figures/07_fourblock_diagnostic.png)。
- 实验参数：[w_config.json](w_config.json)；[逐项指标](RESULTS.csv)；[物理账本](results/PHYSICS_LEDGER.jsonl)；[失败](results/FAILURE_LEDGER.jsonl)。
- [复现命令](REPRODUCE.md)；[原ZIP文件清单](docs/SOURCE_INVENTORY.json)；[原始理论PDF](source_package/Research_W_Report_ZH.pdf)。

物理范围仅二维TM、已声明噪声与材料类。没有posterior概率、神经网络、三维或加速结论。没有哈希/备份流程。
'''
    (root/'START_HERE.md').write_text(start);(root/'README.md').write_text(start)
    gate='''# Research W P0 — Gate decision

|Gate|状态|理由|
|---|---|---|
|0 解析复现|PASS|完整一维参考面、PEC级联、幺正性及带宽负结果独立复现|
|1 二维资格|UNQUALIFIED_SOLVER|128→192修订版多频完整数据白化差0.951–1.229，未达0.1|
|2 Calibration机制|NOT_RUN_UNQUALIFIED_SOLVER|未启动候选搜索，没有合格有限墙体见证|
|3 Target-present joint|NOT_RUN_UNQUALIFIED_SOLVER|没有联合目标/仪器补偿实验|
|4 带宽/微扰鲁棒性|UNRESOLVED_2D|只有一维指定补偿负结果，不能外推二维；未做轻微扰动压力测试|
|5 算法/P1价值|HOLD|未有qualified joint见证；P1没有执行授权|

总体：**P0_B_UNRESOLVED_NUMERICAL_QUALIFICATION**。没有耗尽总预算；是资格门槛阻止搜索。不能把20项测试通过、258个离散状态或四块恒等式误写成物理见证成功。物理存在性与搜索算法价值均没有被否定或确认。

精确数值/统计位置：[results/GATE_DECISION.json](results/GATE_DECISION.json)、[RESULTS.csv](RESULTS.csv)。
'''
    (root/'GATE_DECISION.md').write_text(gate)
    validation=f'''# Physics validation

## 证据分层

1. 标准数学：纤维、四块resolvent、条件性TV界经审查，不作新颖性证明。
2. 一维独立Maxwell：行波矩阵、Fresnel、PEC边界与完整级联核对；[原始结果](results/p0a/results.json)。Gate0 PASS。
3. 二维实际离散全波：两seed、九个不同频率、4源12收、cal及固定目标，完整LU多次散射。[原版](results/qualification/WALL_GRIDS.json)、[匹配源平均修订](results/qualification_source_average/WALL_GRIDS.json)。
4. 求导/伴随：四步长最佳相对差{min(r['relative_error'] for r in dd['finite_differences']):.6g}；real-adjoint缺陷{dd['real_adjoint_relative_error']:.6g}。这只验证离散接口。
5. 互易性：源中心/接收积分原版3.48829e−7；匹配平均后7.20389e−16。[实际两种实现对比](results/qualification/RECIPROCITY_REPAIR.json)。
6. 独立圆柱解析解：修订版nx128仍7.0036%相对散射误差；非充分精度。[数据](results/qualification_source_average/CYLINDER.json)。
7. 独立网格：最细白化差{g['grid_whitened_min']:.6f}–{g['grid_whitened_max']:.6f}>0.1。Gate1 UNQUALIFIED。
8. 原位LU内存修复：同材料数据差0；最大RSS从2.0545GB降到0.5723GB。[parity](results/qualification/ASSEMBLY_PARITY.json)。
9. 8个任务受限四块离散核验：相对误差最大{four_error:.6g}，Born诱导误差{born_min*100:.2f}%–{born_max*100:.2f}%。没有执行全T签名或证明连续精度。[原始数组](results/four_blocks/ARRAYS.npz)。

## 开放边界与范围

出射VIE的外部是无限均匀介质，采用解析Sommerfeld辐射条件；没有人工PML/盒边界，故PML增厚为NOT_APPLICABLE，而非PASS。网格增密检查的是内部系数和场离散。没有连续PDE误差证书；没有把经验网格差当上界。

均匀真空χ=0时，前向直接为解析自由空间线源，没有数值LU；零对比、打包/白化、跨频材料共享、真实离散伴随及完整多次反馈均有对应代码/日志。没有额外收发孔径、增益、频率、目标位移或隐含目标选择。

独立平板检查由P0-A精确Fresnel提供，不冒充有限二维板边缘参考；二维额外采用均匀圆柱Mie作不同物理模型检验。通过同一离散系统的多条解法不能替代这个外部参考。

## 停止

本轮没有任何候选获GO，没有最终见证细网格审核、轻微扰动、未知目标joint补偿或posterior实验。独立细化失败导致这些环节NOT_RUN。20项标准库测试通过，保留一次pytest依赖启动失败。后续必须先获得足够精确的前向平台，不能靠提高噪声或放宽0.1条件继续。
'''
    (root/'PHYSICS_VALIDATION.md').write_text(validation)
    p1='''# P0 → P1 decision

**两条P1路线均HOLD；没有自动执行P1。** 目前失败是二维空间精度，不是被验证的传播唯一性或未知墙体物理NO-GO。

下一步只应先解决原两个世界、原采集下的前向资格。不能因为一维单频见证成立就启动概率系统，也不能因为未运行二维搜索就声称宽带已经消歧。

## P1-A：概率传播等价

只有先得到qualified joint端点、并理解其窄/宽带及扰动适用范围，才考虑一个可独立积分的低维联合prior参考。保留目标/墙/仪器依赖，明确事件，比较direct joint SMC/MH或参考积分。任意candidate权重不是posterior。事件/CDF≤0.02和同准确度费用至少减少20%仍只是未来审核目标，本轮没有取得这些结果。

## P1-B：任务导向递进

在同数据、同prior、同费用下比较joint single-point、progressive single-point、direct joint multi-hypothesis、progressive multi-hypothesis。若只有多解增量没有递进增量，删除progressive独立优势。当前没有资格支持启动这项成像实验。

## 理论与工程预测

四块形式确实比单程/Born保留更多目标反馈，但它不能修复错误或不够精确的前向离散。其完整T端口成本可能高于直接多RHS。优先比较高阶/界面相容的前向表示与本轮真实保存数组；通过资格后再恢复普通Full-J、随机、multistart和task-aware同成本搜索。

若之后只有窄带calibration歧义，宽带joint将其消除，应REDIRECT_TO_ACQUISITION_IDENTIFIABILITY；若存在qualified宽带joint歧义，才值得先做P1-A。现在两种分支都尚未被二维证据支持。
'''
    (root/'P0_TO_P1_DECISION.md').write_text(p1)
    theory=(root/'docs/THEORY_AUDIT_DRAFT.md').read_text().replace('# Research W 数学范围审查草稿','# Research W 理论独立审计')
    theory=theory.replace('本文供主 Codex 审查使用，不作项目最终科学结论或新颖性裁决。','本章已由主线程复核；最终科学状态见 GATE_DECISION.md。数学条件正确不代表二维假设成立或具有论文原创性。')
    theory=theory.replace('没有运行新的物理求解、数值试验、材料搜索或文献检索；没有修改源报告。','本章数学审阅不额外运行物理解；项目实际执行情况另见 PHYSICS_VALIDATION.md。源报告未改写。')
    theory=theory.replace('本次工作核对记录及解析预期，没有再次执行该计算。','本项目P0-A已实际独立执行；本章审阅没有额外重复求解。')
    theory=theory.replace('已有 CPU 回归记录，本次未重跑','本轮中断前实际执行的独立CPU回归记录')
    (root/'THEORY_AUDIT.md').write_text(theory)
    (root/'requirements.txt').write_text('numpy>=2.0\nscipy>=1.14\nmatplotlib>=3.9\n')
    import scipy,matplotlib
    (root/'ENVIRONMENT.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,device='Apple Silicon CPU; BLAS threads=1',existing_runtime='A25 scientific venv; no installs',scope='2D TM only'),indent=2))
    (root/'REPRODUCE.md').write_text('''# 复现

需要Python3.11+、NumPy/SciPy/Matplotlib，无需pytest、CUDA、后台服务。项目实际使用已有A25环境；新环境可安装requirements.txt。所有命令应从本目录运行。设置 `PYTHONPATH=src`；示例的 `python` 指你的科学Python。

运行前阅读w_config.json和PHYSICS_PROTOCOL.md。公开结果保留原样。缓存重画需要从[本轮Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-p0-20261022)下载Research_W_P0_ARRAYS.zip，按内部相对路径解压至本项目根目录。较大的墙体NPZ和逐网格field CSV放Release，其他摘要和四块小数组在Git。

若要独立从头重算，用 `python tools/new_run.py /absolute/new/empty/directory` 生成不含旧结果的新实验目录；它只复制协议/源码/只读输入，不启动求解，不是重置原campaign的账本。之后进入新目录执行下面命令。原公开项目及其负结果保持不变。

```sh
python tools/run_job.py local-analytic-001 python -m research_w.analytic --output results/p0a --source-results source_package/analytic_results/results.json
python tools/run_job.py local-tests-001 python tools/run_tests.py
python tools/run_job.py local-benchmark-001 python -m research_w.qualification benchmark
python tools/run_job.py local-cylinder-001 python -m research_w.qualification cylinder
python tools/run_job.py local-derivatives-001 python -m research_w.qualification derivatives
python tools/run_job.py local-wall-grids-001 python -m research_w.qualification walls
python tools/run_job.py local-parity-001 python -m research_w.qualification assembly-parity
python tools/run_job.py local-reciprocity-001 python -m research_w.qualification reciprocity
python tools/run_job.py local-wall-average-001 python -m research_w.qualification walls-averaged
python tools/run_job.py local-cylinder-average-001 python -m research_w.qualification cylinder-averaged
python tools/run_job.py local-derivatives-average-001 python -m research_w.qualification derivatives-averaged
python tools/run_job.py local-fourblocks-001 python -m research_w.four_blocks
python tools/run_job.py local-figures-001 python tools/plot_results.py
python tools/run_job.py local-summary-001 python tools/summarize.py
python tools/run_job.py local-report-001 python tools/report.py
```

仅重画已公开数组和指标不需要Maxwell：运行plot_results.py、summarize.py、report.py即可。公开 `results/jobs/*/RECEIPT.json` 保存真实执行命令及其scope，包含失败。

P0-A的独立入口：`python -m research_w.analytic --help`；上述命令独立保存在新results/p0a并比较原归档。原ZIP脚本若直接运行会写到它的旁边；只允许在自己的可写副本运行，不能写只读source_package。

本轮没有实现/执行完整四种见证搜索、未知目标joint正确器、P1或3D。它们不是隐藏的未公开成功代码。资格未通过时不能绕过gate去调用未来搜索器。
''')
    body=''.join(f'<figure><img loading="lazy" src="figures/{n}.png"><figcaption>{html.escape(cap)}</figcaption></figure>' for n,cap in [
        ('02_geometry','固定两世界的真实介电率、损耗和空气区域。'),('03_first_source_field','实际第一个源的完整离散总场，非重建或手绘。'),
        ('04_convergence','最细网格白化变化仍大于0.1：资格未通过。'),('01_1d_bandwidth','一维解析回归：完整目标数据会破坏指定补偿。'),
        ('05_reference_pair_bandwidth','两个固定参考seed不是搜索见证，目标质心差为0。'),('07_fourblock_diagnostic','同一离散系统的有限四块恒等式，不是连续见证或加速。')])
    (root/'index.html').write_text('''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research W P0 — 物理核验与未解决精度</title><style>body{max-width:1100px;margin:40px auto;padding:0 24px;font:18px/1.75 system-ui;color:#172631;background:#fafcfc}h1{font-size:2em}aside{background:#fff0d8;border-left:6px solid #b56c18;padding:18px 24px}img{max-width:100%;height:auto;background:white}figure{margin:40px 0}figcaption{color:#4d616c}a{color:#00677d}code{background:#e9f0f1;padding:2px 6px}nav{display:flex;gap:24px;flex-wrap:wrap}</style><h1>未知不规则墙体的 Maxwell 传播歧义</h1><p>Research W · P0 · 2026-10-22</p><aside><b>一维 PASS；二维 UNRESOLVED / UNQUALIFIED_SOLVER。</b><br>二维真实全波已运行，但128→192网格的多频白化差0.951–1.229高于0.1。停止见证搜索与P1；没有宣布歧义存在或不存在。</aside><p>重要区别：离散方程解得准确、四块公式正确，并不能证明空间离散已经足够精确，更不能证明算法优越。</p><nav><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_w_p0_20261022/RESEARCH_REPORT_ZH.md">完整中文报告</a><a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_p0_20261022/RESEARCH_REPORT_ZH.md">匿名raw</a><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_w_p0_20261022">代码与原始证据</a></nav>'''+body+'<p>没有posterior概率、NN、3D或同成本算法优势。全部图片由已保存实验数组生成。</p></html>')
    print('Reports written; no native physical calls')

if __name__=='__main__':main()
