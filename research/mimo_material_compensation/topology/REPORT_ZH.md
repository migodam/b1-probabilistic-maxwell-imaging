# 固定端口拓扑竞争与新读数复拟合

在独立连续背景 fresh_bg_b 上，交错 B 与补偿引导布局追加新数据后，仍各有一对严格相容材料在 0.1484 平方米的区域内产生至少 0.5 的均值差，外区面积加权 L2 变化严格低于预登记的 0.35。常量区域编辑和形状自由均值任务都找到这种见证。两布局首轮少迭代时连 unrestricted 材料也未找到，同数据加深求解后找回了歧义。七对历史材料的 6,435 个切分穷举用于固定对排名，历史对象的全类分离下界仍均为 0。

## 已有材料见证与端口

七个独立对象是 bmtm26_n32_c01/c02/c03/c04/c05/c06/c08_f400_a08_c48。使用发布的 A 噪声数据、两份封存材料和各自的预定 ROI。七对的原 A 联合球残差数值为 0.44–0.81r 一带。对同一离散模型的严格复审中，全部 14 个端点的最大上界为 0.8364r，ROI 均值差为 2.177–2.602。每份材料都符合 [0,4.5]。这为相同对象的全 6,435 个布局各给出一个关于全类最小分离度的安全上界。已证下界都是 0，幅度阈值取 0.5，未剪去任何候选切分。全候选最大上界减所选布局下界的差为 1.552–10.732，单位为 (D/(2r))。[严格回执](fixed_pair/verified_upper_receipts.json) · [逐切分上下界](fixed_pair/class_bounds.npz)

| 原生对象 | A 固定对 (D/(2r)) | B 固定对 | C 固定对 | 最好的八发八收固定对布局 |
|---|---:|---:|---:|---:|
| n32_c01 | 0.350 | 8.024 | 6.122 | 8.889 |
| n32_c02 | 0.366 | 1.419 | 1.018 | 1.545 |
| n32_c03 | 0.436 | 9.139 | 6.897 | 10.372 |
| n32_c04 | 0.434 | 8.894 | 6.527 | 9.743 |
| n32_c05 | 0.424 | 7.803 | 5.890 | 9.613 |
| n32_c06 | 0.377 | 8.014 | 5.870 | 8.940 |
| n32_c08 | 0.357 | 9.115 | 6.957 | 10.723 |

B 的固定对分离能量有 96.44%–99.93% 来自 A 未读取的弧内边。B 保留 A 的 32 条跨弧互易边、移除另外 32 条并补入 32 条弧内边。A 的角色全反转保留全部 64 条无向边。A+B 在 128 个复读数中有 96 个可独立的复坐标，C 对 256 个复读数有 136 个可独立的复坐标。对应实坐标上限分别为 192 和 272。C 的 (D^2) 在七例均大于 B。C 的 16 发各用 1/16 W，A、B 的 8 发各用 1/8 W，总可用源功率之和同为 1 W。C 的联合 Frobenius 半径是 (2r_B)，归一化比在七例均低于 B。这些数值条件与脉冲时长未指定的离散数值仿真相关。[逐边归因和距离](fixed_pair/summary.json) · [图一](figures/fixed_witness_competition.png)

固定仪器为 24 个物理端口，编号 16–23 的八个端口始终在 50 Ω 匹配端接。总散射使用 \(S=I-2Z_{\mathrm{ref}}(Z_A+Z_{\mathrm{ref}}I)^{-1}\)，测量使用 \(S-S_{\mathrm{empty}}\)。单次提取全部候选发射列使用一个物理分解，同一材料的角色切换由投影实现。标定散射矩阵的转置互易误差为 \(3.73\times10^{-15}\)，独立裸 VIE 加加载电路与名义正演的差为 \(9.32\times10^{-15}\)，功率平衡偏差小于 \(10^{-15}\) W。若各方向的复噪声独立同方差且为圆对称 Gaussian，取 \((C_{rt}+C_{tr})/\sqrt{2}\) 得到 B 的均值和噪声方差。4096 次模拟所得方差 0.995（理论 1），该构造在本研究的确定性联合球内只保证投影噪声不超过 (2r_B)。实化伴随内积审计及三个步长的 JVP 差分相对误差另存。[端口、负载及复核](fixed_pair/port_fairness.json) · [导数与相干发射审计](fixed_pair/derivative_and_coherent_receipt.json)

## 逐对象留一的固定对排序

每轮训练以另外六个材料为一个独立集合，并在剩下的整件材料上评分。严格固定对最高值只有事后参考。逐对象 Spearman 排序在 c02 对中心 Jacobian 高达 0.998，其余六件多数接近 0，负值也存在。混合场 (C (I-X_1 G)^{-1}(X_1-X_0)(I-GX_0)^{-1}B) 在每个切分上精确回放有限材料差，使用的两端点已经固定，给出的不是全类下界。中心 Born 指空场对差分，混合近似使用一侧空场与一侧中心背景，中心割线是三次状态算出的代数精确对照。[代理得分](fixed_pair/cut_scores.npz) · [图二](figures/surrogates_and_power.png)

| 方法 | 留出中位 (D/(2r)) | 相对于每件留出材料固定对最高值的中位比例 |
|---|---:|---:|
| A | 0.377 | 0.042 |
| B | 8.024 | 0.896 |
| 平衡随机 | 8.263 | 0.804 |
| 空间频率几何 | 8.096 | 0.842 |
| Born 训练 maximin | 8.604 | 0.917 |
| 中心 Jacobian 训练 maximin | 8.557 | 0.891 |
| 局部 Fisher | 6.309 | 0.656 |
| 旧数据约束外区后的 Fisher | 8.389 | 0.870 |
| 精确 Gaussian 条件 Fisher，四个 ROI 模态 | 8.102 | 0.874 |
| 信号能量 | 7.931 | 0.871 |
| 精确有限对训练中位秩 | 8.316 | 0.896 |
| 精确有限对训练 maximin | 8.557 | 0.891 |
| 留出材料事后固定对最高值 | 9.613 | 1.000 |

相干 oracle 对各布局的已知固定差矩阵取首右奇异向量，同时驱动八个端口。运算使用一个 1 W 发射脉冲、8 个复读数。用于定理的距离在原 64 读数半径下为 3.413–20.042，按 8 读数球则为 9.652–56.687。此上限使用事后对材料、相干源控制和八倍单脉冲可用功率，不能归入开关布局的同硬件逐列实测账。[逐材料得分与相干脉冲](fixed_pair/summary.json) · [留出表](fixed_pair/heldout.csv)

## 主动试验合同

先登记四个新的位置任务和阈值，再只凭带噪 A 得到不同形状、均值、外区变化的候选。从 A 复读、整体反转、B、随机、Born 空间频率、Fisher、考虑旧数据外区自由度的 Fisher、信号能量、稳健 T 判别与补偿引导中在各自规则下选布局。T 判别最大化已生成 A 相容候选对预测距离的最小值。补偿引导在 16 个领先切分加上基线布局中进行候选级比较，允许高端点外区在原宽盒和单独登记的面积加权 L2 半径 0.35 内变化，以旧和新布局的联合线性残差为排序。用于布局排序的是近似优化求得的分数，计算带候选与局部近似性质。同步保存方向型线性下界供复查，全类下界仍为 0。c02 的本轮预登记 ROI 与历史 ROI 不同，冻结前未找到两类相容端点，记录有明确降级为不相容探针的回执。[冻结规则](FROZEN.md) · [对象与活动预算](campaign.json)

每臂费用独立记账，原 A 加至多一轮追加。每轮八个等时脉冲，每发 1/8 W，读 8 路，增加 64 个复读数、一个单位脉冲时长乘 W 的可用源能量。实际阵列切换和接受功率根据驱动端与未驱动端物理角色列出。新读数由独立进程在提案和预测阵列写入后用裸 VIE 加目标加载网络模拟。原与新数据各有一个确定性复 Frobenius 球，反演要求两个残差分别过门，绝对半径在同对象各臂一致，模拟追加读数的噪声方向长度为 (2r/3)。每臂分别用原 A 的数据和该臂自己的新数据重新约束 ROI 零值、指定连通片下限 2、无形状的 ROI 均值 ≤0.5 或 ≥1.5，以及不限制 ROI 的材料。外区从宽盒 [0,4.5] 独立调节。均值跨度达 0.5 的数值解加入若干中间值复拟合，延拓点仍是有限组数值见证。[数据隔离采集入口](simulator.py) · [全波复拟合入口](loop.py)

## 主动试验结果

用每臂自己的 A 与追加读数联合复拟合。在四件材料上，重复 A 与 A 整体反转都保留宽的有限面积均值歧义。表中每个宽度来自两个通过 Arb 独立残差审计的材料。旋转布局列给出不限制 ROI 时两个观测球的最大残差倍数，数值低于 1 才有严格审核的资格。

| 材料 | A 复读严格均值宽度 | A 反向严格均值宽度 | B 无限制最大残差 | T 判别 | 补偿引导 | 同预算旋转布局最小 unrestricted 残差 |
|---|---:|---:|---:|---:|---:|---:|
| n32_c02 | 1.106 | 1.118 | 10.442 | 10.430 | 10.785 | 8.106（random） |
| n32_c07 | 1.616 | 1.587 | 6.000 | 5.230 | 6.235 | 4.844（Fisher） |
| fresh_bg_a | 2.282 | 2.282 | 4.075 | 3.877 | 3.655 | 2.750（random） |
| fresh_bg_b | 1.500 | 1.500 | 3.820 | 3.363 | 2.630 | 2.630（compensation） |

32 个逐材料旋转布局中，0 个在首轮同预算五个材料家族攻击中找到数值相容端点。独立离线复算封存生成材料，四件的旧 A 和所有新布局的严格残差上界均不超过 0.677r。同一数据有物理解。旋转布局下较大的已见残差反映求解支出和初始化，首轮结果不能用于全类排除。全部外宽仍为 4.5。

四件材料的标准反演合计 505.6 秒 CPU，A 数据下的设计合计 119.0 秒 CPU。所有历史、动作及逐臂外区变化、真值材料误差和八个未读取终端的仅离线复算比分开保存。

[分析数据](active_compact.json) · [逐布局及家族表](active_summary.csv) · [首轮残差及见证宽度图](figures/closed_loop_residuals_and_spans.png)

## 对同一读数加深求解后，B 和补偿引导仍然留有歧义

在独立连续背景 fresh_bg_b 上，对原 B 和补偿引导开展同数据加强试验。只读每布局自己的 A 与新 Y。补测为 0。先按每组接收端口独立改写 Dual-CIE 的电流辅助空间，保留发布的三个 Fourier 阶段、β、像素读出及站稳停止规则，再从辅助结果和原全波候选作更长的物理 GN。两次 Dual-CIE 加 GN 约各 51 秒 CPU。严格审计确认 B 和补偿引导的零 ROI、连通斑块、低均值、高均值及 unrestricted 五个家族各有相容端点。

| 数据合同 | 零区材料残差上界，旧 / 新 | 高均值材料残差上界，旧 / 新 | 最大已验证 ROI 均值跨度 | 该对外区面积加权 L2 |
|---|---|---|---:|---:|
| 交错 B | 0.797 / 0.733 | 0.780 / 0.708 | 1.500 | 3.411 |
| 补偿引导 | 0.799 / 0.761 | 0.648 / 0.612 | 2.037 | 3.681 |

增加到 220 步并校正归一化权重的单纯物理 GN 在补偿引导臂也找回零区和高均值端点，严格宽度为 2.168。B 的同预算纯 GN 在 unrestricted 对照仍高于噪声半径。CIE 起点使 B 和补偿臂都达到双球内，首轮的全波失败属于数值求解失败。

以各臂重新求出的零区材料为中心，再用预登记的较紧外区预算 0.35 检验 38 格、0.1484375 平方米的 ROI。两种布局均找到 ROI 常量和形状自由两类备选。各备选的 ROI 均值相对各自中心严格升高至少 0.5000009999998。材料在 [0,4.5]，外区面积加权 L2 严格上界均小于 0.349999965001。以下列出目标均值处的双球残差上界。

| 数据合同 | 固定常量 ROI，旧 / 新 | 形状自由均值，旧 / 新 |
|---|---|---|
| 交错 B | 0.8146 / 0.7437 | 0.7093 / 0.6535 |
| 补偿引导 | 0.8295 / 0.7563 | 0.7987 / 0.7856 |

因此这两个布局对 fresh_bg_b 的同一预登记任务均未迫使所需外区变化超过 0.35。界限对各自从数据构造的零区参考点成立。严格回放验证了有限个材料端点和中间任务幅度，连续路径和最小外区代价仍待证明。其余对象及布局的宽类保持未决。

[双球和预算回执](loop/fresh_bg_b/budget_witness_verified.json) · [Dual-CIE 逐点严格回执](loop/fresh_bg_b/dual_cie_verified.json) · [纯 GN 严格回执](loop/fresh_bg_b/long_gn_verified.json) · [同色标材料图](figures/fresh_bg_b_budget_witnesses.png)

## 96² 原生与更细格点复算

冻结布局只使用七件 32² 历史见证。两件原生 96² 材料的 A、B 固定对分离及用 32² 数据训练的精确 maximin 布局如下。全部为固定对数值结果，复用已经独立在 96² 求出的旧材料。

| 物体 | A | B | 32² 训练 maximin | 事后 96² 固定对最高值 |
|---|---:|---:|---:|---:|
| n96_c01 | 0.182 | 8.557 | 9.719 | 10.269 |
| n96_c02 | 0.393 | 8.026 | 7.416 | 9.406 |

n96_c01 的两份材料还以分片常数细分到 192²，由重新构造的裸空间算子与目标加载电路独立复算。两份细格点响应相对于原 96² A 带噪观测的残差为 0.692r 和 0.728r。这一对材料的跨格差异已测；目标加载的全材料类模型误差上界未取得。此处未做 96² 主动复拟合。
[原生和更细网格数据与支出](native/report.json)

## 复现

从 Inv_SLAM 根目录进入。下列 Python 需要已发布的研究依赖、单 BLAS 线程以及原研究目录的冻结算子。旧结果只读，新文件全部写在本轮 topology 内。已完成的主动试验可直接复用冻结提案、读数和材料。重新生成采集需要在另一份实验目录重跑。c02 的提案使用最早一版完整维度求解和一个自由初值启动，其他提案使用精确零电流主块和额外一个自由初值。c02 的首个复拟合进程在产出逐臂回执前中止，开销另行保存，正式各布局攻击从冻结候选和读数重启。[中止进程开销](interrupted_receipt.json)

```bash
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/rank.py
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/verify.py fixed
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/audit_derivative.py
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/loop.py run
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/native.py
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/verify.py loop n32_c02 n32_c07 fresh_bg_a fresh_bg_b
for arm in B_interleaved compensation; do
  .artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/dual_cie.py "$arm"
  .artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/long_gn.py "$arm"
  .artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/budget_replay.py "$arm"
done
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/verify_supplement.py dual_cie
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/verify_supplement.py long_gn
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/replay_truth.py n32_c02 n32_c07 fresh_bg_a fresh_bg_b
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/budget_check.py
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/figures.py all
.artifacts/z_alg_diag_python/bin/python -B research/mimo_material_compensation/topology/summarize.py
```
