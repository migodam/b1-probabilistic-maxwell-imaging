# R3 独立机制与文献审查

2026-11-03。只读 R1/R2 理论、报告、直接邻近文献与 Atlas 相关限定快照；本任务没有跑反演、读取新隐藏 GT、修改历史文件或发布。主线程保留核心新推导、最终路线、主张及新颖性裁决。

## 判断

主线程拟定的 `m+iγ` 被动损耗 continuation 可以作为有限开发实验，**不能称作原数据的精确 continuation**。人工损耗改善 Maxwell 优化的明确先例已找到：[Yao 等2022](https://doi.org/10.1007/s00158-022-03389-5) §3 Eq18 后使用 Q=5→1000。其任务是光子设计；本轮的未知材料测量拟合不同，仍可能值得比较，但“人工损耗／continuation 首创”不能成立。

冻结 `γ=4,1,.25` 各250步，再 γ=0 原目标corrector1000步，并与额外1750步普通优化比较，能区分**这一损耗路径是否改变有限工作结果**。迭代数相同不等于 LU、RHS 或时间相同；两臂额外工作与失败都需计费。R2失败对象已曝光，所得结论应为 DEV 机制结果；未经新鲜封存测试不能升级为盲测新增覆盖。

## 1. 三条路线的实质与先例

**电流数据可行性与 strong-mode release。** R2 mf12 没有释放固定 DPCS，故不能排除强模态硬锁的影响；但 Atlas10 F-BI-021 已明写将 `V_l` 修正与 APCS 一起优化。该修正是合理的单因子诊断，没有新数学资格。接收数据一致也不意味着共享实材料本构可行。若使用有限罚或数据球，要分别报告 current 的接收残差、本构残差及 material 的原 Maxwell forward 残差。

[CC-CSI2017](https://doi.org/10.1109/TAP.2017.2673758) Eq11–12 已增加 measurement-domain 的交叉本构残差。其 `y−S[m(E+GJ)]` 使用背景传播，不能被称作完整物理 residual 的新证书。原文把大 condition number 与小状态误差可能有大测量误差相连，但“condition number 大”单独并不足以证明绝对算子放大大；R3不应照抄这种缺少尺度的推断。

**非线性 material coordinate／inverse contrast resolvent。** 成熟 CIE 已用 `R=βχ/(1+βχ)`；[Bevacqua–Isernia2021](https://doi.org/10.1109/OJAP.2021.3063248) 比较 CS-EB/NIE/Y0 与组合重写。仅把 m 改写成逆对比、bounded contrast 或一个 Möbius 坐标，没有建立新信息或新 basin 定理。必须保留原先验的拉回；对 inverse contrast 显式处理 m=0，禁止用人为小正数无声删掉合法材料。高对比时有界坐标的逆映射可能病态；重参数化收益必须与普通尺度预条件比较。

**人工损耗／耦合 continuation。** 损耗可以平滑共振，但原观测不随之改变时，stage objective 与真实 lossless 数据生成不一致。大 γ 可能先迫使 real material 补偿不可能的衰减，再在回到零损耗时落回旧 basin。该风险不会被一个低 stage loss 自动排除。[Shin–Cha2009](https://doi.org/10.1111/j.1365-246x.2009.04102.x) 与 [Górszczyk等2017](https://doi.org/10.1002/2016JB013891) 的 time damping/complex-frequency 路线使用时域／频带观测，不能给本轮有限单频数据添加缺失信息。

## 2. data-only transform 命题的审查重点

主线程拟导出 `F_γ=T_γ∘F` 当且仅当 `F_γ` 在每个 F-fiber 上恒定。该命题属于标准“函数经过商集合分解”的结构。建议将 T 定义在 `F(M)` 上：其余 data 空间的任意延拓不由该命题决定。fiber 上恒定只给集合意义存在，**不自动给连续、可测、稳定或可计算**。局部 Lipschitz/微分结论必须单列假设；若有 Lipschitz 常数 C，则沿合法可微材料方向应满足 `||DF_γ v||≤C||DF v||`，从而得到 kernel 必要条件。充分性及邻域正则性不可从该一阶必要条件反推。

用既有 Z-ALG 六实材料／四实观测同响应分支加均匀被动损耗，是比 R1 新源更直接的 altered-medium 检验。**关键数值边界：** 若原 pair 仅浮点近同，而原差δ0非零，新的分离δγ不能单凭比值很大否定任意不连续 exact T；它给出所需 Lipschitz 常数下界，或否定预先声明的稳定 C。若要断言集合意义不存在，必须明确继承哪个精确 finite-model fiber/root 论证，以及新分离怎样超过数值误差；不得把机器小残差升级为连续 PDE 反例。原分支预测与新增 lossy forward 都收费；无需再做已知 resolvent 恒等式的精度测试。

若新的 source/receiver/material prior 在检验时变化，fiber 定义已变，不能称同一数据合同的反例。损耗应按同一 time convention 保持被动，原 real m 与新增 iγ 的作用域、ROI/host 是否一并加损耗必须冻结。

## 3. R1/R2 对新主张的硬边界

R1张量同场条件 `ΔεE_s=0` 说明信息缺失可保持任意 receiver 不变；旧源重组不增加信息。该结论不能直接转为实 scalar-TM 曲条的不可辨识证明。R1六项本构激发谱不含 receiver，不是完整 inverse Jacobian 谱。M 的兼容 H(curl) 和 T 的统一 k<p／完整δ仍未满足，不应拿这两个缺口给新的 scalar continuation 设置一个与机制无关、必失败的先行门。

R2 global 与 native 都2/3；full profile1/3；mf12弯条仍失败并且材料更差。它们排除“这些已执行修正已建立新增覆盖”，没有排除所有 nonlinear search paths。端点未证明驻点，不可称 chaos 或已认证错误局部极小值。损耗或几何路径若有效，应归于探索路径；原数据和最终材料类相同仍不等于新增物理信息。

空间 transport 作为备选有明确 shape-evolution 邻域：[Ferrayé等2003](https://doi.org/10.1163/156939303322235770)。translation／affine resampling 加自由 corrector 可测试几何位移是否绕过 amplitude-only 失败；它不是全新 shape optimization 理论。需记录插值、边界外填充值、box、体积变化与高对比细节平滑；不得用 GT 决定变换方向或 scale。

## 4. 最小决定性判据

保留主线程固定损耗路径和匹配普通延长对照，不扩γ sweep。先问零损耗最终的原数据门是否改善，再单列材料 NRMSE/IoU及总成本。若仅 stage loss好、最终回旧端点，或 ordinary延长同样解释全部改善，损耗机制未显示额外价值。若物理 residual改善而材料恶化，保留这一反证，不换 online ranking 或 joint门救结论。

若 data-only transform 数值检验找到同 fiber 的合法分离，支持的是**该离散合同下不可由这些数据确定 altered response**；人工损耗 heuristic 仍可作为求解策略，这两件事不矛盾。若该 pair未分离，也不证明全域存在 exact T；只能说明选定反例没有否定它。

## 收据与覆盖

见 `../literature/PRIOR_ART.md`、`CLAIM_LEDGER.json` 与 `receipts/`。ScholarQA三查询六操作全部429，final verify也429；Crossref已独立解析所列身份，出版社／作者原文提供A/B证据。Bevacqua完整PDF403，NIE2016本轮只有C级身份。没有穷尽查新或全球原创结论。
