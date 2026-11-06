# MAX-STRUCT-R9：最终研究交接

**本轮科学队列已全部结束。** R1–R8保持不变。本轮没有新增Fresnel实测或NN训练，不将继承证据记为新成果。

[研究网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_max_struct_r9/) · [中文报告](RESEARCH_REPORT_ZH.md) · [理论](THEORY_DERIVATIONS.md) · [英文初稿](PAPER_DRAFT.md) · [完整数组Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)

## 最终判断

**A/A+保留真实的高对比开发正例，但尚未做到“普通问题普遍打平、困难问题稳定击败成熟CIE”。** A改变材料读出的空间残差度量；A+降低其部分费用。Möbius条件更新没有凸化整个Maxwell材料目标。比原A便宜与比普通CIE多恢复材料是不同主张。

实际完成：第一阶段113次反演（47开发＋66确认）、经典Austria的16项Windows比较、小型全矢量三维29项、原3λ三维九臂，以及三起点各八步物理校正。Mac迁移前9个Austria终点及中断费用另存，不混入Windows性能。运行次数不是独立材料样本数。

## 关键证据

|问题|实际结果|判断|
|---|---|---|
|历史b20d弯条|A+ NRMSE约.0406；Dual/Global约.390/.488|有用路径，属于公开开发对象|
|噪声硬DPCS|公开噪声球选秩后，六例single/Dual/Global/A+全部2/6→4/6|共享基线修复，不是A+独有|
|修复后|平均NRMSE .18666/.18882/.18766/.18552；CPU38.84/40.81/42.98/46.37秒|A+未多恢复对象、非最低成本|
|40ac开缝壳|完整current松弛.05807；弱阶段校正.05835，归属CPU19.63 vs28.10秒|CSI/WRI型路径互补；新六例仍2/6，组合仍4/6|
|Austria ε3.5|CIE .18513、A+ .18450、direct .60683、独立NIE .18644、LSIE .61526|CIE优势真实；A+基本打平|
|Austria ε4|CIE .47475、Global .57517、A+ .58300|A+未解决更强散射失败|
|A+工程收益|Austria两条件完整CPU比A少24.5%/35.5%|整任务含不同轨迹；非单独内循环测速|
|小三维29项|共同校正后A/A+/Global终点相近，direct更便宜|无稳定A+增量；真实矢量、小电尺寸|
|原3λ matched|CIE/Global/A+/LSIE χ相对L2为.42727/.42275/.42347/.92953|CIE明显优于LSIE；A+未进一步改善|
|共同八步校正|Global/A+ χ误差降至.24899/.24950；air为.97501|均有有效物理方向，无A+独有收益；不是等总成本比较|

二维NRMSE与三维χ相对L2分母不同。经典NIE/LSIE为R=0参考，其余六臂同prior。全部候选按公开物理目标选择，GT仅封存后评分；公开论文几何不重新称盲测。

## 最值得继续的公式和现象

辅助电流J、e=L(m)J−D_mE、d=SJ−Y满足：
`F(m)−Y=d−S L(m)^{-1}e`。

这是已知恒等式。新证据是原3λ matched辅助data RMS只有.017–.027，完整物理RMS却约.37；实际缺陷方向的传播增益约7–13（不同β范数不可直接排名，上端不是operator-norm上界）。Global/A+同范数结果近乎相同，小系数步没有控制resolvent传播的端口缺陷。

共同八步校正确实改善两个端点，说明它们尚有有效原物理方向。但Austria ε3.5从预存80步继续到300步，材料NRMSE约.136恶化至.185–.190，single CIE的含prior Φ仅由.0043534降到.0043184。不能按GT改选80步，也不能将sqrt(2Φ)叫纯数据残差。

**下一轮优先：利用公开噪声和可信模型误差，区分需要释放的本构缺陷与不该继续拟合的端口残差。** 可研究可观测残差软约束和低罚权阶段的材料方向，不能默认更多A+步、更长corrector或更深网络能同时解决两种相反问题。严格递增loss变换保留局部极小，局部可逆换元保留驻点Hessian惯性。

## 物理与统计边界

- ε9球Mie散射平面误差8/10/16网格为61.6%/55.5%/43.1%；固定16网格积分Green为43.4%。已知界面arithmetic/tensor为27.8%/25.4%，仍不够准确；GT界面未交给反演。
- 原3λ CM+RR全60源30/60网格差53.5%；source59积分后50.8%→36.0%，仍远大于10%噪声。matched是另一个有限模型数据，不是同数据算法提升或连续物理认证。
- 27000实未知、至多21600实观测。维数限制内部常秩点，不能单独证明box边界GT多解。论文相对ε指标受背景影响，故并列χ误差和完整色轴图。
- 经典corrector全部300 CAP；原尺寸九臂达到论文系数停止；短校正全部八步预算停止。均非全局收敛证书。
- 六个确认材料不支持广泛统计优胜；没有新实测、复tensor或NN收益结论。

## 远端与费用

Windows二维四lane固定i9 P核逻辑0/2/4/6；三维4060使用resident complex128 CUDA/cuFFT，host支持P核12/14。队列全部结束，GPU已释放。

同一普通CIE三维current＋material步2.975→1.021秒（2.915倍布局优化）；fine单源RHS CPU99.91→CUDA5.89秒（包含Krylov实现差异）；二维4096材料value＋gradient仅1.400→1.103秒，GPU冷setup34.64秒，所以四P核跑二维、GPU跑三维。不能外推统一全任务加速。

科学/分析CPU已知下界约Mac4.64小时、Windows10.38小时；这是消耗，不是并行历时或GPU时间。最终归档另加入[总账](COST_LEDGER_EXTENSION.md)。LU/RHS/FFT按维度分账；失败、迁移、预计算和拒绝trial保留。未计量CPU、网络、LLM、能源和金额保持null。

## 入口与论文资格

[经典二维](CLASSIC_BENCHMARKS.md) · [小三维](VECTOR3D_EXTENSION.md) · [原3λ九臂](delegated/original3d_feasibility/ORIGINAL_CM_NINE_REPORT_ZH.md) · [短物理校正](delegated/3d_extension/physical_refinement/REPORT_ZH.md) · [物理诊断](VECTOR_DISCRETIZATION_DIAGNOSTIC.md)

[合同](BENCHMARKS.md) · [先例](PRIOR_ART_AND_NOVELTY.md) · [复现](REPRODUCE.md) · [第一阶段审查](delegated/independent_phaseI_review/REVIEW_ZH.md) · [三维审查](delegated/independent_original3d_contract_review/REVIEW_ZH.md) · [经典审查](delegated/independent_classic_review/REVIEW_ZH.md)

目前能写的是有原论文对象与真实矢量实验支撑的机制、基线修复和失效研究。广泛优于成熟CIE的算法论文仍缺：新独立对象上的稳定恢复增量、合格连续物理/实测支撑、等实际成本优势。本轮没有用期望发表的结论替代证据。
