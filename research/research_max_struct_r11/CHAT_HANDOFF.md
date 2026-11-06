# MAX-STRUCT R11 · 执行中的证据交接

**状态：IN_PROGRESS。扩大确认与三维反演正在执行，不能把本页当作最终结论。** R10冻结证据保持。本页先公开已完成的公式、开发正反例和冻结合同；后续会补齐独立结果、全部成本与论文资格判断。

## 已完成的主要发现

1. 分母梯度中的强抵消主要含普通CIE已有的β收缩项。10个真实曝光检查点中，β项范数是Green新项的2.288–5.390倍，早期方向相反；不能把抵消全部归咎于Green。
2. freezeA与保留CIE收缩的contractedT都保住了环，却没有救回开缝壳和桥。固定材料块凸不等于逆问题凸，二者也不是moving分母目标的MM。
3. 新控制BothBlocks把Green分母真正接入current与material两块，保持原FFT/APCS、PR、β、初值、prior和corrector。历史桥的NRMSE由材料-only Green的.337降到.10088，同时环仍约.105；开缝壳仍约.515。桥已在辅助stage2/3改变路径，不是同起点corrector突然变强。
4. CIE在桥上约.097，仍略优于BothBlocks；因此这是修复Green对桥的损害，不能写成桥上打败CIE。环的CIE约.423，Green各版本约.10，是保留的开发优势。
5. 真实同状态current曲率只中等改善，reset方向cos却很低；单一条件数改善不能预测成功。标准resolvent结构给出强Green通道饱和的解释，但尚未建立充分的成功预测规则。

![全部曝光正反对象](FIGURES/exposed_development/MATERIAL_COMPARISON_3x8.png)

## 已冻结并实际启动的扩大研究

|组|独立材料|实际计划反演|硬件与范围|
|---|---:|---:|---|
|主二维确认|48，8拓扑×high/mixed×χ8/12/18|288，6臂|i9 P核，CIE/Dual/Global/fullG19/freezeA/contractedT，480 whole-child CPU与成本曲线|
|BothBlocks后续确认|16，8拓扑×high/mixed|80，5臂|Mac两lane，CIE/Dual/Global/fullG19/BothBlocks；全臂同硬件，不合并Windows原CPU|
|新三维确认|8个规定设计，4拓扑×high/mixed|40 matched +10 cross-grid|4060真3分量dyadic，12³全部1728材料未知量，sourceCIE/Dual/Global/fullGreen/contractedLagged|

二维两组都在各自新对象生成前冻结，尚未读取新GT评分。三维生成端已完成8 matched12及2个18→12观测，反演有限队列会在曝光DEV之后自动接续；这不是已经取得50个三维重建结果。cross-grid两观测不是两个新独立材料，且不与matched成功数合并。

三维分母采用公开固定`span(Einc,G Einc)`的rank120投影。top32最大奇异模式对晚期defect传播遗漏很大；source空间也不完整，不能宣称完整G或连续Maxwell精度。full/freezeA/contractedT已完成一次必要的真实3N共享材料导数接线核对，完整反演正在运行。保持FP64、真实LU/RHS与全部失败成本。

## 可直接阅读

- [公式、证明与边界](THEORY_DERIVATIONS.md)
- [研究如何根据失败转向](RESEARCH_LOOP.md)
- [文献与新颖性定位](PRIOR_ART_AND_NOVELTY.md)
- [真实轨迹图](FIGURES/exposed_development/STAGE_CORRECTOR_TRAJECTORIES.png)
- [主48对象冻结合同](delegated/twod/INDEPENDENT48_CONTRACT_PROPOSAL.md)（文字提案保留；实际冻结状态以configs中的正式JSON为准）
- [BothBlocks16冻结合同](delegated/diagnostics/bothblocks_confirmation/FROZEN_CONTRACT.json)
- [三维正式冻结合同](configs/vector_confirmation_contract.json)

尚不能回答新算法在新48/16或三维上是否改善恢复。本轮不会以三个开发正例或一次gate替代这些正在执行的实验。新确认不得读取GT调参；最终报告会包含全部失败、原物理残差、材料误差、等允许CPU/实际费用/物理LU三种比较及其区别。
