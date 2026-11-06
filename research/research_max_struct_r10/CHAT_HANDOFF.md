# MAX-STRUCT R10 · 执行中交接

本入口先公开已经完成的理论和历史开发证据。新24材料确认及真3D反演仍在执行，不能用下表替代它们。R1–R9未改写。

## 当前明确结论

1. 原A+并非把完整G放入分母：分母只吸收原Gii，随后在局部2×2块传播一次。尚无证据表明Gii数值写错。
2. 新实现 `R_K=β(I+βD_m−D_mK)⁻¹(J−D_mQ)`，K=PGP，材料仍为完整576个实未知；仅小矩阵逆用19/28维。current更新、FFT/APCS、prior和完整原Maxwell corrector继承。
3. Two-fold28与Green强19均恢复旧A+失败的40ac。保留原Gii的控制也成功；暂不需要把Two-fold接收信息当成成功必要条件。
4. 新方法不是全局凸化。固定状态更新中存在辅助目标下降、原物理Phi升高的反例。论文价值应看可重复的适用区域和成本，不要求所有例子获胜。

## 已完成的历史DEV材料结果

|方法|b20d NRMSE|40ac NRMSE|解释|
|---|---:|---:|---|
|成熟CIE旧R6参考|.5531|.6190|旧实现/历史费用单列|
|Dual-CIE（R9）|.3896|.6335|成熟sourcewise对照|
|Global-Dual（R9）|.4880|.6190|内部同current消融|
|A+（R9）|.0406|.6260|局部one-hop材料读出|
|Two-fold28矩阵分母|.0400|.0557|本轮，完整corrector后|
|等秩Green28|.0403|.0563|无接收子空间的控制|
|Green强19|.0409|.0548|更简单候选|
|Green19＋保留原自项|.0404|.0540|只作DEV机制控制|

所有对象已公开，均为EXPOSED/OFFLINE_DIAGNOSTIC。40ac在新分母的corrector前NRMSE约.15–.16，最终好结果包含共同物理corrector的贡献。R9的完整current低罚权路线也曾恢复40ac，因此这里不是项目首次恢复该对象，更不是已胜过所有优化方法。

## 公式驱动的新解释

24个保存状态中，局部块反馈谱半径全部<.168，但投影全域反馈在后期18个状态都>1；19维小分母的minσ仍为.312–.688。局部传播弱不能说明全域反馈可以忽略；全域Neumann级数不收敛与直接矩阵逆近奇异也不是一回事。

Green19残差在24/24状态上比原one-hop更接近完整shifted矩阵分母残差。这是辅助metric证据，不是原Maxwell适定性、GT恢复或全局盆地证书。

## 正在完成的确认

24个独立二维材料覆盖8家族×low/high/mixed；5个原对照加2个冻结新分母候选，共168次运行。新GT仍封存；所有候选只按原Phi在线选择。混合材料支撑评分在GT读取前冻结为共同χ>1，另保留旧Otsu诊断，避免把低介电部分自动抹掉。

真实3D使用0.6λ域、12³全材料、5184矢量current、60外部电偶极源及180接收分量。独立18³生成数据已取得。GPU正在服务其他用户任务，已停止自己的等待进程，改用同物理CPU FFT＋保留LU全RHS后端实际反演；并非将二维堆叠为三维。

## 阅读顺序

先读[理论](THEORY_DERIVATIONS.md)、[研究循环](RESEARCH_LOOP.md)、[历史材料图](delegated/denominator_score/MATERIAL_COMPARISON.png)、[能量审查](delegated/energy_audit/ENERGY_INTERPRETATION_REVIEW.md)、[谱诊断](delegated/denominator_spectrum/DIAGNOSTIC_REVIEW.md)。最终报告将根据真正完成的新确认更新，不能预设成功。
