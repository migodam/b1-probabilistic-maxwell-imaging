# G19实际CPU参考下的强基线续算审计

本对照为POSTHOC_BUDGET_CONTROL，不是新确认算法。全部16×(CIE,Dual)共32对的触发、每对象T=G19实际完整childCPU、原停止原因及目标额度，先在不读GT/材料分数的阶段冻结。仅原CPU<T且不是成功PGTOL停止者触发：7对；21对PGTOL保持原结果，4对原CPU已达到或超过T保持原结果。没有按GT挑对象、扫描阈值、变更prior或改动原80主表/48合同。

7个续算均从原封存minPhi材料冷重启同一r6.refine_observed；L-BFGS历史和物理LU不复用，绝非原轨迹无缝延长。R/box/scale1e6/ftol1e-12/gtol.01/maxls30保持原值，maxiter1000000只作非绑定保护，按剩余whole-child CPU/原停止准则结束。保留原候选，以相同publicPhi选择更好的合格已付状态；全部32选择封存后才GT评分。

结局：CIE材料PASS仍4/16，Dual仍5/16；0门改善、0门回退。G19原7/16未动。均值NRMSE仅有小变化，故这次有规则边界的额外CPU续算未消除已测门差；这不证明额外计算普遍无用或结构的唯一因果优势。

|方法|触发|原PASS→续算视图PASS|mean NRMSE原→视图|新增实际CPU秒|该方法总实际CPU秒|G19参考总CPU秒|
|---|---:|---:|---:|---:|---:|---:|
|cie_single|4|4→4|0.231343→0.231240|33.370744|800.763917|973.533201|
|cie_dual|3|5→5|0.210683→0.210673|29.486995|806.604482|973.533201|

触发额度合计114.111821 CPU秒，实际新增完整childCPU62.857739秒。4个续算触remaining cap、3个以PGTOL自然结束，未空转补时间。触发者正unused合计52.043103秒；atomic/serialization/exit收费超支合计.789021秒。所有新所选状态累计合格paidCPU都≤各自T，迟到状态没有倒填；实际child总花费可以因这些已收费原子尾部略超T，因此不能称自然花费恰好相等。

新增物理成本为7370个576维LU、176880个complex RHS（含已付overrun）；不是新增观测生成。原32个baseline及其G19参考成本早已在80campaign中，只作部署视图，不再次加入研究成本。预算控制controller与postprocess另外计费，GT评分本身零FW。

文件：FROZEN_MANIFEST.json；QUEUE.json；ALL32_SELECTIONS_BEFORE_GT_SCORE.json和ALL32_SELECTED_PREDICTIONS.npz；ALL32_POSTSEAL_SCORES.json/.csv；SUMMARY_BY_METHOD.json；COST_LEDGER.json；POSTPROCESS_EXTERNAL_COSTS.json；各runs/receipt与PAID_STATES。源码快照包含原r6 corrector。
