# 失败与反例

1. 解析 conditional-volume：profile/no-det弱边缘保持O(1)偏差；x²必须双分支且fold厚度sqrtσ；固定弱rank补prior在σ→0失败。
2. V1固定 slope inverse失败且 identity POU 不自动补齐 chart 遗漏质量；原坐标 Gaussian cubature 对 curved/saturated ridge 未收敛。原 code/protocol/results 保留在`v1_failed_preserved/`。
3. V2有界Brent仍假设连续rank1。S1 host正质量rank0层不由该chart表示，324 inverse failures/case不可重归一化。原始`v2/run_v2_as_run.py`和结果保留。
4. V2 `BudgetStop`继承RuntimeError，最初被普通inverse except捕获；未突破硬cap，已显式先re-raise BudgetStop，旧作业中断费用保留。
5. V3初版NumPy view `U=J[:,1]`会修改缓存J；改`.copy()`并重启，从原已付费derivative文件恢复；首个response covariance缺层内latentvariance，保留为不合格，材料矩/质量不受这项遗漏影响。
6. V3 host mode-partition参考仍未收敛；不能引用旧host质量。直接GL readback改善但不构成严格tail/mass证书。
7. rank2额外求积触发32680状态软停止，尚未完成P1额外chart；该parent原坐标+baselines+RTO已完成。
8. 旧参考数值出现TV>1、负KL、q integral≠1，标INVALID_REFERENCE_METRIC，不能转成算法优劣。
9. parent A5六个中心/噪声充分门均FAIL，过保守bound不等于局部chart不可能。seam的三次退化来自decoder，不能归为新Maxwell物理零空间。

所有失败与setup费用已计入共享累计32705states/291784RHS；无退款、无预算重置。
