# 对比度超阈值面积的严格下界

原始观测球内，R10 环形件至少 **2/576 格**的对比度严格高于 4，optfresh 01 开槽块至少 **3/9,216 格**严格高于 0.5。换算面积比例为 **0.3472%** 和 **0.03255%**。0.5% 面积目标分别需要 3 格和 47 格。96² 观测只有十束照明，单次源空间秩论证在最有利情形也至多证明 10 格。 [区间回执](area_receipt.json) · [图 PDF](area_certificate.pdf) · [图 PNG](area_certificate.png)

| 冻结观测 | 阈值 | 原误差半径 | 严格排除的最大高格数 \(q\) | 该 \(q\) 下的残差下界 | 证实高格数 | 离线真值高格数 |
|---|---:|---:|---:|---:|---:|---:|
| `r10_twod_0001` | 4 | 0.010517183214392477，CN99 | 1 | 0.010869825321694438 | ≥2 | 102 |
| `optfresh_01__pw360_10x20__n02` | 0.5 | 0.09096763488916171，原 2% 加生成误差余量 | 2 | 0.2110307543417967 | ≥3 | 1,703 |

R10 沿用原生 24²、12 发 48 收的匹配网格模型及冻结误差统计。optfresh 01 沿用全环 10 发 20 收的原 96² FFT 算子和原始 2% 观测。计数以等面积的原生单元为单位。报告针对已暴露的两件独立仿真材料。 [R10 原始算子](../../z_alg_secant_tasks/highdim/public_vie_inputs.npz) · [R10 原测量](../../z_alg_secant_tasks/highdim/results/r10_twod_0001.npz) · [原半径](../../z_alg_secant_tasks/highdim/observation_radii.json) · [optfresh 生成冻结](../../z_alg_optical_inertia/synthetic/FREEZE.md) · [原生设定和半径审计](../../z_alg_optical_inertia/results/theory/fine_feasibility.json)

## 证明

令 \(M=\operatorname{diag}(m_i)\)，\(E=B+GJ\)，\(J=ME\)，\(D=(SJ-Y)^*(SJ-Y)\)。取统一阈值 \(t\)。对实非负材料，非负扇区权重 \(\mu_i\) 的矩阵和为
\[
 Q_t=\sum_i\mu_i m_i(m_i-t)e_i^*e_i,\qquad e_i=E_{i,:}.
\]
实材料的光学矩阵和恒为零。若最多 \(q\) 格的 \(m_i>t\)，则 \(Q_t\preceq P_t\)，其中 \(P_t\) 只保留这些格子的正项，\(P_t\succeq0\)，\(\operatorname{rank}P_t\le q\)。等于阈值的格子贡献零。每个高格不论幅值大小只贡献至多一个秩。

对固定乘子记 \(Z=\operatorname{diag}((t\mu_i+\mathrm i\lambda_i)/2)\)，\(K=S^*S+\operatorname{diag}(\mu_i)-ZG-G^*Z^*\)，\(L=S^*Y+ZB\)。可以用任何经过验证的可逆下包络 \(\widehat K\preceq K\)，取 \(X=J-\widehat K^{-1}L\)，\(H=Y^*Y-L^*\widehat K^{-1}L\)。代入得
\[
 D \succeq H - \bigl(P_t+X^*\widehat K_-X\bigr).
\]
括号内为半正定，秩至多 \(q+\nu\)，其中 \(\nu=n_-(\widehat K)\)。设 \(h_j=\lambda_j^\downarrow(H)\)。降秩特征值交错和 \(D\succeq0\) 得
\[
 \|SJ-Y\|_F^2\ \ge\ \sum_{j>q+\nu}(h_j)_+.
\]
因此，只要验证右侧严格大于误差半径的平方，就有至少 \(q+1\) 格严格超过阈值。证明覆盖原有限算子的每个物理解，包括系数矩阵奇异而场解存在的材料。两件证书均验证了正的 \(\widehat K\) 绝对谱间隙或紧凑消元间隙。一般遇到奇异 \(K\) 时，可以选择并验证可逆的 \(K-\eta I\preceq K\)，再按新的负惯性计数。

固定一组乘子时，最强可能计数为 \(s-\nu\)。合并多组证书可逐个 \(q\) 取最大保守下界，仍受 \(s\) 限制。秩约束允许高对比度值任意大，也允许不同高格的源行相依。这个计算只给个数，不能定位像素或给出高格的幅值分布。等面积换算属于有限网格的材料模型。

### 等价的剪截视角

设 \(M_t=\operatorname{diag}(\min(m_i,t))\)。若剪截材料的原生 VIE 可逆，则
\[
 (I-M_tG)(J-J_t)=(M-M_t)E,\qquad
 \operatorname{rank}(SJ-SJ_t)\le q .
\]
也就是说，任意幅值的少数高格在整组读数中也只增加低秩扰动。对 \(D_t=(SJ_t-Y)^*(SJ_t-Y)\) 的残差奇异值做秩交错，再对低于阈值的完整材料盒用原证书，同样只移除 \(q\) 个源谱项。剪截推导要求该剪截点的物理响应有定义。上面的扇区证明直接处理原材料的存在场，因而无需剪截点的正演逆，也无需全局对比度—场 Lipschitz 界。

## 原球与灵敏度

R10 的逐格秩门在 \(q=1\) 严格超过原 CN99 半径，\(q=2\) 的最强已验下界为 0.00534598791243215。数值提案采用四乘四分块，即 16 个乘子区。用于原球结论的 \(q=1\) 空间乘子有 \(\nu=0\)。以同一含噪中心另作灵敏度演示时，假设误差半径为 \(0.005\|Y\|_F\approx0.005012485466\) 可条件性得到 3 格。保存真材料在同一原生模型中的严格残差约为 0.009931082542，超过这个更小的半径。原球科学结论维持为至少 2 格。 [分块权重提案](ring_spatial_proposals.json) · [真材料离线范围](truth_offline.json)

optfresh 01 的严格第三项单独给出 0.2110307543417967 的残差下界。使用原 2% 半径和加过生成误差余量的半径，计数均为 3。原生网格实现的下包络保留 625 个傅里叶模式和 20 个接收分量，完整保留全部 9,216 个材料自由度。真值的 1,703 格用于离线尺度对照。 [区间回执](area_receipt.json) · [真材料离线范围](truth_offline.json)

## 重放与费用

从项目根运行，Python 环境需 NumPy、SciPy、python-flint、threadpoolctl、Matplotlib。第一条重放数值提案，也可直接使用冻结 JSON。第三条仅作离线评分，公开复现两项证书不依赖它。

```sh
python -B research/z_alg_xudong_cie/area/search_ring.py
python -B research/z_alg_xudong_cie/area/area_certificate.py \
  --extra-ring-profiles research/z_alg_xudong_cie/area/ring_spatial_proposals.json
python -B research/z_alg_xudong_cie/area/score_area_truth.py
python -B research/z_alg_xudong_cie/area/figure.py
```

最终区间重放用 Apple M4、120 位 Arb、2 个区间线程和 1 个数值线程，墙钟约 11.85 秒。R10 五组 576 阶状态证书耗时 5.79 秒，optfresh 的 645 阶紧凑惯性和源谱证书耗时 6.06 秒。R10 分块提案耗时 2.01 秒，293 次浮点 LDL、3,504 列驻点右端和 292 次报告的源谱分解。此前的标量粗扫和两档加密探索合计约 9.86 秒，均只使用公开数据。区间认证与提案阶段均未运行材料正演。离线计数复用两个已存的真场误差区间。 [区间证书源码](area_certificate.py) · [区间回执](area_receipt.json) · [数值提案回执](ring_spatial_proposals.json)
