# 历史强 CIE 的电流—体素交替

MAX-STRUCT R10 中的成熟 CIE 是固定材料更新电流，再固定电流更新 **全部 576 个共享实材料体素**。R10 比较继承的 CIE-FFT 用快速算子完成电流步，材料步调用 FBG-R3 的共享实体素读出。所有方法共同使用原 Maxwell 物理目标作阶段选点和最后校正。

## 电流更新

对接收算子 $S$ 做 SVD，记保留的奇异向量和奇异值为 $U_r,\Sigma_r,V_r$。数据确定电流与剩余电流写作

$$
J_D=V_r\Sigma_r^{-1}U_r^*Y,\qquad
J=J_D+Ua,\qquad U=(I-V_rV_r^*)F_{\mathrm{mf}}.
$$

$F_{\mathrm{mf}}$ 是库中所选模式的逆 DFT。固定 $m$，令 $D_m=\operatorname{diag}(m)$，$P_{\beta_j}(m)=\operatorname{diag}[\beta_j/(1+\beta_jm_i)]$。对系数 $a$ 的辅助目标是

$$
\tfrac12\lVert S J-Y\rVert_{\Omega_d}^2
+\tfrac{\rho}{2}\lVert P_{\beta_j}(m)[(I-D_mG)J-D_mE_0]\rVert_{\Omega_s}^2,
\qquad \lVert X\rVert_{\Omega}^2=\operatorname{Re}\operatorname{tr}(X^*X\Omega).
$$

R10 的 $\rho=1$。源权重为固定对角矩阵，数据列权为 $1/\|Y_s\|^2$，状态列权为 $1/\|E_{0,s}\|^2$。CIE-single 和 Dual-CIE 每轮各源执行一次 Polak–Ribière+ 共轭梯度。令 $g$ 为当前复系数的实内积梯度，$d=-g+\gamma d_{\rm old}$，其中 $\gamma=\max(0,\operatorname{Re}\langle g,g-g_{\rm old}\rangle/\|g_{\rm old}\|^2)$。每列若 $\operatorname{Re}\langle g,d\rangle\ge0$ 则重启为 $-g$，再用该固定材料二次式的精确一维步长

$$
\alpha=\frac{-\operatorname{Re}\langle g,d\rangle}{\|SUd\|_{\Omega_d}^2
+\rho\|P_{\beta_j}(m)(I-D_mG)Ud\|_{\Omega_s}^2},\qquad
a\leftarrow a+\alpha d.
$$

上式对单列使用对应源权，GlobalDual 对整矩阵用同一标量步长。嵌套 Fourier 掩模之间按模式零填充旧系数。R10 的 DPCS 秩取数值可支持范围内使数据投影尾范数不超过原始 CN99 噪声半径的最小值。[R10 选秩](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r10/src/twod_adapter.py#L110) · [实际电流更新和诊断](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r8/src/cie_fft.py#L228)

## 体素更新

固定 $J$，总场 $Q=E_0+GJ$ 随之固定。对每个体素定义 $V_i=Q_i+\beta_mJ_i$、$T_i=\beta_mJ_i$ 与 $r_i=\beta_m m_i/(1+\beta_m m_i)$。变换后的状态残差是

$$
\frac{\beta_m}{1+\beta_m m_i}(J_i-m_iQ_i)=T_i-r_iV_i.
$$

一个实标量 $r_i$ 同时拟合该体素的所有源。无先验最小二乘解是

$$
r_i=\operatorname{clip}\left(
\frac{\operatorname{Re}\sum_s w_{s} \overline{V_{is}}T_{is}}
{\sum_s w_{s}|V_{is}|^2},\ 
\frac{\beta_m\ell}{1+\beta_m\ell},\
\frac{\beta_m u}{1+\beta_m u}
\right),\qquad
m_i=\frac{r_i}{\beta_m(1-r_i)}.
$$

该盒范围内的变换单调。分母数值为零时保留旧 $r_i$。有共同先验时，此点仅作为种子。实际代码比较它与旧材料的固定电流目标，用较优者启动盒约束 L-BFGS-B，目标为 $\tfrac{\rho}{2}\sum_{i,s}w_s|T_{is}-r_i(m_i)V_{is}|^2+R(m)$，梯度用 $-\beta_mV_{is}/(1+\beta_m m_i)^2$。内循环最多 60 次，ftol 为 $10^{-12}$，gtol 为 $10^{-8}$。若内循环候选非有限或比两个启动点中较优者还差，则保留较优启动点。先验 $R(m)$ 为差分平方均值和平滑 TV 均值，各乘 $10^{-4}$。差分按 19.5 缩放，TV 平滑常数 0.01。[原材料函数 `_readout`](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_fbg_r3_experiments/src/cie_tsom_baseline.py#L61)

## 历史停止和诊断

R10 主组从全体素 0.1 起步，Fourier MF 为 3/6/10/12，各执行 500/500/500/1500 步。材料 $\beta_m$ 依次为 1/.2/.2/.2。Single 的 $\beta_j=\beta_m$，Dual 的 $\beta_j$ 恒为 1。代码记录相对系数变化平方、当前步辅助目标改变量、两种状态项、材料变化和内循环终止。冻结配置把系数早停阈值设为 0，因此这些阶段用迭代配额结束。单 $\beta$ 的两个子步对应同一个阶段辅助目标，双 $\beta$ 对应不同状态权，不能由单步下降推出共享目标收敛。

初始材料与各阶段端点各自做精确 Maxwell 正演，按 $\Phi(m)=\tfrac12\|Y_{\mathrm{pred}}(m)-Y\|_F^2/\|Y\|_F^2+R(m)$ 选启动点。最后用同一盒约束原目标 L-BFGS-B 最多 3000 迭代，保存最低已付 $\Phi$ 点。校正器以 $10^6\Phi$ 优化，gtol 0.01 对应原目标投影梯度无穷范数 $10^{-8}$，还有相对目标、迭代及求值上限。系数变化、辅助目标、原 $\Phi$、原数据残差和离线材料误差分别解释。[冻结阶段](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r6/src/benchmark_solver.py#L34) · [原物理校正](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r6/src/benchmark_solver.py#L164)

## 强基线与具体失败对象

- **FBG-R3 的已披露 576 维 E3 组。** Dual-CIE 加共同原物理校正联合恢复 10/12，成功对象为数字 1–9 和双岛，环与弯条失败。Single 为 8/12。Dual 的 11/12 条轨迹被运行时间截断，结果对应当时预算。[历史结果与停止状态](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_fbg_r3_experiments/CIE_TSOM_BASELINE_AUDIT.md#L96)
- **MAX-STRUCT R3 的旧弯条。** `obs_b20d5a55677d4eca` 从 MAX-STRUCT R2 保存。它有自己的源、接收和体素坐标，名义噪声 0.1%。R2 存档的 native Dual-CIE 终点为 NRMSE .510695、Otsu IoU .245902、相对物理残差 .100735，残差为 $2417.64\sigma$，CN99 半径为 $25.165\sigma$。R3 的 11 个控制沿用这个开发对象。去噪 native 和放开强模态仍失败；人工损耗路线把目标相对残差降至 .051017，同一例材料 NRMSE 升至 .817495。该组历史 NRMSE 分母固定 12，IoU 对比预测独立 Otsu 与存档面积支撑。`cohorts.load_r3_curved_bar()` 读原观测，`load_r3_curved_bar_truth()` 和 `legacy_r3_material_metrics()` 供离线复评。[R2 逐对象原表](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r2/results/BLIND_SCORES.csv) · [R3 十一控制](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r3/results/DEVELOPMENT_SCORES.md)
- **Z-ALG R3 的同轮编号。** `research_z_alg_r3` 是另一组 $12\times12$ Full-J、GlobalDual 与坐标路线实验。其 V6 案例和成功定义不对应上述 CIE 弯条或 R10 的 576 维组。[Z-ALG R3 原表](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_z_alg_r3/results/HOLDOUT_V6_SUMMARY.md)
- **MAX-STRUCT R10 的原主 24 案。** CIE-single 材料恢复 11/24，平均幅度 NRMSE .162871，每对象实际完整 CPU 均值 157.705 秒。Dual-CIE 为 10/24、.163460 和 153.982 秒。CIE-single 有两个数据不合格对象，另有 11 个对象的数据误差及真值均在各自 CN99 球内，但材料评分未合格。三份历史坏端点在 R12 被用于续跑开发，来源分别是 CIE 环 `0001`、G19 壳 `0004` 和 G19 桥 `0022`。[R10 冻结主表和失败标识](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r10/delegated/twod_score/scored/FRESH_TWOD_REVIEW.md) · [R10 三例轨迹](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/ec3a37da6f3ee90beb61e753067d4f639f4118f5/research/research_max_struct_r10/delegated/twod_score/scored/THREE_CASE_INTERPRETATION.md)

|R10 案例|原方法终点|幅度 NRMSE|固定支撑 IoU|相对数据残差|原诊断|
|---|---|---:|---:|---:|---|
|`0001` 高对比环|CIE-single|.42288|.35333|.13684|超过 CN99。保存原目标投影梯度约 $9.68\times10^{-9}$。G19 在同观测上为 .09889 NRMSE。|
|`0016` 高对比包覆包体|CIE-single|.46502|.07299|.28755|超过 CN99。表中 Dual-CIE、A+ 和 G19 也到达类似坏端点。|
|`0004` 高对比开缝壳|G19|.51514|.08696|.58980|超过 CN99。CIE-single 恢复，幅度 NRMSE .11109。|
|`0022` 高对比桥|G19|.33700|.33036|.16850|超过 CN99。以相对目标变化终止，保存原目标梯度约 $3.18\times10^{-5}$。CIE-single 恢复，幅度 NRMSE .09741。|

`0000` 和 `0021` 在 .15 材料阈值附近发生二值变化，量级有别于表中大分叉。其余 11 个 CIE-single 材料未合格但数据已合格对象为 `0003`、`0007`、`0008`、`0009`、`0010`、`0011`、`0012`、`0013`、`0014`、`0018`、`0020`。原始记录支持逐类比较当前更新、原物理目标和材料恢复。
