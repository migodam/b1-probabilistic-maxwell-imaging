# 原生 96² 局部补偿

在 400 MHz、2 m 见方、32 发 32 收、2% 实现噪声的六个冻结观测上，固定上限 1 和 v3 的最强重放权重 `scalar_grid_601`。只要逐像素非负框架的支持至多 52 像素，原方案定理 2.4 的全盒证书在六个观测上均无法过噪声半径。这里对像素选择、非负权重和带正负特征值的中心矩阵 $M_Y$ 都取了统一上界。上限 1 的材料可行性仍未知。

按独立审读再次对全部 9,216 格排序后，通用范围已扩至 **任意 59 格**，对 52 格的上限也进一步收紧，并检验了显式弱正模。下文保留最初得到的端点。[修订后逐件单侧结果](deflated_sparse_tightening.json) · [正模分拆与严格适用范围](deflated.md)

## 严格结果

|上限与权重|负惯性|负、正侧谱隙下界 $K/\mu$|旧距离/半径最大值|旧排除|
|---|---:|---:|---:|---:|
|0.75，stable scalar optical|4|0.00981327、0.03125|2.20360443|5 / 6|
|1，scalar_grid_601|13|0.03017404、0.001953125|0.30105522|0 / 6|
|1，stable scalar optical|13|0.05512557、0.015625|0.26099909|0 / 6|
|1.25，scalar_grid_566|15|0.13006815、0.015625|0.11331174|0 / 6|

上限 1 的基线最大值对应 `optfresh_04__pw360_32x32__pair2pct`。上限 1.25 对应 `optfresh_03__pw360_32x32__pair2pct`。在中心扣分矩阵要求 $M_Y\succeq0$ 时，去掉中心扣分的乐观控制要求至少补偿 7、5、5、5、5、7 个方向。以 416 个已选像素能够补偿的上界 4 代入，六件最好的乐观距离/半径上界为 0.9934882。带任意符号的中心矩阵须另算中心项。[完整源谱](u100_best.json) · [乐观控制逐件数组](OUTWARD_OPTIMISTIC_CONTROL.json)

|上限 1 的支持范围|补偿维数 $c$ 的严格上界|依据|
|---|---:|---|
|任意至多 26 像素|0|13 维公共负子空间上的补偿迹上界 0.88113140|
|任意至多 52 像素|2|11 维公共负子空间上的迹上界 0.95401746|
|冻结 QR 像素前缀 8、26、52 行|0|在分别固定的全均匀材料下，完整有符号行 Gram 为正定|
|冻结扩展 208 行|3|12 维迹上界 2.43088371|
|冻结扩展 416 行|4|12 维迹上界 3.10889072|
|全域 9,216 像素内的任意支持|8|9 维迹上界 4.60606637|

冻结 8、26、52 行的有符号 Gram 最小特征值严格下界为 0.97906988、0.94712391 和 0.83351776，对应全均匀材料 0.40833、0.39167 和 0.41667。对这些固定支持，单条被补偿负方向也无法通过要求全盒成立的框架。这个结论允许在所选联合行上使用任意半正定行权矩阵。全域还有一个独立的一维公共负方向见证，实际加权响应至多 0.71325214，所需负能量至少 $1-2\times10^{-15}$。[稀疏上界](u100_best_strict.json) · [联合行 Gram](u100_best_joint_cut.json) · [任意支持全域见证](u100_best_all_pixel_strict.json) · [密集维数上界](u100_best_dense_outward.json)

当支持至多 52，允许任意符号的 $M_Y$ 后，六件的距离/半径统一上界为 0.799968、0.768009、0.978116、0.878793、0.870330、0.598073。上限 0.75 时，任意至多 16 像素都只能取 $c=0$，原来未排除的第六件在支持至多 4 像素时的上界为 0.820343。上限 1.25 时，任意至多 60 像素都只能取 $c=0$。支持至多 30 像素时，六件上界均低于半径，最大为 0.898574。支持至多 60 像素时，第 01、05 件在该上界下仍开放。超过 52 像素时，上限 1 的任意符号中心矩阵路线仍有余量。所选 208、416 行的 $c$ 结论覆盖 $M_Y$ 的任意符号，排除距离上限则使用 $M_Y\succeq0$ 的乐观控制。[0.75 外包](u075_strict.json) · [1.25 外包](u125_strict.json)

## 达到噪声球所需的负方向电流

上限 1、同一固定权重下，每件材料若要达到其噪声球，都需要非零的**实际电流负谱分量** $j=\|\Pi_-J\|_F$。这里电流 $J$ 保持原始定义。用正空间消元得到 $H_+=Y^*Y-L^*K_+^\dagger L$。对数据选出的固定源矩阵 $0\preceq Z\preceq I$，令 $h=\operatorname{tr}(ZH_+)$、$b=\|V_-^*LZ\|_F$、$a=\|(-K)_+\|_2$，则残差平方至少为 $h-2bj-aj^2$。将三个量外包后，六件各自得到以下必要下界。

|观测|$j$ 下界，原生电流单位|同一源矩阵上 $\sqrt h/r$ 下界|
|---|---:|---:|
|optfresh_01|15.75970|1.28509|
|optfresh_02|9.05648|1.19750|
|optfresh_03|14.42828|1.38886|
|optfresh_04|22.00850|1.66512|
|optfresh_05|19.76867|1.42410|
|optfresh_06|7.65458|1.27032|

这些界适用于所有达到噪声球的物理材料。逐个条件同时复用已验证的真实驻点误差、负子空间泄漏和原始出入射矩阵。六件能否达到噪声球仍待判定。[逐件谱与不等式外包](required_currents.json)

## 计算怎样同时保留完整中心和原生网格

驻点使用完整的 $L=S^*Y+ZB$，一次解 32 个接收端和 32 个入射端，组合成六组 $J_0$、$E_0=B+GJ_0$ 及 $H$。每个成员的线性求解残差都由 Arb 复核。下方的正框架仅处理 $K/\mu$ 以改善数值尺度，中心项和补偿期间的行矩阵包含对应的 $\sqrt{\mu}$ 因子。读取已存的 192² Green 乘子、B32、S32、各个 Y 和已有半径，全程保留 9,216 个独立实材料单元和原接收残差 Frobenius 球。

正空间端点 Gram 用投影 CG 求出，随后把投影误差、FFT 误差、解残差及正谱隙合为矩阵上下包。主例的 104×104 Gram 上包为记载的中点加 $4.74487\times10^{-9} I$。首轮发现复数 ARPACK 在重根簇上的 Ritz 基 Gram 最小值只有 0.740。直接用于投影 CG 会放大端点误差。重新做 QR 和小矩阵 Rayleigh 对角化后，Ritz Gram 下界达到 0.9999999999999859，严格的正方向泄漏不超过 $2.20\times10^{-11}$。两次原始记录均保留。驻点使用独立线性求解，首轮驻点的 Arb 回执继续有效。[重正交](u100_best_conditioned.json) · [初次端点误差](u100_best_strict_initial_nonorthogonal_full.json)

稀疏上界以 $w_i\le 1/q_i(0)$ 为起点。对正投影 $\Pi_+$，Cauchy–Schwarz 给 $q_i(0)\ge \|\Pi_+e_i\|^4/(e_i^* K_+ e_i)$。负投影的逐像素杠杆量和对角负能量由外包 Ritz 残差控制。由此得到每个像素允许的权重上界。选一组固定负方向和一套合法逐像素材料，取最坏的前 $k$ 个行容量相加。若整个 $d$ 维负子空间上的加权迹小于 1，所有可补偿子空间与它相交的向量都会违反定理 2.4 的左上块条件，所以 $c\le\nu-d$。

任意符号的中心矩阵由同一定理右下块给出 $-M_Y\preceq d^*d$。选源空间尾部的固定压缩，再选择一个使该压缩上中心预算较小的合法逐像素材料。交错关系和正部迹的次可加性给距离上界。驻点及完整场的 Arb 误差传播进入每个上界。压缩和负空间两步可选择不同的材料点，因为框架及中心矩阵的条件须在全盒成立。

联合行检查在固定均匀材料点直接验证 $A K^{-1} A^*=\Gamma-(AN)(AN)^*\succ0$。等价地，$(AN)^*\Gamma^{-1}(AN)\prec I$。正框架限制的任意半正定行权满足 $W\preceq\Gamma^{-1}$，从而原定理左上块在该固定支持无法补偿任何方向。密集检查仍使用同一个原生量。对几个正空间见证，逐像素选更强的框架端点，得到相应权重的线性不等式。用 LP 只寻找对偶系数。逐系数 Arb 下包和剩余权重上包构成最终的严格对偶。若某个 $d$ 维空间的压缩迹至多 $\beta$，在该压缩上至多 `floor(beta)` 个特征值可达 1。交错关系给 $c\le \nu-d+\lfloor\beta\rfloor$。

## 中心感知 SDP

先对两端负行用带主元 QR 在其列秩以内选独立行，尾部补齐固定支持。主例的所选 52 像素用严格端点 Gram 控制正框架，对各行 $[r_0,r_1]$ 用二阶系数空间构建完整区间 SOS，再用源矩阵目标，不把中心压成标量。主例固定方向最大容量预筛远小于 1，联合 Gram 独立证实该支持的 $c=0$ 边界。保留带符号中心修正的 $c=0$ SDP，观察 04 件。首轮 Clarabel 失败，SCS 得到 `optimal_inaccurate`，其可用结果再由 Arb 核验。该方案比值下界 0.30104574，冻结旧下界 0.30105522。两者取最大值后沿用旧下界。附加两个同合同检查也低于各自旧下界。上限 0.75 的第 06 件得到 0.33023727，上限 1.25 的第 03 件得到 0.11330250。[0.75 SDP](u075_center_sdp_6_4.json) · [1.25 SDP](u125_center_sdp_3_4.json)

密集的一维固定方向搜索允许 208、416、9,216 行的必要 LP 容量上界分别达到 3.46、4.25、4.38。它们只是放宽框架条件的上界，最终允许的维数见上表。[SDP 回执](u100_best_center_sdp_4_4.json) · [密集搜索记录](u100_best_dense_frame.json)

## 复算和成本

[合同](CONTRACT.md) · [细分费用](COST_LEDGER.tsv)

从 Inv_SLAM 根目录运行。依赖 NumPy、SciPy、python-flint、threadpoolctl 和 CVXPY/SCS。实际环境为 Apple M4，单 BLAS 线程和两 Arb 线程，112 位精度。线性求解、失败尝试、数学验证和 solver 求值分别保留。并行任务的墙钟各自计量。以下顺序保留独立复算所需的原始数组。

```zsh
native_dir=research/z_alg_optical_inertia/ip_unified/numerics/native
native_python=.artifacts/z_alg_diag_python/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
for native_case in u100_best u075 u125 u100_stable
do
  "$native_python" "$native_dir/replay.py" "$native_case"
done
for native_case in u100_best u075 u125
do
  "$native_python" "$native_dir/condition.py" "$native_case"
  "$native_python" "$native_dir/outward.py" "$native_case"
done
"$native_python" "$native_dir/global_frame_probe.py"
"$native_python" "$native_dir/global_outward.py"
"$native_python" "$native_dir/joint_cut.py" u100_best
"$native_python" "$native_dir/dense_frame.py"
"$native_python" "$native_dir/dense_outward.py"
"$native_python" "$native_dir/summarize.py"
"$native_python" "$native_dir/required_currents.py"
"$native_python" "$native_dir/center_sdp.py" u100_best --object 4 --rows-factor 4
"$native_python" "$native_dir/center_sdp.py" u075 --object 6 --rows-factor 4
"$native_python" "$native_dir/center_sdp.py" u125 --object 3 --rows-factor 4
```

每个 cap 的 `.npz` 保存六件 $J_0,E_0,H$ 及线性响应，`*_conditioned.npz` 保存重正交负基、像素、正 Gram 和投影解，`*_strict.npz` 保存逐行上界和中心评估所用的材料点。密集与全域的 `.npz` 保存正见证、LP 对偶和逐像素严格密度。冻结数据位于 [历史阵列准备](../../../physical_revision/array/FREEZE.md) 和 [合成观测说明](../../../synthetic/README.md)，没有新增真值或观测。
