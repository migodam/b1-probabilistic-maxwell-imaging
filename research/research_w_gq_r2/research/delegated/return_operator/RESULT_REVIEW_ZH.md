# E3 返回表示与一次复对称投影修复：实测交接

最终机器状态为 **COMPUTED_UNCERTIFIED**。第一轮与主线程明确授权的一次机制修复均完整完成，各 64 行；没有物理求解失败。保留原始失败机制、全部 paid actions、每秩误差与两个独立 wrapper 回执。后续没有增加墙、目标、秩、PDE 或冷 Bayesian 阶段。

科学裁决留给主 Codex。以下是可直接复核的有限公共 FE 空间结果，不是连续 Maxwell、全目标族、概率或算子范数证书。

## 固定范围

`f=1`，S0 两源、16 接收器；公开矩形墙 `public_wall(0)` 与 stepped/cavity 墙 `public_wall(2)`；合同中的 2 train + 前 2 heldout 目标。共同自由 FE 维数 **23,575**，公共目标 ROI 端口 **2,257**。

所有目标使用解析 incident 的 cut-polygon 弱载荷 `bX`；目标自作用由精确同空间 `(L0-VX)` sparse FE 因子处理。固定 seed 28102931、32 个复 Gaussian range 列、零 oversampling、嵌套 QR 前缀秩 4/8/16/32。全部 action/core 数值谱采用 Euclidean FE coefficient coordinates；ROI 场误差使用独立保存的固定 P3 L2 Gram 矩阵。

## 第一轮负面证据完整保留

原始模型为 `H_r=Q_r(Q_r* H Q_r)Q_r*`。它是声明过的同子空间双侧近似，但在此复对称、非 Hermitian 场景，主线程指出其输出 range 与输入空间的复共轭关系需要处理。

第一轮秩 32 的返回模型最坏白化响应 norm 为 **0.2309404**，目标诱导 ROI 场 L2 相对误差仍为 **0.1034679**，真实闭合 current 相对缺陷为 **0.0817206**。同一 paid Gaussian 输入上的 H action 相对误差仍为 **0.967–0.972**，尽管 `H Q` 投影残差约 1e-8。因此不能用小 core 谱或对 Q 的 action 残差声称完整 H 已收敛。

第一轮输出：`results/e3_return_operator/`，186 个文件、142,930,266 bytes；修复只读这些缓存，不覆盖该目录。第一轮机器汇总在本目录 `RESULT_REVIEW.json`。

## 唯一允许的机制修复

主线程授权保持同一 Q、目标模板、Bw、完整参考与全部秩，只改右投影：

```
H_r = Q_r [Q_r* H conj(Q_r)] Q_r^T
j_r = j0 + KXQ_r core_r [I-Q_r^T KXQ_r core_r]^-1 Q_r^T j0
```

只新增 `Gw conj(Q)` 与 `G0 conj(Q)` 的 paid actions。精确目标模板 `q/j0/KXQ/VQ`、Bw、uw、L0/Lw/V 和 full reference 全部重放已有 raw cache；没有重做目标模板或 full reference。脚本执行前逐项确认 mesh、free/port DOF、源/收、直接场、L0 与解析 forcing RHS 和原缓存一致，没有哈希检查。

以下每行取同一固定 8 个 wall-target 组合中的最大误差；都是 full-feedback 模型与直接 full FE 参考的比较：

|秩|修复后返回：白化响应 norm|修复后返回：目标诱导 ROI L2 相对误差|修复后返回：真实 current 闭合相对缺陷|同 H 基底 whole-C：白化响应 norm|
|---:|---:|---:|---:|---:|
|4|0.0734857|0.0340701|0.0267400|5.79516|
|8|2.04769e-4|1.69250e-4|1.41556e-4|6.08755|
|16|4.04446e-8|2.58322e-8|2.21071e-8|6.08290|
|32|4.46426e-11|2.56824e-11|2.18696e-11|5.92853|

修复后秩 32 的完整 FE 方程 RHS-relative residual 最大为 **4.35238e-12**；小闭环 2-condition 最大 **1.334705**。H 的原 paid Gaussian action 相对误差降到 **1.53e-9–3.94e-9**，小 core transpose-symmetry 相对缺陷约 **1e-14**。这些仍是已测方向/已测目标的数值，不是独立 operator-norm 认证，也不是连续场误差界。

whole-C 对照采用完全相同、仅拟合 H 的 Q，而非为 Cw 独立训练的最优子空间。它表现差不能推广为所有 whole-operator/普通 ROM 都差。C0+H 是标准散射分解，不能据本实现主张分解本身原创。

## 预约、CPU 与核验

|运行|实际 RHS|实际 FE|wrapper CPU 秒|墙钟秒|采样峰值 RSS bytes|
|---|---:|---:|---:|---:|---:|
|原始同 Q* 投影|708|19|13.647190|13.856988|911,818,752|
|一次 transpose 修复|256|3|7.498273|7.829529|666,910,720|
|本返回表示任务合计|**964**|**22**|**21.145463**|**21.686517**|最大 911,818,752|

两个 job 各自的 ledger reservation 与计划完全一致。原始运行 100 次已计费线性 action 的最大 backward 为 2.446e-18；修复 68 次的最大为 2.014e-18。接收伴随 Bw 和 paid true-field reconstruction 的最大绝对响应差分别为 1.403e-16、2.293e-16。两个固定纯矩阵实现核验分别为 3.52e-16 与 1.60e-16，只核对代码顺序，不计为物理实验。

第一轮的 708 RHS 包含 256 range/core、4 background primal、32 receiver adjoint、272 exact target-self template、16 direct full reference、128 true-feedback validation。修复的 256 RHS 为 128 conjugate-basis actions + 128 full-feedback validation。所有 RHS/LU 在 action 前由 `r2.physics` 自动预约，失败不会退费。收费 CPU 仅取 wrapper completion；文件中的内部 CPU 字段为信息计时。

这里没有等成本加速结论。第一次使用的 8 个 direct full 参考只需 16 RHS/8 FE，已清楚单列；返回表示的 setup、模板、原始不足和修复成本不能隐藏。E3 cold Bayesian 为 **NOT_RUN_SEPARATE_GATE**。

## 交付定位

- 原始代码：`src/r2/return_model.py`、`tools/return_operator.py`。
- 修复代码：`tools/return_transpose_repair.py`，默认只输出 PLAN_ONLY。
- 原始与修复 raw：`results/e3_return_operator/`、`results/e3_return_transpose_repair/`。
- 每秩逐例数据：两目录各自的 `PER_RANK_ERRORS.json` 与 `.csv`；每行保存 current、field、data、closure 与 paid full-feedback reconstruction。
- 原始账本核对：本目录 `RESULT_REVIEW.json`；修复账本与聚合：`TRANSPOSE_REPAIR_REVIEW.json`。
- 计费回执：`results/jobs/return-operator-frozen/RECEIPT.json` 与 `results/jobs/return-transpose-repair-frozen/RECEIPT.json`。

主线程需要介入的仅为最终科学解释、工程门映射和整份 R2 报告融合；本子任务不作新 novelty、最优 ROM、速度或 probability acceptance 裁决。
