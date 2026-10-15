# 本轮理论澄清

## 噪声球与竞争似然不相同

在固定 sigma 的 512 维标准实高斯噪声下，令 δ=W(F_b−F_a)、d=||δ||、y=F_a+noise，则 Δell=d²/2−ξ·δ。两端同时位于 tau 噪声球，只给出一个宽松的集合条件；它不保证 |Δell|≤log(10)。点似然比为 exp(−Δell)，仍不是 posterior 质量比。

## 原先验改变了可补偿方向

自由 complex box 中允许实虚部分开变化，以及大面积 DC/低频补偿。原 P0 的端点必须由真实参数生成，经同一个非负实 pattern 和固定 v 映射，包含原峰值归一化。P0 允许方向比完整 complex box 小，phase 测试通过仍不足以证明生成器支持。历史自由补偿成功不能推出 P0 内存在同样的补偿。

原 generator 的 rasterization 可能对小几何变化返回完全相同材料。这里复用真实 rasterizer 的有限参数 secants，没有把零几何导数当作物理可辨识性证明，也没有悄悄平滑后替代 target。多个 latent 对应同一材料时按精确材料数组复用已付费端点，cache provenance 单列。

## A17 的适用边界

普通 full material Jacobian 已使用真实 resolvent。E1 的 Jacobian 排序只提出原生成器端点；所有接受判断都用 full solver。没有 uniform nonlinear remainder bound，不能把线性预测或谱排序当成不可补偿证书。没有调用新的 A17 学习/ROM，也没有在新增原先验 active constraints 后套用无约束 paired-GN 保证。

上述为已知概率与物理恒等式在本问题中的合同澄清，不作为原创定理。当前没有强 IEEE TAP 理论／算法贡献证据。
