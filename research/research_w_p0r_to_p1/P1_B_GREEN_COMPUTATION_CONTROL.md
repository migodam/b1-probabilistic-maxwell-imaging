# Green计算复用：16个新contrast的真实工程对照

这是固定墙、固定圆形目标几何、单一频率f=1、四源/十二收上的局部工程控制，不是五频posterior或一般墙体更新算法。目标contrast在[2.94,3.06]的16个预定值变化，σ=.02固定；warm表示新contrast使用同一base LU，不是重复返回同一个响应缓存。

全部16×7=112行完成。圆边界64/512的独立Full资格端点最大白化差0.007284，小于0.05。exact Woodbury使用全部1228个结构活跃行，没有SVD截断或伪逆；单标量contrast不意味着该离散材料更新只有一个端口。

| 求值器 | 最大白化响应误差 | 平均warm CPU/秒 | 平均warm wall/秒 | warm新RHS合计 |
|---|---:|---:|---:|---:|
| Optimized Full | 0（比较参考） | 7.5273 | 9.2687 | 64 |
| Exact all-row Woodbury | 1.45e−12 | 4.0768 | 4.5532 | 0 |
| Ordinary field goal-ROM V | 1.02e−6 | 2.2821 | 2.4597 | 0 |
| Ordinary field ROM span(V,R0V) | 8.42e−7 | 2.2954 | 2.5986 | 0 |
| Current Galerkin V | 0.0007161 | 2.2597 | 2.3337 | 0 |
| Same-current A17 algebra control | 0.0007161 | 2.2611 | 2.3492 | 0 |
| Baseline-LU preconditioned GMRES | 1.47e−12 | 81.9030 | 92.3811 | 416 |

warm表不含前置setup，也不是已测定的standalone端到端速度比。共享ROM setup中32 RHS包含primal4、dual12、R0V16，不得再把该聚合与三个子项相加。Exact额外setup为1228个base-LU RHS、约75.02 CPU秒；还需base、投影与响应准备。整个实际工程scope支付1760 RHS、21次LU、root wrapper CPU2189.726秒、wall2442.255秒；峰值RSS约4.049GB。scope内部测量CPU2174.920与root wrapper有范围差，不能二次相加。另有80次reduced dense-core LU（16次1228维exact core、48次16维ROM/current core、16次32维enriched core），CPU已在wrapper中；它们与21次full FE材料–频率LU分开记录，不能以warm表的Full-LU=0隐去小矩阵因子分解。

## 不能据这张表得出的结论

Krylov使用Full已经付费组装的A1；各非Full求值也使用其系数采样ε。源码记录了当前矩阵/输出准备的共享来源，但没有无重叠的standalone冷部署费用。另一方面，Krylov416个base action的实际LU kernel合计CPU只有9.844秒，远低于其warm全管线CPU；后者包含逐action预算读取、残差核验、事件记录、持久化和输出。因此不能把表中慢值解释为Krylov线性代数本身必然慢，也不能用只计kernel的值宣布部署胜出。此处暴露了instrumentation/工程开销，需要独立新benchmark，不能从现有时钟反推一个未实测的优化后费用。

没有实际运行response精度筛选→失败Full fallback的部署policy；所有参考Full本来就支付了。差的近似不能借这个pool免费fallback，单频小响应误差也不能替代完整五频|ΔΦ|或MH决策资格。没有strict likelihood区间，不称认证早拒绝。

## A17的独立增量

同一个current trial/test V内，A17与Galerkin矩阵仅相差QR的VᴴV与I舍入；最大场相对差与白化响应差见tables/green_control/RESULTS.json。普通field ROM享有相同已付费R0V增广权利，而且本例响应误差更小。这是有限state/current algebra控制，不是原A17论文全部算法的重现实验，不影响原论文其他局部GN结论。

裁决：局部response复用有准确性证据，普通ROM是强对照；本轮没有A17独立增量、没有standalone全成本或完整posterior节省20%的合格主张。保留Full/ordinary-ROM与全部setup记录，停止把Schur改写当成加速贡献。
