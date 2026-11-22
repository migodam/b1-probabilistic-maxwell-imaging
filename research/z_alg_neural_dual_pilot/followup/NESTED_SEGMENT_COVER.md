# 固定权重线段的严格覆盖

上限为 0.75。权重线段连接已证明正定的 Paper 1 16 区端点 0 和 2。端点具有同一个物理尺度。交接协议记录在 `SEGMENT_COVER_FROZEN.json`，先冻结节点、等分方式和判据，再读各批观测。本试验只使用带噪端口读数及原半径。

## 核对一个观测的整条线段

记 $H(t)$ 为权重线段上的原生 10×10 源矩阵。正定的状态矩阵允许直接消去电流。由于消去前的二次式是关于权重的仿射函数，对电流取下确界给出 $H(t)$ 的 Loewner 凹性。

令 $a<b$ 为两个已取样权重点。对 $x\le a$，把 $a$ 写成 $x$ 和 $b$ 的凸组合，就有

$$
H(x)\preceq H(a)+\frac{x-a}{b-a}\bigl(H(b)-H(a)\bigr).
$$

$x\ge b$ 由对称推导得到同样的外推上包络。相邻节点之间用左侧最近弦和右侧最近弦。最边缘的区间只有一条外弦。每个相邻节点区间再等分成四段，只在这些细分点计算仿射外推的源矩阵，不新增原生求解。

定义 $\phi(H)=\operatorname{tr}(H_+)$，为正特征值之和。该函数对 Loewner 序单调且凸。因此，每条仿射上包络在一个细分区间内的 $\phi$ 都以该段两个端点的最大值为上界。对可用的外弦取最小上界，再对全部细分区间取最大值，覆盖所有实数 $t\in[0,1]$。原生证书距离是 $\sqrt{\phi(H(t))}$。所有谱端点都通过向外舍入复核。重复的外弦与二进制权重点共用区间结果。

节点的最大严格下端和整段严格上端夹住连续最优距离，其差是最优性间隙的上界。节点严格超过半径时记录 `PASS`。整段严格上端低于半径时记录 `BELOW`。其余记录 `UNRESOLVED`。`BELOW` 是这对固定权重及其线段的结论，材料可行性需另行判断。

比较对象包括最佳端点、固定精确中点，以及对每个相同预计算预算取用的最佳节点表。内部节点的下端严格超过两个端点上端，才记录优势。仅查节点表可证明过线。连续外推另能证明整段在半径以内。

## 数据、费用和复现

- 3 个节点复用前一轮已保存的 90 个原生右端。5 个节点新增 60 个右端，9 个节点相对 3 点新增 180 个。每个新增节点使用同一个 96² 算子和原区间驻点核验，Arb 内对冻结系数作精确二进制插值。其正定性来自原端点证明的凸组合。
- 逐观测记录同预算节点表与加上覆盖证明的墙钟及 CPU。历史验证数据、独立 32 件材料、28 件历史压力测试材料及独立 8 件旧形态重新采集的读数各自报告。噪声视图始终按同一材料配对。

从工作区根执行下面的入口，其中 `prepare` 可单独保存 9 份区间响应。缓存存在时，计算各批次仅需源空间运算。

~~~bash
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/segment_cover.py prepare
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/segment_cover.py validation
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/segment_cover.py fresh
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/segment_cover.py historical
PYTHONDONTWRITEBYTECODE=1 .artifacts/z_alg_diag_python/bin/python -B \
  research/z_alg_neural_dual_pilot/followup/segment_cover.py external8
~~~

单个节点保存的 $Q(t)$ 可用于任意同合同的观测。`results/segment_enclosures.json` 保存严格复区间及预计算费用，四个 `results/*_segment_cover.json` 保存逐视图结果、严格改善的对象与三种预算的对照。

## 结果

| 数据批次 | 3 点 过线／整段低于／未定 | 5 点 | 9 点 |
|---|---:|---:|---:|
| 12 件验证材料，12 个视图 | 3／0／9 | 3／9／0 | 3／9／0 |
| 32 件独立材料，96 个配对视图 | 4／0／92 | 4／90／2 | 4／92／0 |
| 28 件历史材料，84 个配对视图 | 5／0／79 | 5／78／1 | 5／78／1 |
| 8 件重新采集材料，24 个配对视图 | 0／0／24 | 0／24／0 | 0／24／0 |

5 点已证明外部 8 件材料的 24 个视图在整条线段上证书恒为零。9 点在独立批次额外判定了 `ip_port_17` 的 1% 和 2% 视图，两份同样是整段严格零。唯一仍未定的历史视图是 `s2d_004` 的 2% 观测，9 点把连续最优距离夹在约 $[0.07036,0.08494]$，误差半径约为 0.07993。

`ip_port_09` 的 9 点最佳权重为 $t=0.625$。在 1%、2%、10% 噪声下，严格距离下端分别为 0.13129、0.13248、0.13996，同视图的最佳端点上端为 0.10442、0.10537、0.11061。该材料在 1%、2% 的两个视图已由端点过线。四个批次都没有新增过线。各层节点下端与整段上端的排序逐项核对。192 个已有视图的端点和中点合计 576 份区间与之前的计算均相交。

3、5、9 点分别支付 90、150、270 个原生右端。累计墙钟为 27.26、84.49、197.53 秒，进程 CPU 为 29.05、66.66、141.87 秒。新增节点准备时存在并发作业。独立批次逐视图中位 CPU，节点表分别为 4.24、6.87、12.26 毫秒，加上连续覆盖分别为 9.94、28.47、70.98 毫秒。
