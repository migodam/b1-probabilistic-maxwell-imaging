# WQ VIE 离散选项：实现和负结果

本子任务只改动 `src/wq/vie_upgrade.py`，并在本目录留下非求解积分审计。没有由本 worker 运行 Maxwell 状态、LU 或 RHS。父线程统一运行并计费的圆柱实验已经表明：下面的双重平均选项在已测网格上比旧方案更差，不能作为通过资格的修复。

## 已实现的选项

`VIEUpgrade` 沿用 `p0.vie.VIE` 的单位线源、`exp(-iωt)`、`g=iH₀⁽¹⁾/4`、`k²` 源积分、解析出射边界、材料范围、源优先数据顺序和原残差阈值。默认 `field_average_axes=(False,False)`、`field_mode='full'`，完整场和 Jacobian 直接复用旧求解路径。`quadrature` 和 `source_cell_average` 保持原含义。

- `field_mode='active'` 只保留有对比单元的总场，跳过全域空气场重建。完整接收数据与散射数据仍会计算。未重建的 `state.field` 项明确写为复 NaN；`state.field_mode` 和元数据会说明这一点。此模式不允许全域材料 Jacobian。
- `field_average_axes=(True,True)` 将检验函数改为单元面积平均，源积分仍采用常数 pulse 基函数。也提供仅 x 或仅 y 平均。必须明确启用 `source_cell_average=True`，使入射场采用相同检验平均。平面波入射若由调用者覆盖，仍须乘相应方向的 sinc 因子。
- `field_quadrature` 控制非自项的检验平均。双轴平均的自项独立用源/观测方格重叠的三角权重积分；不会把粗观测采样的对角误差误认成物理误差。非自项建立成本约正比于 `quadrature² × field_quadrature²`，不建议无根据地全局提高到很高阶。
- `solve_subcells(material)` 显式积分 `χ(y)g(x−y)`，总场在每个母单元内保持常数。`frozen_subcell_material` 调用旧 `wall_values`、`in_query`、相同圆盘 probe 覆盖和 `ε∞+iσ/f`，并保留旧 `material_grid` 的单元均值数组。微单元自身或边角的奇异积分用真实矩形积分替换。该路径未运行物理求解，Jacobian 明确不提供。

典型调用：

```python
from wq.vie_upgrade import VIEUpgrade, frozen_subcell_material

p = VIEUpgrade(96, 1., tx, rx, book,
               source_cell_average=True, field_mode='active')
old_discrete_state = p.solve(epsilon)

g = VIEUpgrade(96, 1., tx, rx, book,
               source_cell_average=True, field_mode='active',
               field_average_axes=(True, True), field_quadrature=4)
# 下面的调用属于物理状态，必须由统一预算的父线程安排：
# averaged_state = g.solve(epsilon)
# material = frozen_subcell_material(p, seed, subsamples=4, probe=probe)
# weighted_state = p.solve_subcells(material)
```

## 为什么这仍不是高阶场离散

平均检验改变了 pulse/test 配对，但没有增加母单元内的场自由度。混合单元仍把材料均值与常数场结合；它无法恢复界面附近的真实场变化。微单元加权方案纠正源积分中的材料位置，却同样保留常数总场。因此，两者都不能称为高阶点 Nyström 或局部双线性场求解。

本轮没有实现每母单元 2×2 独立场点的高阶 Nyström，也没有实现跨单元双线性场重建。前者会增加未知量和直接 LU 成本；后者需要明确的空气侧场自由度、边界插值和相应导数合同。不能把仅提高 Green 积分阶数当作这些方法。

## 由父线程运行的双重平均负结果

原始证据在 `results/galerkin-cylinder/RESULTS.json` 和同目录复数据数组。以下为相同圆柱、入射平均和单位归一化下的相对全数据误差。

| nx | 原方案（父线程提供） | 双重平均 | 双重平均白化误差 |
|---:|---:|---:|---:|
| 64 | 0.2261 | 0.4025081 | 40.26586 |
| 96 | 0.0984519 | 0.1839315 | 18.40003 |
| 128 | 0.07004 | 0.1182396 | 11.82838 |

双重平均的线性残差约 `8.05e-16` 至 `1.55e-15`；这说明离散系统求解精确，不能解释为对连续 Maxwell 参考准确。该选项在这些实验上变差，保留为负证据，不作为资格通过依据。

父线程在 nx=96 上比较旧离散的源积分阶数 q4/q8/q16，圆柱相对误差约为 `0.09845190 / 0.09843743 / 0.098437431`。这一受控圆柱实验中，提高非自项积分阶数没有消除主要误差；它不构成对所有墙场和所有频率的误差分解。

## 非求解积分及几何审计

`NON_SOLVE_AUDIT.json` 保存所有数值；`audit_integration.py` 可复现，不会创建 VIE 后端或调用 solve。

| 积分检查 | 相对差 |
|---|---:|
| 角点积分 = 双边长方格自项的四分之一 | 3.70e-16 |
| 非中心自项 vs 独立 Duffy192 | 1.27e-10 |
| 4×4 微单元划分之和 = 母方格自项 | 9.31e-16 |
| 双方格平均自项 vs 独立 Duffy192 | 5.35e-10 |
| 左右不同对比的加权积分 vs 分半矩形参考 | 9.01e-17 |

最后一个检查中，`mean(χ)×∫g` 捷径相对参考差约 `0.05859`；实际加权实现没有采用该捷径。Duffy192 是数值积分参考，不是 PDE 误差界。

几何审计完整使用旧 nx=96 网格、两个旧 seed、cal 和旧 `probes()[8]` 的 task，逐项检查 subsamples 4/8/16/32。下表是相邻采样等级之间的面积加权系数 L2 差，未进行数据求解：

| seed / 情形 | 4→8 | 8→16 | 16→32 |
|---|---:|---:|---:|
| 2026102201 / cal | 0.1023104 | 0.0370344 | 0.0125005 |
| 2026102201 / task | 0.1031645 | 0.0374601 | 0.0128763 |
| 2026102202 / cal | 0.0801909 | 0.0313178 | 0.0132601 |
| 2026102202 / task | 0.0812777 | 0.0318201 | 0.0136150 |

这些系数差不能代替白化接收数据差。额外在 nx=32 上检查 `frozen_subcell_material` 的 4/8/16/32 四级采样；其 `epsilon` 数组与旧 `material_grid` 全部逐元素精确相同。子样本对比均值的最大舍入差为 `1.03e-14`。

## 资源和执行边界

非求解审计区域共测得 **3.200507 CPU 秒**：积分 0.085503 秒、nx=96 几何 2.750420 秒、子样本恒等检查 0.364574 秒。首次运行未记录解释器启动和导入 CPU，JSON 已注明计时范围；runner 已改为在今后运行时额外记录原始进程 CPU。该 worker 使用物理 RHS=0、物理 LU=0。

加权微单元方案只有一个真实 n×n 直接矩阵，并按行小块组装。缓存的平移核有独立 192 MiB 上限，超过就失败，不会默默减少采样。nx=192、subsamples=8 的核缓存约 107.3 MiB；subsamples=16 超限。直接 LU 保留旧 4500 活跃单元硬上限。实际峰值内存尚未测量，因此 1.5 GB 目标仍需由正式状态测量确认。

所有实际 LU 和 RHS 均先预约同一 Book。此文件不提高预算、放松 `1e-10` 残差或 `0.1` 白化网格阈值，不启动候选搜索，不修改旧世界，也不声称完成求解器资格。后续是否选择更高阶场表示及是否通过科学 gate，由父 Codex 线程判断。
