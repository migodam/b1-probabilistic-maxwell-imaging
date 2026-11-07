# 三个大轨迹分叉：0001、0004、0022

所有数值来自冻结168结果、原生保存stage、原Phi所选corrector起点及缓存GT current。没有新增inverse、Maxwell forward、阈值或预算。完整七臂、三时点及保存梯度见`THREE_CASE_NATIVE_CORRECTOR_DETAILS.csv`；全24数据仍全部保留。Fresh Twofold实际r=30。

## 0001 high annulus：矩阵分母的真实大正例

| 方法 | 原生stage3 NRMSE | corrector前NRMSE | 最终NRMSE / fixedIoU | 最终物理current误差平方能量 | 最终Phi |
|---|---:|---:|---|---:|---:|
| CIE | .30746 | .30746 | .42288 / .35333 | 1795.94 | .00938538 |
| A+global | .26957 | .26957 | .42275 / .35570 | 1797.37 | .00938778 |
| Twofold30 | .13610 | .13610 | .10248 / .71329 | 146.59 | .00005250 |
| G19 | .18953 | .18953 | .09889 / .75556 | 134.99 | .00005255 |

各臂均从原Phi选中的stage3开始corrector。矩阵路线在原生末stage已经较接近GT，继承corrector再改善；CIE/A+的corrector则降低其自身Phi但材料误差增大。CIE以投影梯度准则停止，保存的projected-gradient∞约9.68e−9，仍未拟合到CN99（whitened residual328.26，对应门25.17）；G19/Twofold通过CN99。这是已有可行更低Phi端点的**wrong-data轨迹/停止点**证据，不能以“CIE已经数值converged”当作恢复成功。

G19/Twofold的corrector都达到3000迭代上限，故不能宣称它们完全converged；其改进也包含corrector贡献。DPCS rank13、末APCS秩563、共同GT affine floor .07933，无法解释CIE1795.94对G19134.99的差异为额外表示维度。G19末joint物理误差515.29、aux518.69、gap3.60、cross−.195；这不是将aux误差误当physical误差后制造的正例。

## 0004 high slotted shell：CIE成功，全部A系走向坏数据拟合

| 方法 | 原生stage3 NRMSE | 公开Phi选择的corrector起点 | 最终NRMSE / fixedIoU | 最终物理current误差平方能量 | 最终Phi |
|---|---:|---|---|---:|---:|
| CIE | .25694 | stage3，.25694 | .11109 / .75943 | 151.12 | .00006020 |
| A+global | .55343 | 较早stage，.51219 | .51514 / .08696 | 2670.06 | .17393422 |
| Twofold30 | .54232 | stage1，.50788 | .51514 / .08696 | 2670.06 | .17393422 |
| G19 | .56771 | stage1，.50778 | .51514 / .08696 | 2670.05 | .17393422 |

这里“原生末stage”不是matrix corrector的起点；原Phi保留了较早stage1。G19/Twofold以projected-gradient准则在482/580迭代停止，保存PG∞分别约8.91e−9/6.89e−9；允许的3000上限没有用尽。它们的whitened data residual约1415.67，远大于25.17，而CIE为23.24。因此不能把失败简单归为共同预算太短，或noise-compatible材料偏差；存在CIE的合法、低Phi端点，支持不同原生路线进入不同坏数据停止域的解释。它仍不是全局盆地边界/唯一性证明。

共同DPCS rank13、APCS563、floor .02469。G19末joint aux3286.73、physical4291.65、gap1910.15、cross−2915.07；大负cross说明aux看起来较小不代表物理current误差小。末族oracle容量足以把GT current近似得远好于这些误差，故没有“静态APCS floor迫使它失败”的证据；这不等于保证材料/物理一致性/硬lock共同约束下能找到该oracle。

## 0022 high bridge：末stage偏离，corrector也不能补救

| 方法 | 原生stage3 NRMSE | 公开Phi选择的corrector起点 | 最终NRMSE / fixedIoU | 最终物理current误差平方能量 | 最终Phi |
|---|---:|---|---|---:|---:|
| CIE | .36196 | stage3，.36196 | .09741 / .71774 | 153.00 | .00005462 |
| A+global | .15956 | stage3，.15956 | .09608 / .71200 | 150.12 | .00005461 |
| Twofold30 | .39765 | stage1，.26189 | .43720 / .11364 | 3199.95 | .01038257 |
| G19 | .39511 | stage1，.26006 | .33700 / .33036 | 3830.03 | .01420888 |

CIE/A+最终数据合格；Twofold/G19均数据不合格（whitened residual345.42/404.20）。Twofold在3000迭代上限终止，保存PG∞约1.04e−6；因此只能报告给定预算失败，不能称已到稳定局部minimum。G19在236迭代以relative-F准则停止，PG∞仍约3.18e−5；这个success flag也不是一阶收敛到1e−8的证明。两种停止方式都必须保留，不能一律归因“已converged但错盆地”或一律归因“预算已耗尽”。

DPCS rank12、APCS564、floor .03695。G19末joint physical2897.33、aux2216.09、gap746.12、cross−1427.35，最终physical升至3830.03；当前和末端材料目标都不能用单个aux能量替代。两个矩阵候选在同current表示/预算合同下未达到CIE/A+已经展示的低Phi支路，是强反例；没有证据要求重新增加rank、换阈值或加预算补救本轮结论。

## 必须与两种其它失败分开

0000/0021是NRMSE在.15门附近约.0038–.0061差异造成的二值门变化，IoU仍通过；不是上述三个大分叉。另有11个CIE对象虽数据通过CN99却材料失败，GT与估计都在同一个单case noise球中，是有限噪声允许的材料偏差实例；不能称wrong-data basin。具体见`CIE_GT_CN99_FEASIBLE_PAIRS.csv`。

原生材料规则、有限inner优化、Phi选择及共同corrector共同产生这些观察轨迹。现有保存状态足以反对“普遍改善”“只有额外current维度”“所有失败只是多迭代即可”之类强归因，但不足以给每例一个唯一机制解释或全局landscape定理。root之后新增right修正和conditional实验未混入此冻结比较。
