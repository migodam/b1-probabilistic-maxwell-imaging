# R12 高对比度 true3D 补充

两个独立 seeded 非球形三维对象：chi8，以及混合 chi2/8；再增加第一个对象的8³正演→6³反演配对。配对不是第三个独立对象。seed2026111213，在随机几何和新观察前冻结。

几何是带空腔的随机长方体并集/优先级材料分区；6³和8³均通过连续分区的解析体积交计算真实 cell 平均材料，未使用 cell 中心标签近似。12个最终预测及 CPU前缀选择全部封存后才离线读取truth。

每个逆问题216个共享实材料、648个矢量电流分量、60个原 R8 电偶极照明、180接收通道。统一初始chi=.25，box[0,8]，保留原source-normalized物理Phi和原quadraticprior。

CIE/Dual保留原R8 current_PR、material、physical与prepare，120cycles、inner12、mf1→2；Dual晚期currentbeta3。AP/DR均3000步、audit50、scale6/omega1、groupedports；projection320CPU秒上限含setup。共同corrector至多1000迭代/60CPU秒。各方案的proposal费用不同，不声称总费用相等。

R8 projected-IFFT basis默认经过SVD range白化；这里完整披露该坐标选择，不能称为逐字2019论文复现或原始gamma坐标。物理值缓存只复用已付费LU得到的场，不改变方程。共享R8[0,4]优化器未改写，本补充独立wrapper向原优化步骤提供[0,8]。

|对象|方法|材料相对L2|proposal材料相对L2|原Phi|source-normalized residual|CPU秒|corrector停止|
|---|---|---:|---:|---:|---:|---:|---|
|random_high8|CIE_sourcewise|0.129785|0.457843|0.00072734924|0.0116519|4.27|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8|Dual_sourcewise|0.129783|0.454820|0.00072734924|0.0116518|4.40|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8|DR|0.129783|0.061631|0.00072734924|0.0116519|38.06|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8|AP|0.129786|0.103885|0.00072734924|0.0116519|44.01|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_mixed2_8|CIE_sourcewise|0.227237|0.618636|0.00037537091|0.0114161|7.61|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_mixed2_8|Dual_sourcewise|0.227232|0.608627|0.00037537091|0.011416|7.54|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_mixed2_8|DR|0.227239|0.078315|0.00037537091|0.011416|45.80|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_mixed2_8|AP|0.227260|0.216416|0.00037537091|0.0114161|42.86|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8_crossgrid8|CIE_sourcewise|0.432218|0.397404|0.0034159219|0.0464793|7.78|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8_crossgrid8|Dual_sourcewise|0.432220|0.395368|0.0034159219|0.0464793|8.25|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8_crossgrid8|DR|0.434481|0.491405|0.0034150704|0.0464126|37.88|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|random_high8_crossgrid8|AP|0.432217|0.424022|0.0034159219|0.0464794|28.39|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|

两独立 matched 对象的平均材料相对L2（不混入cross-grid配对）：

- CIE_sourcewise: 0.178511
- Dual_sourcewise: 0.178508
- DR: 0.178511
- AP: 0.178523

本补充的明确负结果：两个 matched 对象经过共同 corrector 后，四臂最终材料误差与原Phi均非常接近，DR没有最终材料收益且总CPU更高。DR在proposal阶段的材料相对L2为0.061631/0.078315，共同corrector进一步降低原Phi后却变为0.129783/0.227239。cross-grid中DR原Phi略低，材料相对L2反而略高（0.434481 vs 约0.43222）。这三个观察均保留，不能只报告proposal材料优势或更低Phi。

该小规模高对比度离散挑战不是连续 Maxwell 精度证明。chi8对应epsilon_r9；最大内部 cell 相位约1.885（6³）/1.414（8³），网格很粗。cross-grid实际 clean数据差异与原始CN99噪声球分别保留；只用sigma定义的球没有包含模型离散误差。

|对象|8³/6³ clean相对差异|coarse truth最大chi|forward最大chi|
|---|---:|---:|---:|
|random_high8|0|8|8|
|random_mixed2_8|0|7.80691|7.80691|
|random_high8_crossgrid8|0.261567|8|8|

OFFLINE_PREFIX_SCORES 使用先前已按原Phi冻结的10/30/60/120CPU秒候选，没有用GT排名。相对函数容差终止不是局部最优或唯一性证书；低物理残差与材料恢复单独列出。全部负结果、原始完整projection audits、成本和corrector停止原因保留。

Solver CPU合计 276.860 秒，physicalLU 3913，sourceRHS 234780，adjointRHS 169260；projection H分解 6，setupRHS 1080，iterationRHS 1080000。详见COST_LEDGER，嵌套timers/counters不重复相加。

复现：generate_high.py（仅新目录）；run_high.py run；整组seal后score_high.py。无训练、无在线truth、无GPU占用、无push。科学主张由主线程最终评审。
