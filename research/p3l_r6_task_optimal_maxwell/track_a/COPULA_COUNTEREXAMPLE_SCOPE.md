# Prior-fibre copula control：标签与执行范围

实现状态：**PREPARED_NOT_EXECUTED_BY_IMPLEMENTING_WORKER**。入口为`track_a/copula_counterexample.py`；仅做了语法编译检查，没有创建Book/Target或调用Full。由root在当前D物理job结束、singleton资源锁释放后显式执行。max528 RHS是本supplement上限，不是实际支出；原A实际157696、剩余2304的账保持未改，root执行后按实际receipt更新。

两个材料点只读`raw/REGULAR_material_GH21.npz:mode`与`raw/PURE_OVERLAP_material_GH21.npz:mode`，不读GT或D truth。joint标签是**prior fibre conditional on a fixed material**，transition用GL128离散求积条件，pure质量单独保留；它们不是完整posterior samples。所有其余field latent固定0，属性保留原Gaussianized[-4,4]截断对应的合法native区间。factorized control只置乱field/Q相对于D，仍处于原prior支持，但刻意改变材料copula。

输出将保存native q、joint/factorized材料、实际factorized Full prediction、一个共享joint reference、reference current/field、原始及σ归一化响应差、置乱和求积标签、fresh-query/cache-hit标记及实际RHS。reference cache只按decoded complex128材料精确字节匹配复用；保存的32个joint材料不被虚构为32次Full查询。没有训练、没有新模型、没有修改原A主阈值、原积分代码或既有分析；最后科学解释由root审查。
