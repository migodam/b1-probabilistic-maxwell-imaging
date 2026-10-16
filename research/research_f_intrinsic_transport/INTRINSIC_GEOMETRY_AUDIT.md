# F0 原先验材料几何审计

本次只读取旧 CHART/JACOBIAN 与 generator；decoder-only 求导及矩阵代数。新增 Maxwell states、LU、complex RHS、Mie、GPU 训练均为 **0**。旧 45 states、轨迹、预算账本及 RUN_CLOSED 均未修改。

## 度量与定义

材料坐标为 m(chi)=[Re chi;Im chi]/24，故欧氏范数等于 24×24 材料 RMS。白化数据按每源 Re32/Im32，乘 sqrt(2)/sigma。缓存 J_white_mass 正是这两种度量之间的线性算子。

对缓存方向 v（原 Gaussian-CDF 参数 ||v||=1），w=Dg v；SVD 去除共同 Gaussian 权重缩放 gauge 后取材料正交 Q。d2_parallel=QQ^T D²g，d2_perp=(I-QQ^T)D²g。分别计算 C_coordinate=J d2_parallel、C_prior=J d2_perp。二阶量按 ||w||² 归一化，不把参数单位与材料单位混用。

这里所称 intrinsic 项是 prior 在固定欧氏材料度量中的外在第二基本形式，不是 Riemann 内蕴曲率张量。去 gauge 不改变原始概率测度，也不创建均匀七维 prior。

## 已完成的分解

| anchor | ||w|| | ||C_coordinate|| | ||C_prior|| | ||raw|| | cosine(coord,prior) | ||Qout C_prior|| |
|---|---:|---:|---:|---:|---:|---:|
| 250033_weak | 0.000500363315 | 3.64871074 | 0.00839118281 | 3.6523145 | 0.4285317 | 0.000648559834 |
| 250034_ring | 0.01526316 | 9.2525104 | 3.50531456 | 7.52726016 | -0.6357196 | 1.43599371 |

表中数据响应使用原参数单位；完整逐步和单位材料速度平方结果见 CSV。范数不能加减代替向量分解：coordinate/prior 的 cosine 与 cancellation_ratio 已保存。

## 缺失与可用边界

**两个 anchor 的 exact Maxwell 方向 Hessian 向量均 MISSING。** 原实际执行代码将 bmax/q 留在局部变量，只向 METRICS 写 Maxwell 曲率范数。逐个检查旧 results 全部 NPZ 的 key，未找到该向量或 q；curvature endpoint 只存有限步 state，不能替代方向 Hessian。本次不调用缓存 LU、不解 RHS，也不从粗有限步拟合“exact Hessian”。

因而 C_Maxwell 与 C_coordinate/C_prior 的 cosine、三项总和、总 intrinsic 响应和总余核响应均缺失。旧 Maxwell scalar norm 仅提供三角不等式界：max(0,||Qout C_prior||-||C_Maxwell||) ≤ ||Qout(C_Maxwell+C_prior)|| ≤ ||Qout C_prior||+||C_Maxwell||；该界依赖旧 norm 的正确性，不是独立重算。

| anchor | ||C_Maxwell+C_prior|| 范数界 | ||Qout(C_Maxwell+C_prior)|| 范数界 |
|---|---:|---:|
| 250033_weak | [0.00829597379, 0.00848639182] | [0.00055335082, 0.000743768849] |
| 250034_ring | [3.32763722, 3.6829919] | [1.25831637, 1.61367106] |

上述界为参数单位下的保守区间；CSV 同时保存除以 ||w||² 的材料单位结果。这些区间并不恢复缺失向量，也不验证旧 bmax norm 的物理计算。

## 差分、gauge 与重参数化

执行前固定 h=1e-3,5e-4,2.5e-4，全部报告，不筛选最漂亮一步；另独立推导 fixed active-max 的解析 D²g 交叉验证。JSON 给出每步解析误差、相邻变化、减法消去估计、CDF 饱和、active-max cell、侧点 rank。有限点一致不构成整个邻域无边界或恒秩的严格证明。

三个步长的误差并非严格单调下降：最小步长开始受到减法消去影响，因此最终摘要使用独立解析 fixed-max D²，三步 FD 全量保留。特别是 weak anchor 的 prior normal 项很小，其误差应相对该小项评估，而不能只看 raw D² 的误差。JSON 同时保存 normal、J normal 和其余核投影相对解析结果的误差。

同 image 局部坐标变换 u=u0+z+0.5*(0.25e0)*(v^Tz)^2 在原点导数为 I；固定相同材料速度。raw D²g 增加 Dg(0.25e0)，法向项和其数据余核投影不变。已同时执行解析和三步 FD 对照；映射 Jacobian 行列式在所检查点非零。不声称概率目标在新坐标下自动均匀。

两个 prior 材料 rank 都是 7，JQ 也满列秩；当前不存在非零 exact prior tangent kernel。所查缓存 v 是 approximate weak direction，余核二阶响应只能称 fixed-velocity finite-noise 几何诊断，不能称已证实 exact fiber 障碍实例。

平坦 image 弯曲坐标控制 g(u,v)=(u+1000v²,v) 给 raw curvature=2000、prior normal curvature=0；geodesic 参数加速度消去材料伪曲率。

## 历史边界与交付

anchor 是公开且已暴露的 TRAIN 已知材料，自己的 GT 就是起点；本F0子阶段没有独立anchor或新反演；之后F1的实际结果单独见PHYSICS_ANCHOR_RESULTS。保留旧 sign-tracking bug 与已发布的轨迹重放说明：此次没有重跑旧轨迹。

geometry numerical checks 与完整科学 gate 分开：本脚本只给数值一致性结果，父线程负责 Gate1 与下一阶段裁决。

新增时间 0.076896 秒；计数 {"npz_header_reads": 153, "npz_files_read": 6, "decoder_evaluations": 40, "analytic_decoder_first_derivatives": 14, "svds": 16, "matrix_products": 240, "analytic_decoder_second_derivatives": 2}。

来源路径、完整缓存 key inventory、数值误差、条件数和所有 FAIL/MISSING 均在 GEOMETRY_TESTS.json；衍生向量在 GEOMETRY_VECTORS.npz。
