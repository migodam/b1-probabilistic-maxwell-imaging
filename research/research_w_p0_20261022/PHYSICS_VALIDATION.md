# Physics validation

## 证据分层

1. 标准数学：纤维、四块resolvent、条件性TV界经审查，不作新颖性证明。
2. 一维独立Maxwell：行波矩阵、Fresnel、PEC边界与完整级联核对；[原始结果](results/p0a/results.json)。Gate0 PASS。
3. 二维实际离散全波：两seed、九个不同频率、4源12收、cal及固定目标，完整LU多次散射。[原版](results/qualification/WALL_GRIDS.json)、[匹配源平均修订](results/qualification_source_average/WALL_GRIDS.json)。
4. 求导/伴随：四步长最佳相对差3.50982e-10；real-adjoint缺陷1.91186e-15。这只验证离散接口。
5. 互易性：源中心/接收积分原版3.48829e−7；匹配平均后7.20389e−16。[实际两种实现对比](results/qualification/RECIPROCITY_REPAIR.json)。
6. 独立圆柱解析解：修订版nx128仍7.0036%相对散射误差；非充分精度。[数据](results/qualification_source_average/CYLINDER.json)。
7. 独立网格：最细白化差0.950860–1.229102>0.1。Gate1 UNQUALIFIED。
8. 原位LU内存修复：同材料数据差0；最大RSS从2.0545GB降到0.5723GB。[parity](results/qualification/ASSEMBLY_PARITY.json)。
9. 8个任务受限四块离散核验：相对误差最大8.6675e-16，Born诱导误差100.55%–337.49%。没有执行全T签名或证明连续精度。[原始数组](results/four_blocks/ARRAYS.npz)。

## 开放边界与范围

出射VIE的外部是无限均匀介质，采用解析Sommerfeld辐射条件；没有人工PML/盒边界，故PML增厚为NOT_APPLICABLE，而非PASS。网格增密检查的是内部系数和场离散。没有连续PDE误差证书；没有把经验网格差当上界。

均匀真空χ=0时，前向直接为解析自由空间线源，没有数值LU；零对比、打包/白化、跨频材料共享、真实离散伴随及完整多次反馈均有对应代码/日志。没有额外收发孔径、增益、频率、目标位移或隐含目标选择。

独立平板检查由P0-A精确Fresnel提供，不冒充有限二维板边缘参考；二维额外采用均匀圆柱Mie作不同物理模型检验。通过同一离散系统的多条解法不能替代这个外部参考。

## 停止

本轮没有任何候选获GO，没有最终见证细网格审核、轻微扰动、未知目标joint补偿或posterior实验。独立细化失败导致这些环节NOT_RUN。20项标准库测试通过，保留一次pytest依赖启动失败。后续必须先获得足够精确的前向平台，不能靠提高噪声或放宽0.1条件继续。
