# R7 数学与实验合同

材料距离为完整576格复场实Hilbert RMS，`||chi1-chi0||/24`；所有 latent 为Gaussianized66，两个原uniform属性用精确normal CDF恢复。原prior推前、sigma、source-major packing与原R4/A30代码一致。Maxwell导数包含sqrt2/sigma白化，material Gram先去数值decoder-null，再求广义谱；数值cutoff不是exact rank证明。

E0/E1 0新增RHS；E2硬上限50000。当前用户任务优先于附件旧建议：只比较原8照明的两互斥子组，不增加旋转采集；新locked的GT仅在封存后评分，在线中心由测量恢复。低维gauge不能外推Full66，三次decoder seam不能称为非零材料一阶速度的Maxwell三阶补偿。

端点noise-ball、DeltaPhi<=2、pair whitened距离<=1、合理prior启发式分别报告；不将任何一项称posterior质量。全部精确/条件理论与有限离散数值资格分开。全66维posterior不认证。

完整冻结配置见 configs/E2_FREEZE.json。
