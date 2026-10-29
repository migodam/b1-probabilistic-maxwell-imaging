# 命题、实现和证据状态

|命题|数学条件|本轮检查|裁决/限制|
|---|---|---|---|
|零阶介电微结构趋向均值|有界对比度、零均值、统一 resolvent、合适载荷和观测|两个真实宏观墙、三尺度，界面被 P3 网格解析|PROVED_UNDER_ASSUMPTIONS 为附件条件推导；NUMERICALLY_SUPPORTED 为有限采样。未估计统一常数，未证明渐近率|
|目标稳定界|必须另有宏观参数/目标的有限分离裕度|known-wall 和 physical finite reference 目标恢复；宏观模型在线失败|普遍有限下界 NOT_ESTABLISHED，不能把假设 c_X 当测量结果|
|镜像横向保护|合法墙扰动和仪器 nuisance 偶对称，居中目标，镜像白化|局部信息保留 .999996，破坏后下降；FD通过|NUMERICALLY_SUPPORTED_LOCAL。网格非完全镜像，泄漏不为机器零；不保护深度/材料|
|满射墙导数诱导局部补偿|合法内点上真实导数满射，去除冗余观测|A0/A1 浮点满行；合法有限全波曲线，独立网格|IFT条件数学成立；未证明连续满射。实际有限离散见证 NUMERICALLY_SUPPORTED|
|新增嵌套观测信息不减|同 nuisance 和旧数据、相容噪声/先验|同噪声切片，数据Schur/finite reference；旧补偿在新数据下失败但可重优化|不保证每个噪声个例或 misspecified 模型算法误差下降|
|二阶补偿|必须保留 b_i/b_ij、几何导数和全部合法核方向|dense mixed-RHS 核对；两个一阶补偿的真实二阶 secant|不是核空间全优化，也不是全局排除证书|
|C0+H / T0|一致离散坐标、正确 bX、保留 d/a/B|2墙4目标，exact targetself + reduced return，完整方程核验|NUMERICALLY_SUPPORTED，同源标准 Schur/T-matrix；不是新代数理论|
|平面端口指数尾|共同迹范数、已知均匀间隔、反射算子界；有限窗尾单独控制|未用周期计数认证开放 ROI；仅计算 FEM actions|条件理论，不是本轮全 ROI 秩定理|
|返回误差到 likelihood|全路径稳定性及二阶余项|cotangent dense FD；实际 field/response error|物理统一误差界 NOT_ESTABLISHED，不由小 core condition 自动推出|
|目标 contrast|保留增量与背景残差投影、正确实白化|独立 dense 恒等式与已付物理表精确差分|公式接线 NUMERICALLY_VERIFIED；不是免费可用的分组证书|
|条件混合收缩|所有背景的条件分布直径/组内界，相关prior正确|旧118墙直径近1；新微组小但有5个>.02尾部|普遍丢弃质量 REJECTED；有限物理组经验线索保留|
|概率早停|未查背景/目标质量或条件分歧有合法便宜包络|默认区间保守，缺少新物理 bound|CERTIFICATE_COST_FAILURE/NOT_ESTABLISHED，不启动 cold Bayes|
|实际成本优势|同质量、同数据、同模板权限、冷启动真实计费|全部setup/修正/核验收费；direct8例16RHS明确|NET_SPEEDUP_NOT_ESTABLISHED|

纯代数记录：[THEORY_CHECKS.json](results/THEORY_CHECKS.json)；关键物理导数记录：[INFORMATION.json](results/compensation/INFORMATION.json)；物理有限见证：[最终第三层数组](results/final_pair_refinement/)。

本轮没有做严格连续 PDE 误差上界、全参数盒上的稳定常数、统一随机 action 认证或一般全 prior posterior 认证。没有用随机抽样的最大值替代数学 supremum。
