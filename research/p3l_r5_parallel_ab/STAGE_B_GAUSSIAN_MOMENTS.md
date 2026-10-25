# B：确定性Gaussian材料矩

实现属性—场相关时的一维外积分和CDF二阶矩，核验实际provider恢复到原x后场—属性covariance为零，才使用独立特例。相关、极端tail与退化例对独立积分误差小于4e-15；高K MC分块复核。

K32 ddof=0与population variance的31/32差异保留。确定性矩+冻结head MSE=.00704106，匹配重训=.00703415；这是已有VAL256结果。生产Q64稳定，Q16快路径只在已测分布核验，不宣称任意Gaussian通用固定节点精确。

[数值与失败细节](stage_b_material_risk/REPORT_ZH.md) · [主数值验证](stage_b_material_risk/MOMENT_VALIDATION.json) · [相关CDF验证](stage_b_material_risk/JOINT_CDF_VALIDATION.json)
