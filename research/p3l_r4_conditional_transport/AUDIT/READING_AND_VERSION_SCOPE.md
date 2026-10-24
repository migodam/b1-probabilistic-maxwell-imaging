# 输入阅读、版本与权限边界

研究输入是用户2026-10-25明确授权的Research P3L-R4任务，以及附件中2026-10-23正式R3报告。主线程阅读完整主报告、执行方案、P1源码与R2主要报告；条件概率代理完整核对附录B，物理代理完整核对附录A/C及有限反馈源代码，材料风险代理核对条件矩及历史统计。PDF与LaTeX正文均真实存在，不将理论包视为已经运行的新实验。

新任务明确取代附件“尚未授权、只能Stage0”的执行限制。附件的数学和数值要求继续有效。其SHA字段仅作为原始文字保留，本轮不做SHA检查。

新工作区从远程main的c7594d4f095bb52a657bcaa08055d16bcaaafcfe建立。与用户固定6eadb4690c962aaa95821fc4056b677f0f924a6b对照，src/p1_material、src/p3r2和原TARGET配置没有内容差异。外部A30工作树HEAD为717d3d26a290248bd2b9aeb55f3d03694b407b84，检查时干净；adapter内部旧SOURCE_COMMIT指向其vendored来源，两者含义分列。

历史R2数据只读复用，新Book、授权与所有计算写入R4目录。历史R2剩余额度不属于本轮。合成矩阵测试、旧缓存方程回放、新Full Maxwell端点、监督风险动作与FULL posterior资格分别记录。

原P1是24x24二维scalar-TM/VIE，绝不外推为三维vector Maxwell。完整概率定义在Gaussianized66维latent参考测度上，经原饱和decoder推送；材料空间可能有退化层/原子，不补乘材料Gram行列式。

原P1/F3/LXZ/P2/P3-L/R2负结果保留。R3中“full feedback代数成立”“网络population目标正确”和“下一轮性能建议”均不等于R4实测性能。
