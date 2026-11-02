# 冻结、授权和运行回执解释

R2在R1发布提交88f19920的独立detached worktree开发；没有重设旧主工作目录，也没有修改R1文件。根据本轮用户明确要求，新增结果提交现有 `migodam/b1-probabilistic-maxwell-imaging`；GitHub发布不改变科学判定。用户要求不做SHA验证，因此code_hash/config_hash未知，不用新计算哈希替代源码与参数记录。

团队复用本机CPU环境及现有operator/data；没有CUDA、driver、daemon、新训练或外部付费实验。顶层负责核心推导、路线选择与最终解释，三个既有子任务分别处理实现、文献和独立评分。模型/provider未改动，运行时没有提供可信的计费token/货币meter；不填猜测费用或model身份。

原始回执入口：

- 二维12次：`results/blind/<opaque_id>/<arm>/RECEIPT.json`，含源码、时间、prior、初始化、schedule、stage事件、停止、物理目标、残差、成本与truth_access。
- 生成：`GENERATION_RECEIPT.json`；新profile导数：`PROFILE_GRADIENT_WIRING.json`。
- 3D：`results/active/PLAN.json`、`ACQUISITION.json`、`FIT.json`、`SCORE.json`、`COVERAGE.json`；每次完整DDA另有solver.json、各进程PROCESS.json。
- 实现接线/端点：`results/diagnostics/SOURCE_WIRING.json`、`ENDPOINTS.json`、`SUMMARY.json`。
- 单次事后capacity：`results/capacity_control/<id>/RECEIPT.json`，独立score在`results/POSTHOC_CAPACITY_SCORE.json`。
- 封存与评分：`PREDICTIONS_SEALED.json`、`results/BLIND_SCORES.json`；包含首次GT访问时间。

`SEALED_PREDICTION`是运行/保存状态，不是科学PASS；成功由冻结噪声/材料/IoU合同的独立评分判定。重复嵌套的cost/history是可读冗余，不加倍收费。unknown计时、线性代数硬件操作数和历史组装费用不补零。没有把辅助state Cholesky记作物理LU，没有把线程CPU求和叫并行wall。

本轮阶段数据完全保留；sourcewise原生端点通过仅转发的_readout wrapper捕获，原算法numerics未变。第三方全文留本地，不打包；近邻论文身份/DOI与检索失败回执公开。
