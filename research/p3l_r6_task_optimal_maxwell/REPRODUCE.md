# R6 的复现入口与原始证据

本目录是独立 R6 实验；历史 R4、R5、R5-B 的结果及账本只读。报告中的统计值应从已保存的预测和收费收据重算，不必重新求解 Maxwell。关闭后的中央控制器拒绝新增物理工作，也不会自动删除目录、重置额度或覆盖失败任务。

## 代码与必要输入

物理与 A/C/D 入口在本研究目录。B 的可执行 Python 包位于本仓库根目录 `track_b/`，输出写入本目录 `track_b/`；二者均为本轮新增资产，证据包必须同时保留。这样保留实际运行代码，避免事后搬迁路径而声称原实验使用了另一份代码。

原始输入由 P0 资产清单记录：R4 的 P1 prior、A30 背景与原优化器；R5 冻结 Gaussian provider、三 seed local heads、214k provider；32768 joint TRAIN、Gaussian teacher 特征、1024 expensive TRAIN context 与256 VAL。它们不是本轮新增 Maxwell 状态。依赖包只归档实际使用的输入，不包含整个历史仓库、Python虚拟环境或历史 LU 缓存。

原目录关系为同一 `Gaussian/B1/` 下的 `b1-p3l-r4-conditional-transport`、`b1-p3l-r5-parallel-ab`、`b1-p3l-r5-b-confirm`、`b1-p3l-r6-task-optimal`，以及 `Gaussian/A30/a30-maxwell-material-manifold`。源码使用相对 sibling 定位和少量冻结的绝对身份记录。异机恢复应先按 release manifest 恢复这些关系，再明确重定位身份清单；不能把路径改变误报为同一原始身份验证通过。

## 已存证据的复查

在仓库根使用既有 `../../.venv_nn/bin/python`。原环境版本见 `audit/ENVIRONMENT.json`。所有 GPU 训练均为本机 MPS；物理为 CPU、单 worker、BLAS单线程。

```sh
# A：只读已经付费的求积数组与复响应；不会新增 RHS
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/track_a/analysis/posthoc_latent_q.py

# B：已有模型/预测的诊断与报告；不是新训练或新数据
../../.venv_nn/bin/python track_b/diagnostics.py
../../.venv_nn/bin/python track_b/report.py
```

Track D 的主 scorer 为 `track_d/score.py`：先核对全512预测封存身份和独立 GT_UNLOCK，再打开真值。若评分目录已存在，不应擅自覆盖；可在隔离副本中按同一 bootstrap seed 重算，原收据保留。补充 direct-prior 对照有独立封存及探索性统计，不能混入原全局8项比较族。

## 新物理计算与训练

实际运行入口包括 `track_a/stratified.py`、`track_a/local_material_reference.py`、`track_a/compensation_vector_audit.py`、`track_a/copula_counterexample.py`；`track_b/run.py`；`track_c/anchors.py`、匹配训练与 direct-prior 控制；`track_d/online.py` 的 generate/anchors/predict 分进程步骤。每项配置与真实命令在相应协议、日志和原始 job 中保留。

这些入口是研究源码，不是“忽略旧文件并从头重跑”的按钮。新 campaign 需要新授权、新目录、新父对象或明确重放身份，并保留旧账本。禁止删掉 `CLOSURE.json` 后在原实验上重新收费，也禁止将读出的历史缓存重新算成本轮新求解。

## 精度与资格

原始数组采用无损 NPZ/zstd 归档，保留 complex128 Maxwell 输出、拒绝/失败、模型权重、预测、原始计时和 bootstrap 数据。逐字节解压比对只证明归档完整，不证明 posterior 资格、连续 PDE 误差或算法创新。所有 PASS/HOLD/FAIL/NOT_RUN 和局部 Q1/Q2 限制均应随证据保留。
