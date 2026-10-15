# 复核与复现

本轮已封账。不要删除 RUN_CLOSED.json、COST_LEDGER.jsonl 或换目录伪装成未消费额度。

1. 阅读 START_HERE.md、PROTOCOL_FROZEN.json、SOURCE_MANIFEST.json 和 EXECUTION_DEVIATIONS.json。
2. 无新增物理复核：使用同仓库 `src/prior_constrained/prior.py` 重放保留的 TRAIN 或 fresh generator 参数，读取现有 `E0_AUDIT.csv`、`CACHED_PAIR_SCREEN.csv`、`E1_SEARCH.csv`、`PAIR_NOISE_REPLAY.csv` 和对应 NPZ。`ell=norm(pred-y)^2/sigma^2`，实白化等价于 sqrt(2)/sigma。
3. 无物理重画：在仓库根运行 `PYTHONPATH=src python tools/plot_prior_constrained.py`。依赖 NumPy、SciPy、Matplotlib，使用既有 Gaussian/.venv_nn 环境，无新安装、无 GPU。
4. 实际物理入口为 `PYTHONPATH=src python -m prior_constrained.witness_search`，其完成裁决及封账检查会拒绝再次执行。需要从新明确授权和新冻结协议开始才能复现实验，不以删除记录绕过。原 A30 adapter、background S/G/incident 与 N3/历史缓存来源见 manifest。

失败的两次调用日志与数组原样保留：E1_execution.log、E1_completion.log、offline_only/E1/、offline_only/E1_completion/；完成执行位于 E1_completed_execution.log、offline_only/E1_execution_2/。各执行共享同一成本账本。实际目录名以 E1_SEARCH.csv 的精确路径为准，不依靠手工编号猜测。

E1原始端点与参数属于 offline-only benchmark construction；未来 online 只可获得 y、sigma、acquisition、公开 background、N3 与独立 TRAIN 材料 reservoir。不能给它读取这里的隐藏端点、生成参数、truth identity 或 clean means。
