# 复核与重现

主配置和模型在新物理计算前冻结：`configs/FREEZE.json`、`FROZEN_CONFIRMATION_PROTOCOL.md`、`configs/PHYSICS_RECONSTRUCTION.json`、`configs/STATISTICAL_PROTOCOL.json`。模型权重保持原 R4/R5 来源，不在本轮训练。

本地解释器为现有 `Gaussian/.venv_nn/bin/python`；NumPy 2.5.3、SciPy 1.18.1、PyTorch 2.14.0，详见 `audit/ENVIRONMENT.json`。单线程CPU是共同评估后端，不能用另一后端重算图后替代已评分数组。

## 不新增物理的复核

1. 恢复本轮证据包到隔离目录；使用原 R5 Release 中的必要复用资产和 BLAKE2b checkpoint 清单。保留 B1 下 `b1-p3l-r4-conditional-transport`、`b1-p3l-r5-parallel-ab` 和本确认工作区的相对名称。清单中的绝对路径标明原始来源；跨机器恢复须在恢复副本中明确做路径前缀映射，不能修改原封账文件或模型字节。
2. 从已封存 `raw/offline_scores/PARENT_LOSSES.npz` 重新运行 `src/statistics_confirm.py`，输出到新路径；比较风险、20,000次bootstrap与四比较校正。统计程序不读取图像GT，不进行求解或训练。
3. 所有正式预测在 `raw/predictions/`，相应内容身份和封存时间见各批及全局 seal。评分程序 `src/score_predictions.py` 先核验两批全部预测，才打开 `raw/offline_only`。如需重算评分，指定新的输出路径，保留原评分。
4. 图由 `src/plot_confirmation.py` 从实际评分包生成；固定案例和事后误差选例的身份不同，图注保留。

```bash
python -B research/p3l_r5_b_independent_confirmation/src/statistics_confirm.py \
  research/p3l_r5_b_independent_confirmation/raw/offline_scores/PARENT_LOSSES.npz \
  --output /private/tmp/r5b-confirm-statistics-replay.json
```

## 从零物理重跑

完整重跑是新的付费物理实验，不能复用已关闭账本。需要在新的授权目录中保留相同模型、parent清单和配置，建立独立账本。`generate_observations.py`、`run_anchors.py`、`predict_frozen.py`、`score_predictions.py` 分属生成、在线反演、预测、离线评分进程；不可用一次合并的GT读写流程替代它们。

原 `run_anchors.py` 直接调用冻结 R4 Native96 / Probit64，不是廉价代理。每次8源forward计8 RHS，receiver-J额外32 RHS，失败尝试不退款。新parent、line search、优化停止、噪声、边界、坐标目标与原配置一致。输出保留Native与Probit诊断，不把迭代上限称为已收敛。

原始端点、轨迹、材料与预测采用无损 `.npz`，不逐endpoint落盘密集LU。所有复杂材料风险对已存储的复数数组计算；seed平均损失与三seed图像平均的风险分别保存。

## 费用重放

`RESOURCE_LEDGER.jsonl` 是追加式实际费用账本；`RESOURCE_LEDGER_SNAPSHOT.json` 为可核对快照。原始账本与全部运行日志随证据包提供。`audit/timing/CONFIRMATION_TIMING.json` 保存独立计时重复、p10/p50/p90与费用范围；checkpoint冷加载的对象级冷启动不等于清空OS缓存的整机冷启动。历史teacher、训练数据和强context费用在 `audit/HISTORICAL_OFFLINE_COST_SCOPE.json` 单列，不与本轮新物理混算。
