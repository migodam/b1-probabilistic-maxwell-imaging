# 复现七件 R14 最小峰值对照

从本目录使用同一工作树运行。需要 Python、NumPy、SciPy、python-flint、threadpoolctl 和 matplotlib。复跑请用原结果的单独副本，保留本轮的壁钟费用和封存记录。

```sh
python -B spatial.py
python -B minpeak.py
python -B summarize.py
```

- 输入来自 `../../z_alg_optical_inertia/validation/max_struct/public_inputs.npz`，包含同一个 24² 盘脉冲算子和七份完整的复端口观测。核对使用原 R14 工作树 `../../delegated/max_struct_r1/r14_worktree/research/research_max_struct_r14` 中的 `data_public/independent24`。选历史初值还读取原 `RESULTS/independent24_v1` 的 `primary.npz`、`WRAPPER_PROVENANCE.json` 和 `RECEIPT.json`。
- `minpeak.py` 从同一观测的历史端点选择初值。在原 [0.02,20] 材料盒内逐步压低峰值上限并用精确切向的 Gauss–Newton 求解，每次接受检查调用完整离散 Maxwell 正演。最后原生区间审计通过才登记为可行点。回执含全部峰值尝试、停止原因、费用和工作量。可行数组含 576 格材料、实际 48×12 复端口及初值。
- `spatial.py` 读取已发表的 `verified.json` 和 `numerical_grid.json`，用完全相同的 G、B、S、Y 和噪声半径做 8×8 更新。每件数值择权约一秒。按区间计算接受提案，独立记录复核时间。
- `summarize.py` 先核对同批最强标量提案，封存在线输出，再读取原 `offline_scalar_scores.json` 中的区间标量评分。工作流程不需要读取真值材料数组。
- 局部求解失败只能调整后续提案。严格下界由独立的惯性证书给出。两端及生成峰值的机器可读字段见 [summary.json](summary.json)。
