# 复现入口与文件恢复

环境使用现有 Apple Silicon CPU 的 Gaussian/.venv_nn；从本工作树根运行解释器 `../../.venv_nn/bin/python`，不要把其符号链接解析成系统Python。依赖numpy、scipy、matplotlib；未新装CUDA或框架。

完整原始证据按固定路径保存在Release的五个tar.zst分片，分别为P1、N1及三份N2。下载到隔离副本，逐片解压到仓库根；归档内部已带research/p3l_r12_bilateral_fixedpoint相对路径。reading kit另含必要算子、原prior、复用锚点/真值、源码、紧凑CLOUD、评分与图。阅读包不是全部原始Jacobians。

默认复现只读取已经发生的结果，不支付新Maxwell：

```
../../.venv_nn/bin/python -m research.p3l_r12_bilateral_fixedpoint.shared.portable_score --dataset P1 --out /tmp/r12_p1_offline
../../.venv_nn/bin/python -m research.p3l_r12_bilateral_fixedpoint.shared.portable_score --dataset N2 --out /tmp/r12_n2_offline
```

脚本在指定新目录写结果，不更改冻结预测。N原盲读回执保留；重复评分不是新的盲测或独立确认。

物理入口为shared/prepare_campaign.py、run_development.py、run_digits.py（N1）与run_digits_n2.py。它们的中央预算已STOP，不能通过重置账本续用本轮额度。若下一次另行授权重跑，应在新的隔离实验目录重新冻结时限与额度，并使用inputs/reuse中的相同观测/锚点；历史兄弟目录的绝对路径在本地manifest中保留为溯源，portable_score将其解析到复用输入副本。

所有实际状态在results/*/states，完整J与双边L/R都随付费状态保存。root COST_LEDGER是原始追加中央账本；RESOURCE_SUMMARY是最终聚合，standalone COMMON分摊不能再累加为实际物理费用。

shared/finalize_results.py、agent_c/offline_score.py、agent_c_max/plot_n2_results.py只做已有数据的统计和绘图。site builder只渲染文档，不调用求解器。
