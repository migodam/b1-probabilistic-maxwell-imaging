# 复现与费用

代码在src/lxz，冻结协议在configs。旧P2的COMMON_RESPONSE_BANK与DEVELOPMENT_CASES、ADAPTIVE_PHYSICAL_REFERENCE、PILOT_REFERENCE在既有仓库。C1权重使用c1-joint-p2@83056881的p2/results/models/full_distill-26101321.pt，原Release joint-p2-20261014可取；不得改成P1二维模型。

零物理重新出报告：PYTHONPATH=src <numpy/scipy环境python> tools/report_lxz.py；图库：tools/plot_lxz.py。只读已有rawchains/results与旧参考，不做Maxwell。

proxy.py controls为本轮必要的finite-state详细平衡检查。proxy.py run/--round2会调用proxy，完成目录禁止静默重跑；repair参数只用于已记录的索引修正。它不运行fullMie。

local2d.py、derivative_refine.py运行真实Maxwell，需要原A30/a26 source和背景缓存。用户本轮96states已全部消费，RUN_CLOSED与ledger阻止任何新增物理；勿删除/换runID重置额度。复核结果优先读NPZ，需重跑应先取得新授权。

首次失败35state原始endpoint数组缺失不补造。历史执行源码926f847/87c7680及全部异常、修复revision保留在Git和日志，rawgap不影响费用计数，不把这些状态纳入成功候选比较。
