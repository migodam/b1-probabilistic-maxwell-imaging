# 复现与审计

本轮原始账本已封闭，运行器默认拒绝重跑已完成阶段。下载证据Release到现有B1根目录；不要删CLOSED或ledger来重置本轮额度。

依赖：Python3.13、NumPy2.5.3、SciPy1.18.1、Torch2.14.0、Matplotlib3.11.1。物理使用A30 commit717d3d26a290248bd2b9aeb55f3d03694b407b84；可用B1_P2_A30_ROOT指定已存在checkout。背景数组和selected-anchor caches在Release。当前源码入口tools/run_material_fiber.py，阶段A/B/C/D/E由统一Book计费。

tools/analyze_material_fiber.py为隔离GT评价；它在封账后拒绝覆盖既有主分析。报告重建入口tools/report_material_fiber.py只读已有数组，不产生Maxwell。tools/check_material_fiber_evidence.py检查已保存物理方程和测度代数，不支付新RHS。完整从零物理复算需要独立的新执行目录、明确新预算和保留原账本，而不是把本证据包当免费solver bank。

主要代码：src/material_fiber/core.py、optimize.py、explore.py、global_search.py、learn.py、acceleration.py、baseline_audit.py。所有随机种子与已暴露补充协议在configs。模型checkpoint小于1MB/个；生成过程候选/失败日志与raw endpoints保留。

Release保留所有raw材料/current/field/prediction、训练weights、原prior与背景、选中anchor的LU/receiver transfer；未交接的其余LU属于可重建计算缓存，仍在本地保存，未作为免费部署成本。任何重新分解、RHS或额外频率都应按新预算计费。
