# MAX-STRUCT-R8

本轮计算、封存和评分已完成：42个独立材料、750项主实验/初始化消融、126项固定偏置扩展；另完成原论文物理合同适配、25项真实Fresnel反演及52个小型矢量3D端点。端点数不是独立样本数。原尺度3D尝试为截断，3D神经网络与3D实测未运行。

主要结果：轻量NN接原物理优化改善了独立二维材料误差；直接接sourcewise/global A+往往损失这个收益。固定电流偏置只部分修复，最低原物理目标的双支路选择也可能使材料更差。没有证明A+在经典、实测或3D上普遍胜出。

- [中文研究报告](RESEARCH_REPORT_ZH.md) · [Chat交接](CHAT_HANDOFF.md) · [NN完整分析](NEURAL_WARM_STARTS.md)
- [数学推导](THEORY_DERIVATIONS.md) · [新颖性与先例](PRIOR_ART_AND_NOVELTY.md) · [英文论文初稿](PAPER_DRAFT.md)
- [论文对象](PAPER_BENCHMARKS.md) · [Fresnel实测](FRESNEL_REPORT.md) · [三维报告](VECTOR3D_REPORT.md)
- [原始数字与配对统计](delegated/warm_analysis/CONFIRMATION_REPORT.md) · [所有42对象图](FIGURES/warm_material_all_42.pdf)
- [成本](COST_LEDGER.md) · [复现](REPRODUCE.md) · [失败](FAILURES.md) · [完整数组分包索引](delegated/cost_package/ARCHIVE_INDEX.md)

起点为 R7 发布提交 6eb37585，独立工作树 codex/max-struct-r8-20261106。R1–R7 只读复用。Chat 工具返回该引用的全部五轮、没有更早 cursor；已将每条消息的上限增加至 20000 字符以完整读取。历史 R4–R6 解释同时从仓库原始报告补全。Chat 和 Pro 文本视为待审查来源，不视为实验事实。

主线程负责理论、NN、公共比较与最终判断；三个原生 Sol/high 子任务分别处理论文基准、实测适配、三维，机械统计/成本工作采用Sol/medium。训练实际在Apple CPU完成；物理比较使用Mac与既有i9有限CPU队列，每进程BLAS1。本轮没有使用MPS/CUDA，也未安装驱动、服务或新调度框架。

数据生成与训练允许读取自己的训练真值；新测试反演只读取公开观测及冻结权重。训练/验证/测试按独立材料划分。新测试预测封存后才评分；原 R7 对象此后仅属公开开发数据。
