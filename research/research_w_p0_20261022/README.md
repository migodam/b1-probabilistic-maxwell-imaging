# Research W P0：执行入口

**结论：P0-A PASS；P0-B UNRESOLVED / UNQUALIFIED_SOLVER。**

一维完整Maxwell反例及其带宽负结果正确复现。二维已经运行两个有限、不规则、非均匀、有损墙体的真实多频全波状态，但最细网格白化差0.951–1.229，大于冻结要求0.1。故停止四种见证搜索；calibration/joint没有合格反例，P1没有授权执行。这不是歧义不存在的证明。

- 完整中文解释：[RESEARCH_REPORT_ZH.md](RESEARCH_REPORT_ZH.md)。
- 状态与原因：[GATE_DECISION.md](GATE_DECISION.md)。
- 最重要图片：[真实几何](figures/02_geometry.png)、[复场](figures/03_first_source_field.png)、[未通过的收敛](figures/04_convergence.png)、[四块有限诊断](figures/07_fourblock_diagnostic.png)。
- 实验参数：[w_config.json](w_config.json)；[逐项指标](RESULTS.csv)；[物理账本](results/PHYSICS_LEDGER.jsonl)；[失败](results/FAILURE_LEDGER.jsonl)。
- [复现命令](REPRODUCE.md)；[原ZIP文件清单](docs/SOURCE_INVENTORY.json)；[原始理论PDF](source_package/Research_W_Report_ZH.pdf)。

物理范围仅二维TM、已声明噪声与材料类。没有posterior概率、神经网络、三维或加速结论。没有哈希/备份流程。
