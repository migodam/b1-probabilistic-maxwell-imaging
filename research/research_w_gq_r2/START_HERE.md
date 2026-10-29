**发布说明（2026-10-29）：** 此为用户授权后的公开版本。[网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r2/) · [完整原始证据 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-gq-r2-20261029) · [匿名 raw](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_gq_r2/REPORT_ZH.md)。执行结束时的本地未发布说明及打包回执保留为历史状态；当前发布记录单列在 `publication/`，科研数组、判决和费用不变。

# Research W-GQ-R2：真实全波执行入口

**结论：物理机制得到新证据；一般未知墙在线识别与新的 Bayesian 净加速未成立。**

- 两类普通非均匀墙，两纹理材料 RMS 差 **1.13137**。最细尺度的 12 个 held-out 目标、四发三频最大响应白化距离 **0.159971**；独立网格支持其有限任务传播近等效。
- 明确 **4 墙×12 目标 finite prior** 的目标位置 MAP 全选对，均值 RMSE **0.079224 λ₀**，墙纹理熵仍 **0.998388 bit**。这是受控 reference，不是一般 unknown-support 在线成功。
- 更宽合法墙材料中，目标横移 **0.05 λ₀** 后，A2 重新补偿距离在三离散层为 **0.889152 / 0.912855 / 0.913910**。增加通道可破坏旧 pair，但不保证排除重新补偿。
- 保留精确目标自散射，只压缩墙返回：rank 8 最坏白化误差 **2.05×10⁻⁴**。普通 Schur/T-matrix 同享此结构；没有独立加速主张。
- 未知宏观墙 finite 在线 A2 位置 RMSE **0.68343**；连续普通 joint FWI 最低 objective 选错分支。不得以真值较好分支冒充算法结果。
- 原 118 墙的条件目标共识未通过，冷 Bayesian **NOT_RUN**。本轮 **722 FE / 5212 complex RHS**，阶段费用与失败全部保留。

优先阅读 [完整中文报告](REPORT_ZH.md)、[机器裁决](GO_NO_GO.json)、[GPT Pro handoff](NEXT_PRO_HANDOFF.md)。复现见 [REPRODUCE.md](REPRODUCE.md)；最终费用以 [PHYSICS_COST.json](PHYSICS_COST.json) 和各 job receipt 为准。

重点图片：

1. [真实微结构与尺度](figures/F02_MICRO_SCALE_AND_ALL_PROBES.png)
2. [目标准确而墙纹理仍未知](figures/F13_CONTROLLED_MICRO_RISK_TEXTURE_AND_TV_TAILS.png)
3. [全波补偿与三层网格](figures/F14_FINAL_FINITE_COMPENSATION_AND_THIRD_GRID.png)
4. [全部 FWI 起点和错误的低 objective 分支](figures/F15_CONTINUOUS_ALL_STARTS_OBJECTIVE_AND_TRUE_ERROR.png)
5. [返回表示、修正和费用](figures/F08_RETURN_RANK_ERROR_SPECTRUM_AND_PAID_COST.png)

本轮只完成本地研究归档，未自动 push/发布。固定历史网页：[BGGA 原报告](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/73ba6eec3746548bf21259f627b635678709df00/research/research_w_bgga/RESEARCH_W_BGGA_REPORT_ZH.md)。新结果不要引用成该旧 commit 的内容。
