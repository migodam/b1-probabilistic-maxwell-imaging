# Research F：材料几何与独立物理反演

**已执行并封账。裁决：几何审计 GO；独立反演取得两份合法且噪声相容的材料，但未达到事前冻结的驻点标准，F2 NOT_RUN；learned Explorer NO-GO。** 这不是不存在其他材料的证明。

[中文独立报告](RESEARCH_REPORT_ZH.md) · [网页版](report.html) · [材料与数值图集](gallery.html) · [全部指标](METRICS.json) · [成本](RESOURCE_SUMMARY.json) · [完整证据 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-f-intrinsic-20261017)

- 过去巨大 decoder 曲率中，weak 例几乎都是可消除的坐标加速度；ring 例还保留明显真实先验曲率与抵消。完整 Maxwell 二阶向量缺失，保持 MISSING。
- 新增24个二维 full states、960个 complex RHS，物理执行5.534秒；未重跑旧45-state诊断。
- 独立 full-physics TRF 反演不用N3/GT初始化；两条件是同一个已暴露 smooth parent 的两个背景，已知全局几何条件，不是一般未知形状反演。
- 0个经过F2验证的有限传输样本，0个新增后验模式；概率字段为null，不启动NN。

建议顺序：报告 → INTRINSIC_GEOMETRY_AUDIT → PHYSICS_ANCHOR_RESULTS → GO_NO_GO_DECISION → REPRODUCE。
