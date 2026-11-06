# MAX-STRUCT R9：A/A+高对比机制、噪声硬锁与非凸搜索

[中文研究报告](RESEARCH_REPORT_ZH.md) · [公式与条件](THEORY_DERIVATIONS.md) · [经典论文对象](CLASSIC_BENCHMARKS.md) · [三维扩展合同](EXTENSION_PROTOCOL.md) · [远端与CUDA](REMOTE_COMPUTE.md) · [成本](COST_LEDGER.md) · [复现](REPRODUCE.md) · [英文论文初稿](PAPER_DRAFT.md)

**目前的科学判断：**A-material有高对比开发正例；A+降低其材料读出费用，但尚无普遍优于成熟CIE的证据。第一阶段113次反演完成：公开噪声决定DPCS约束的修正，让独立六例中的普通CIE、Dual、Global及A+均从2/6提高到4/6；这主要是共享基线修复。完整current松弛路径恢复了一个历史失败壳体，简化后更便宜，但没有增加新六例的恢复数。

用户要求扩大到原NIE/CIE案例及真正三维。本目录第二阶段正执行64² Austria、全矢量小型对象和原3λ四立方体；只将有实际回执的结果作为已完成证据。未完成的对象不提前写进成功表。具体实施差异、原文停止公式纠正、欠分辨控制与Mie物理误差均保留。

R1–R8冻结证据不覆盖。代码复用原全波算子与FFT-TSOM，公开观测、预测与封存后评分分开。普通/高对比、同离散/跨网格、二维标量/三维矢量、公开开发/新确认不混合分母。低残差、运行正常与材料恢复是不同判断。
