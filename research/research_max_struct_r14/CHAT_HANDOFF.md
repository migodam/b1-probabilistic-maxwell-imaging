# MAX-STRUCT R14：一页交接

**不是A+被否定，也不是“一致metric必胜”。A+的部分二维增量得到独立验证；一致非局部11也产生了原栈全部失手的新正例。根处闭环谱不足以预测真实有限程序。**

- 实际576次反演，另有同完整state干预与根/非根机制诊断；所有运行和失败计费。
- 曝光48：00/01(A+)/10/11为29/35/30/31，01相对11多0017/0024/0026/0045且无丢失。旧R13 A+36冻结不改；不同停止/计时只作历史参考。
- 新独立24：A+12，CIE9，Dual9，G19 11，DR14，11为11；A+对CIE/Dual各多3且无丢失，成对共同CPU与共同LU/RHS池仍12比9。
- **新r14_new_0021只有11恢复**：NRMSE.08696、IoU.9252、residual/δ.9283；A+、CIE、Dual、G19、DR、direct和晚切11全差。反过来0012/0015是01好、11差。
- 01→11切换已实际独立验证，仍同12集合，没有新增恢复；不宣传成突破。公开01→G19→DR回放14/24，平均160.28CPU秒，DR单方法平均116.57秒，费用不能省略。
- self控制与共享χ-LBFGS显示0017/0045需要offdiag Green，0017不是Arc专属收益；0024也可由global-Dual解释。
- ideal01ρ1.002712；真实finite PR/guard使根附近有限响应不同。非根metric导数力在56预定first-step上常与构成力同阶，const梯度上坡38/56而full梯度0/56。数组/PR/weights/stage齐全供后续理论。
- 五旧严重失败仍未解决，有已知更低Φ合法GT；19/20端点仍有微小物理下降，不能叫严格局部极小证书。
- 三维24 runs是真dyadic Maxwell、3几何4条件；未建立普遍新增恢复。新Fresnel/Austria/NN arms没有运行。

不要只读联合标签：新0002的A+与DR连续NRMSE都约.053，IoU .667/.710的差不能当成大basin救援。主要结构/部分恢复与旧联合标准分开解释，不能GT改门槛。

首先读[中文完整报告](RESEARCH_REPORT_ZH.md)、[理论](THEORY_DERIVATIONS.md)、[轨迹](TRAJECTORY_ANALYSIS.md)、[栈关系](TECH_STACK_MAP.md)、[Pro材料指南](PRO_HANDOFF.md)。新24全图/配对/费用在[独立分析](analysis/independent24_v1/report.html)；完整原档、source/public与离线truth见DOWNLOADS。英文稿为[PAPER_DRAFT.md](PAPER_DRAFT.md)。

本轮定位是有真实二维算法增量的完整机制研究稿。无需所有case全胜；强TAP/强SIAM仍缺更有预测力的有限程序理论与更强外部/三维证据。下一步优先解释“为什么0012/0015与0021需要相反的早期metric”，而非继续堆阶数或按GT选算法。
