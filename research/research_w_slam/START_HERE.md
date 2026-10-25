# Research W-SLAM 实验入口

2026-10-25 · 实际计算完成 · 二维TM、开发研究

**结论：联合墙–目标补偿及回访纠错有真实Maxwell证据；普通可回访联合推断已经同样有效。本轮新重开规则没有计算优势，一般未知墙体在线恢复仍未解决。**

首先阅读[完整中文研究报告](RESEARCH_W_SLAM_REPORT_ZH.md)与[科学裁决](SCIENTIFIC_DECISION.md)。机器判决见[GATE_SUMMARY.json](GATE_SUMMARY.json)。

## 实际完成

- E0：FEM/PML、Mie、伴随/差分/互易及独立网格；QUALIFIED_DISCRETE_FULLWAVE，非连续PDE认证。
- E1：墙、位置、内禀材料三组Jacobian，多观测/运动公式和有限补偿；完整合同仍INCONCLUSIVE。
- A0：72条measurement-only连续主拟合，未知内禀/3×3扩展及标准伴随补优化；A0_ONLINE_NOT_REACHED，保留未收敛。
- E2/E3：2950个合法世界的全波参考表，12个parent×2噪声、11规则实际概率地图和回访。普通smoother位置误差.346λ0，hard freeze为.531λ0；新调度无净收益。
- E4：384项有限posterior比较，13个Gaussian未定义保留分母；边际正确不保证联合相关正确。
- E5：24案例、7策略、独立未来噪声；TaskEVSI实际优势区间跨零，不宣布新测量策略胜出。

有限概率prior已知目标内禀与仪器，2×2支撑不能表示孔洞。A0和有限实验不同prior，不能横向当作算法胜利。未开启locked、NN、三维或新A17 ROM。

## 建议查看

![回访精度与费用](figures/pilot_actual/E3_error_vs_CPU.png)

[物理协议](PHYSICS_PROTOCOL.md) · [理论条件](DYNAMIC_COMPENSATION_THEORY.md) · [E1物理](E1_MAXWELL_OBSERVABILITY.md) · [A0](A0_ONLINE_RECONSTRUCTION.md) · [概率belief](E2_BELIEF_MAP.md) · [回访消融](E3_REVISITATION_ABLATION.md) · [后验近似](E4_POSTERIOR_VALIDATION.md) · [主动测量](E5_ACTIVE_SENSING.md)

数据与费用：`results/pilot_summary_tables/`、[成本](COST_SUMMARY.json)、[物理账本](PHYSICS_LEDGER.jsonl)、[失败账本](FAILURE_LEDGER.jsonl)。完整NPZ/复场/联合权重按[复现说明](REPRODUCE.md)的相对目录放在本项目Release，源码及图表在Git。共享表实际准备约9288 CPU秒、11800 RHS/2950FE，不能拿cached毫秒数声称完整成像加速。

原Research W目录和旧结论未修改。下一步建议见[SCIENTIFIC_DECISION.md](SCIENTIFIC_DECISION.md)，本轮不自动扩大后续campaign。

公开数组：[Research W-SLAM Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-slam-20261025)。六个ZIP解包到本stage根目录，成员不多套一层目录，详见[RELEASE_ASSETS.json](RELEASE_ASSETS.json)。
