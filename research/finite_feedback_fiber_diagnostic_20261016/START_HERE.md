# Maxwell有限反馈与GT差量诊断 · 2026-10-16

**局部弱能量不等于可连通的材料解族。** 同一公开development对象上，N3–GT差量的98.09%–99.06%属于预设阈值下的弱材料方向；但N3和physics-only R基线均未成为原P0、原噪声下的合格anchor。没有建立通向GT的合法等数据路径。

真实full-Maxwell的有限反馈已核验；高阶显著提高预测精度。两个合法TRAIN chart只得到小幅连续材料变化、已采样端点noise-compatible，0个显著新结构。主要有限步障碍是材料decoder曲率，不是此处的Maxwell共振。

- [完整中文报告](DIAGNOSTIC_REPORT_ZH.md) · [网页报告](report.html) · [数值图与材料图](gallery.html)
- [GT强弱分解](STRONG_WEAK_GT_DECOMPOSITION.md) · [有限反馈](FINITE_FEEDBACK_ORDER_COMPARISON.md)
- [局部路径与曲率](LOCAL_PATH_AND_CURVATURE.md) · [anchor/先验审计](PRIOR_SUPPORT_AND_ANCHOR_AUDIT.md)
- [下一轮Pro应回答的问题](THEORY_QUESTIONS_FOR_PRO.md) · [机器可读数据](METRICS.json)
- [独立资源账本](RESOURCE_LEDGER.jsonl) · [复现及证据数组](REPRODUCE.md)

45/48新增二维states，1880/2048 complex RHS，受控wall 11.146秒；GPU训练/Mie/3D均0。逐endpoint立即保存，没有LXZ原始数组丢失的重演。已完成预设诊断，剩余3 states/168 RHS不追加寻找正结果。旧研究与预算未改。

TRAIN已知材料不是未知成像成功；GT方向只属ORACLE_GEOMETRY_DIAGNOSTIC。无posterior/mode概率声明。
