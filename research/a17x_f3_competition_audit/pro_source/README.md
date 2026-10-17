# READ THIS FIRST: F3 evidence update (2026-10-17)

The original 1,340-line Chinese theory report, proofs and algebra tests are preserved. Research F3 was completed AFTER the original A17-X report. New pinned F3 evidence closes original P0 local Explorer tests; it does NOT test A17 current-space exchange. FIRST read `protocol/F3_FIRST_GATE.md`, `integration/F3_A17X_DECISION_20261017.md` and pinned F3 sources listed in `evidence/f3_pinned/F3_EVIDENCE_INDEX.md`. Do not launch archived P0 full-state budget automatically. Prefer a zero-new-physics readiness audit, followed only by a separately frozen and authorized P1 study. Use `protocol/PASTE_TO_CODEX.txt` for the shortest instruction.

---

# A17-X Research Package

2026-10-17。中文理论研究、证据审计与下一轮 Codex 交接。

**首先阅读** `report/RESEARCH_A17X_FULL_REPORT_ZH.md`。完整报告包含原 A17 定理边界、有限双侧 exchange 及尾项界、多源反馈可实现材料图册、显式 source-chart 密度、原创性审计、成本、P1 分支与 Q1–Q10。

**执行交接** `protocol/CODEX_NEXT_ROUND_ZH.md`。这是下一轮实验提示词，不是已执行的物理实验。建议新预算 128 full states / 4096 complex RHS，禁止直接启动训练。

**代数核验**：

```bash
python checks/check_finite_source_transport.py
```

需要 NumPy 和 SciPy。结果写入 `checks/algebra_results.json`。固定随机种子为 20261017。运行只处理小型合成复矩阵，不访问网络、不运行真实 Maxwell 后端、不训练 GPU、不读取或假设 P1。

**最重要的判定**：有限材料传输机制有条件成立；A17 的真实增量仍需与精确局部 T-matrix / Woodbury 和普通 goal-oriented ROM 对照。不能把 DCTMC、CSI、Schur、隐函数、Davis–Kahan 或 MH 重新命名为本项目新理论。

单源标量正例不代表多源三维成功；人工 H/O 正例不代表真实天线场景。所有物理实验和 posterior 性能仍为未运行/未证实。

来源与访问范围见报告第 0、21 节及 `sources/SOURCE_AUDIT.md`。`CLAIM_LEDGER.json` 区分已推导代数、合成核验、待实验结论与明确非主张。