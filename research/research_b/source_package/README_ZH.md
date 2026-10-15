# Research B 研究包

日期：2026-10-14。完整中文理论报告与可执行的无 Maxwell 概率/代数核心。

## 阅读入口
- `RESEARCH_B_REPORT_ZH.md`：完整模型、证明/反例、算法、E0–E3、接口、成本/停止规则、新颖性与 Go/No-Go。
- `SOURCE_AND_NOVELTY_LEDGER.md`：项目证据、文献链接、阅读层级和访问限制。
- `CODEX_PROMPT_ZH.md`：可直接交给 Codex；先做零新增 Maxwell 调用 E0。
- `schemas/`：两个共享接口的 JSON Schema。
- `configs/EXPERIMENT_PROTOCOL.json`：E0 硬 cap 0；E1 4096 为尚未授权的独立预算建议。
- `code/`：有限 bank 风险/EIG、共享 gain、nuisance Schur、排名区间。
- `checks/THEORY_CHECKS.json`：11 项已经运行且通过的代数/概率检查，Maxwell calls=0。

## 运行
需要 Python、NumPy、SciPy；Schema 验证另需 jsonschema。
```bash
python code/verify_theory.py
```
代码不包含 Maxwell solver，也不自动下载数据或启动训练。`acquisition_core.py` 是可复用的有限-bank参考实现，不是已接通 Joint P2 的完整 E0 runner；实际缓存 adapter、E0 排名及物理 E1/E2/E3 尚未执行。

## 主要边界
有限 bank 与连续 posterior 的误差不同于外层 MC error；未知材料质量保持 null。单频 C1 原始 proposal 不自动适用于新的多频/history。full likelihood correction 必须有合法 proposal 密度或另一个明确采样机制。任何无法稳定区分的世界可以继续保留。

没有新最优频率、主动成像 MAE 改善、ROM 加速或强 TAP 结果已经成立的主张。
