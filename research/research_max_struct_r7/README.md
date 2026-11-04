# MAX-STRUCT-R7

主要入口：[中文完整报告](RESEARCH_REPORT_ZH.md)、[给Chat交接](CHAT_HANDOFF.md)、[论文初稿](PAPER_DRAFT.md)。

本轮已有真实算法与实验，不是待执行计划。A+为Möbius二次状态驻点候选+原prior-MM接受；固定电流loss非增不保证原物理成像收敛。主30新对象×6方法、完整披露开发和回归、精确有理有限动作、真实Fresnel先导与真正矢量DDA反演都保存正负结果。根据主30信号追加的家族确认单列，不能合并成一个事前注册50对象分母。

- `THEORY_DERIVATIONS.md` / `METHODS_AND_ABLATIONS.md`：实际公式、证明和实现区别。
- `PRIOR_ART_AND_NOVELTY.md`：已有工作、证据分级与原创性限制。
- `BENCHMARKS.md` / `BENCHMARK_3D.md` / `FRESNEL_AUDIT.md`：物理、噪声、材料和评分合同。
- `RESULTS/scoring/`：180主预测的独立评分、配对区间、CPU/LU/RHS检查点匹配、GT可行见证；`extra/`、`rational/`分母独立。
- `RESULTS/targeted/`及目标家族合同：后续冻结确认，不覆写原30。
- `src/` / `configs/` / `experiments/`：复用R1–R6算子的可复现入口；不改历史代码。
- `audit/`：来源冻结、配置/预测封存、内部审查、费用和发布回执。
- `FAILURES.md` / `RESEARCH_LOOP.md` / `COST_LEDGER.md`：所有否证、修正、转向及费用缺口。

大型逐步轨迹、已支付状态和预算数组以GitHub Release分卷提供，所有分卷解压到仓库根。Git保存方法、配置、最终预测/指标/图及审查；第三方原始测量、论文全文和附件ZIP仅本地保留，官方获取脚本/来源链接公开。详见`REPRODUCE.md`与`audit/PUBLICATION_SCOPE.json`。

所有GT在预测封存后离线评分；评分后公开，不再称未来盲测。无NN训练，无GT在线选根，无SHA校验，未把运行完成等同科学成功。
