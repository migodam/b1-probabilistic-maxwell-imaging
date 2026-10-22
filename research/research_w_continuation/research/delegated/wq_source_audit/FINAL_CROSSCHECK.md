# 发布前机械交叉核对：差异清单

2026-10-23。只读报告、标量结果和账本；未运行物理、载入大NPZ、做哈希或修改科学判决。除主线程明确授权的37处来源链接修复外，没有改写报告。以下仅列实际发现的差异。

1. **最小特征约束尚未实现核验。** 冻结配置 `materials.minimum_feature=0.1`，但逆表示为48×32独立像素，dx=4/48≈0.08333、dy=3/32=0.09375；现有平滑提案及box/top-area投影没有强制该最小特征。报告“同一box/面积/先验”不能被读成全部冻结材料合法性已核实。主线程已指定标为 **minimum-feature legality UNVERIFIED / protocol implementation gap**：候选不能用于正面见证声明；本轮最小dcal仍为19.912258993，费用保留，这一缺口也不能转成不存在结论。证据：[冻结配置](../../../docs/P0_FROZEN_CONFIG.json)第127行；[PixelGrid/投影](../../../src/wq/inverse_support.py)第20–55、108–141行；[平滑提案](../../../src/wq/search.py)第32–34行；[报告](../../../RESEARCH_W_CONTINUATION_REPORT_ZH.md)第32行。此项不设计新准则、不更改gate。

2. **37个本机来源链接已修复。** `docs/OLD_EVIDENCE_AUDIT.md` 第156–179行原指向 `/Volumes/.../research-w-p0/`，公开网页无法解析这些本机路径。按主线程授权，已逐个改为旧来源固定commit `db7146cae2d74ba231577d42f225dbfc9a389268` 下 `research/research_w_p0_20261022/...` 的GitHub链接，并将`:行号`改为`#L行号`；保留全部正文与范围。未改源文件或旧gate。见[修复后的来源表](../../../docs/OLD_EVIDENCE_AUDIT.md)。

3. **CPU字段名称比实际计时范围更强。** [COST_SUMMARY.json](../../../COST_SUMMARY.json)第10行将Q的335.773484秒命名为`measured_job_process_cpu_seconds`；它包含[Q账本](../../../results/PHYSICS_LEDGER.jsonl)第1185行worker的3.200507秒，而该行明确仅测数值区域，排除了imports/startup。算术与charged总数一致，不应删费用；但字段/注释应保留“包含局部区域计时，非全部完整进程实测”的限制。完整作业CPU部分为332.572977秒，另加该区域计时才得到335.773484秒；60秒保守准备分配仍另外列入charged。

4. **一个来源链接的类型标签不符。** [P0_B_PHYSICS_VALIDATION.md](../../../P0_B_PHYSICS_VALIDATION.md)第27行的“原复核数组”实际指向`results/p0b/finer-endpoints/RESULTS.json`，该文件是标量结果表，不是NPZ数组。应只改标签为“原复核结果”，或将数组入口指向Release清单；不表示原始NPZ缺失，也不要求重新求解。

`ARRAY_MANIFEST.json`在检查时尚未生成；主线程明确由`tools/prepare_publication.py`稍后生成并随Release打包，故不将这一已知出版尾项列为物理数据缺失。
