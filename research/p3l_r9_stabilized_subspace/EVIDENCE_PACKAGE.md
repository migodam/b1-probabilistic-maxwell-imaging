# R9完整证据与复现

本轮复用现有B1仓库及Pages。完整数组、逐状态账本、预测封存、原始合成GT、所有失败和历史依赖以无损Release分片保存；Git保存源码、配置、主要CSV/JSON、推导和图片。

[Release入口](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr9-stabilized-volume-20261030) · [复现说明](REPRODUCE.md) · [分片身份清单](EVIDENCE_MANIFEST.json)

## 两层布局

1. `R9_fresh_000_007`、`008_015`、`016_023`、`024_031`四个tar.zst，均为repository-relative：32新parent的完整端点、场/J、visited候选、拒绝和原始job日志。
2. `R9_development_and_probability_arrays.tar.zst`：旧对象开发对照、独立概率积分、精确ROI后端的完整数组。
3. `R9_complete_metadata_scoring_and_sources.tar.zst`：repo-relative源码、冻结配置、封存、原始评分、全部账本与合成GT。`private/`表示实验揭盲前的GT隔离区，完整预测封存与获准评分后作为OFFLINE_ONLY科研证据公开；在线代码仍不可读取。
4. `R9_historical_dependencies.tar.zst`：单独Gaussian-relative依赖树，保留R4/R6/R7及A30的原相对结构，只读历史输入不重新计费。
5. `Research_P3LR9S_Reading_Kit.zip`：轻量阅读包，含报告、公式、代码、主要统计与图，不代表完整原始数组。

同一层的分片解压到同一个新目录。两层不能混放；按照REPRODUCE中的relocation helper组装独立运行树。原始归档不改写；路径迁移生成明确标记的派生副本，保存原身份及迁移记录。

## 科学身份与成本

- 新32parent生成前科学冻结提交：`75184170`。
- 最终报告与代码提交：由本Release关联tag及`PUBLICATION_IDENTITY.json`记录。
- 源身份、代码与配置修订、预揭盲封存时间均保留。
- 新物理分别计22776个576阶RHS与1747416个256阶RHS，不能混算等效费用。
- 压缩、传输、离线重评分与网页发布不记为新增Maxwell实验。
- 每片给出文件数、字节数及BLAKE2b身份。完整下载与公开可访问性由PUBLIC_READBACK.json记录；没有把网址存在当作科学资格。

## 结果资格

32parent只支持本次小规模机制确认。候选数和GT覆盖不等于posterior mass；启发式权重不是已认证posterior。精确ROI消元是标准算子机制的实际实现，58.65%只属于匹配物理核；固定对动态的实际整链差约2.75%。高概率GT覆盖、20%固定空间加速及强TAP独立贡献均未成立。

无需重新生成GT或Maxwell状态即可重评分。历史资料、第三方论文全文、私人对话和虚拟环境不夹带进入科研证据包。全部科学负结果保留。

## 已完成的独立迁移复现

在新的隔离目录中，按公开归档布局恢复179份依赖和97份OFFLINE_ONLY合成文件，先验证32parent/96方法组的封存，再运行一次冻结评分。七份主要评分与成本JSON的47013个数值完全一致；原始文件不变，新增物理为0。该检查仅证明保存证据可搬迁并重评分，不是重新生成32对象或新的科学重复。

[实际复现审计](audit/reproducibility/RELOCATION_REPLAY_AUDIT.md)
