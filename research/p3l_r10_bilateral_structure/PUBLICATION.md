# 同一B1仓库中的公开交付

本页给出发布目标与无损归档布局；最终远端状态以 `PUBLIC_READBACK.json` 为准。

- [阅读网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_r10_bilateral_structure/)
- [源码与报告目录](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/p3l_r10_bilateral_structure)
- [原始中文Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/p3l_r10_bilateral_structure/RESEARCH_REPORT_ZH.md)
- [同一B1 Release与下载](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr10-bilateral-20261030)
- [持续公式集](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/maxwell_formula_atlas/)

## 归档内容

| 资产 | 用途 |
|---|---|
| Research_P3LR10BI_Reading_Kit.zip | 报告、公式、源码、冻结配置、紧凑原始指标、独立QA、图库与论文草稿 |
| R10_DEV8_full_evidence.tar.zst | DEV8全部付费状态、Jacobian、端点、候选、封存记录与逐调用账本 |
| R10_Fresh10_full_evidence.tar.zst | 新10parent的观测、私有生成GT（评分后公开）、全部状态、失败、封存与评分 |
| R10_physics_probability_full_evidence.tar.zst | 传递/积分及失败、两次压缩实验、836节点重放与汇总逐调用账本 |
| R10_dependency_inputs.tar.zst | 423个只读历史依赖，R4 prior与原配置、R9精确ROI缓存、DEV8上下文及历史选中端点 |
| RELEASE_MANIFEST.json | 文件大小、BLAKE2b身份及恢复说明 |

精简阅读包不是全部训练或历史证据。依赖包按B1 sibling布局恢复，不覆盖历史文件；当前新证据包按仓库相对路径恢复。历史原始结果不改写。具体恢复与只读重评分见 [REPRODUCE](REPRODUCE.md)。
