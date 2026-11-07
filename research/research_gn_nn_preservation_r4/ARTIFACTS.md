# 数据、权重与完整复现归档

已发布到现有B1仓库的 [GN–NN R4 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-preservation-r4-20261107)。共17个附件，8.393 GB；没有新建仓库。

[全部20套训练权重，13.59 MB](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-models.part001.zip) · [完整归档索引](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-release-index.json)

普通Git保留源代码、中文报告、紧凑结果表、科学PNG及1,296页无损WebP图册。大数组、完整评分JSON、分组能量表、全部候选和权重通过Release发布；阅读网页无需下载8.4 GB归档。

## 附件清单

|附件|实际字节|内容说明|
|---|---:|---|
|[DEPENDENCY_CORRECTION.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/DEPENDENCY_CORRECTION.json)|650|归档索引／包装回执／依赖更正|
|[gn-nn-r4-analysis-tables.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-analysis-tables.part001.zip)|592,562,389|完整评分JSON、E0与新TEST能量明细|
|[gn-nn-r4-anchors.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-anchors.part001.zip)|1,226,563,383|严格GN与全部传统控制端点及失败准备证据|
|[gn-nn-r4-canonical-receipts.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-canonical-receipts.part001.zip)|27,531,634|实际执行、失败、成本与冻结回执|
|[gn-nn-r4-data-evaluator.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-data-evaluator.part001.zip)|76,154,459|GT和离线评分清单（候选封存后揭晓）|
|[gn-nn-r4-data-public.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-data-public.part001.zip)|19,603,162|公开观测、固定物理算子及参考Jacobian|
|[gn-nn-r4-evaluation.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-evaluation.part001.zip)|1,499,939,818|全部真实物理候选、在线排名与原始seal|
|[gn-nn-r4-evaluation.part002.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-evaluation.part002.zip)|1,499,721,016|全部真实物理候选、在线排名与原始seal|
|[gn-nn-r4-evaluation.part003.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-evaluation.part003.zip)|346,636,417|全部真实物理候选、在线排名与原始seal|
|[gn-nn-r4-fixed-projectors.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-fixed-projectors.part001.zip)|13,743,089|固定完整Green子空间及接收投影|
|[gn-nn-r4-frozen-python-dependencies.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-frozen-python-dependencies.part001.zip)|3,039,266|实际导入的冻结旧研究Python依赖，含相对导入及初始化器|
|[gn-nn-r4-models.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-models.part001.zip)|13,587,789|全部20套权重、训练配置及曲线|
|[gn-nn-r4-predictions.part001.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-predictions.part001.zip)|1,499,621,262|网络原始提案及阶段电流|
|[gn-nn-r4-predictions.part002.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-predictions.part002.zip)|1,499,580,191|网络原始提案及阶段电流|
|[gn-nn-r4-predictions.part003.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-predictions.part003.zip)|49,306,267|网络原始提案及阶段电流|
|[gn-nn-r4-release-index.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/gn-nn-r4-release-index.json)|25,227,929|归档索引／包装回执／依赖更正|
|[PACKAGING_RECEIPT.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-preservation-r4-20261107/PACKAGING_RECEIPT.json)|937|归档索引／包装回执／依赖更正|

## 恢复与复算

需已有Python依赖，版本见 `audit/FINAL_ENVIRONMENT.json`。以下为复现指令，不宣称本轮另外重复训练或重跑物理。单纯查看结果和图像无需执行。

```bash
git clone --filter=blob:none --sparse --branch gn-nn-preservation-r4-20261107 https://github.com/migodam/b1-probabilistic-maxwell-imaging.git b1-r4
git -C b1-r4 sparse-checkout set research/research_gn_nn_preservation_r4
gh release download gn-nn-preservation-r4-20261107 --repo migodam/b1-probabilistic-maxwell-imaging --dir r4-downloads
python b1-r4/research/research_gn_nn_preservation_r4/experiments/package_release.py restore --index r4-downloads/gn-nn-r4-release-index.json --target r4-restored
rsync -a r4-restored/research/ b1-r4/research/
```

恢复工具保留记录的mtime。原始seal含原机路径，作为历史证据保留；在新目录生成另外的本地seal，不能修改原seal冒充原机运行。具体评分命令见 [恢复说明](delegated/data/RELEASE_PACKAGING.md#restoring-bytes-and-preserving-provenance)。重评分和重画图不需要重新生成GN或重新训练。

新R4完整复算不需要旧R3训练数据；历史E0原始审计另依赖 [已发布R3数据／评估归档](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-completion-r3-20261107)。该部分为已披露诊断，不冒充新测试。

包装时记录的科学源码基准为 `e27a413e2d8c41d9b5c88b4e68c800a8c06dfd4b`；Release标签目标为 `feee27654527eb1ba2dceb583ff6edd8682c4cf8`，差异为报告补充。后续网页／下载回执提交不改变冻结权重或数值结果。

PACKAGING_RECEIPT的“未上传”字段记录包装阶段本身，后续公开发布由本页及 `audit/PUBLICATION_RELEASE.json` 记录。


发布后的依赖封装更正：原AST导入扫描未递归相对导入及包初始化器，补入space.py、solver.py和cuda_physics.py三个原有冻结Python文件，并同步更换小型依赖ZIP与索引。全部观测、权重、候选、评分和科学结论不变；更正记录为Release附件 DEPENDENCY_CORRECTION.json。没有为此重新运行实验。
