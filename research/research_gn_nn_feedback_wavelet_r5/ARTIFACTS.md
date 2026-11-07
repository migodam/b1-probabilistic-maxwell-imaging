# 数据、权重与完整数组

已发布到现有B1仓库的[GN–NN R5 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/gn-nn-feedback-wavelet-r5-20261108)。11套权重、两阶段公共预测与候选、新噪声感知GN以及精确的旧TRAIN/VAL依赖均已上传。普通Git保存报告、代码、紧凑逐案例表及全部科学图；没有新建仓库。

所有对象仍是DEVELOPMENT；下载包包含TRAIN监督和离线VAL评分材料，公共推理不读取GT。新鲜TEST为NOT_RUN。

|附件|实际字节|内容|
|---|---:|---|
|[gn-nn-r5-models.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-feedback-wavelet-r5-20261108/gn-nn-r5-models.zip)|2,339,480|11套checkpoint、config、history、原物理费用；小型便捷副本|
|[gn-nn-r5-reproducibility.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-feedback-wavelet-r5-20261108/gn-nn-r5-reproducibility.zip)|769,992,162|冻结R5源码、数学、权重、全部原始预测／阶段状态、候选、在线排名、离线评分、费用、新GN输入和图像|
|[gn-nn-r5-frozen-inputs.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-feedback-wavelet-r5-20261108/gn-nn-r5-frozen-inputs.zip)|231,002,029|精确复用的TRAIN128／VAL78观测与材料、严格GN、固定G/E/S／谱／参考J，以及完整相对导入所需旧Python包|
|[gn-nn-r5-release-index.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/gn-nn-feedback-wavelet-r5-20261108/gn-nn-r5-release-index.json)|846,414|逐文件路径与字节、来源commit、包装费用与恢复说明；不包含哈希检查|

两个完整ZIP合计约1.001GB，小型权重ZIP与其中权重重复。若只需查看结果，网页已包含全部主图和逐案例表，不必下载完整数组。

## 在另一个目录恢复

```bash
git clone --filter=blob:none --sparse --branch gn-nn-feedback-wavelet-r5-20261108 https://github.com/migodam/b1-probabilistic-maxwell-imaging.git b1-r5
git -C b1-r5 sparse-checkout set research/research_gn_nn_feedback_wavelet_r5
gh release download gn-nn-feedback-wavelet-r5-20261108 --repo migodam/b1-probabilistic-maxwell-imaging --pattern 'gn-nn-r5-reproducibility.zip' --pattern 'gn-nn-r5-frozen-inputs.zip' --dir r5-downloads
python -m zipfile -e r5-downloads/gn-nn-r5-reproducibility.zip b1-r5
python -m zipfile -e r5-downloads/gn-nn-r5-frozen-inputs.zip b1-r5
```

归档保留`research/...`层级；先恢复上述两个完整ZIP即可，不需要另下载整个旧R4的8.4GB归档。小型models.zip只是便捷下载。依赖包包含原冻结实现，不是重新编写的物理后端。

数值源码基准为`5930a6bde26541f50402d21146cd87ca84634ee7`。归档内ARTIFACTS的“发布整理中”是该包装时点的历史状态；本网页及`audit/PUBLICATION_RELEASE.json`记录之后的实际公开发布。后续网页／下载说明不改变冻结权重、候选或科学结论。

报告与复现指令见[REPRODUCE.md](REPRODUCE.md)。移机复算生成新回执，不能改写原机封存路径冒充原实验。新环境计时属于新复算，不冒充原机性能。
