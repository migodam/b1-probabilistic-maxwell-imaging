# 完整证据下载

本轮沿用现有B1仓库与网页，不新建repository。Git保存源码、公开输入、封存后真值、图、最终汇总、费用与审计；大量逐运行材料／预测NPZ和压缩轨迹放在完整Release ZIP，本地研究目录仍保留全部原件。完整下载包用于一次性复算；它不额外提供未执行实验。

- [完整复现ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-fbg-r3-e-20261101/fbg-r3-e-full-reproducibility.zip)
- [轻量阅读包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-fbg-r3-e-20261101/fbg-r3-e-reading-pack.zip)
- [Release说明](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-fbg-r3-e-20261101)

完整包在research/research_fbg_r3_experiments/下保留相对路径，并携带本轮直接调用的R4物理、Q6/Q20和R2/CUDA依赖源码与所需历史公开观测。Windows CUDA需要已有兼容PyTorch/CUDA环境；没有捆绑驱动或安装服务。CPU小网格入口见REPRODUCE.md。已有输出不要覆写，用新replay tag；GT仍只供封存后的评分。

阅读包只含报告、数学、源代码、代表图、结构化最终摘要与manifest，体积小；它不冒充完整原始数据集。工具缓存不入包；已测失败及未知成本不被清理掉。
