# 公开证据与本地原始场的边界

[同一现有仓库的W证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-p0r-to-p1-20261024)提供两份独立ZIP分卷，以及单独的库存与恢复说明。不新建仓库。先下载两个ZIP，解压到同一个空目录；若REASSEMBLY.json列出单文件byte chunks，按offset顺序拼接并核对记录字节数。

主要文件：

- `w-evidence-part001.zip`、`w-evidence-part002.zip`：固定源码、vendor、配置、完整现存复响应、材料/latent、Jacobian诊断、全部SMC轨迹、评分、失败与cost元数据；每卷不超过1.8GB。
- `INVENTORY.json/INVENTORY.csv`：逐文件来源、字节数、included状态和原因。
- `LOCAL_FULL_FIELDS_INVENTORY.json`：未上传的原始full fields仍在原本地路径；不把库存当成云端备份。
- `PACKAGE_MANIFEST.json`、`REASSEMBLY.json`、`RESTORE.md`：格式、分片和恢复合同。

明确公开的完整复场子集为第一joint见证cal/joint×两世界×五频，共20个完整频率archive及端点伴随文件。它包含该子集的材料/网格/复场；其余物理端点的完整场保留本地，公开响应、材料参数、decoder与费用足以检查报告中的指标，不宣称全部本地全场均已上传。

P1A/P1B的OFFLINE_TRUTH与生成计划只在CAMPAIGN_CLOSED后作为不解码的完整文件流拷贝进证据包；它们是已结束实验的离线重现证据，不能成为在线输入。当前小型MaterialHypothesisPackage完全不含GT，并可单独下载重新加权。

Python环境、系统runtime、.deps和旧压缩包不混入新分卷；版本要求与冻结源保留。没有SHA-256检查、内容哈希认证或新Maxwell求解。打包仅进行存在性、字节数与流拷贝前后size/mtime/inode核对；公开发布和匿名读回的记录另列，不充当科学Gate。
