# R12 原始证据归档（已公开并验收）

所有源文件保持原位。代码、JSON/CSV/MD、图片与小预测留Git；巨大数组及运行bundle保存在分包资产。单包严格小于1.8GB。未计算内容hash；核对ZIP中央目录的路径、条数、字节和实际包体大小。

|资产相对路径|字节|文件条目|内容索引|
|---|---:|---:|---|
|`release_assets/max-struct-r12-local-01.zip.zst`|1197758756|11|[max-struct-r12-local-01.index.json](audit/release/max-struct-r12-local-01.index.json)|
|`release_assets/max-struct-r12-local-02.zip.zst`|1072155131|5|[max-struct-r12-local-02.index.json](audit/release/max-struct-r12-local-02.index.json)|
|`release_assets/max-struct-r12-local-03.zip.zst`|1101604121|10|[max-struct-r12-local-03.index.json](audit/release/max-struct-r12-local-03.index.json)|
|`release_assets/max-struct-r12-local-04.zip.zst`|1207851369|4|[max-struct-r12-local-04.index.json](audit/release/max-struct-r12-local-04.index.json)|
|`release_assets/max-struct-r12-local-05.zip.zst`|1150850375|4|[max-struct-r12-local-05.index.json](audit/release/max-struct-r12-local-05.index.json)|
|`release_assets/max-struct-r12-local-06.zip.zst`|586229526|2|[max-struct-r12-local-06.index.json](audit/release/max-struct-r12-local-06.index.json)|
|`release_assets/max-struct-r12-local-supplement-01.zip`|199543411|124|[max-struct-r12-local-supplement-01.index.json](audit/release/max-struct-r12-local-supplement-01.index.json)|
|`release_assets/max-struct-r12-remote-01.zip.zst`|1207041413|312|[max-struct-r12-remote-01.index.json](audit/release/max-struct-r12-remote-01.index.json)|
|`release_assets/max-struct-r12-remote-02.zip.zst`|931818541|2460|[max-struct-r12-remote-02.index.json](audit/release/max-struct-r12-remote-02.index.json)|
|`release_assets/max-struct-r12-remote-03.zip.zst`|989774680|2396|[max-struct-r12-remote-03.index.json](audit/release/max-struct-r12-remote-03.index.json)|
|`release_assets/max-struct-r12-remote-04.zip.zst`|980505288|2442|[max-struct-r12-remote-04.index.json](audit/release/max-struct-r12-remote-04.index.json)|
|`release_assets/max-struct-r12-remote-05.zip.zst`|791050745|3616|[max-struct-r12-remote-05.index.json](audit/release/max-struct-r12-remote-05.index.json)|

逻辑包总计11416183356字节，11386个ZIP条目、11371个去重科学路径。完整映射在[ARTIFACTS.json](ARTIFACTS.json)，最终覆盖复核见[FINAL_ARCHIVAL_COVERAGE_REVIEW](audit/release/FINAL_ARCHIVAL_COVERAGE_REVIEW.json)。

远端basin原始states恢复到RESULTS/basin_remote；confirmation完整RESULTS/queues/engineering在delegated/confirmation/remote_complete；实际部署的R12源与共享legacy源在remote_runtime中。共享源按canonical路径计一次；15个已存历史basin复制文件作为逐字节相同的别名保留，见[别名索引](audit/release/CANONICAL_PATH_ALIASES.json)，不重新打包、不计第二次科学运行。历史输入bundle亦不是新实验。没有私人聊天、密钥文件、原始Fresnel测量或仓库全history。

发布采用128MiB运输parts，另见[PUBLICATION_PARTS_MANIFEST.json](PUBLICATION_PARTS_MANIFEST.json)。这些parts只切分逻辑包字节、不重算/重压缩科学证据；先按运输索引连接parts恢复原wrapper，再进行以下解码。不要把逻辑包和parts当成两份独立科学证据或重复收费。

恢复：.zip.zst先用`zstd -d <asset>.zip.zst -o <asset>.zip`无损解码；然后在空目录运行`python -m zipfile -e <asset>.zip <empty-destination>`。普通.zip直接第二步。归档路径均相对仓库，不要解到正在写入的科学工作树。Git最新科学文稿、图和代码优先于大包中的历史文稿快照；大包专门保留原始数值/执行证据，不为小文稿更新重打。

费用和运输见audit/release/{LOCAL_PACKAGING_RECEIPT,REMOTE_PACKAGING_RECEIPT,REMOTE_ZSTD_RECEIPT}.json及TRANSPORT_COSTS.jsonl。远端压缩只读取已结束任务，写入专用export目录，不触碰其他会话。准备脚本没有发布、push或删除源；实际draft/release上传由root另行记账。

## 最终发布验收

12个逻辑归档以89个运输分片公开；逐个匿名Range回读确认可下载与完整远端尺寸。见[公开回执](audit/PUBLICATION_READBACK.json)、[Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r12-20261111)及[运输清单](PUBLICATION_PARTS_MANIFEST.json)。使用restore_release_assets.py按顺序恢复；旧归档内文稿以当前Git文稿为准。
