# R8归档已完成并公开上传

实际打包完成时间：2026-11-06T00:11:17.019786+00:00。52个ZIP、29,823个唯一原成员、1,126份独立compact已生成；压缩后共5,308,028,340 bytes，最大单包173,504,125 bytes。无SHA检查，无源文件删除或改写。全部52包随后已上传至[公开Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r8-20261106)，匿名访问见[发布回执](../../audit/PUBLICATION_RECEIPT.json)。ARCHIVE_INDEX.json的PACKAGED_NOT_UPLOADED为打包完成时的原快照，不替代后续发布回执。

实际依据：[打包执行回执](PACKAGING_EXECUTION_RECEIPT.json)、[52包成员与尺寸核对](ARCHIVE_VERIFICATION.json)、[全部归档绝对路径](ARCHIVE_PATHS.txt)。最终费用已把打包过程单独计一次，见根[COST_LEDGER.md](../../COST_LEDGER.md)。

以下为打包前stat计划的历史快照，文件源字节不是压缩下载量；不代表尚未执行。实际ZIP以执行回执和核对文件为准。

# R8发布档案计划

Git候选 18187 文件 / 87628690 bytes；release 29823 文件 / 16334522276 bytes。

默认计划只stat原文件与读取公开manifest；未读数组/GT。执行结果以PACKAGING_EXECUTION_RECEIPT为准；所有数值都是当前快照，最终完成后重跑刷新。

|包|成员数|非压缩字节|
|---|---:|---:|
|max-struct-r8-audit-part001.zip|2|22647|
|max-struct-r8-delegated_cost_package-part001.zip|3|27671648|
|max-struct-r8-delegated_warm_analysis-part001.zip|4|8913841|
|max-struct-r8-delegated_warm_anchor_analysis-part001.zip|8|19004215|
|max-struct-r8-delegated_warm_direct_budget_analysis-part001.zip|2|1217497|
|max-struct-r8-experiments-part001.zip|2|17273|
|max-struct-r8-fresnel-part001.zip|98|11946615|
|max-struct-r8-history_warm-part001.zip|615|368414076|
|max-struct-r8-paper_baselines-part001.zip|4|52265|
|max-struct-r8-paper_benchmarks-part001.zip|137|50994697|
|max-struct-r8-physics_operators-part001.zip|7|19732380|
|max-struct-r8-scored_truth-part001.zip|54|370073|
|max-struct-r8-training-part001.zip|3|71363|
|max-struct-r8-training_data-part001.zip|4|30672525|
|max-struct-r8-vector3d-part001.zip|64|79009718|
|max-struct-r8-warm_anchor_confirmation-part001.zip|766|472330846|
|max-struct-r8-warm_anchor_confirmation-part002.zip|780|475198905|
|max-struct-r8-warm_anchor_confirmation-part003.zip|769|467415063|
|max-struct-r8-warm_anchor_confirmation-part004.zip|754|479374722|
|max-struct-r8-warm_anchor_confirmation-part005.zip|765|477404414|
|max-struct-r8-warm_anchor_confirmation-part006.zip|224|122221098|
|max-struct-r8-warm_anchor_development-part001.zip|250|154992955|
|max-struct-r8-warm_anchor_development_extended-part001.zip|746|471543017|
|max-struct-r8-warm_development-part001.zip|917|477161173|
|max-struct-r8-warm_development-part002.zip|553|308362797|
|max-struct-r8-warm_development_analysis-part001.zip|8|3473679|
|max-struct-r8-warm_primary-part001.zip|765|468677098|
|max-struct-r8-warm_primary-part002.zip|952|475797045|
|max-struct-r8-warm_primary-part003.zip|836|477728048|
|max-struct-r8-warm_primary-part004.zip|908|476344261|
|max-struct-r8-warm_primary-part005.zip|863|479995497|
|max-struct-r8-warm_primary-part006.zip|988|479083738|
|max-struct-r8-warm_primary-part007.zip|961|479999880|
|max-struct-r8-warm_primary-part008.zip|903|479477346|
|max-struct-r8-warm_primary-part009.zip|934|479572673|
|max-struct-r8-warm_primary-part010.zip|903|472241592|
|max-struct-r8-warm_primary-part011.zip|1024|476385034|
|max-struct-r8-warm_primary-part012.zip|900|472597367|
|max-struct-r8-warm_primary-part013.zip|854|477316610|
|max-struct-r8-warm_primary-part014.zip|900|479189304|
|max-struct-r8-warm_primary-part015.zip|847|465332997|
|max-struct-r8-warm_primary-part016.zip|757|479294882|
|max-struct-r8-warm_primary-part017.zip|894|478026247|
|max-struct-r8-warm_primary-part018.zip|956|470177239|
|max-struct-r8-warm_primary-part019.zip|856|479779712|
|max-struct-r8-warm_primary-part020.zip|969|479878380|
|max-struct-r8-warm_primary-part021.zip|938|479479776|
|max-struct-r8-warm_primary-part022.zip|819|471087312|
|max-struct-r8-warm_primary-part023.zip|893|470931126|
|max-struct-r8-warm_primary-part024.zip|936|467114693|
|max-struct-r8-warm_primary-part025.zip|726|318698404|
|max-struct-r8-warm_scoring-part001.zip|2|72726513|

完整原prediction(j/Q等)、轨迹/预算/checkpoints、训练数据及physics算子保留于release原路径。Git的compact为独立新文件，不改写源文件，不参与封存。默认所有private/truth路径HOLD且不读取。显式--include-scored-truth仅在750/42与126/42两套完整seal、各自评分回执和750/126行SCORES都通过后把自有GT单独分包；原Chat和第三方raw永远排除。Fresnel原Y只从官方重取代码取得，不分发原测量文件。

最终执行：`python delegated/cost_package/package_release.py --execute --finalized --archive-dir /absolute/outside/R8/release`。
执行器拒绝覆盖原档案/compact；zip生成后只核对成员路径、文件数、尺寸和CRC元数据，不执行SHA。恢复时在新目录按repo相对路径解压。
IGNORE_RULES.txt为本轮明确文件路径建议；不是已修改的gitignore。新实验出现后必须刷新。
