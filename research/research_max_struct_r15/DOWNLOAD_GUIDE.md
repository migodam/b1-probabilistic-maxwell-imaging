# 下载与恢复

下载清单在网页及`DOWNLOADS.json`。只有实际发布且匿名回读成功的附件才进入该清单；不把草稿release或计划文件当公开结果。

## 全量科学原始文件

原始归档V2共16个ZIP：`r15_main264_raw_part001.zip`至`part014.zip`，以及fork56和source/public两个分组的实际附件名。合计7,364,794,746 bytes，对应10,384个原始文件、8,723,094,459 uncompressed bytes。每卷是独立标准ZIP，不是需专用软件连接的分割压缩流。

将16卷全部解压到同一个**新的R15根目录**，保留内部相对路径。分组各有文件清单与恢复说明；不要只下载main的一卷后就认定完整。完整raw保留264任务的材料/field/current/PR/checkpoint、候选池、56fork和公开输入/源码。old-endpoint GT衍生诊断、后来离线分析和最终报告由其他实际附件及Git提供。

同一release另有实际继承库快照、全状态分析、队列/封存/归档metadata。具体文件名、字节数以`DOWNLOADS.json`为准，不用输入source ZIP冒充原始实验结果。ZIP本身的常规CRC属于无损解压，不要求SHA校验。

## 快捷Pro研究包与关键状态选集

快捷Pro包提供最终报告、准确公式、全部摘要/逐对象表、图、公开数据、旧端点离线诊断及选定完整状态，方便下一轮数学分析。它不是全量候选/周期raw的替代。实际成员和排除项在包内manifest、PRO_PACKET_CONTENTS及构建回执。

关键状态选集的原独立归档为`R15_PRO_KEY_FULLSTATES_V1.zip`，188,513,757 bytes、80成员：4对象×5base fork的完整初末状态及12个背景cache，4对象×G19/S19/邻域的主final＋固定cache，以及4个真实checkpoint512源。全部68个非null缓存引用齐全；8个hop/reset不需背景cache属于正常缺失。

独立关键状态选集的路径相对R15根；快捷Pro包路径相对仓库根。恢复到对应层级，勿重复嵌套。Windows fork元数据中的反斜杠需要POSIX读取器解析时转换，原状态字节不改。完整运行依赖public输入、runtime及实际继承库，独立选集本身不声称包含这些；快捷Pro包则列出它纳入的依赖。

## GT与保存范围

`data_public`只含公开观测和合法先验。旧端点材料误差诊断及R14已有GT仅用于封存后的离线分析，并在Pro包标明OFFLINE_SCORING_ONLY；不传入在线算法。任何没有secondary候选的对象仍保留primary回退和缺失原因，不从24分母删除。

原封存时间早于fullstate恢复，旧seal是历史证据，不能直接充当新绝对路径树的seal。需要重新评分时生成新的离线派生seal，保留原seal与当前报告。详见[复现说明](REPRODUCE.md)。
