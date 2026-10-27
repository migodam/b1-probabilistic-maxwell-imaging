# R6 原始证据、紧凑阅读包与恢复

报告、图和源码复用现有B1。大型数值证据以独立无损分卷提供，不塞入Git历史。数据都是本轮合成研究数据；历史正式预算未被覆盖。

## 阅读与下载

- [完整原始证据Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr6-task-optimal-20261027)
- [GPT Pro阅读ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3lr6-task-optimal-20261027/Research_P3LR6_Results_and_GPT_Pro_Handoff.zip)
- [复现入口](REPRODUCE.md)
- [512对象紧凑数值manifest](GPT_PRO_HANDOFF/COMPACT_NUMERICAL_MANIFEST.json)

阅读ZIP供下一轮推理，包含最终报告、关键代码、图、表与两条件紧凑材料/预测；不是全部原始物理轨迹或冷训练数据。完整raw、失败记录、field/current/J、模型、bootstrap和必要历史输入在Release分卷中。

## 原始文件的定位

本地稳定数据/代码/必要输入快照共410091个文件、13,142,518,072原始bytes；12个独立tar.zst压缩为8,370,745,338bytes，全部逐成员逐字节还原比对通过，原件保留。最终报告/图库/交付脚本另作小包；以上是稳定批次口径，不把它称为所有最终文件总量。

Release的INVENTORY与每卷manifest保留相对Inv_SLAM根目录的原路径。请按manifest与restore_lossless.py恢复到全新目录；不要覆盖历史工作树。每卷严格小于2GB。恢复校验不运行SHA-256；有原始目录时可以逐字节比较，没有原始目录时不得声称做过不存在的字节对照。

下列大数组只在证据卷中，不是网页丢失的数据：

- research/p3l_r6_task_optimal_maxwell/track_a/raw/COMPENSATION_VECTORS.npz
- research/p3l_r6_task_optimal_maxwell/track_a/raw/COMPENSATION_BASE_CACHE.npz
- research/p3l_r6_task_optimal_maxwell/track_a/analysis/POSTHOC_LATENT_Q_NODE_MOMENTS.npz
- 各track的raw/contexts/models、根RESOURCE_LEDGER.jsonl和完整预测/GT
- R6实际B执行代码原位置为仓库根track_b；science分支和归档保留原样，main的R6/repro_source/repo_root/track_b保留逐字节相同副本

最终源码提交、发布main提交、实际资产字节、上传及匿名回读状态由交付receipt记录。无损或上传成功不会升级概率资格和科研GO；所有HOLD/FAIL/NOT_RUN随包保留。
