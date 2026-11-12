# R13 原始附件：发布工具与恢复合同

当前仅完成工具与本地清单规划，**没有由本工具发布R13 Release，也没有归档新的封存输入/真值**。真实上传、终态验收和最终公开入口由父线程完成；在此之前不列假下载地址。

## 包装入口与已做检查

`src/package_artifacts.py`默认是metadata-only `plan`：只读终态queue/ZIP中央目录，列出既有campaign archive与runtime包，不重压缩、不切分、不读取NPZ数组。当前规划回执位于 `delegated/evidence/release_tooling_plan_v1/PACKAGE_PLAN.json`：10个实际已有资产、2个尚缺的已完成本地队列归档、0个被排除策略阻止的已有archive。缺的是development_v3_local与development_v4_local完整包；这不代表其他远端队列已经完成或全部结果已下载。

临时非研究字节夹具验证了分片顺序、4/4/3字节连接、原文件保留和未授权拒绝逻辑。没有用研究数据执行prepare、下载或上传。

## 包的边界

1. 已结束队列的 `archives/*_full_results.zip` 与各自ARCHIVE_RECEIPT：完整RESULTS、queue、必要失败/恢复证据。
2. 实际冻结 `runtime*.zip`：保留原ZIP及其来源记录。不同runtime有不同运行版本，不能覆盖合并成一个“统一实际源码”。
3. 获授权后才生成 `max-struct-r13-postseal-inputs-sources.zip`：本研究公开输入、封存真值、R13源码/config/自有推导/原Pro研究包与来源元数据。内嵌 `RELEASE_REPRODUCIBILITY_SOURCE_INDEX.json`，保留旧R8 NN权重路径及继承代码依赖；不宣称本轮重训或独立重写Maxwell后端。
4. 获授权后才生成 `max-struct-r13-postseal-analysis-tables.zip`：真实analysis的JSON/CSV/JSONL/TSV/NPZ大表。代码、图和文稿继续留Git；大RESULTS、公开NPZ和巨大候选表走Release。

明确排除：私人chat readback原始载荷、敏感key/env/credential文件、`.git`、venv、cache、虚拟scorer夹具、原始第三方`.exp`测量、第三方全文PDF。自有FIGURES与原Pro研究报告PDF保留；source-access元数据、自有证明和用户研究任务保留。符号链接不跟随，不把旧transport ZIP再次塞入source包。工具不会写入R12目录。

## 终态与授权后准备

先检查实际已结束的campaign是否已用现有 `src/archive_campaign.py`归档。不要为了打包重跑科学，也不要在正在写入的队列上归档。

```text
python src/package_artifacts.py --root R13_ROOT --out NEW_PLAN_DIRECTORY --mode plan
```

新的input/truth/source或analysis包需要父线程授权文件。该文件须明确 `authorized_release_packaging=true`、列出**全部实际configs声明的confirmation_*与physics_v1 campaign**，并提供对应seal映射（或使用默认audit/CAMPAIGN_seal.json）。脚本检查每个seal的parent授权、终态queue、原文件字节/mtime快照；漏掉任何声明campaign会拒绝。

```text
python src/package_artifacts.py --root R13_ROOT --out NEW_EXTERNAL_PACKAGE_DIRECTORY --mode prepare --include-postseal --authorization PARENT_PACKAGING_AUTHORIZATION.json
```

prepare拒绝已有输出目录，保留所有原始ZIP，不删除来源。大于1,000,000,000字节的逻辑包按原始二进制顺序拆成不超过该上限的 `.part001`、`.part002`等；较小包本身就是一个运输part。不重压缩原campaign ZIP，不计算SHA。记录每part精确字节、顺序和总和。

`TRANSPORT_MANIFEST.json`采用既有R12的 `logical_assets → name/bytes/parts → name/bytes`合同，未公开时 `published=false`且不写URL。每个新ZIP有完整成员/字节index；`PACKAGING_RECEIPT.json`记录实际prepare CPU/wall，内部逐包和split计时不应再加到whole过程上。

## 真实上传后才创建DOWNLOADS.json

该脚本**不上传资产**。父线程实际上传后，提供逐part的真实URL、精确bytes与 `verified=true`回执。再显式调用：

```text
python src/package_artifacts.py --root R13_ROOT --mode downloads --transport-manifest PREPARED/TRANSPORT_MANIFEST.json --upload-verification ACTUAL_UPLOAD_VERIFICATION.json --published
```

脚本仅接受既有项目的稳定公开GitHub Release下载URL，不接受猜测/临时凭据URL。它对每个URL做小Range公开尺寸回读，精确匹配part bytes后才在R13根创建 `DOWNLOADS.json`。没有 `--published`或逐part真实回执时拒绝；不会把计划URL放进网站。此尺寸检查不冒充SHA或完整内容哈希认证。

`src/build_site.py`会在真实DOWNLOADS.json出现后收录它；此前只有当前说明文档，不产生“已公开”卡片。

## 恢复到新的外盘目录

R13入口 `src/restore_release_assets.py`包装原有 [R12恢复实现](../research_max_struct_r12/restore_release_assets.py)，默认manifest是**R13的DOWNLOADS.json**，显式传入R12 helper，因此不使用/改写R12输出路径。不安装依赖、不跑实验、不检查SHA。

```text
python src/restore_release_assets.py --asset LOGICAL_ARCHIVE_NAME.zip --output-dir NEW_EXTERNAL_DOWNLOAD_DIRECTORY --extract-dir NEW_EXTERNAL_EXTRACT_DIRECTORY
```

`--asset`可重复，用实际manifest里的逻辑包名。恢复器按manifest顺序连接运输parts，逐part与最终逻辑包检查精确字节数，保留完整下载ZIP。可选extract把每个逻辑包解到独立子目录，保留不同runtime版本，避免相互覆盖。输出必须为新/空目录，extract必须为新目录；拒绝绝对/traversal ZIP路径和symlink成员。

先保留各包独立来源，再在新的独立R13工作目录中按实际运行source/runtime版本组织需要的campaign与input/analysis。不能把多个不同版本runtime的src直接无说明覆盖。继承的R11/R8代码与原NN权重仍依赖已有项目来源；嵌入source index及runtime冻结记录说明这一关系。

**原seal中的绝对路径和mtime不具可移植性。** 恢复后必须在新的独立目录重新用seal_results记录当前预测文件和新root，父线程核对恢复来源后另外授权；不能直接把旧absolute seal视为新目录已通过。重新封存不要求重跑科学，也不引入SHA步骤。
