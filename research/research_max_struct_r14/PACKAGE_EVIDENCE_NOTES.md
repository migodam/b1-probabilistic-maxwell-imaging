# 完整证据归档与恢复：机械交付说明

新增工具`src/package_evidence.py`；本次只执行本地pilot的读取和小型manifest写出，未生成大归档、未上传、未读GT、未重跑数值或评分、未计算SHA。现有pilot结果1175个文件共1,372,041,636bytes，最大单文件NATIVE_HISTORY.json为13,189,983bytes；没有发生单文件超限，不开发假设性切块重组。

## 已审查的两个旧入口

`archive_campaign.py`保留RESULTS、queues和config原件，NPZ/NPY/PNG/PDF/ZIP用ZIP_STORED，JSON压缩不改变文本内容；它不删除原件。但原650MB只是按原文件字节累计的近似分卷，不保证含ZIP开销后的实际大小，也不处理单文件超限；没有排除private/cache/token流程。新工具不改此旧脚本，而采用600MB保守目标，写完stat要求实际每卷≤650MB。

`export_scoring.py`明确把initial/raw/corrected/primary/secondary中的NPZ改成仅m数组，供已封存预测的离线评分运输；其CANDIDATE_POOL等按原样复制。这是衍生material-only包，不能代替完整J/Q/辅助场、PR历史、全checkpoint和付费候选证据。新工具完全不调用该脚本，不载入NPZ，也不cast/重算/补造未保存状态。

## 计划模式与最终执行

默认`plan`只扫描公开允许树及终态campaign，保存PACKAGE_PLAN.json、每包MANIFEST、FILE_INVENTORY.csv、RESTORE.md；每文件记录相对路径、完整字节数、mtime_ns及所在卷，不记录内容哈希。原始NPZ/JSON不打开为科学数组、不筛字段。所有入包文件统一ZIP_STORED，绝不重新压NPZ。

```text
python src/package_evidence.py plan --root /external/research_max_struct_r14 --campaign dualmetric_pilot12_v1 --campaign other_terminal_campaign --include-diagnostics --output /external/research_max_struct_r14/audit/packaging_plan/final_v1
```

包含：source/configs/runtime/公开NPZ包；每个终态campaign全RESULTS与queues（包括全部JSON/JSONL、日志、付费candidate pool、JOINT_FINAL、所有full checkpoint及已保存全字段endpoint）；可选analysis/audit/FIGURES完整诊断包，包括局部根和完整状态干预数组。source包带现有根文档；第三方全文、cached Chat、sealed/private目录绝不入包。最终请对**全部实际campaign**重建计划，包括失败/未成功task，不用pilot计划冒充完整最终结果。

```text
python src/package_evidence.py execute --root /external/research_max_struct_r14 --plan /external/.../PACKAGE_PLAN.json --output /external/release_shards
```

execute先核终态campaign及完整文件路径集合，逐文件核计划byte size/mtime，再逐字节写ZIP_STORED；源文件变化就拒绝陈旧计划。每卷目标600,000,000bytes，实际stat硬上限650,000,000bytes。单个原文件若超过600MB直接明确报错并保留计划错误回执；不悄悄丢文件。执行写PACKAGE_RECEIPT或ERROR，保留已付费失败尝试，不删除科学原件。本机macOS的大归档/恢复输出必须在/Volumes外盘，默认小计划也位于研究外盘，不用内盘/tmp缓存大包。

归档成员目录和大小核对属于字节机械核对，不是科学或发布验收。不运行ZIP testzip全包重复读取；恢复通过正常zipfile流式读触发ZIP自带CRC，无SHA。source和raw的config重复时恢复只接受逐块字节相等的原文件，不覆盖差异内容。

## 恢复及Git范围

RESTORE.md给出完整恢复命令。工具以新目标目录恢复原相对路径，拒绝越界成员、缺卷、额外成员或原件大小差异；无需numpy/优化器，只有标准库。下载URLs、上传、匿名可达性核验、最终data closure和科学验收由主线程完成，不预填虚假发布URL。

R14局部`.gitignore`已补：RESULTS/queues/predictions/exports/archives、NPZ/NPY/PT/ZIP/log、private_operational/sealed_truth/CHAT_READBACK和第三方fulltext。公开数据metadata、configs、紧凑分析CSV/JSON、图、作者文档和源代码仍可进Git；原始公开NPZ通过release源包提供。git check-ignore已确认关键raw/private路径忽略而package_evidence.py不忽略。当前R14无已跟踪文件；ignore只防新git add，若主线程之后强制跟踪大文件，需另行检查，不能认为ignore会移除已跟踪内容。

本工具的token检查只扫描候选文本中的保守credential字面模式，命中时拒绝打包并仅输出文件名，不打印secret，也不静默改写/删证据。它不替代父线程最终权限/隐私审查；未触碰私有树内容。

## 文本压缩修订（以后以此为准）

JSON/CSV/MD和其余文本后缀现用ZIP_DEFLATED level1，NPZ/NPY/PNG/PDF及其他二进制仍STORED。前文“统一STORED”是旧实现说明，已被此修订替代。每file manifest记录archive_compression/compression_level；执行核目录compression，恢复也逐file核对。解压原始文件bytes不变，不读科学数组，不转换数值。仍按raw600MB目标分卷；如果单个文本本身超过600MB（实际已有1.8GBJSON），单独放一卷、实际压缩大小运行后必须≤650MB，否则明确失败保留原件，不开发片段重组或进一步优化。

Git额外忽略已实际发现的大表NATIVE_STAGE_HISTORY、NATIVE_STEPS、PAID_CANDIDATE_HISTORY、PAID_INCUMBENT_CURVES.json、SOURCE_CONFIGURATIONS.json和TYPED_STATE_FIELDS.json；原件保留外盘并由include-diagnostics计划进入完整release。
