# 主线程交接

已写 `experiments/build_cost_ledger.py`，标准库脚本；仅读取明确列出的回执/manifest和路径stat，无GT数组读取、评分或物理调用。未改其他研究报告、gitignore、R1–R7，未打包、未上传。

刷新命令（从R8研究目录）：

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" experiments/build_cost_ledger.py
```

输出：
- `COST_LEDGER_IN_PROGRESS.json`：明确入口的非重叠加和行、嵌套非加和分解、失败下界、缺口与队列历时context。始终IN_PROGRESS，主线程最终审查才可改完成结论。
- `COST_GAPS_AND_PACKAGING_ZH.md`：人读汇总、成本缺口、失败路径、发布清单提案。
- `PUBLICATION_INVENTORY.json` / `.tsv`：全部文件大小/类型、去向、分包具体成员。扫描当时仍有运行任务，必须完成后刷新。

机械验证：重复(source,record)不存在；warm setup从outer CPU/wall减去后另列，逐行重建原总量一致；所有已知CPU非负；local/raw/chat边界没有进入Git/release候选。未做SHA检查。

成本边界：
1. warm外层总CPU先扣setup，只将setup分成另一个加和行，aux/current/corrector/validation都是online内部成本，不能再次加。训练含每epoch验证；训练history累计timer不加。实际initializer批次和独立部署benchmark分别收费。
2. paper现有COST_LEDGER逐条导入一次；它已减去原64grid复用辅助的CPU/RHS。原DATA/FAILURE/INTERRUPTED不再重复导入。32grid共享装配从`setup_shared`仅按case/timer唯一值加一次，因源码在每臂timer之前装配。**主线程若新增paper实验，应先刷新paper_cost_ledger.py再运行本脚本。**
3. vector52臂只用ONLINE_COST_LEDGER一次；COMPLETE/PR_COMPLETE/共享SIMULATION都只作context。PR48最初已完成五臂后失败，resume总计不含这五臂，不能用进程总计减52臂总计。8³未缓存中断OS CPU344.46秒已含完成LSIE227.39秒，后者不再单加。原3λ整过程379.59CPU秒单列，不能伪拆成精确阶段成本。
4. Fresnel只计最终每臂RECEIPT；辅助回执/PROGRESS不是额外费用。四个工程中断按独立下界计，第一起精确CPU缺失。原3D序列化失败CPU未测，初次A+数组NOT_RETAINED。
5. 全部task-wall和只能称工作量，不是项目历时。队列batch/elapsed不加进CPU；不归一化不同硬件/维度的LU；金额和能耗始终未知。
6. 主test scorer只在`audit/PRIMARY_PREDICTION_SEAL.json`已存在时考虑白名单评分回执；当前未评分。若root最终采用其他真实seal/scoring路径，应明确添加它们到注册表，不能递归捞所有SCORES以冒充评分成本。

发布提案：Git保留代码/配置/紧凑最终结果/回执/失败证据/报告/图，history、CORRECTOR完整轨迹、stage/corrector中间数组、budget以及大文件转独立Release evidence分包；training data/weights另一个分包。采用R7约600MiB非压缩上限分包，不执行压缩。`sealed_truth`和`data/paper_benchmarks_private`统一HOLD到主线程最终封存与独立评分后审核。原Chat inputs、data_local原测量、第三方PDF/HTML/download及*.local.*、audit临时zip、SSH身份文件、本地vendor/cache严格exclude。

公开前还需主线程检查日志/回执中可能包含的绝对路径和连接元信息；本任务的发布扫描按授权只stat文件，不宣称做过内容脱敏审计。全部实际自有数据/失败/历史须被Git或Release覆盖，HOLD和exclude需最终回执列明。

## 精简Git与可执行分包（更新）

**发布去向以新 `ARCHIVE_INDEX.json` / `.md` 为准，取代早先PUBLICATION_INVENTORY的约220MB Git提案。** 新脚本默认只计划：

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" delegated/cost_package/package_release.py
```

主线程确认全部实验结束、审查seal后才执行（本worker未执行真实打包）：

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" delegated/cost_package/package_release.py --execute --finalized --archive-dir /absolute/outside/R8/release
```

需要同时公开本轮GT用于复现时，加 `--include-scored-truth`。此flag严格读取**公开JSON回执**核对：
- `audit/WARM_ALL_PUBLIC_SEAL.json` status `ALL_FROZEN_WARM_TASKS_SEALED`、tasks750、objects42、entries长度750；
- `RESULTS/warm_scoring/RECEIPT.json` status `SCORED_ALL_FROZEN_PREDICTIONS`、tasks750、objects42；
- `RESULTS/warm_scoring/SCORES.json` rows长度750。

全部满足才将自有truth文件放独立 `scored_truth` 分包；默认仍HOLD，计划及compact从不打开私有数组。GT分包发布后42对象成为公开开发对象，不再称新盲测。原Chat/第三方raw永远排除。成本脚本当前seal路径已修正为准确的WARM_ALL_PUBLIC_SEAL，history_warm单臂回执已加动态注册。

新计划将Git候选降到约26MB（不含执行后新compact与随后新增结果，确切数值见实时index）。完整RESULTS数组（含prediction j/Q等）、完整history/budget/corrector/checkpoints、data_training及大physics算子保留原成员路径至release；Git留代码、配置、全部receipts、scalar scores、图、权重、public小观测，以及**独立新文件** compact_endpoints。compact仅复制已有m/Y/prediction/gradient白名单，原值和精度不改；源文件无gradient就不捏造，不调用物理。所有原文件不改写、不删除。

按研究分支/相对路径排序，每包源字节最多480,000,000，生成后ZIP强制≤500,000,000。单一源文件>480MB会明确拒绝，不能悄悄超限或删文件。无hash检查。归档后检查member路径与尺寸；在新目录解压即可按repository-relative路径恢复原结果。生成IGNORE_RULES.txt只是建议，不修改gitignore；已tracked文件是否转release由root审查Git状态处理，不能删除磁盘原件。

已用临时合成fixture验证执行器：compact保留m/Y/gradient、剔除J；原源文件bytes/mtime完全不变；ZIP原成员可恢复相同bytes；默认truthHOLD；缺seal拒绝truth路径；严格750/42公开JSON gate允许计划。fixture不是实际实验，未打开本轮GT，未进行真实研究文件压缩/上传，临时目录已自动移除。

## 新开发队列/双频更新

已登记 `configs/warm_anchor_development_manifest.json` 对应8个 `RESULTS/warm_anchor_development/*/RECEIPT.json`，与750主任务分开；`RESULTS/warm_anchor_development_analysis/RECEIPT.json` 是独立validation分析费用。未出现的receipt保留expected gap，不编造零。开发anchor的warm外层setup仍扣出单列，其他nested费用不再加。

history公共manifest动态读取21任务；新增 `data_public/history_warm/PREPROCESS_RECEIPT.json`（一次features+load+inference）、21单臂外层、`RESULTS/history_warm/SCORING_RECEIPT.json`。`ALL_PREDICTIONS_SEALED.json`及该分支重建COST_LEDGER只作context，不重加已读单臂或跨R6/R7旧baseline费用。

`experiments/fresnel_multifrequency.py` 约定两个mode `control_8` / `joint_6_8`，输出 `RESULTS/fresnel/FoamTwinDielTM_multifrequency_<mode>/RECEIPT.json`，按独立outer记录；setup/optimization/两频heldout validation都包含其中，不重复计时。physical_cost_optimizer和physical_cost_validation按frequency保留，绝不把6/8GHz的动作表加成另一组timer。anchor使用原8GHz已支付预测，不再加旧run费用。实际FAILED outer可单独收费；同目录INTERRUPTED下界有FAILED时仅context。MULTIFREQUENCY_BATCH_RECEIPT包含等待前序任务与两个孩子，不加和。

新packaging计划自动覆盖上述manifest和全部明确RESULTS分支；multifrequency compact白名单额外复制已有prediction6/prediction8（都是模型预测，不是原测量Y）。仍默认truthHOLD、未读GT、未生成真实compact/ZIP、未发布。当前机器快照与核对见 `ADDITIONAL_REGISTRY_UPDATE.json`，最终完成后再次刷新。

## 冻结126确认扩展与24验证补充

新增明确回执路径 `RESULTS/warm_anchor_confirmation/*/RECEIPT.json`（126任务，与原750分开）及 `RESULTS/warm_anchor_development_extended/*/RECEIPT.json`（24任务，与pilot8分开）。各任务完整外层费用收费，setup仍从同一outer扣出单列；失败也保留。补充状态从对应manifest动态算，不写死本地完成数。分析器真实路径登记为 `delegated/warm_anchor_analysis/{PILOT,DEVELOPMENT,CONFIRMATION}_RECEIPT.json`；若没有计时字段，按未测缺口保留。

远程独立checkout的 `audit/ANCHOR_XINAN_PREPARED.json` 单独记录已测transfer wall（CPU为空），不把将来8lane队列工作再加一次。其他ANCHOR_XINAN控制回执作为非加和context；完整单任务结果从原RESULTS根接收后自动计费。`transfer_anchor_xinan.zip`及其他audit zip仍严格排除；不能以临时传输包替代保留正式Release成员。

发布GT guard已升级：只有原750/42 seal+scoring+750行SCORES，**以及** `audit/WARM_ANCHOR_PUBLIC_SEAL.json` status `ALL_ANCHOR_TASKS_SEALED`、126任务/42对象/126entries，`delegated/warm_anchor_analysis/CONFIRMATION_RECEIPT.json` status `SEPARATE126_ANCHOR_OFFLINE_SCORED`、126/42，以及CONFIRMATION_SCORES的126rows全部满足，显式 `--include-scored-truth` 才能放行。缺126 seal/scoring时即使750已完成也拒绝。默认继续HOLD。

126扩展全部原prediction/state/history/budget数组均按原相对路径进入Release，不以compact代替；所有公开JSON/代码/失败回执留Git或完整Release，不删原件。新增两个RESULTS分支已加至packaging显式roots。未实际读testGT、运行分析器、生成真实compact、压缩或上传；已用临时合成公开JSON验证双seal门控。状态见 `ANCHOR_EXTENSION_REGISTRY_UPDATE.json`。
