# R13 Git/Pages 发布体积规划

只按实际路径、字节和Git index状态分类；没有加载NPZ/真值、查看新结果质量、staging、commit或push。当前root主分支状态不由本规划改变。

## 当前候选体积

|类别|文件数|字节|
|---|---:|---:|
|Prospective Git allowlist|612|30484790|
|其中：实际图/自有报告PDF|210|15278097|
|其中：代码/数学/文稿/license|141|2437313|
|其中：compact summary/audit metadata|261|12769380|
|Release only|4396|2467104407|
|Git与assets均排除|120|2123296|
|父线程另行判断|1|9316|

当前index已经跟踪、但不在prospective allowlist中的文件：0。完整路径在JSON中；本工具没有unstage/untrack任何文件。

## 规则与边界

- allowlist保留R13代码/config/自有数学和文稿、实际PNG/PDF、compact summary CSV、费用/audit/source receipts及license/attribution。保留原Pro自有报告和算法代码；没有把第三方全文当作自己的资产。
- 所有RESULTS/raw states、public/sealed NPZ、runtime/原始archive、巨大candidate/history/paid-prefix以及逐对象大表走Release；不让Pages承载原始数组。
- 私人CHAT_READBACK/PRO_CHAT_READBACK载荷、fake fixtures、cache、env/keys、第三方全文PDF与原始第三方.exp均从Git和assets排除。SOURCE_ACCESS等来源元数据保留。
- DOWNLOADS.json只有实际上传验收后published=true才入Git；当前规划不制造Release链接。
- 单文件超过50,000,000B不自动入Git。compact JSON上限2,000,000B、CSV上限1,000,000B是发布体积审阅规则，不是科学成功门。
- 当前文件名/大小不证明科学结论、图像验收或独立对象数量。父线程最终选择精确文件集合。

## 超过50MB的实际文件

- `research/research_max_struct_r13/analysis/development_v1/ALL_PAID_AND_OVERRUN_CANDIDATES.json` — 57780498B；release_only；untracked_by_index
- `research/research_max_struct_r13/archives/confirmation_local_reference16_full_results.zip` — 603921051B；release_only；untracked_by_index
- `research/research_max_struct_r13/archives/development_v1_full_results.zip` — 114088735B；release_only；untracked_by_index
- `research/research_max_struct_r13/archives/representation_v1_full_results.zip` — 78621625B；release_only；untracked_by_index
- `research/research_max_struct_r13/runtime_confirmation_C2.zip` — 68607341B；release_only；untracked_by_index
- `research/research_max_struct_r13/runtime_confirmation_C2_final.zip` — 68607527B；release_only；untracked_by_index

## 交接

机器可读精确路径集合在GIT_PUBLICATION_PLAN.json的git_allowlist/release_only/exclude_from_git_and_assets中。不要对整个R13目录、RESULTS、analysis或data_public执行宽泛git add。先按此计划审阅，再由父线程显式选择；本规划没有改变任何科学文件或Git index。

规划文件自身未计入候选体积，以避免自引用；新增结果/图生成后应按相同规则复核最终文件集合。详细重复cost inventory表改走Release，compact COST_SUMMARY及权威queue/audit回执继续候选。

当前唯一另行判断项是原Pro包的历史SHA256SUMS文本：本轮没有执行或验证它。原始启动日志归Release，活动handoff lock从Git/assets排除。

## 公开选择索引与操作隐私补充

公开的 `predictions/min_change_paid_pool/**/*.json`（≤2,000,000B）可进入Git精确allowlist；所有预测NPZ仍仅走Release。该索引记录预GT的付费候选保留选择，不是新搜索、同目标优化提升或节省corrector费用的证据。

`private_operational/**`、`audit/transfers/**`、`audit/REMOTE_TRANSFER_PROGRESS.json`、`audit/remote_metadata/**`以及两类CHAT_READBACK载荷一律从Git与所有assets排除；活动文件没有移动或修改。此次仅对原快照及这些新增目标路径做元数据更新，不冒充当前全目录重新盘点。

父线程已接受compact `delegated/physics/PHYSICS_CONDITIONS.csv`作为Git科学条件表，解除先前Release分类；未执行Git index修改。
