# R8最终构建与发布交接（已打包，尚未上传）

实际compact/ZIP已于2026-11-06T00:11:17.019786+00:00完成：52包压缩后共5,308,028,340 bytes，1,126份compact；逐包成员/尺寸核对通过。未上传。依据：[打包执行回执](PACKAGING_EXECUTION_RECEIPT.json)、[归档核对](ARCHIVE_VERIFICATION.json)、[全部路径](ARCHIVE_PATHS.txt)。主线程拥有科学审查、Git提交、Release上传与匿名链接核验。

以下命令与运行前提醒保留为历史/复现说明，已执行的实际状态以上述回执为准；不得把历史“尚未执行”提醒当当前状态。

## 当前计费规则

- 主750与anchor126逐份实际run_case outer收费；setup从outer扣出另列，nested auxiliary/corrector/validation不得加第二次。迁移后的真实平台优先来自receipt平台/解释器，旧对象分配fallback扣除迁移7对象。最终分配目标为主Mac405/XINAN345、anchorMac66/XINAN60；远程465/81队列acknowledgments包含141个迁移记录，不能称实际物理求解。
- 8份缺threadpoolctl的实际失败进程独立计费，即使后续同task_id成功；失败队列总计不另加。
- projection_transport公开公式诊断outer只计一次，不重新读取状态或跑物理；Green/FFT/SVD动作保留。
- FETCH_HISTORY每行是一次真实远程压缩操作；remote_pack CPU/wall加一次。控制总wall减nested remote_pack wall得到运输/解压/SSH余项，localCPU为空。最新FETCH.json是镜像context，queue elapsed不是传输或项目历时。HARDWARE_REALLOCATION control只记已测4.125秒wall，CPU未知，绝不计物理。
- Fresnel当前23个成功outer自动登记；之后两份corrected-air Dual/A+全预算DEV结果若落RESULTS/fresnel/<run>/RECEIPT.json会自动增至25。所有old/new outer都保留；optimizer/current/setup/两频validation嵌套不重复收费。AIR0_FINALIZATION独立元数据QA CPU另外记录；含等候的批次wrapper只作context。
- 仍不能声称完整精确总CPU：原3D序列化失败未测、部分报告/调度/控制CPU未测，工程中断有下界；能耗/金额未知。NULL不是0。

## 最终刷新命令

在repo根运行（下面路径均为实际项目路径）：

```sh
cd "/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r8_worktree"
R8_PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
"$R8_PY" -m research.research_max_struct_r8.experiments.paper_cost_ledger
"$R8_PY" research/research_max_struct_r8/experiments/build_cost_ledger.py
"$R8_PY" research/research_max_struct_r8/delegated/cost_package/package_release.py
```

这是plan-only。不要仅凭队列dispatch结束触发GT评分。主线程核对750+126真实terminal、单独公开seal及未来25个Fresnelouter后，按自己的授权顺序运行评分；本交接不触发任何scorer。

只有两套750/126公开seal及各自完整评分后，root可用下列命令产生最终成本记录和实际Release包。`--finalized`代表root确认全部授权工作已完成，账本状态为FINAL_RECORDED_COSTS_WITH_EXPLICIT_GAPS，不把费用齐全等同科学成功：

```sh
"$R8_PY" research/research_max_struct_r8/experiments/build_cost_ledger.py --finalized
"$R8_PY" research/research_max_struct_r8/delegated/cost_package/package_release.py --include-scored-truth
"$R8_PY" research/research_max_struct_r8/delegated/cost_package/package_release.py --execute --finalized --include-scored-truth --archive-dir "/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r8_release_20261106"
```

最后一条已于2026-11-06T00:11:17.019786+00:00实际执行完成，原目录已有52包与compact，请勿对相同目录重复执行。复现时目的目录须为空/新建，不覆盖既有包或compact。真值flag要求：
1. audit/WARM_ALL_PUBLIC_SEAL.json: ALL_FROZEN_WARM_TASKS_SEALED、tasks750、objects42、entries750；warm_scoring/RECEIPT: SCORED_ALL_FROZEN_PREDICTIONS、750/42；SCORES750rows。
2. audit/WARM_ANCHOR_PUBLIC_SEAL.json: ALL_ANCHOR_TASKS_SEALED、126/42、entries126；delegated/warm_anchor_analysis/CONFIRMATION_RECEIPT: SEPARATE126_ANCHOR_OFFLINE_SCORED、126/42；CONFIRMATION_SCORES126rows。

该flag未满足就拒绝，默认全部truthHOLD。GT发布后42对象是公开开发集，不能再宣称盲测。

## 包与Git审查

ARCHIVE_INDEX.json为权威成员索引，IGNORE_RULES.txt只是建议，实际gitignore与已tracked处理由root执行。完整原prediction(j/Q等)、history/budget/checkpoints、training data和大operators按原相对路径放Release；独立compact只复制已有m/Y/gradient或预测6/8GHz值，不改原文件。每包源最多480,000,000 bytes，ZIP最终≤500,000,000 bytes，成员路径/尺寸校验，无SHA。下载全部包后在新repo副本根解压即可恢复原结果。

自有FIGURES/*.pdf与PNG留Git；第三方PDF/全文仍exclude。inputs原Chat、data_local原测量/原Y、*.local.*下载、audit临时ZIP（含transfer_anchor_xinan与incremental public_results）、SSH身份/.ssh/.env/credential文件、本地vendor/cache永远不发布。官方重取代码留存。当前defaultHOLD的truth需最终flag后独立scored_truth分包，不能被静默遗漏或宣称已公开。

已核对路径去向；已知credential marker文本检查范围限小型代码/回执/日志，不读取GT或完整状态，不构成未知凭据格式的全面保证。公开回执中的绝对路径/机器名为来源元信息，原件不可为美化而改写。发布前root确认Git/Release覆盖全部自有证据，所有exclude/HOLD列入最终回执；root负责真实上传链接和匿名访问验证，不用计划URL代替已发布结果。

## 最终费用Markdown与打包后更新

新增根COST_LEDGER.md由下列命令根据最终账与公开标量部署表生成（不读原状态/GT）：

```sh
"$R8_PY" research/research_max_struct_r8/delegated/cost_package/write_cost_report.py
```

最终真实pack执行完成后，必须补记实际打包CPU/wall，按次序重跑：

```sh
"$R8_PY" research/research_max_struct_r8/experiments/build_cost_ledger.py --finalized
"$R8_PY" research/research_max_struct_r8/delegated/cost_package/write_cost_report.py
```

PACKAGING_EXECUTION_RECEIPT只计一次outer；其成员/文件数不是额外成本。最终COST_LEDGER_FINAL.json现在明确留Git，允许这次打包后费用刷新；原已封存科学源文件不改动，实际ZIP成员不递归包含自身打包回执。ARCHIVE_INDEX和PACKAGING_EXECUTION_RECEIPT为真实包的来源依据，打包后不要用一份新的plan覆盖实际包索引来声称同一批ZIP已重建。

late audit/remote_imports仅≤256KB的小JSON/log/source文本留Git；是复制来源镜像，不是新增计算，也不进入Release科学数组。其他mirror payload排除，canonical RESULTS完整原件不删、不因内容重复而省略。目录量与唯一成员检查见PACKAGE_VOLUME_AUDIT.json，无内容hash或GT读取。
