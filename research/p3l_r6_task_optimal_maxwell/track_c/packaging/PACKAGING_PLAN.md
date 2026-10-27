# 关闭后分阶段打包计划

**当前状态：PLAN ONLY。锁定 D 正在运行；尚未启动归档、压缩或 raw 递归扫描。**

目标是逐阶段保存本轮已经产生的原始文件、失败尝试、独立账本、图、模型和报告，保持原始字节、路径、符号链接与文件元数据。打包不是科学验收，也不会改变 PASS/FAIL/HOLD/NOT_RUN、D 的状态或 GO/TAP 判断。

## 执行顺序

1. 根线程完成 D 物理、锁定评分和最终报告；写入正式 `CLOSURE.json`，确认中央没有活动 job/action。
2. 根线程显式执行 metadata inventory，先 A，再 B，再 C；D 可在其最终证据关闭后同样处理。脚本只收集名称、类型、大小、mtime、权限与链接目标，不计算 SHA，也不读取模型/真值内容作新分析。
3. 根线程检查 inventory 的范围、大小和外部依赖声明，再显式执行单阶段 tar.gz。每阶段顺序执行，避免多个压缩任务相互竞争资源。
4. 保存完整 inventory、归档和独立资源收据。失败的 `.partial` 留存，不覆盖已有正式包。

## 每阶段范围

- 对应 `track_a` / `track_b` / `track_c` / `track_d` 的全部本轮文件；C 的 packaging 输出目录本身排除，避免递归包含归档。
- `raw/jobs` 中 `START.json.stage` 与该阶段匹配的真实 job 全部文件。
- 共同的 `configs`、`shared`、`src`、中央账本与关闭收据；这些是原始配置/成本/源码关系所必需的上下文。
- C 附带 CDEV64 的 online 与 offline_truth 原始 NPZ，以及生成收据；历史 TRAIN512 分离的 online/full/labels、所有匹配模型和失败收据已在 track_c 范围内。
- D 按 `parents_D*.json` 的 parent 名单包含相应 online/offline_truth，只有正式关闭后才允许打包。

原始 R4/R5 的外部不可变源与缓存不会被伪装成本轮新数据。包中保留已有引用/来源清单；默认包**不宣称完全自包含的冷重建环境**。是否另打历史依赖包由根线程根据最终公开范围决定，脚本不会擅自复制整个历史仓库或历史 LU。

## 未来执行命令

以下命令只是计划，没有在当前 D 计时期间执行。

```bash
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/track_c/packaging/package_stages.py --stage A --mode inventory
../../.venv_nn/bin/python research/p3l_r6_task_optimal_maxwell/track_c/packaging/package_stages.py --stage A --mode pack
```

将 A 改为 B/C/D 可选择其他阶段。默认输出在本目录 `output/`。没有 `CLOSURE.json` 或中央仍有活动工作时，inventory 与 pack 都会停止；`--mode plan` 只打印范围，不递归扫描。

## 明确限制

- inventory 是文件元数据快照，不是物理、统计、模型或内容审查。
- 打包只允许在关闭后运行；本次准备不消耗 D 计时窗口的压缩 CPU 或大规模 raw I/O。
- 不执行 SHA-256 或额外内容哈希。压缩由标准 tar/gzip 格式完成，不修改源文件。
- 源文件在 inventory 与写包之间发生大小或 mtime 变化会停止该包，保留失败收据，不能悄悄覆盖或重新定义证据。
