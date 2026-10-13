# A32 研究包

研究日期：2026-10-12。

## 阅读入口

先读 `A32_NEXT_DECISION_ZH.md` 获得下一阶段裁决，再读 `A32_THEORY_REPORT_ZH.md`。完整报告含目录、原始证据审计、材料概率空间、T1–T10、证明与反例、文献原创性地图、四条算法路线、实际小规模核验及最小实验门槛。

核心判断：推进有限尺度 topology-stratified material posterior；研究 A17 双侧材料因子对后验 likelihood ratio 的精确作用。先建立可信参考概率，不启动大型生成网络，不以 A30 已失败的全域证书作为必须成功的主流程。

## 文件

- `A32_THEORY_REPORT_ZH.md`：完整理论报告。
- `A32_NEXT_DECISION_ZH.md`：下一阶段的主线、停止条件与最先验证。
- `a32_checks.py`：独立 NumPy/SciPy 小矩阵与解析二维 TM 检查。
- `A32_CHECKS.json`：当前14项检查；区分代数、有限 prior、解析 scalar TM 与真实成像。
- `A32_CHECKS_INITIAL.json`：初版结果，保留其过易双频例与概率 bin 浮点边界问题，不作为校准通过证据。
- `A32_RESOURCE_LEDGER.json`：两轮解析调用共74次；80次累计上限。
- `A31_RECHECK.json`：原 A31 22项检查的独立复核。
- `upstream_a31/a31_matrix_checks.py`：用户上传 A31 原脚本的未修改副本，仅用于复核；不是 A32 新算法。
- `MANIFEST.json`：包内文件哈希、原始附件哈希及 A30 固定版本。

## 本轮真正完成了什么

A31 原22项检查重新运行，22/22；A32 新14项检查通过。A30 四项父样本 gain 使用公开 CSV 的已保存评分重新聚合，不是重跑其优化器或原 Maxwell 实验。A32 的物理演示是新建的二维无限介质圆柱 TM 全波解析解，不是三维向量成像。

新解析材料–频率调用累计74次，包括初版、重复和控制。A31 复核另有一个9维微型矢量偶极状态/2 RHS、两个解析 slab 案例，详见其 JSON。没有 GPU、训练、历史 truth-array 访问、locked-set 访问或 current-selection 扩展。最后一次 support-bound 检查加强仅重跑代数，新增物理调用0。

## 独立复现

依赖：Python、NumPy、SciPy。无需网络、CUDA、用户数据或外部仓库。

当前结果目录保留累计账本。在这里直接重跑整个 A32 脚本会超过80次累计上限并被拒绝，不要删除账本伪装成免费重跑。

独立复现请只把 `a32_checks.py` 复制到一个明确声明的新空目录，在该目录执行：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python a32_checks.py
```

新的运行生成新的 `A32_CHECKS.json` 与 `A32_RESOURCE_LEDGER.json`；当前单频脚本一次含28次解析调用。结果的主要数值应可复现，墙钟因设备变化不要求一致。

复核 A31 同样将其原脚本复制到单独新空目录再运行，避免覆盖已有记录。两个脚本的物理定义和计费口径不同，不能将 A31 的小矩阵状态与 A32 解析频率调用混称大规模 Maxwell 性能。

## 证据边界

证明针对文中明确假设；浮点核验不是区间算术认证。18-state 的有限 prior 概率检查不是任意连续材料后验，也没有证明真实世界校准。当前未建立三维复杂拓扑概率、线上 feedback 加速或同成本算法优势。

对原始项目与文献的阅读范围逐项列于报告。A22 独立档案未取得；DiffusionEMIS 和部分其他论文只有摘要/元数据。没有把这些未读内容写成不存在。候选理论的发表优先权未确认，不声明 first-ever。