# 下一轮 Maxwell 实验从这里接手

**以数据族、完整观测 ID 和父运行定位问题。** 稳定入口是仓库与外置工作区根的 `AGENTS.md`。本页用于重新进入研究，耗时作业的最新状态看各轮本地 `TASK_RESUME.md`。[实验手册](README.md)

|已有判断|下次引用时带上的条件|
|---|---|
|A+-MM / Arc 在 R13 48 件上各恢复 36 件，Dual-CIE 32、CIE 30|共同时间前缀材料收益仍在。R13 的 48 件与 R11 的 48 件有不同身份；其 16 个半环材料全部同网格。跨网格半环的表现须另取观测|
|标准半环上 SRAP 相对 Dual 的全图误差有配对收益|30 秒 raw 21/1/2 胜平负；raw 内部误差中位数从 0.554 升到 0.754，背景从 1.768 降到 0.523。轻校正后的 `s2d_010` 双环和半环 `s2d_007` 过高对比同时值得保留。[完整材料图册](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/srap_cie_array_study_20261121/figures/confirmation_atlas_30s/)|
|三维矢量与实测已有真实证据|三维 12 件、两阵列、双极化、24³→12³；Mie 球从 12³ 到 24³ 误差约 7.95%→6.15%。8 GHz 实测 Twin 有 733 个训练复通道，另有 1,205 个未见源通道及 390 个未见接收通道，须按两个留出集分别报告|
|2019 原文停机阈值存在历史误读|原文 Eq.16 用逐源相对系数变化平方均值的平方根，阈值 0.001。R8 三维小型 Fourier-PR 两臂实际固定 240 次；R8 二维 Fresnel 早期和最终代码分开。不要把元数据修正写成重新证明了材料排名|

继续做算法比较时，先看[对照错题本](BASELINES_AND_ERRATA_ZH.md)和[设备与算子](OPERATORS_ZH.md)。需要新数据的依据可以具体到“已有数据缺少哪个配对条件”；若已存数组足以裁决，就使用[可运行的复算](WALKTHROUGHS_ZH.md)。这一轮建立了目录、走查和复算，新增反演为零。

维护者的更新顺序是先落本轮证据，再登记数据或算子，写下影响科学结论的修正，结清存储，最后提交。需要恢复旧工作时沿[项目系列索引](PROJECT_FAMILIES_ZH.md)和[上下文交接规则](DESIGN_AND_CONTEXT_ZH.md)回到原任务、原附件、运行目录和未完成作业。磁盘清单的本地入口为 `communication/EXPERIMENT_STORAGE_INDEX.md`；旧数据删除范围见该轮批准记录。[存储办法](STORAGE_ZH.md)

网页由本目录 Markdown、CSV 和 `render_site.py` 生成。维护者用安装了 Python-Markdown 的解释器运行 `python render_site.py`，将源文件及对应网页一起提交到现有 B1 站点。正文、登记表和公开复算结果均可从 GitHub 单独读取。
