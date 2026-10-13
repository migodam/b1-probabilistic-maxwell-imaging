# B1 可视化交付审计

已在既有脚本上完成布局修复和证据范围修正。输出为 8 张 200 dpi 科学图与自包含中文 `dashboard.html`（10,532,253 bytes）。仅使用保存结果和派生数组，没有正演、推断重跑、训练、安装或锁定/历史数组读取。没有读取或修改共享成本账本。

## 数据与概率范围

- 主竞争情形保持：H1 **50.6418%**、H2 **40.4265%**、H0 **8.9317%**。默认 Top-2；H0 是实际剩余建模质量，Top-2 没有重新归一化。
- 全部质量条件于声明的三个模板，三项合计为 1。该先验按定义在模板外赋 0；更广泛物理几何的后验质量**未定义、未量化**，不能当作第四分支，也不能解释成物理不可能。
- 圆心与边界固定；图为三维同心球的二维截面，未进行位置推断。各分支的代表场使用条件均值材料参数；不把事件概率说成精确图像概率。
- 代表图重建只做确定性几何映射；12 个情形的加权代表图与已保存 `posterior_mean`、`void_probability` 的逐像素最大差均为 **0**。
- 原始全网格 5/7 的四个 `REFERENCE_UNRESOLVED` 保留。细化图单独显示 `EMPIRICALLY_REFINED_REFERENCE`，没有统一插值余项界或尾部误差证书。经验细化不改写原始失败。

## 静态图与可读性

1. `fig01_competitive_representatives.png`：三个模板及事件边际均值的 Re/Im ε 截面；实部背景 1、虚部背景 0；同一数量共用色标；显示四位小数百分比。
2. `fig02_top2_remaining.png`：竞争 Top-2 与 H0 余量；六主情形精确排序表避免极小堆叠段上的字溢出。
3. `fig03_competitive_cavity_uncertainty.png`：空腔占据与复 ε 边际标准差；修正了原先错误的实部 SD 色标；固定边界包含 0.24、0.46、0.5、0.6。
4. `fig04_cases_prior_negative.png`：六情形、先验敏感性、模板外真值负对照；极小非零质量以科学记数法显示。
5. `fig05_calibration.png`：24 独立母体 / 720 噪声；使用 native reviewer 的母体聚类 95% 区间。材料边际覆盖与已知真 H 的 oracle 条件覆盖分开展示；红外圈为分箱观测减报告差值区间不含 0（描述性、未校正多重比较）。不声明“校准通过”。
6. `fig06_reference_status.png`：原始全网格 3/5 与 5/7 状态；重要性采样概率差 0.02565 超过 0.02 容差；48 有限样本自助区间不认证遗漏模态。
7. `fig07_refined_reference_comparison.png`：原始与经验细化的概率并排，积分17/25与插值poly5/7诊断；原/细标签间距已修复。
8. `fig08_independent_forward_checks.png`：已保存的 24 点插值检查与 4 点独立 treams 实现对照。有限点检不推出全先验误差界；所有检查已在总预算384/384中。

**视觉检查完成**：通过 `view_image` 逐张检查全部 8 张图；针对原图的标题/字幕、极小质量标注、图例/柱值、页脚与色标错误作了修复；最后重看图7，没有发现文字遮挡。图内的物理与科学解释仍由父级 Codex 最终评审。自动图像检查不是用户审美或科学验收。

## HTML 范围与验证

- 6 主情形与 6 经验细化情形分组；默认主竞争情形 + Top-2 代表图；可查看全部3模板与事件边际均值。
- 空腔占据/材料不确定度随情形切换；第三模板余量单独列出，更广泛几何不画成概率柱。
- 所有图像和数据内联，44 个 PNG 数据块均成功解码；没有远程脚本、图像或样式引用。
- `check_dashboard_logic.js` 在纯 Node 模拟文档元素上检查了12情形、2代表图模式、7放大按钮、默认概率和零外部资源；结果 `PASS_CODE_LEVEL`。它是代码级逻辑验证，不是浏览器渲染。
- **浏览器实机渲染与移动端视觉验收：NOT_RUN**。CUA安全策略拒绝 `file://`，并明确禁止替代浏览器/间接路径绕过；未尝试localhost服务或任何绕过。响应式规则存在，但不声称已实机验收。

## 执行成本

两个完整制图进程（含解释器启动/导入）：wall **19.764659625 s**，CPU **18.824556 s**。
最后一次制图：wall **9.023235500 s**，CPU **8.852960 s**。
代码级HTML逻辑核验：wall **0.080419334 s**，CPU **0.071934 s**。
上述实测阶段之和：wall **19.845078959 s**，CPU **18.896490 s**（不包含人工查看工具往返时间）。

**physics_calls=0, budget_consumed=0, nn_training=0, locked_or_historical_arrays_read=0**。BLAS/OMP/MKL/Accelerate线程数均设为1。成本回执位于本目录的 `receipt.json` 与 `render_runs.json`；既有工人的旧回执仅备存为 `receipt_before_visual_review.json`，未合并到本次测量。

## 复现

使用现有 Gaussian/.venv_nn 运行 `research/delegated/visualization/make_figures.py`，并设置四个线程环境变量为1。`run_reviewed_render.py` 已执行该设置和进程计时。脚本/审计/回执均仅写在本目录；成品仅写 `figures/` 和 `dashboard.html`。未修改 src/config/results/reports。
