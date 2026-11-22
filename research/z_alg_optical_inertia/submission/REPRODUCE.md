# 论文、数据与复现入口

[英文论文 PDF](manuscript_ip.pdf) · [发布页 PDF](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-submission-20261122/manuscript_ip.pdf) · [同目录 LaTeX 源包](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-submission-20261122/Z_ALG_Optical_Inertia_LaTeX_2026-11-22.zip) · [完整科研 ZIP](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-optical-inertia-submission-20261122/Z_ALG_Optical_Inertia_Submission_2026-11-22.zip) · [上一轮复现说明](../REPRODUCE.md)

完整包把论文、区间验证器、冻结复观测、误差球、结果和历史依赖放回原有的 `research/` 路径。另一个 ZIP 装的是单层目录 LaTeX 源。英文学术稿的 PDF 可直接阅读。使用 LaTeX 源时，从解压目录执行 `tectonic -X compile manuscript_ip.tex` 或 `latexmk -pdf manuscript_ip.tex`。 [排版与审稿选项](typeset/README.md)

## 从选定结果开始复算

在科研附件根安装 NumPy、SciPy、python-flint、threadpoolctl、matplotlib。进入 `research/z_alg_optical_inertia/submission` 后可分别运行下列命令。原始 JSON 记载原机费用。重算会更新部分同名回执，如需保留原始计时请使用一份独立解压目录。

| 想核对的结果 | 输入、输出与命令 |
|---|---|
| 全局双约束单源 SDP 的严格可行性 | `python -B global/verify_global.py`。输入是已存观测和原生 96² 几何，输出逐观测的严格余量。 [区间回执](global/certificates.json) |
| 2×2／4×4 空间分块下的单源和联合界 | 从 `single/` 执行报告列出的冻结提案回放。四宫格逐一复核局部矩约束。十六宫格回放正谱及匹配联合结果。 [中文对照与复现](single/REPORT_ZH.md) · [严格正性与残差回执](single/coherent_general_verified.json) |
| 30 份配对全环观测与 150 个噪声／上限单元 | 从 `curve/` 按报告运行。可从已封存的观察数组直接验证乘子。只有重生配对噪声与可行性审计时会读取原干净仿真场。 [噪声曲线与命令](curve/README_ZH.md) · [观测清单](curve/observations/manifest.json) · [逐项 CSV](curve/noise_curve.csv) |
| 原生高阈值四宫格的负惯性障碍 | 按顺序运行 `python -B highcaps/screen.py` 和 `python -B highcaps/verify.py`。几何本身给出严格结果，原始读数无需载入。数值扩展见 `highcaps/skew.py`。 [证明与适用范围](highcaps/REPORT.md) · [区间回执](highcaps/certificates.json) |
| 普通 Gauss–Newton 重建与材料峰值证书 | 用 `diagnostic/run.py` 保存每次付费重建和复端口预测，再运行 `diagnostic/audit.py` 区间校核。两份终点和对应上限证书均可从公开输入复算。R10 真材料沿用封存约定，评分标量和已成图公开。 [复现细节](diagnostic/REPRODUCE.md) · [英文结果](diagnostic/RESULTS_EN.md) · [区间数据](diagnostic/interval_and_offline.json) |
| 旧半环、区域定位、双频与 Fresnel | 这些模块沿用旧实验证据。完整附件包含旧合成算子、冻结观测和输出，以及 Fresnel 重算代码。原始 `.exp` 测量文件从采集团队提供的档案取回。 [旧实验复现说明](../REPRODUCE.md) · [Fresnel 误差账本](../fresnel/REPORT_ZH.md) |

## 重建研究网页与发布文件

网页从 Markdown 生成，沿用上一轮样式和 MathJax 排版。用以下命令只预览已经写入的页面。

```sh
uv run --with markdown-it-py==4.2.0 python -B research/z_alg_optical_inertia/submission/publish_pages.py --preview
```

在项目工作区中，使用定稿 PDF、上一轮完整附件与已有 B1 检出目录，可执行下列命令生成网页、两个 ZIP 和单独的论文 PDF 资产，并把网页和导航写入 B1 检出目录。命令不会发布 Git 提交。参数支持指定旧 ZIP、发布文件目录、B1 检出目录和 release tag。完整科研 ZIP 已包含复现实验所需的上轮输入，无需重新生成旧 ZIP。

```sh
uv run --with markdown-it-py==4.2.0 --with markdown==3.11.3 python -B research/z_alg_optical_inertia/submission/build_release.py
```
