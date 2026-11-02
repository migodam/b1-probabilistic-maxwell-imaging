# 复现：小模型、CPU、保留全部输出

冻结主仓库基准：7317601b43b914525c21f09e5a60020552ea677b。新代码独立存于本目录；numerics 需要的 A9 dyadic Green 与 D1/D3 冻结适配、4×4 公开算子已附带，不需要访问原工作目录或 private_truth。来源在 `audit/D_INPUT_PROVENANCE.json`、`research/delegated/implementation/BACKEND_AUDIT.md`。不计算 SHA；回执对应字段 null 并解释。

已用环境：Python3.13.13、NumPy2.5.3、SciPy1.18.1、Matplotlib3.11.1、threadpoolctl3.7.0；macOS CPU、complex128。版本是实际环境记录，不宣称只有这些版本可用。环境依赖见 `requirements.txt`。本轮复用了既有环境，未安装 CUDA、daemon 或新框架。

## 只读结果

直接打开 `REPORT.html` 或 `RESEARCH_REPORT_ZH.md`。JSON/NPZ 为原始结果。`predictions/` 与 `sealed_truth/` 的说明明确本轮没有盲成像任务；`data_public/dev_anchors/` 的 ε 是公开开发锚点。

## 重新执行数值研究

先把交付包解压到一个**新目录**，从其中 `research_max_struct_r1/` 运行。脚本会写该副本的结果；不要在封存原目录重跑覆盖本轮回执。使用自己的 Python 环境，或本机既有解释器：

```sh
PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
"$PYTHON" experiments/qualify_backend.py
"$PYTHON" src/tensor_families.py
"$PYTHON" src/dynamics_audit.py
"$PYTHON" experiments/plot_backend.py
"$PYTHON" experiments/plot_families.py
```

第一项仅两个球网格，对解析 Mie；第二项四个锚点/分支组合加一次固定粗方向插值和一个噪声射线；第三项16材料真实 Hankel–Green 局部动力学。没有大扫或额外数字。`tensor_families.py` 显式单 BLAS 线程；其他程序按各自实现记录。复跑 wall/CPU 应变化，不能期待逐位复现时间。

`experiments/review_f_metadata.py` 是对首轮回执的非物理审查与已存场秩分类；可在新副本运行补充 `RANK_CLASSES.json`，它不新增正演。`audit/consolidate_receipts.py` 仅汇总已有回执，但其中 D2 历史提取若涉及历史来源应按 `D_INPUT_PROVENANCE`/审查说明只读，不能凭缺失源码重造旧历史。完整本轮审计回执已经随包给出。

## 计数和随机性

F 使用 seed20261102，仅噪声方向用随机数；初始球/材料、两档网格、源库、.12 cell 扰动、.02弱源等在首次 F runner 写成后一次性执行，没有结果驱动调参。`configs/initial_plan.json` 是执行前总体计划；它不是覆盖全部实现细节的独立预注册系统。

计数：F 自身11 LU/56 forward RHS＋208额外三角 RHS；球资格2 LU/4 RHS；动力学 D1/D3 及其共享 setup 合计2 LU/12 RHS。总计15 LU/72 forward RHS。D2 历史记录审计本身没有新增正演。

回执审查后修正了 DOF、同场激发量与 constitutive residual 的字段区分、支持布尔值和计时边界；没有重算数值数组。当前源码已包含这些输出修正。历史细化插值程序对粗 q 作线性延拓再乘解析 bump，而非直接插值 ε；必须保留这一语义。

## 未运行项不是复现失败

M 没有可用的兼容 H(curl) 表示 runner；T 未冻结合法 k/δ；ROM 没有时序数据。因此不存在可声称已经可运行的 M/T 成像命令。不存在隐藏的后续 NN、GPU 或大规模图任务。

打包未包含在线下载的第三方全文 PDF/txt。已核查文献的 DOI、定位、claim ledger 和错误记录保留；如需重读全文，使用报告列出的作者/出版社链接。原始本地审查目录仍保留下载证据。
