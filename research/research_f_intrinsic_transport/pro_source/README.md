# Research F 研究包

**主题：非线性 Maxwell 解集、有限反馈、先验几何与可学习材料探索。**
日期：2026-10-16。

## 阅读顺序
先读 `report/RESEARCH_REPORT_ZH.md` 的执行摘要及第6节。
完整报告覆盖用户全部十个核心问题、高阶推导、条件命题、反例、
N3源码审计、五篇近邻、算法选择、E0-E5与论文边界。

接下来读 `CODEX_NEXT_PROMPT.md`；它是可以直接交给Codex的下一轮任务。
默认只读缓存代数，新增Maxwell预算为0。
新独立anchor/finite-exploration预算只是后续建议，不在本次自动授权。

## 文件
- `report/RESEARCH_REPORT_ZH.md`：完整中文报告，含公式、证明步骤、来源。
- `CODEX_NEXT_PROMPT.md`：可直接交接的默认audit-only Prompt。
- `EXPERIMENT_PROTOCOL.json`：机器可读实验边界与后续门槛。
- `sources/SOURCE_AUDIT.md`：实际访问范围、缺失、文献读取深度与证据分级。
- `checks/algebra_checks.py`：小矩阵/两单元代数检查。
- `checks/algebra_checks_results.json`：本次实际运行结果。
- `MANIFEST.json`：打包文件的大小与SHA256，便于传递时核对完整性。

## 核心裁决
当前证据更支持在合法prior内研究有限噪声posterior，而不是假设N3附近已有精确fiber。
原decoder二阶项需要拆成坐标加速度和真实prior第二基本形式。
NN应在有效物理传输已经存在后学习有限补偿；N3保留为单图决策基线。
A17的候选增量是有限双侧残差/传输计算，不是自动的多解或采样保证。

## 本次计算范围
没有新增Maxwell网格仿真、Mie、GPU训练或历史重跑。
`checks` 的小矩阵计算不代表复杂电磁成像、泛化或posterior benchmark。

运行轻量检查：
```bash
python -m pip install numpy scipy
python checks/algebra_checks.py
```
脚本会覆盖同目录的代数检查结果JSON；不会访问网络、仓库或GPU。
报告为UTF-8 Markdown，建议用支持LaTeX数学的Markdown阅读器打开。
本包不包含第三方论文全文、用户原始二进制证据或字体文件。
