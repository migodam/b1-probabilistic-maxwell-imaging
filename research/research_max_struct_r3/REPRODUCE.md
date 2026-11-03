# R3复现与依赖

从完整 `migodam/b1-probabilistic-maxwell-imaging` 仓库运行，保留同仓库R1/R2、FBG-R3、R4-DIVERSITY及Z-ALG-E1目录。ZIP是本轮源码/证据包，**不是不依赖仓库的独立solver**。实际输入路径和算法入口都在本仓库中；没有旧session私有算子或生成数据依赖。

执行环境为Apple Silicon CPU、complex128、BLAS 1线程；复用已有venv，未安装新环境/GPU。依赖numpy/scipy/threadpoolctl/matplotlib及其版本见audit/ENVIRONMENT.json。原工作区使用 `Gaussian/.venv_nn/bin/python`；外部复现将下面的 `python` 指向包含这些包的解释器。

## 只读结果

优先读results/DEVELOPMENT_SCORES.md、audit/FINAL_COST_TOTALS.json及所有17个外层数值回执。独立评分使用保存的m/Y，不再次跑forward。

## 重算的安全工作副本

生产runner会拒绝覆盖已封存文件。另复制整个仓库到专用重算目录，先在**该副本内**将R3的 `results/`、`predictions/` 和原 `audit/PASSIVE_GRADIENT_WIRING.json` 移到 `reproduction_original/` 保存，建空 `results/`、`predictions/`、`audit/` 输出目录。保留R1/R2全部依赖与R3源码；不要在发布证据目录中删结果。所有这些重算仍是公开开发例，不会重新成为blind。

在仓库根按依赖顺序运行（命令确有对应入口）：

```sh
python research/research_max_struct_r3/experiments/fixed_material_diagnostics.py
python research/research_max_struct_r3/experiments/development_controls.py
python research/research_max_struct_r3/experiments/passive_path.py passive_path
python research/research_max_struct_r3/experiments/passive_path.py extra_original
python research/research_max_struct_r3/experiments/known_controls.py
python research/research_max_struct_r3/experiments/finite_coordinate_paths.py
python research/research_max_struct_r3/experiments/all_pair_swaps.py
python research/research_max_struct_r3/experiments/paired_material_rays.py
python research/research_max_struct_r3/experiments/dressed_neighbor_rays.py
python research/research_max_struct_r3/experiments/noise_floor_distribution.py
python research/research_max_struct_r3/experiments/data_fiber.py
python research/research_max_struct_r3/experiments/physical_onecell_barrier.py
python research/research_max_struct_r3/experiments/vector_scope.py
python research/research_max_struct_r3/experiments/score_development.py
python research/research_max_struct_r3/experiments/plot_mechanisms.py
```

known_controls会先运行一次新passive adapter方向接线，再产出GlobalDual与scaled-original两个控制。有限路径依赖scaled-original保存的端点。各run JSON记录实际退出原因、候选与费用；不同SciPy/BLAS版本的ftol和近实root可产生少量数值差，不能逐bit要求同轨迹或改成功门适配结果。

独立重建图入口在 `audit/` 的审查脚本中；发布脚本 `experiments/export_report.py` 只生成网页/ZIP，不运行物理。源导航脚本仅给仓库总览增添本轮入口，不编辑R1/R2报告。

`experiments/generate_blind.py`及DRAFT数据/成本/评分配置是未执行后续模板；不要把它纳入上述已完成复现链，也不能把未冻结算法直接用于新确认。后续新盲测须先确定具体机制与对照，再一次生成，独立持有GT，所有预测封存后才评分。

代码与配置hash保持null（用户禁止hash核验）。git基版本只作来源记录；不重跑历史公式到机器精度，不重复所有历史tests。当前所有重要预测/根/候选表均已随本轮发布。
