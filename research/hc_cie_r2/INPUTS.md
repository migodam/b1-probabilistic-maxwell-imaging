# 复现输入附件

`HC_CIE_R2_Replay_Inputs.zip` 保存本轮程序引用的历史冻结观测、算子、端点、离线评分键和圆盘小模型数值。ZIP 为 **90,741,531 字节**，含 770 个输入条目及一份逐文件索引。每个 NPZ 保留所需数组的原数值、形状与 dtype。附件由文件及数组键白名单生成，发布的 HC-CIE-R2 结果另在当前研究目录。

## 解压到已有项目

从已包含 HC-CIE-R2 源码的主仓克隆根运行。历史代码由原仓库提供；附件路径兼容主仓的平铺目录结构。`unzip -n` 保留已有文件，建议使用新克隆或新工作树，避免已有的不同实验输入抢先被读取。

~~~bash
unzip -n /path/to/HC_CIE_R2_Replay_Inputs.zip
export HC_CIE_R2_INPUT_ROOT="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
~~~

复算环境使用 Python、NumPy、SciPy、Matplotlib、threadpoolctl、SymPy、mpmath。三维需要支持复双精度的 PyTorch。各实验可按对应报告选择运行入口。以下命令可独立复跑强凸投影、圆盘坏点与两种严格材料证书，输出写往新的工作目录。

~~~bash
TOY_SOURCE=research/hc_cie_r2/replay_inputs/toy
TOY_OUT=hc_cie_r2_replay/toy
python research/hc_cie_r2/src/toy_replay.py --source "$TOY_SOURCE" --out "$TOY_OUT"
python research/hc_cie_r2/src/toy_certificates.py --source "$TOY_SOURCE" --out "$TOY_OUT"
python research/hc_cie_r2/src/toy_material_interval.py --source "$TOY_SOURCE" --out "$TOY_OUT"
~~~

原 ZIP 和两份报告 PDF 无需参与数值复算。若希望复做 PDF 审核，可用自己持有的两份 PDF 运行对应脚本。 [小模型完整复现说明](TOY_AND_CERTIFICATES.md)

## 输入范围

| 原始系列 | 附件内容和使用位置 |
|---|---|
| R13 | `data_public/fresh48` 全 48 份带噪观测。CIE、Dual-CIE、A+-Arc、G19、DR 的 240 个 `primary.npz` 各保留 `m` 和 `Y_prediction`。两份 `corrected.npz` 仅保留实权检查用 `m`，另留原端点计时字段。 |
| R6 / FBG / R2 | 互相可分的两条已评分弯条、四个历史端点、复算所需的公开算子和观测字段、成本记录。公开精简任务索引仅含两条观测。完整主仓也有其原历史配置，解压可直接复用。 |
| R14 / R15、R16 | 位于历史 R14 工作树的三份 `r14_new_*` 观测按实际 R15 目录归档。R16 所调用的 R10 24 份 `twod` 和 R11 两个视图各 8 份均保留，原因是原 R16 加载器在取本轮子集前先验证完整子群。 |
| SRAP 二维 | 外积谱与模型差距线使用的 22 份合成观测和对应端点。另含两份原 SRAP 已处理的 8 GHz 训练掩码观测及原端点，供原物理边界重算。 |
| SRAP 三维 | 半球及面阵各 12 个原带噪端口和误差半径，所有 Dual `light_180` 材料。继承的 `vector/physics.py` 在主仓直接重建传播，无需另存三维 Green 矩阵。 |
| 外部复算 | 三份原 Austria 合成观测，四份观测对应的 NIE 与 Dual-CIE 历史端点和选定成本字段。原计分表和正演代码保留在主仓。 |
| 小模型 | 观测几何、固定复矩阵、驻点种子和复核用的精确数值放在 `research/hc_cie_r2/replay_inputs/toy/research_work`。 |

R13 离线标签放在 `sealed_truth/fresh48`。R6 评分只使用原来的 `sealed_truth` 或 `truth_after_freeze` 目录及公开精简映射 `offline_replay_mapping.json`。SRAP 的 `data/truth` 仅含模型差距线准备阶段使用的材料/家族键。端点材料是历史公开预测。脚本要求封存后才执行离线评分，监测改变量的真值和观测由各自冻结生成脚本重建。附件没有新实验预测、完整运行会话、提示、凭据、原 PDF 或历史私有映射文件。 [原 R13 归档](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r13-20261112) · [原 SRAP 归档](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/srap-cie-array-study-20261121)

## 外部 Fresnel 复算

外部复算的另一套 FoamTwinDielTM 8 GHz、13 个训练源分割，需读者自备原始 `FoamTwinDielTM.exp`。它和附件中的旧 SRAP 两份已处理训练掩码是不同观测版本。附件保留 13 源分割的历史预测端点，原始测量文件、该分割的训练观测与留出场按原研究政策由读者在本地准备。指定目录后，原代码重建所需的训练与留出数组，两份都写入原研究的 `.runs/` 下。评分也使用同一 `FRESNEL_2005_DIR`。 [Institute Fresnel 原测量与出处](https://doi.org/10.1088/0266-5611/21/6/S09) · [SRAP 取得与适配说明](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/srap_cie_array_study_20261121/REPRODUCE.md)

~~~bash
export FRESNEL_2005_DIR=/path/to/authorised/Fresnel2005
export MAXWELL_KILLTEST_WORKSPACE="$PWD"
python research/top_venue_killtest_20261120/run_killtest.py \
  prepare --workers 1 --case FoamTwinDielTM_8GHz_13tx
python research/hc_cie_r2/src/classical_run.py \
  --config research/hc_cie_r2/configs/classical_v1.json \
  --out hc_cie_r2_replay/classical --workers 1
~~~

## 白名单制作与平铺验证

制作脚本从已有工作区历史资产读取，由 `--toy-source` 指向原附件解压目录，输出平铺 ZIP。附件内 `research/hc_cie_r2/replay_inputs/manifest.json` 记录每个条目的原始相对来源、读取角色、数组 dtype 和形状。 [制作入口](src/collect_replay_inputs.py)

在隔离平铺目录里，从原仓提交 `5f806c9ac` 取所需物理代码，再解压最终附件。21 份谱研究观测全部按正常加载器打开。R13 A+-Arc 端点完整正演的相对端口差为 \(1.45\times10^{-15}\)，两个 R6 端点端口一致且评分材料可读。旧 SRAP 两套已处理 Fresnel 掩码完成四条原物理记录。三维半球 `hemi_00` 的审计得到 1.0440δ，原方程相对缺陷为 9.66×10⁻¹¹。外部 13 源分割从另行提供的原件重建，四份外部奥地利/Foam 历史端点的重新正演端口与保存值相同。数值小模型、精确有理证书及两合同材料区间证书均在独立目录通过。原始工作区保持只读。
