# 复现入口

从现有B1仓库根运行，Python需NumPy/SciPy/PyTorch/Matplotlib/Pillow/threadpoolctl。无需安装CUDA、驱动、服务或新调度系统。当前主机使用既有`Gaussian/.venv_nn`；训练MPS、最终物理complex128/float64。

以下路径均相对仓库根，新实验目录简称`research/research_gn_nn_preservation_r4`。先从本轮Release恢复公开的数据/锚点/权重/预测归档，保留research目录层级。图像生成使用系统字体；没有同名字体的平台请使用冻结数据，不要替换字体后声称复现同一数据。

## 数据与严格GN

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python research/research_gn_nn_preservation_r4/experiments/prepare_data.py
python research/research_gn_nn_preservation_r4/experiments/run_anchors.py --split train --split val --split test --workers 2
python research/research_gn_nn_preservation_r4/experiments/release_anchors.py
```

生成器在新目录运行；不要覆盖已冻结的已发布数据。严格主锚点是`anchors/<obs_id>.npz`，discrepancy/weakreg是另外的后缀。远端SSH工具依赖用户已有A9 transport，仅用于本次资源接入；本地复现无需该私人传输环境。

## 从头训练

完整TRAIN/VAL严格锚点准备后，使用实际准备回执发布的数据就绪状态。不要用部分数据或旧R3权重替代：

```bash
python research/research_gn_nn_preservation_r4/experiments/train_models.py --data-ready --arms v_physics j_full j_null g_weak g_dark --seeds 7 --device mps --predict-val
```

alpha、额外种子和场约束的实际冻结选择见最终方法清单。`--loss const`或`--loss internal`分别运行单项消融，均为新模型。已存在的训练模型不应被覆盖。

## 公开输入推理、验收与离线评分

`experiments/predict_models.py`接受明确checkpoint和只含公共锚点路径的JSON列表；不读取形状标签或GT。预测输出仍须完整物理验收。

单个新观测的端到端入口为`experiments/infer_observation.py --observation OBS.npz --operators OPS.npz --checkpoint MODEL/best.pt --output OUTPUT --device mps`。它直接生成严格GN、运行网络和原物理校正，输出材料与完整费用，不需要数据清单、类别或GT。CPU/CUDA可改device，最终物理保持complex128。

```bash
python research/research_gn_nn_preservation_r4/experiments/run_evaluation.py --runs j_full_seed7_cascade2 --splits val
python research/research_gn_nn_preservation_r4/experiments/score_predictions.py seal --split val --runs j_full_seed7_cascade2
python research/research_gn_nn_preservation_r4/experiments/score_predictions.py score --split val
```

TEST在最终`audit/METHOD_FREEZE.json`后，按冻结run列表执行同样流程并将split改为test。不能使用TEST评分重新选择checkpoint、alpha、损失或seed。完整物理能量使用同一封存与冻结run清单：

```bash
python research/research_gn_nn_preservation_r4/experiments/test_energy.py --seal research/research_gn_nn_preservation_r4/evaluation/TEST_SEALED.json --runs-json research/research_gn_nn_preservation_r4/audit/FINAL_RUNS.json
python research/research_gn_nn_preservation_r4/experiments/plot_results.py
```

## 成熟结构CIE参考

默认CLI仅列公开计划，数值运行需明确`--execute`。CIE从公共均匀材料.1冷启动，不能免费使用NN的GN锚点。主参考为固定13母体×两种噪声；扩大覆盖时必须保留所有预先定义对象。

```bash
python research/research_gn_nn_preservation_r4/experiments/run_cie_reference.py --split test --scope figure13 --methods cie_dual --workers 2 --execute --test-methods-frozen
```

所有figure身份标记由清单生成；TRAIN/VAL/TEST图不可混用。raw、物理校正top1、离线best@K分别命名。最终结果与实际运行覆盖以结果表及回执为准，代码可运行不等于相应实验已经运行。
