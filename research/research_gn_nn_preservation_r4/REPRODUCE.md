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

实际soft-alpha选择为1，因此复用J-full；受约束C*为G-weak。Voxel/J-full/G-weak均补17/27种子。`--loss const`或`--loss internal`分别运行J-full/G-weak的单项消融，均为新模型。完整20个NN及6个传统／继续方法的路径与种子见 `audit/FINAL_RUNS.json`，不能事后按TEST删方法。已存在训练权重不覆盖。

追加RNN的实际训练入口为：

```bash
python research/research_gn_nn_preservation_r4/experiments/train_recurrent.py --authorized --rho 1 --modes refreshed_q feedback1 feedback_hinge --device mps --predict-val
```

可选完整VIE微调已通过预声明VAL门，实际运行两个各320次更新的匹配臂。它们从G-weak-const的同一checkpoint开始，本阶段分别为纯材料监督与材料监督加完整VIE hinge；没有在控制臂额外保留const项。实际合同与更正记录位于 `audit/OPTIONAL_FT_*.json`。

## 公开输入推理、验收与离线评分

`experiments/predict_models.py`接受明确checkpoint和只含公共锚点路径的JSON列表；不读取形状标签或GT。预测输出仍须完整物理验收。

单个新观测的端到端入口为`experiments/infer_observation.py --observation OBS.npz --operators OPS.npz --checkpoint MODEL/best.pt --output OUTPUT --device mps`。它直接生成严格GN、运行网络和原物理校正，输出材料与完整费用，不需要数据清单、类别或GT。CPU/CUDA可改device，最终物理保持complex128。

```bash
python research/research_gn_nn_preservation_r4/experiments/run_evaluation.py --runs j_full_seed7_cascade2 --splits val
python research/research_gn_nn_preservation_r4/experiments/score_predictions.py seal --split val --runs j_full_seed7_cascade2
python research/research_gn_nn_preservation_r4/experiments/score_predictions.py score --split val
```

TEST方法与checkpoint在08:03:14UTC冻结。全量公开推理使用 `experiments/predict_frozen_test.py`，只读online清单与公共锚点，一个固定算子上传复用20个模型；它不会读取GT。原执行队列按固定M4/i9分配，准确命令与覆盖见 [双机执行说明](delegated/data/TEST_EXECUTION.md)。不能使用TEST评分重新选择checkpoint、alpha、损失或seed。完整物理能量使用同一封存与冻结run清单：

```bash
python research/research_gn_nn_preservation_r4/experiments/test_energy.py --seal research/research_gn_nn_preservation_r4/evaluation/TEST_SEALED.json --runs-json research/research_gn_nn_preservation_r4/audit/FINAL_RUNS.json
python research/research_gn_nn_preservation_r4/experiments/plot_results.py --tag main --runs discrepancy_allscales45_anchor cie_dual_cold55 v_physics_seed7_cascade2 j_full_seed7_cascade2 g_weak_seed7_cascade2_const recurrent_feedback1_seed7_cascade2
python research/research_gn_nn_preservation_r4/experiments/plot_test_energy.py
python research/research_gn_nn_preservation_r4/experiments/plot_generalization.py
python research/research_gn_nn_preservation_r4/experiments/plot_completion_costs.py
```

## 成熟结构CIE参考

默认CLI仅列公开计划，数值运行需明确`--execute`。CIE从公共均匀材料.1冷启动，不能免费使用NN的GN锚点。最终范围为全部108母体／216观测，按两个主机分配；13母体只是主图展示清单，不是总实验分母。

```bash
python research/research_gn_nn_preservation_r4/experiments/run_cie_reference.py --split test --scope all108 --methods cie_dual --workers 2 --execution-partition local --execute --test-methods-frozen
```

该命令只执行固定local子集72份观测；其余144份由原i9队列执行。它不是单机全216的伪命令。独立平台上运行单个未知观测可直接用前述 `infer_observation.py`，不需要原双机主机配置、类别、清单或GT；批量NN前向也不受原主机分配限制。重新执行的硬件／耗时属于新复算，不冒充原实验计时。

## 归档结果的直接重评分

Release保留权重、公开输入、全部候选、真实物理输出和离线数据。恢复入口见 [归档说明](delegated/data/RELEASE_PACKAGING.md)。原seal记录原机绝对路径，不应编辑成新机器路径；在恢复目录另建本地seal，再用同一评分器得到复算表。无需为重画图重复运行GN或训练。

主科学图为PNG；全部对象／全部方法图册使用**无损WebP**以控制静态网站体积。两者均由真实数值数组绘制，固定材料色阶，不使用生成图像。

所有figure身份标记由清单生成；TRAIN/VAL/TEST图不可混用。raw、物理校正top1、离线best@K分别命名。最终结果与实际运行覆盖以结果表及回执为准，代码可运行不等于相应实验已经运行。
