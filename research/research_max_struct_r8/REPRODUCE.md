# R8 复现

所有命令从 **仓库根目录** 执行。R8直接复用R1–R7、FBG及Z-ALG的现有模块，不能只复制R8一个目录就运行。旧证据只读；复跑使用新输出目录或独立checkout。任务回执中的原绝对路径是执行来源记录，不要求复现者具有相同用户名或磁盘。

## 环境与数据

实际Mac环境：Python3.13.13、NumPy2.5.3、SciPy1.18.1、PyTorch2.14.0、Matplotlib3.11.1、threadpoolctl3.7.0；准确回执为 `audit/REPRODUCTION_ENVIRONMENT.json`。原训练CPU执行，不需要CUDA或MPS。i9采用现有Windows Python3.10.9、NumPy2.2.6、SciPy1.15.3、Torch2.6.0+cu124，**本轮物理队列未用GPU**。threadpoolctl在任务私有依赖目录提供；首次缺依赖的失败保留。依赖/BLAS差异可能导致非线性长期轨迹的数值差异，不承诺跨平台位级相同。

公共观测、冻结初值、模型权重、配置和紧凑结果随本轮发布。完整电流/场、预算轨迹和训练数据按发布索引分包，解压到仓库根目录保留相对路径。评分完成后发布的合成真值仅用于独立复核，不能反馈到反演。原Fresnel测量及第三方全文需从合法来源自行取得；不把缺失实测替换为合成。

每个物理进程限制BLAS为1线程。任务本身已设置并用threadpoolctl约束；批次并行度按本机资源选择。复现无需安装新调度器、驱动或服务。

## 最小完整单例

下面是真实入口，不是只计算代理目标。第一个例子运行冻结NN材料初值接global A+、完整FFT-TSOM辅助阶段及原Maxwell校正，第二个接sourcewise A+。输出目录必须尚不存在。

```bash
python research/research_max_struct_r8/experiments/run_case.py \
  --manifest research/research_max_struct_r8/configs/warm_all_manifest.json \
  --task-id r8_test_0000__a_plus__nn \
  --output-dir research/research_max_struct_r8/reproductions/case0000_nn_global

python research/research_max_struct_r8/experiments/run_case.py \
  --manifest research/research_max_struct_r8/configs/warm_all_manifest.json \
  --task-id r8_test_0000__a_plus_sourcewise__nn \
  --output-dir research/research_max_struct_r8/reproductions/case0000_nn_sourcewise

python research/research_max_struct_r8/experiments/run_case.py \
  --manifest research/research_max_struct_r8/configs/warm_anchor_confirmation_manifest.json \
  --task-id r8_test_0000__a_plus__nn_affine_anchor \
  --output-dir research/research_max_struct_r8/reproductions/case0000_nn_anchor
```

将task-id换成同对象的 `cie_dual__nn`、`physical_lbfgs__nn`、`a_plus__uniform` 等即可直接使用相应冻结比较。`PUBLIC_TASK.json`、`RECEIPT.json`、`prediction.npz`、预算快照和auxiliary/current状态自动保存。默认预测选择为已支付原Phi最小值，代码不打开GT。

## 原始训练与测试生成

只在没有相应生成回执/输出的**独立重建工作副本**中执行以下命令。现有发布目录保存原冻结回执，不能覆盖。种子和对象数写在生成器中，训练按独立验证误差选择epoch。

```bash
python research/research_max_struct_r8/experiments/generate_warm_data.py --split train
python research/research_max_struct_r8/experiments/generate_warm_data.py --split validation
python research/research_max_struct_r8/experiments/train_warm_start.py --mode physics --epochs 120 --device cpu
python research/research_max_struct_r8/experiments/train_warm_start.py --mode born --epochs 120 --device cpu
python research/research_max_struct_r8/experiments/generate_warm_data.py --split test
python research/research_max_struct_r8/experiments/prepare_warm_initializers.py \
  --cases research/research_max_struct_r8/data_public/warm_test/test_manifest.json \
  --weights research/research_max_struct_r8/RESULTS/training/physics/weights.pt \
  --name primary_physics
```

Born-only控制把weights换成`RESULTS/training/born/weights.pt`、name换成`primary_born_network`。测试生成器写独立`sealed_truth/warm_test`；反演器只读清单中的public_files和冻结initialization_file。历史执行曾先将评分分母从公共索引移到私有评分元数据，审计保留；当前生成器已直接采用分离合同。

论文采集域修复是另一份网络：依次执行 `generate_domain_training.py --split train`、`--split validation`、`train_domain_warm_start.py`。它复用原网络并加入256/64个32²平面波材料，未用Austria标签，不能代替主确认原权重或省略前序训练费。

## 批次和评分顺序

在清洁输出副本中，`run_queue.py --manifest ... --name UNIQUE --workers N`运行有限清单。原750清单是`warm_all_manifest.json`；新偏置126是`warm_anchor_confirmation_manifest.json`。部署时把`HARDWARE_ASSIGNMENT.json`与`HARDWARE_REALLOCATION.json`合读：后者将7个全部尚未执行对象整体移到空闲Mac，每对象所有方法仍在同机。实际主比较Mac405/i9345、偏置Mac66/i960；远程dispatcher的迁移确认不是物理运行。本轮已支付的Mac/i9部署特征计时另存，不把混合CPU秒当同一速度基准。

```bash
python research/research_max_struct_r8/experiments/run_queue.py \
  --manifest research/research_max_struct_r8/configs/warm_all_manifest.json \
  --name replay_original750 --workers 4
python research/research_max_struct_r8/experiments/run_queue.py \
  --manifest research/research_max_struct_r8/configs/warm_anchor_confirmation_manifest.json \
  --name replay_anchor126 --workers 4
```

这些命令的清单输出路径指向正式RESULTS，因此仅在清洁复跑副本运行；不应在已含原预测的发布副本直接执行。单例复现优先使用上面的`--output-dir`。

原评分入口 `score_warm_starts.py seal`只封存公开预测及不含GT的portfolio选择；随后另封存126偏置预测，再执行原750评分与偏置评分。两个集合都终结前不读取测试标签。封存使用文件路径、大小和时间戳，不做SHA检查。分析脚本的具体CLI可用`--help`查看；确认分析要求对应公开seal和完整评分回执。已有结果的再次分析也不得改变冻结预测；重新生成同名seal/score应在新的审计副本进行。

封存时间戳只检查**原执行期间**文件有没有再改，不是跨机器内容认证。Git checkout和ZIP解压不会保留原纳秒mtime，因此不能要求下载副本直接通过原运行的mtime检查。只做发布后评分复核时，在新的副本中把原两份seal、原portfolio选择、warm_scoring及相应anchor分析产物移入独立`released_original_audit/`保留，再对不变的已发布预测重新执行seal→score顺序。原审计文件仍可查看；这次是**已公开数据的审计回放，不是新的盲确认**。不得在本轮原证据工作树重新封存以抹除修改，也不应把恢复时间戳当作证明科学结果的办法。

## 原论文、历史、实测与3D

- **历史四对象**：`history_warm_prepare.py`、`history_warm_run.py`及补充清单见`experiments/history_warm_manifest.json`；执行21项后才由`history_warm_score.py`评分，`history_warm_interpret.py`生成解释。公开旧GT使用不会让它们重新成为盲测。
- **Austria**：`paper_benchmark.py build/solve/score`是32²探索入口；`paper_original64.py`运行原150/64的eps3.5，`paper_compare64.py --arm METHOD`是eps4五臂。具体arm及长corrector/NN域修复命令、失败恢复范围见`PAPER_BENCHMARKS.md`和`PAPER_FOLLOWUPS.md`。恢复运行会计入已复用auxiliary成本，不能按“全新一次”重复收费。
- **Fresnel**：先运行`fresnel_validate.py`、`fresnel_reference_poses.py`、`fresnel_foam_validate.py`，然后`fresnel_invert.py --case dielTM_dec8f --arm cie_dual`及相应twin/其他arm；2005入口`fresnel_foam_invert.py --side 60 --arm cie_dual`。NN回归添加`--initialization neural`。独立NIE为`fresnel_foam_nie.py`，多频为`fresnel_multifrequency.py`。这些入口实际读测量；获取和相位合同见`FRESNEL_REPORT.md`，不能只拿名义几何生成Y后称实测。
- **真3D**：`vector_benchmark.py`生成小型dyadic任务和exact-current探索；`vector_pr_benchmark.py`/`vector_pr240_benchmark.py`保留相应历史坐标比较；正式原Fourier系数PR为`vector_pr240_fourier.py`。`vector_mie_geometry.py`、`vector_refinement.py`、`vector_hollow_cube.py`分别提供Mie几何修复、跨网格和真正空心立方体。`vector_original_size.py`是明确截断的原电尺寸尝试。`vector_report.py`汇总，不把探索型exact-current称原CIE。

## 审计与发布

`build_cost_ledger.py`只读外层回执，不能启动solver或评分；同一次任务的setup/auxiliary/corrector嵌套时间不重复相加。`delegated/cost_package/package_release.py`生成分包计划；完整原数组保持不变，另建紧凑派生文件。失败、被打断的尝试、预算截断、训练及独立物理验证都保留；无法恢复的计时明确为null或下界。

精确原命令、过程日志、配置快照、冻结时间与硬件分配在audit和各delegated目录。历史日志中的controller PID只用于记录自己的当次队列，不是可在其他机器照抄的进程操作指令。

远程最终回传只把新的预测放回原RESULTS路径；远端打包中携带的旧controller审计写入`audit/remote_imports/队列/时间/`，不覆盖本地已经封存的配置或seal。较早、首次封存前的回传使用原路径，回执记录其时序。迁移确认只代表调度，不产生额外物理求解。

