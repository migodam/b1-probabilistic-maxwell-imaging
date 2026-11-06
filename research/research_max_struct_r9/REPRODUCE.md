# R9复现与完整证据

## 已实际执行

**最终版本：**第一阶段113次、小三维29项、经典Windows16项、原尺寸三维九臂及独立三起点物理校正均已完成。代码、紧凑记录、报告与图在Git，完整大数组在下述同仓库Release。各分支实际命令、环境与复现输出隔离规则见后文。

47个公开开发逆问题、66个新六材料逆问题、2个公开曲率诊断。全部113个科学反演正常结束；初次解释器错误导致的28次工程失败另存。R1–R8未改写。第一阶段环境见[audit/ENVIRONMENT.json](audit/ENVIRONMENT.json)：Mac arm64、Python3.13.13、NumPy2.5.3、SciPy1.18.1、Matplotlib3.11.1、threadpoolctl3.7.0；每进程1线程，最多4进程，无GPU或训练。后续远端Windows为NumPy2.2.6、SciPy1.15.3、PyTorch2.6.0+cu124；真实P核/CUDA范围见REMOTE_COMPUTE.md，不能把第一阶段CPU说明套到整轮。

原始执行清单保留当时绝对路径；不要据此创建同名历史目录。便携脚本将公开输入映射到本目录随附NPZ，并在新目录生成清单。基线使用当前仓库内已冻结的R8源文件，核心物理仍复用FBG-R3/R4。

## 下载

代码、报告、公开观测、离线标签、紧凑评分和主要图在Git仓库。完整逐步current/material、预算轨迹、全部回执/日志/失败和绘图数组位于[同仓库R9 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)。分包内部路径从`research/research_max_struct_r9/`起，解压到仓库根目录即可。分包清单提供路径和字节数，不要求SHA检查。第三方论文全文缓存不再分发，文献URL与提取证据账本保留。

## 第一阶段公共反演复现

先进入完整仓库根目录，使用已有兼容Python环境；基础依赖为numpy、scipy、matplotlib、threadpoolctl。不需要CUDA、torch、NN权重或新服务。下面的`python`指上述环境的解释器，`--python`保留venv入口，**不要将其软链接解析为基础解释器**。

```sh
python research/research_max_struct_r9/experiments/prepare_reproduction.py \
  --cohort development \
  --destination research/research_max_struct_r9/reproduction/development

python research/research_max_struct_r9/experiments/run_queue.py \
  --r8-repo "$PWD" --python "$(command -v python)" --workers 4 \
  --receipt-dir research/research_max_struct_r9/reproduction/development_queue \
  --manifests \
  research/research_max_struct_r9/reproduction/development/baseline_development.json \
  research/research_max_struct_r9/reproduction/development/baseline_historical_fixed.json \
  research/research_max_struct_r9/reproduction/development/baseline_cost_controls.json \
  research/research_max_struct_r9/reproduction/development/landscape_development.json
```

这重放最初43项，不覆盖原结果。`--cohort confirmation`准备初始60项清单：`confirmation_baselines.json`、`confirmation_landscape.json`、`confirmation_single_noise.json`；同样用新receipt目录运行。所有输出目录须在R9内以保护历史结果。

六个弱阶段补充是从已付费profile第一阶段保存的m继续，原清单在`configs/confirmation_weak_handoff.json`，公共初始化也随包提供。若重放该recipe，应将`reused_initialization_provenance`指向本次profile的`profile_0_0.01_physical.npz`，导出其`m`为`m_init`，使用同一3000步原物理corrector。**不能把提供的暖初值当免费：**独立成本需加第一阶段原预算点及setup；总实验账本共享阶段只收费一次。原执行与成本计算均已保存，无须重生成材料来复算评分。

若只复现一项，用`experiments/run_baseline.py --manifest ... --task-id ... --r8-repo "$PWD"`或`experiments/run_landscape.py --manifest ... --task-id ...`。算法合同来自清单而非文件名推测。

## 封存后评分

便携评分器检查所列全部预测已终止并保存公开清单，之后才读随附`sealed_truth/development/`或`sealed_truth/confirmation/`。例如对上述43项：

```sh
python research/research_max_struct_r9/experiments/score_reproduction.py \
  --cohort development \
  --output research/research_max_struct_r9/reproduction/development_scores \
  --manifests \
  research/research_max_struct_r9/reproduction/development/baseline_development.json \
  research/research_max_struct_r9/reproduction/development/baseline_historical_fixed.json \
  research/research_max_struct_r9/reproduction/development/baseline_cost_controls.json \
  research/research_max_struct_r9/reproduction/development/landscape_development.json
```

原机完整证据使用`score_development.py seal/score`与`score_confirmation.py seal/score/figures`。后者要求完整66项及冻结的两组portfolio，先按原Phi决定portfolio候选，再读GT。历史依赖/绝对路径在原封存中保持不变，便携评分入口不依赖这些旧路径。

`plot_development.py`绘制开发图；`summarize_results.py`只读取已评分行，生成逐方法总表及事后30/60 CPU秒切片。不得用这些GT曲线反改本轮在线停止或选根。

## 理论与诊断入口

- `src/baseline_adapter.py`：公开噪声球秩选择＋原R8求解器，没有新forward。
- `src/landscape_solver.py`：完整current profile、固定/延拓/重启和数据目标同伦。
- `experiments/curvature_diagnostic.py`：两个冻结端点的有限模型自由坐标Hessian和有限probe。
- `experiments/generate_confirmation.py`：六例一次生成，独立truth目录；重生成需新的工作目录，不能覆盖原证据。
- `audit/gradient_wiring/RECEIPT.json`：唯一必要的新梯度接线检查；不是重复Atlas恒等式。
- `delegated/R9_ALL_EXECUTION_HANDOFF.json`：全部实际队列和失败索引。

物理值和最终材料会对BLAS/平台的有限精度轨迹敏感，不承诺跨硬件逐位一致。原R6 dense与R8 nativeFFT已明确区分，Windows旧CPU不当本机等成本对照。第一阶段是有限二维同离散模型。新增经典二维与全矢量三维使用下列独立入口；本轮没有新实测或全局凸化资格。

## 扩展的实际入口

- [经典二维的远端队列及完整分析命令](delegated/classic_extension/REPRODUCE.md)：冻结16项、i9四P核；历史队列与新复现输出必须分开，不能覆盖已有评分。
- [原3λ CM九臂与IGT诊断](delegated/original3d_feasibility/REPRODUCE_ORIGINAL_CM.md)：完整repo、已有Windows/PyTorch CUDA，显式新run-dir。已执行的物理入口与只做CLI核验的路径包装明确区分。
- [小三维29项和物理资格](VECTOR3D_EXTENSION.md)：不能将小电尺寸模型当成原3λ任务，也不将数值低残差当Mie准确性。
- `experiments/analyze_original_residual_transport.py`只读九臂已保存范数，计算已发生残差方向的传播增益区间；无GT数组或新forward。
- `experiments/compile_extension_costs.py --final`汇总完整旧账与独立增补。存在在途任务时仍报告IN_PROGRESS，不能靠传入--final把它变成完成。

原尺寸共同短物理校正另在`configs/original_physical_refinement.json`和`RESULTS/original_vector3d/physical_refinement/`记录，不覆盖九臂stage2或修改原论文停止。实际运行入口为`run_original_physical_refinement.py`（Windows P核12/14、已有CUDA环境）；历史入口固定输出且拒绝覆盖。需要重跑时使用如下新增薄包装器，复用同一配置和求解函数，只改输出目录：

```sh
python research/research_max_struct_r9/experiments/reproduce_physical_refinement.py run \
  --output-dir reproduction/physical_refinement
# 三臂完成并完整取回新输出目录后，单独离线评分：
python research/research_max_struct_r9/experiments/reproduce_physical_refinement.py score \
  --output-dir reproduction/physical_refinement
```

该包装器只做CLI/语法核验，没有为了验证包装器重新付费跑三臂；物理算法和原始三臂已实际执行。评分先检查三个完整封存，再读取已有truth；输出REPRODUCTION_SCORES.json，保护旧图、报告和分数。最初输入为已付九臂Global/A+端点，复现完整独立成本需加它们的前缀，不把暖状态当免费。

归档脚本`package_evidence.py`无损保存全部本轮科学文件，排除第三方论文全文缓存、Python缓存及生成包自身。`audit/PUBLICATION_MANIFEST.json`按分包记录路径/大小，`PACKAGING_RECEIPT.json`记录归档费用。分包是归档时快照，后续发布回执和最终归档成本以Git中的同名文档为准；无需重复打包科学数组以加入自指的发布元数据。
