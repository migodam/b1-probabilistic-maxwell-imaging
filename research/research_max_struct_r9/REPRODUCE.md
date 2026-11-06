# R9复现与完整证据

## 已实际执行

**中间发布说明：**下面的113次为已完成第一阶段。三维29项另见VECTOR3D_EXTENSION.md，经典二维16项与原尺寸三维仍在执行。当前只发布已经完成的代码、报告、紧凑结果和图；下面的Release地址是最终归档目的地，尚未创建，不应视为已上传全部大数组。最终队列结束后统一补充远端命令及归档核验。

47个公开开发逆问题、66个新六材料逆问题、2个公开曲率诊断。全部113个科学反演正常结束；初次解释器错误导致的28次工程失败另存。R1–R8未改写。运行环境见[audit/ENVIRONMENT.json](audit/ENVIRONMENT.json)：Mac arm64、Python3.13.13、NumPy2.5.3、SciPy1.18.1、Matplotlib3.11.1、threadpoolctl3.7.0；每进程1线程，最多4进程，无GPU或训练。

原始执行清单保留当时绝对路径；不要据此创建同名历史目录。便携脚本将公开输入映射到本目录随附NPZ，并在新目录生成清单。基线使用当前仓库内已冻结的R8源文件，核心物理仍复用FBG-R3/R4。

## 下载

代码、报告、公开观测、离线标签、紧凑评分和主要图在Git仓库。完整逐步current/material、预算轨迹、全部回执/日志/失败和绘图数组位于[同仓库R9 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-max-struct-r9-20261106)。分包内部路径从`research/research_max_struct_r9/`起，解压到仓库根目录即可。分包清单提供路径和字节数，不要求SHA检查。第三方论文全文缓存不再分发，文献URL与提取证据账本保留。

## 公共反演复现

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

物理值和最终材料会对BLAS/平台的有限精度轨迹敏感，不承诺跨硬件逐位一致。原R6 dense与R8 nativeFFT已明确区分，Windows旧CPU不当本机等成本对照。所有新测试是有限二维同离散模型；本包不提供新3D、实测或全局凸化资格。
