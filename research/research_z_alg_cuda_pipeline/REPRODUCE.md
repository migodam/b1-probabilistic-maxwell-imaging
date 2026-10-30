# 复现与其他 session 接入

本目录是性能工程研究，不是新的盲成像实验。只使用 R4 已公开的匿名观测及 G/E0/S；不读取 `private_truth`、`truth_after_freeze`、数字标签或生成 seed。R2/R3/R4 冻结文件没有被修改。

## 已实际使用的机器

- XINAN：i9-14900HX，Windows/WDDM；本轮新进程仅使用自身 P 核 logical mask `21845`（0/2/4/6/8/10/12/14），没有修改其他进程。
- RTX 4060 Laptop，24 SM，约 8 GiB；驱动 592.82；已有 Python `D:\python\python.exe`，PyTorch 2.6.0+cu124、CUDA build 12.4、NumPy/SciPy。
- `threadpoolctl` 使用原 R4 项目已有的 `vendor/`，未重新安装。无新 CUDA/驱动/daemon/调度框架。
- 所有实验保持 complex128/float64。没有 TF32/FP16 替换。

本机实际路径：

```text
新代码：D:\AI\ZALG_CUDA_PIPELINE_20261030
只读公开数据：D:\AI\ZALG_R4_DIVERSITY_20261030\project\public_data
只读原实现：D:\AI\ZALG_R4_DIVERSITY_20261030\project\src
```

## 固定工作量吞吐实验

从新代码目录执行：

```powershell
D:\python\python.exe benchmark.py --legacy-root D:\AI\ZALG_R4_DIVERSITY_20261030\project --evaluations 64 --output results/throughput.json
D:\python\python.exe tune_backend.py --legacy-root D:\AI\ZALG_R4_DIVERSITY_20261030\project
D:\python\python.exe select_svd.py --legacy-root D:\AI\ZALG_R4_DIVERSITY_20261030\project
```

`throughput_v1.json` 是实际第一轮完整结果的保存文件名；脚本默认写 `throughput.json`。不要覆盖已有结果。`benchmark.py` 故意在最后进行一次 Graph 捕获尝试，原默认库的 `CAPTURE_FAILED` 是结果的一部分；库选择实验单独记录强制 cuSOLVER 后的成功捕获。两者不是重复精度校验。

Graph 不是当前默认生产路径。`select_svd.py` 在一个公开状态的实际 regularized Jacobian 上比较现成 SciPy gesdd 与 cuSOLVER gesvd/gesvdj，计入 GPU 输入/输出。它不比较 GT。

## 12 个完整优化任务

```powershell
D:\python\python.exe -m src.executor run manifest_1.json
D:\python\python.exe -m src.executor run manifest_2.json
D:\python\python.exe -m src.executor run manifest_4.json
D:\python\python.exe -m src.executor run manifest_4_hybrid.json
```

前三份 manifest 只有 lane 数和输出目录不同。每份为三个匿名旧公开观测 × 四个方法，8→12→24 的共同材料插值，最终使用完整 576 阶非线性物理；每阶段 maxiter=12、maxfun=24。时间上限仅作有限运行边界，本次主要停止来自既定迭代/函数调用数，不通过给某个方法更多时间制造性能收益。

第四份显式增加 `svd_backend: "hybrid"`。它在本进程内把 SciPy TRF 的 SVD 函数引用接到成熟 cuSOLVER/SciPy 实现：大矩阵且 GPU SVD slot 空闲时使用 CUDA，否则立即执行 CPU SVD。最多一个 GPU SVD，避免其他任务为它排队。它保持 TRF 信赖域规则和 FP64 economy SVD；采用了不同线性代数库后，不宣称任意非线性路径逐位相同。结果必须计入实际函数/RHS 次数。

### 不同项目如何接入

1. 提供与本目录一致的 scalar-TM `G`、`E0`（576×12）、`S`（48×576）和匿名 `Y`（48×12），保存为两份 NPZ。该入口没有被声称为任意 3D/decoder 的现成后端。
2. 复制一份 manifest，修改 `operators_path`、各 task 的 `public_path`、`output_dir`；全部路径可为绝对路径。`vendor_path` 仅在本机依赖放在原 vendor 时需要。
3. `method` 为 `native_trf`、`reciprocal_trf`、`data_lbfgs`、`globaldual`。每个任务独立持有 stream/LU，所有任务共用一次 G/E/S 上传。
4. 默认使用现成的库选择，不要一律设置 `TORCH_LINALG_PREFER_CUSOLVER=1`；本次 batch=16 在默认库下更快。是否采用 hybrid SVD 以端到端结果为准。
5. 将多目标、多初值等独立任务一次放进有限任务列表。小观测集合预载入；worker 完成后立即取下一项，主线程独立写回执。所有任务结束后进程退出，没有常驻服务。
6. 多照明在同一个材料下应当作为多 RHS，共用 LU。不要为 12 个照明分别开 12 个 GPU 进程；更不要把 stream 个数解释成固定 SM 分区。

若从 B1 仓库根导入：

```python
from research.research_z_alg_cuda_pipeline.src.device import SharedOperators, ObjectiveLane
shared = SharedOperators(operator_arrays)
lane = ObjectiveLane(shared, observation_arrays, stage=24, objective="data")
loss, gradient = lane.evaluate(initial_material)
```

`shared.score_batch(materials_cuda, targets_cuda)` 是 GPU 常驻候选评分入口，返回 full nonlinear relative residual 和 LU info。批大小由调用者决定；输入已经在 GPU，吞吐数字不含候选生成或完整反演费用。

## 汇总、绘图

下载任务 JSON、runner summary 和三个吞吐/算子选择文件到本目录 `results/` 后执行：

```sh
python research/research_z_alg_cuda_pipeline/analyze.py
python research/research_z_alg_cuda_pipeline/build_site.py
```

脚本生成 `RESULTS_SUMMARY.json`、`COST_LEDGER.csv` 和三张真实计时图。汇总使用实际回执，不重新运行物理、不读取 GT。并行 task 的 wall 不相加代替总 wall；共享 setup 每个 run 记一次。Python/Torch 导入与 SSH 网络传输不在 executor 计时窗内，不能把其 wall 声称为完整远程冷启动时间。

## 已记录失败和范围

- 首次 benchmark 缺少原 R4 的 vendor 路径，因 `threadpoolctl` import 中断；后续显式复用已有 vendor，未安装依赖。
- 默认库 CUDA Graph 捕获失败；cuSOLVER 后捕获成功，但没有成为本轮最快的候选评分路径。
- 现有 Torch profiler 只支持 CPU activity，nsys 不在 PATH；SM occupancy、硬件 kernel 重叠比例均为 **NOT_MEASURED**。图中的 task 时间线不是 GPU kernel/SM 时间线。
- 高频数值条件诊断没有搬进热循环；返回 LU/Gram 状态并在输出时记录完整物理残差。病态 source Gram 的旧 QR 回退不在本 fast path，不能把它描述为原冻结 solver 的逐位等价替换。
- 速度结果不意味着更高 GT recovery、更多解分支、全局收敛或概率 posterior。
