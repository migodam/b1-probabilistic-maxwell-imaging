# 可以复用的 Maxwell 算子与成熟 CUDA 实现

本目录回答两个问题：现有工作是否已经有更好的算子，以及如何让其他 session 直接使用。下面区分“本轮已运行”“源码可移植”和“网络中已存在但本机未安装”。阅读源码不等于复现了它的实验。

## 当前默认路径：复用成熟库，压缩调用之间的空档

本轮没有手写 LU、Cholesky、BLAS 或 FFT 内核。`src/device.py` 使用已经安装的 PyTorch 2.6.0+cu124，调用其现成的 CUDA 线性代数与张量算子；库选择继续使用 PyTorch 的默认策略。新实现集中处理算子常驻、保留 LU、一次返回、流依赖和队列周转。

实际接口：

```python
from research.research_z_alg_cuda_pipeline.src.device import SharedOperators, ObjectiveLane

shared = SharedOperators(public_operators)  # 一次上传 G/E0/S
lane = ObjectiveLane(shared, public_observation, stage=24, objective="data")
value, gradient = lane.evaluate(material_coefficients)
endpoint = lane.close(material_coefficients)
```

`objective="dual"` 使用 R3/R4 的完整材料相关总场 Gram；`residual`、`jacobian` 可供 SciPy 的 Native/Reciprocal TRF 调用。CPU 优化器仍需每次目标/梯度返回，但只等待自己的 completion event，不强迫所有 lane 停顿。全部源在一次 LU 下同时求解。

已运行合同：24×24 独立实材料、12 源、48 复接收，完整 scalar-TM disk-pulse VIE，complex128。材料 box、正则尺度与原 R4 一致。接口没有被验证为任意 latent decoder、非方形网格、被动复材料或 3D vector-DDA 的通用后端。

## 项目内可复用资产

路径相对于 `/Volumes/migodam's-external-brain/Research/Inv_SLAM`，其他 worktree 应先定位同名资产。不要修改冻结源文件。

| 来源及精确入口 | 已优化的部分 | 可以怎样复用 | 本轮采用情况与限制 |
|---|---|---|---|
| `research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity/src/cuda_physics.py` | complex128 CUDA LU、共享源 RHS、接收转置 Full-J、源伴随、small Gram | 本轮物理公式、导数和旧执行对照的直接来源 | **已采用/已运行**；旧版每个 operation 做全设备同步且多次 `.item()`/`.cpu()`，本轮移除这些调用边界开销 |
| 同仓库 `research/research_z_alg_r4_diversity/remote_runner.py` | Windows P 核绑定、长驻 worker、任务与成本回执 | 复用远程连接和只约束自身进程的资源原则 | **已参考**；新 wrapper 用单 CUDA 进程和私有 stream，未改旧 4 CPU+2 CUDA 进程实验 |
| `Gaussian/A9/code/a9_engine.py::_Solver` | 3D DenseDDA 的 CUDA LU 和多 RHS/伴随接口 | 向量 Maxwell 模型应优先复用其模型约定和 RHS 计数，而非另写物理 | **只读审计**；原版 solve 会返回 NumPy，若移植驻留层应把连续张量运算留在设备；不能把 3N 场和 N 个实材料直接混同 |
| `Gaussian/A17/SHAPES_3D_BASELINES_20261030/research/delegated/endpoint/gpu_refined_endpoint_batch.py::GpuRefinedCoreEndpointBatch` | 压缩 endpoint 批量 QR、Hermitian SPD Cholesky、direct/Woodbury 分支、一次打包返回；原方法有至多三次同系统 FP64 修正 | 适用于原电流/材料压缩子问题；保留原 SPD 正则、q、RHS 与有效性语义 | **只读审计，未移植/未复跑**；这是可优先引用的现成批处理实现，不是新的 VIE 正演器 |
| `Gaussian/A17/SHAPES_3D_BASELINES_20261030/remote_fast_gpu_parallel.py` 及 `remote_gpu_parallel.py` | 一个 CUDA 进程、两线程作业、完成一个立即补一个、预算合并 | 采用同样的无批次屏障任务周转思想 | **机制已采用**；原文件没有显式给每个任务私有 stream，本轮补齐流所有权；不采用其全局提高优先级设置 |
| `Gaussian/A18_C/route-a/shared_backend/fft_maxwell.py`；Route B 同名文件 | 均匀 3D 向量网格的 padded FFT 卷积、伴随传播、CUDA FFT 路径 | 大网格可优先复用传播形式；把匹配的 GMRES/Krylov 向量也留在 GPU 才能避免迭代往返 | **只读审计**；当前 CUDA `L` 每次同步并回 CPU，SciPy GMRES 仍在 CPU；不是已经完整驻留的 GPU 线性求解。不得不经推导替代本轮 scalar-TM G |
| 同仓库 `research/research_z_alg_r2/src/physics.py::ReducedKernel` | 固定变化材料集合的精确有理核和 Schur 切片复用 | 对确实固定的补集消元，再决定小核是否值得上传 CUDA | **源码可用，未重跑**；4/6 阶求根不应为了“用 GPU”强搬；补集变化使缓存条件改变 |
| `Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/stable_roi_operator.py` | 固定宿主 ROI 的完整非线性 Schur 消元与 Jacobian | 有固定外部背景时省掉重复全域求解 | **源码可用，未移植**；当前全部 576 材料未知，不能把外部材料偷固定来获得加速 |
| `Gaussian/B1/b1-p3l-r10-bilateral-structure/research/p3l_r10_bilateral_structure/src/reduced_response.py` | 固定形状的 Schur/三角求解和 block-Krylov 响应复用 | 同形状大量幅度查询时优先重用；先计冷 setup，再计查询 | **源码可用，未重跑**；形状改变后旧分解不能无条件精确复用；该文件还包含原实验审批/预算/误差控制，不是随处运行的通用入口 |

## 网上已有的成熟算子：采用与暂不采用的判断

| 已有实现 | 真正解决的计算 | 本轮决定 |
|---|---|---|
| [NVIDIA cuBLAS / cuSOLVER](https://docs.nvidia.com/cuda/cusolver/index.html) | 复矩阵乘法、带 pivot 的 LU、QR、Cholesky，以及多 RHS 解 | 已通过 PyTorch 使用；不重写基础线性代数。复用同一 LU 的 primal/adjoint 比单独调用多个 `solve` 更合适 |
| [PyTorch `lu_factor_ex`](https://docs.pytorch.org/docs/stable/generated/torch.linalg.lu_factor_ex.html) | 返回设备 info，`check_errors=False` 避免立即为了检查状态而同步 | 已采用；info 随真正需要的目标和梯度一起返回。不会丢弃失败状态 |
| [PyTorch CUDA streams/events](https://docs.pytorch.org/docs/stable/notes/cuda.html) | 独立任务排队、精确依赖、pinned memory 的异步拷贝 | 已采用；主要工程收益来自去掉不必要的全设备屏障和细碎 host 往返 |
| [PyTorch cuSOLVER/MAGMA backend selection](https://docs.pytorch.org/docs/stable/backends.html) | 在成熟库之间选择，不改变线性代数问题 | 已实测 default 与 cusolver。强制 cuSOLVER 使本次 graph 捕获成功，但 batch=16 评分反而更慢；不设置全局强制默认 |
| [cuFFT / CuPy FFT](https://docs.cupy.dev/en/stable/user_guide/fft.html) | 固定均匀网格传播的卷积 FFT | 项目已存在对应向量实现；大规模优先评估结构匹配的 FFT，而非无限放大 dense LU。本轮没有新装 CuPy |
| [CuPy GMRES / LinearOperator](https://docs.cupy.dev/en/stable/reference/generated/cupyx.scipy.sparse.linalg.gmres.html) | GPU 数组接口上的迭代线性求解，可连接 FFT matvec 和预条件 | 可作为更大网格的备选；**NOT_RUN/未安装**，不把它列为当前吞吐结果，也不假设高对比 GMRES 必定快 |
| [NVIDIA CUDA Best Practices](https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html) | 数据驻留、合并传输、pinned buffers、独立 stream | 作为实现依据；先解决数据流和任务周转，再讨论手写 kernel/warp 优化 |
| [NVIDIA MPS](https://docs.nvidia.com/deploy/mps/when-to-use-mps.html) | Linux/QNX 多进程共享 CUDA 的资源仲裁 | 当前 Windows/WDDM 不适用；本轮不安装服务 |
| [CUDA 12.4 Green Contexts](https://docs.nvidia.com/cuda/archive/12.4.0/cuda-driver-api/group__CUDA__GREEN__CONTEXTS.html) | 驱动层的 SM 资源划分 | **NOT_RUN**；标准 PyTorch stream 不等于 SM 分区。为吞吐固定切碎 24 SM 未必有利，且分区本身也不保证并发进展 |

## 给其他 session 的选择规则

1. **先选物理上等价的算子。** 场/电流变量、材料共享、实化顺序、ROI 和背景条件必须匹配。已有传播 SVD 不意味着一般非均匀 `diag(m)` 与传播矩阵交换。
2. **多源共享 LU，多目标共享 G/E/S，多候选分批。** 独立多初值、不同照明子集或目标可以作为队列任务；同一材料的多个源应是同一个 LU 的 RHS，不能当成多个独立 factorization。
3. **采用本轮实际测过的入口。** `python -m src.executor run manifest.json --lanes 4`；输入为公开观测，不需要 GT。完整复现命令见 `REPRODUCE.md`。
4. **性能设置是按 workload 的。** 梯度型任务优先评估 4 lane；候选正演评分优先评估 batch=16。TRF 的 CPU SVD/优化仍有成本；不要声称只给它 CUDA forward 就已经全 GPU。
   对本轮混合 TRF/L-BFGS/GlobalDual 任务，推荐已测 `manifest_4_hybrid.json`：4 lane 加 `src/hybrid_svd.py`。大 SVD 在一个空闲 GPU slot 上调用成熟 cuSOLVER gesvd，否则立即调用 CPU gesdd；64/144 维 SVD 留在 CPU。12 个完整任务在相同 631 评估/635 LU/26,880 RHS 下由 99.5 s 降到 56.7 s。该适配器使用当前 SciPy private TRF 的 SVD 引用；其他版本需明确接口适配。
5. **复用现有测量，不给每个 session 加一套重复精度检查。** 新的模型或矩阵尺寸只补一个会改变实现选择的吞吐/端到端实验。不会替代实际物理验收。

本机根目录 `AGENTS.md` 和本 Git 仓库根 `AGENTS.md` 已加入这些入口。未向其他 session 发送消息，也未改变它们正在运行的任务。
