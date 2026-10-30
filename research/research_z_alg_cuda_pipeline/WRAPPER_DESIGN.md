# 单进程 CUDA FIFO wrapper

`src/executor.py` 在一个进程内只创建一次 `SharedOperators`，并启动 1、2 或 4 个线程。每个线程立即从 FIFO 领取下一个任务；任务结束后把结果投递至内存队列，随即继续领任务。主线程独立把 JSON 原子写盘，因此慢写盘不阻塞 worker 领取后续任务。线程池在整份 manifest 执行期间长驻，没有批次屏障、重试、后台服务或额外框架。

## 运行

在 repository 根目录运行：

```sh
python -m research.research_z_alg_cuda_pipeline.src.executor run /absolute/path/manifest.json --lanes 2
```

也可在本研究目录运行：

```sh
python src/executor.py run /absolute/path/manifest.json --lanes 4
```

`--lanes` 可省略，此时使用 manifest 的 `lanes`，默认 1。每个相对文件路径均相对于 manifest 所在目录解析。JSON 示例（路径需换成实际 public 文件）：

```json
{
  "operators_path": "public/operators.npz",
  "output_dir": "results/pipeline_lanes2",
  "lanes": 2,
  "device": 0,
  "tasks": [
    {
      "task_id": "case001_data",
      "public_path": "public/case001.npz",
      "method": "data",
      "config": {
        "stages": [8, 12, 24],
        "wall_seconds": 45,
        "maxiter": 2000,
        "maxfun": 4000,
        "lower": 0.02,
        "upper": 20,
        "initial_material": 0.1,
        "smooth_weight": 0.0001,
        "tv_weight": 0.0001,
        "tv_epsilon": 0.01
      }
    }
  ]
}
```

任务可额外给 `result_path` 覆盖默认的 `output_dir/task_id.json`。`device_module` 可指定导入模块，默认是同目录的 `device`。共享输入只读取一次，小型 public NPZ 集合在 GPU setup 时预载入；task 领取不再等待读盘，不再次复制固定算子。

## Device 接口

```python
shared = SharedOperators(public_operators)
lane = ObjectiveLane(shared, public_case, stage=24, objective="data", config={...})
lane.stream                 # 独立 torch.cuda.Stream
lane.n_coeffs               # 当前 stage 的系数数目
value, gradient = lane.evaluate(coefficients)
residual = lane.residual(coefficients)
jacobian = lane.jacobian(coefficients)
physical_receipt = lane.close(coefficients)
```

`close` 返回 JSON 可序列化字典，包含最终材料、物理残差、成本等，并保证返回的 host 数据已可读。residual/Jacobian 均针对材料系数；残差含相同正则项，使 `.5 * residual @ residual` 等于原 data loss 加共同正则。物理算子由 shared 只读常驻，各 lane 的目标数据、LU、状态、workspace、缓存及 CUDA stream 独立持有。每个 stage 创建一个 lane；同一 stage 的所有优化回调复用它。任务结束即释放该 lane。

共享构造完成后主线程记录一个 CUDA Event，所有 lane stream 用 `wait_event` 依赖它。所有回调与 `close` 均在对应 stream context 内执行。设备层若使用独立 copy stream，应以 copy-complete event 建立 copy → compute 依赖，再由 lane event 同步所需 host 返回值。wrapper 不调用 `torch.cuda.synchronize`。host 回传所需同步、CUDA 运算与显存生命周期由 device 实现管理。ObjectiveLane 构造期间若读取 shared 上的 device 数据，构造器自身应在第一次读取前等待 shared 的初始化事件；常规构造仅创建 stream/复制案例输入，因此 wrapper 在构造后建立的共享依赖足够覆盖优化调用。

## 优化及多尺度语义

- `method="data"`：scipy L-BFGS-B，归一化到 `[0,1]` 的材料系数，显式解析梯度。
- `method="dual"`：相同 L-BFGS-B，目标由 lane 的 dual 分支给出；最终 `close` 仍须提供原物理 data 指标。此 wrapper 暂不自动增加 dual 后 data 校正阶段。
- `native_trf`：scipy `least_squares(method="trf", tr_solver="exact", x_scale="jac")`，解析 Full-J。
- `reciprocal_trf`：相同 TRF，变量为材料系数倒数；wrapper 精确施加 Jacobian 列链式因子 `-m**2`。该参数化没有改变物理材料网格。
- `data_lbfgs` 为 `data` 的别名。`globaldual` 在粗阶段采用 data L-BFGS-B，最后阶段先 dual proposal 再 data correction；由 `dual_fraction` 划分最后阶段的迭代/函数调用分配。它是本轮明确的执行合同，不等同旧冻结流程的时间分配。

默认只运行 stage 24。可选 `config.stages=[8,12,24]` 对所有方法开放相同 stage 权限。升尺度采用端点对齐双线性系数插值，与旧实现的线性插值基矩阵等价。`maxiter/maxfun` 对每个 stage 独立生效，`wall_seconds` 是各 stage 共享的任务优化 wall 上限；不分配旧研究 `.20/.25/.55` 阶段份额，不自动扩展预算。wall 在优化回调边界检查，正在执行的 CUDA 操作及最终物理核算不能被强制中断，收据保留实际耗时。初始化、public 文件读取和输出另有 worker/task/总 wall 记录。wrapper 不设置 LU/RHS 门槛。

各回调同时带回原物理目标；在线返回按完整 data loss 加共同正则选择的材料，包含各阶段/校正候选。dual proposal 值和物理 data 值分开记录，选择不访问 GT。若 Full-J 方法尚未提供 `residual/jacobian`，任务返回 `UNSUPPORTED_METHOD`，不会改用另一个算法。其它异常返回 `FAIL` 和 traceback，成本为 `INCOMPLETE_EXCEPTION`；没有自动重试。

可选 `svd_backend="hybrid"` 在本进程内把 SciPy TRF 的 SVD 调用适配至现成 CPU gesdd / cuSOLVER gesvd：小矩阵使用 CPU；大矩阵仅在非阻塞取得一个 GPU slot 后使用 CUDA，否则直接 CPU。它没有改变 TRF 的信赖域规则。源码及适用版本见 `src/hybrid_svd.py` 和 `REPRODUCE.md`。

## 输出与测量边界

每个任务的 JSON 包含全部 stage 的设备物理收据、优化状态、目标轨迹、累计 optimizer evaluation 数、worker wall、排队时间和完成顺序。`runner_summary.json` 包含总 wall（含共享初始化与最终 JSON 写盘）、共享导入次数、lane 数、成功/失败数量、各任务耗时与累计 optimizer evaluations。optimizer evaluation 数只计 optimizer 的目标/残差调用；Full-J 调用、close 的额外物理核算和精确 RHS/LU 成本以设备层收据为准。summary 写文件自身的最后几毫秒未计入保存的 wall。

FIFO 保证取任务的顺序；线程调度可能使开始事件和完成事件顺序不同。JSON writer 使用主线程，worker 不执行压缩/文件输出。线程间共享 BLAS 线程上限统一设为 1，避免每任务进入/退出全局线程配置所造成的竞争。

本 wrapper 已由主代理在原 4060 上完成 1/2/4 lane 及 4 lane hybrid 共 48 个任务；每一配置为相同 12 个任务，631 目标评估、635 LU、26,880 RHS。真实结果在 `RESULTS_SUMMARY.json`。它不修改旧冻结研究，不声称等价复现旧调度，也不据单纯 CUDA 并发推断科学收益。
