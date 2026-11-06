# Original3λ CM九臂复现入口

需要**完整项目仓库**，因为适配器复用R1 `a9_green.py`/`vector_backend.py`和R8辅助模块；只复制R9单目录不足。已有GPU环境须支持PyTorch complex128 CUDA/cuFFT/cuBLAS/cuSOLVER。本次实际环境为Windows、RTX4060Laptop8GB、PyTorch2.6.0+cu124、NumPy2.2.6、SciPy1.15.3、psutil7.2.2；已有驱动/环境直接复用，不安装CUDA/CuPy或新框架。

实际数值程序使用已核实的i9 P核logical12/14、BLAS/OpenMP1线程、Torch CPU线程1。数值入口依赖Windows的psutil affinity/peak_wset；其他主机需要明确记录的**仅亲和性/计时兼容适配**，不能静默改变算子、精度、RTOL、预算或模型。CPU与GPU/库版本不保证逐位相同，完整状态验收使用原矩阵相对残差≤2e−6。

`experiments/reproduce_original_cm.py`只是对实际已执行入口作输出路径路由，不调度任务。`--run-dir`必须显式指定；没有历史D盘/SSH目录默认值，不会覆盖已封存或已记录输出。以下PowerShell命令从完整仓库根目录执行，`python`指向已经存在的上述环境：

```powershell
$originalRun = Join-Path (Get-Location) "runs/original_cm_reproduction"
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" benchmark
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" generate
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" qualify
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" nine
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" score
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" plot
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" cost
```

`benchmark`是同一完整source59 CPU/GPU GMRES对照；合格GPU状态可复用。`generate`生成CM60³数据，逐源检查点可续算，仅对严格匹配材料/源RHS并通过原矩阵残差验证的状态复用。`qualify`是已经预声明的单次CM30³真值响应，供跨网格诊断和matched clean数据复用，不另增forward generation。`nine`实际进入`run_original_cm_nine.py`，一个驻留G/E/S/Vr、两个Y单独Jd及5+4固定臂，严格拒绝覆盖已有cohort。`score`在全部9seal存在前拒绝材料GT读取；`plot`输出两个cohort图与完整阶段/系数轨迹；`cost`生成明确additive/nonadditive的ACCOUNTING_FINAL。

本轮实际直接入口依次为：CPU `vector_original_extension.py generate`（两个工程中断保留费用）；Windows `benchmark_original_cuda.py`；`vector_original_cuda.py generate`（复用52Mac+1GPU检查点，新增7源）；`vector_original_cuda.py qualify-discretization`；`run_original_cm_nine.py`。离线评分、绘图与最终账使用`score_original_cm_nine.py`、`plot_original_cm_nine.py`、`finalize_original_cm_accounting.py`。薄路径入口已做CLI/语法核验，没有为了包装层额外重跑整套数值实验。新目录重算复现最终科学协议，不要求重复历史废弃工程时间。

若要单独重现已经完成的**两个RHS IGT负结果诊断**，仅在CM clean数据与coarse reference已存在时运行：

```powershell
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" igt-prepare
python research/research_max_struct_r9/experiments/reproduce_original_cm.py --run-dir "$originalRun" igt-solve
```

IGT不是九臂算子，没有60源IGT包或IGT反演；不能因source59差异有局部下降就称作者物理复现。九臂原尺寸Φ=mean_s(||pred_s−Y_s||²/||Y_s||²)，没有1/2且R=0，不能数值横比PhaseI Φ。CN10%是期望相对范数，本轮realized≈.099775414；matched保持同方向和相同比例规则，不重新归一化旧fine观测。

全扩展成本编译由主线程在本地完整记录到齐后执行：

```text
python research/research_max_struct_r9/experiments/compile_extension_costs.py --final
```

总账仅加不重叠进程行；每臂摊销/缓存归属估计为nonadditive，不能当独立冷启动或重复收费。历史Mac/Windows混合CPU总和是费用记录，不是单平台速度比较；每进程wall可能重叠或含暂停，不能相加称项目makespan。
