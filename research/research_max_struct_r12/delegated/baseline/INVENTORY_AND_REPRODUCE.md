# R12 Pro toy：成熟基线适配记录

本目录只新增适配器与结果，不修改 R6–R11 源码及 Pro 原文件。科学解释、创新判断归 root。全部反演只读取 `public/` 的 G/E0/S/Y 和公开噪声 σ；`score` 是读取原 truth 的独立离线步骤，必须先封存整组20预测。

## 实际实现与继承链

|功能|实际入口|本适配器如何使用|
|---|---|---|
|成熟 CIE 来源|`research/research_fbg_r3_experiments/src/cie_tsom_baseline.py::run_cie/_readout`|R6 observed copy 与 R8 native action 集成继承共享实材料 modified-contrast 读出|
|R6 原配置|`research/research_max_struct_r6/src/benchmark_solver.py::arm_config`，第34行|作为 β、PR、源权重、内层与原 corrector 参数依据；其硬编码24²及 [.02,20] 校验不适用于 toy，因此不调用 runner|
|R8 实际辅助求解|`research/research_max_struct_r8/src/cie_fft.py::run_cie_fft/_run_cie`，第21/71行|直接调用原函数，保留 PR+、精确固定材料二次线搜索、频率 transport、DPCS 仿射偏移与阶段端点原物理选择|
|CIE/Dual/Global/A+ dispatch|`research/research_max_struct_r8/src/benchmark_solver.py::run_auxiliary`|CIE 两臂 sourcewise；Global/A+ 用 global Frobenius PR。A+ 自动连接原 Möbius 读出|
|原 A+ 读出|`research/research_max_struct_r7/src/aplus_readout.py::mobius_material_readout`，第103行|原 four-color Möbius/MM、arc stationary candidates、两 sweeps、inherited gtol 与 fallback，未简化为短 alternating prototype|
|原 G19/fullGreen|`research/research_max_struct_r10/src/lowrank_denominator.py::prepare_green_strong/matrix_material_readout`，第128/86行|相同公开 G 右奇异谱 .1 阈值、Galerkin PGP、完整材料导数、每次8步 bounded L-BFGS 及原 fallback|
|R11 fullG19 继承接口|`research/research_max_struct_r11/src/material_readouts.py::prepare_full/full_material_readout`|该接口转发到上一行；未采用 freezeA/contractedT/BothBlocks 变体|
|R9 公开 noise DPCS|`research/research_max_struct_r9/src/baseline_adapter.py::noise_ball_rank`，第111行|按完整99%圆复高斯噪声球挑选最小满足 tail 条件的受支持 receiver rank，σ直接取 Pro 已公开记录|
|原完整物理模型|`research/research_fbg_r3_experiments/src/physics.py::MaxwellModel`，第16行|导入精确 Pro G/E/S；原 retained-LU、多 source RHS、analytic adjoint、原物理 Φ；无新正演近似|
|原物理 corrector|`research/research_max_struct_r6/src/benchmark_solver.py::refine_observed`，第164行|复制其数值 body，仅把 bounds 改成 toy [0,12]。scale1e6，ftol1e−12，gtol.01，maxls30，maxiter3000，maxfun12000；保留 paid minimum-Phi incumbent|

## 明确的适配边界

- 精确100个独立实材料未知量，6 sources，24 receivers，初值全部 .1，实 box[0,12]；不施加 support、GT 或形状先验。保留原 Regularizer 接口但 smooth=tv=0，因 Pro 合同为无先验。
- Pro 的 equal-area disk 半径 dx/√π 大于 dx/2，成熟 FFT 类的 nonoverlap guard 会拒绝该几何。适配器只在运行时把 R8 的 Green-action 构造器替换成**提供的精确 dense G 乘法及伴随**。没有重造 Green，没有改变给定 G；不能将这些运行叫做 FFT 复杂度证据。继承 receipt 中 “FFT Green” 的固定文字不适用，实际操作看 `dense_G_*` 计数。
- 原24²的 mf3/6/10/12 在10²网格不全合法。本目录保留两套预先固定且不按 GT 选择的映射：`scaled` 为按10/24缩放并 half-up rounding 得 mf1/3/4/5；root 另授权的 `clipped` 为把原 mf 仅截到 Nyquist：mf3/5/5/5。二者各20运行，不能用表现较好者冒充唯一历史原合同。
- 两套各阶段迭代500/500/500/1500，合计3000，无新的 wall 截断。material β=1/.2/.2/.2；CIE 的 current β同 material；Dual、Global、A+、Green 的 current β恒1。source-wise/global 规则未改。coefficient_tol=0；阶段按 MAXITER自然终止。
- mf5含恰好100个不重复 Fourier frequency，投影后自由 APCS rank 小于100，`CAPACITY.json`明确列出。若 receiver rank12，APCS为88 complex/source，即1056 real free current DOF；DPCS offset固定。材料仍是100 real DOF，不能拿这些量互相替代。
- Green 谱阈值 .1 在本 Pro G 上**自然给出rank19**，不是强行设定。Green denominator 的 full 是完整材料导数，不是 K=G。在线标签保守写为 `green_strong_full`，本 toy 可注明实际 rank19。
- 辅助 current 的 DPCS 受限仿射族仍不是 Pro projection 方法的全部自由 current。material DOF相同不代表 admissible current families相同。
- custom Green 的 history 字段 `state_term_material` 仍是 inherited onehop diagnostic；真正新材料目标看 `material_final_block_loss`。final 原物理 Φ是统一的。
- 正演与伴随、current/denominator操作均计数；CPU是实际耗时，固定迭代预算不能推出等耗时性能。`PAID_STATES.npz`及评分文件额外给共同CPU前缀下最低已付费Φ快照；初始点时间记0，是保守起点标签，阶段及corrector时间实测。快照只覆盖阶段端点与corrector评估，不涵盖未做正演的每个辅助步。
- `OFFLINE_SCORES*.json` 读取 truth 只在对应20结果全部 `SEALED_PUBLIC_PREDICTION` 后。材料指标不参与 online selection/stop。

## 可执行命令

工作目录为现有 `r12_worktree`。解释器固定为 `/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`，脚本强制所有 BLAS/OMP thread 为1，并用 threadpoolctl 限制1。没有额外安装依赖。

```bash
PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
SCRIPT=research/research_max_struct_r12/delegated/baseline/toy_baselines.py
"$PY" "$SCRIPT" prepare
"$PY" "$SCRIPT" run --case ring_high --arm cie_single
"$PY" "$SCRIPT" run --suite clipped --case ring_high --arm cie_single
```

对四个 case（ring_high、split_shell_high、bridge_high、disk_low）与五个 arm（cie_single、cie_dual、global_dual、a_plus、green_strong_full）各运行一次。不得在已有输出路径重跑覆盖。两套各20预测完成后执行：

```bash
"$PY" "$SCRIPT" seal
"$PY" "$SCRIPT" seal --suite clipped
"$PY" "$SCRIPT" score
"$PY" "$SCRIPT" score --suite clipped
```

`runs/` 与 `runs_clipped/` 分别保留两套结果；每个 task 包含原配置、公开DPCS、current capacity、完整辅助history、corrector stop/history、操作成本与精确prediction。`SEALED_COHORT*.json`为 offline truth 前的整组索引。无 push，无 hash 检查，无在线GT。
