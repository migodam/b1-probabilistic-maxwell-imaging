# Observable19 接线说明

仅实现 supplement 内的新 factory、wrapper、配置与验证记录。理论与最终科学解释由父线程负责。原 R14/C2、R4–R13 源码均不编辑，不启动远端或144任务，不读取GT。

固定集合为已曝光24案例 × 六臂：四个 primary G/S × Galerkin/right，加两个事先声明的 anchor-observation × Galerkin/right control。所有 task 的 `algorithm='g19'`，`method/label` 保留六臂身份；576个真实材料未知数保留。protocol 与原 `R14/configs/independent24_v1.json` 完全相等，公共输入只引用 `../data_public/independent24`。完整清单见 `configs/observable19_exposed24_v1.json`，最终冻结核验见 `audit/MANIFEST_FINAL_BEFORE_RUNS.json`。这批数据明确是 **EXPOSED24_FOLLOWUP**。

## 原代码 seam

- `src/run_independent.py` 只用原 C2 的完整 main：同进程 `runpy` 调用父目录 `runtime/equal_budget_C2_runner.py`。对其与 `runtime/frozen_C2_runner.py` 的既有单行 `deadline=cap` 差异进行直接源码比较，不复制/重写主程序。
- 远端 legacy 固定为 `D:/AI/MAX_STRUCT_R11_DEV_20261107`，不假设它含 R13 package。`runtime/frozen_full_current.py`、`field_current.py`、`material_events.py` 是实际 R13 原文件的逐字节副本；来源和直接 byte equality 记录在 `audit/RUNTIME_SOURCE_MAP.json`。
- factory 用 importlib 载入 owned `frozen_full_current.py`。其原 `run_full_current` 通过原 `_local_function` 绑定私有 `prepare_full/full_material_readout`；原 package/module globals 不更改。C2 的顶层 `from full_current import run_full_current` 只在该任务进程中使用一个私有 import shim，`finally` 恢复此前 import slot。
- 实际导入路径依次为 supplement src/runtime、原 R14/runtime、legacy root/vendor。没有 R13 package 依赖；R8 `cie_fft`、R6 Fourier/DPCS、R11 `material_readouts` 和 R10 `lowrank_denominator` 来自 legacy 的真实路径。
- 原 `run_full_current(...mode='restricted',g19=True)` 保留付费 physical-current warm start、原 real PR/line search、stage/β/source metrics/prior/停止与审计。原 R10 `matrix_material_readout` 保留 moving denominator 完整解析导数、LBFGS8、原 scalar seed/old 比较与接受规则。原 C2 minPhi 候选、secondary、prefix、proposal≤72及绝对120余额 corrector 均不变。

## 六种 basis/作用

`g19_internal_galerkin` 使用原 R11 `prepare_full` 代码及其完整原准备链，要求原 G strong(.1) rank 恰好19。通过私有 SVD observer 捕获已经返回的 G/S 谱，不新增 SVD、不改变任何返回值。其他臂付费计算 G/S SVD；取指定算子的前19个右奇异向量，使用原 `LowRankDenominator` 的 `projection='galerkin'/'right'`。不通过 gap 选维，也不对 S/O0 奇异值取逆。

- Galerkin: 原 class 构造 `B=U(U*GU)`，作用 `KX=B(U*X)`。
- Right: 原 class 构造 `B=GU`，作用 `KX=B(U*X)`。
- Raw observation basis 从未scaled的 S 计算。原CN whitening 是共同scalar，不改变右子空间；factory无需noise σ，实际 source/noise 合同仍由原 C2 消费。
- Anchor observation 在相同公开 `m0=.1`、原C2已经付费的 `model._current` 上构造 `O0=S(I-D_m0G)^-1`。只使用 `state['lu'] is model._physics.factor` 的匹配缓存，在原 `_physics.operation` 内做一次 `lu_solve(lu,S.conj().T,trans=2).conj().T`，按 nr 个 adjoint RHS 收费，额外LU=0、额外forward=0。不形成完整逆，不调用 `model.state` 补anchor，不更换anchor。它是 current defect→port transport，不标为 material Jacobian。

每个任务保存完整 `DENOMINATOR_BASIS.npz` 的 U/B，保存已经计算的 G/S/O0 全部奇异值、σ19/σ20 与 relative gap=(σ19−σ20)/σ1；不形成 O0 的臂明确给null原因，不为审计另算SVD。谱与gap只描述固定19的切点，不改变rank。初始缓存/支持rank不足时保存失败，不降rank或悄悄新增forward。

## 固定 readout 快照

六臂同一个 readout observer，在累计调用1/128/512/1500保存输入 J/Q/m_before、实际 beta/rho/Os/Ls/U/B 和已返回的 m_after。原 helper 直接调用，没有新优化、候选或在线GT诊断。快照明确写 `MATERIAL_READOUT_INPUT_OUTPUT_NOT_FULL_PR_CHECKPOINT`，不包含完整PR memory。若原调用被 BaseException/BudgetEnd 中断，只保存已有输入并标 `INTERRUPTED`，不伪造返回材料。copy/I/O全部处于原整进程计时之内；原小LU/triangular RHS/helper计数保留。

R13 adapter 原 receipt 硬编码的 setup SVD 数目与这些 factory 不完全相同。原报告保留为 `inherited_denominator_setup_cost_before_factory_accounting`，`denominator_setup_cost` 改写成实际 factory 计数，数值算法不改：internal Galerkin为原 G/S/residual 三SVD；其余 raw G/S 控制为两个SVD；anchor再加一个O0 SVD和nr个原物理伴随RHS。rank19 dense G×U 建立费、所有小系统solve、source product及IO均按实际记录，physical RHS与small/block RHS不混加。

## 已做有限验证及限制

`audit/wiring_validation_v2` 是最终 unscaled S/O0、portable runtime 与谱observer版本。只用公开 r14_new_0000，两侧各三个固定原native更新，不做corrector恢复、不读GT。m/J/coeff/U/B/K(E0)全部bitwise相等、最大差0；两侧均3次material readout、各36个19维smallLU、864个small triangular RHS。另一个prepare-only anchor复用同一初始缓存，新增48 adjoint RHS、0LU/forward；O0 rank19资格通过。

24案例纯公开 G/S header均通过：原G强rank均19，G数值支持rank576，S数值支持rank24/48。仅首例O0已做实际cached-LU资格检查，其余23例会在各自真实setup中检查，不能把raw S的数值rank当作O0已完成验证。GP理论性质由父线程文件说明，这里不额外设质量gate。

初版固定同一案例验证也保留在 `wiring_validation_v1`，没有隐藏或退款。两版有限验证及一次真实 import-only preflight总OS子进程CPU4.550838s，全部验证共8LU、192原物理complex RHS；没有远端/144任务启动。header/小系统/IO费用和各source见每版receipt与 `audit/VALIDATION_COST_TOTALS.json`。`audit/import_preflight` 已确认 owned runtime三文件被实际resolve；本机仅复用已有ABI兼容psutil，未安装环境。远端应用其现有psutil。

## 远端真实入口（此处未执行）

先在原有数值Python环境运行 import-only preflight，再由父线程审查/启动：

```powershell
$r14Root = 'D:/AI/MAX_STRUCT_R14_DUALMETRIC_20261115'
$supplementRoot = "$r14Root/supplement_observable19"
$legacyRoot = 'D:/AI/MAX_STRUCT_R11_DEV_20261107'
python "$supplementRoot/src/run_independent.py" --root $supplementRoot --legacy-repo $legacyRoot --manifest "$supplementRoot/configs/observable19_exposed24_v1.json" --task-id r14_new_0000__g19_anchor_observation_right --check-imports
```

完整有限队列可直接复用现有 R14 `queue_independent.py`，因为它按给定root查找本supplement的 `src/run_independent.py`；同case六臂仍在同核顺序执行，外部整子进程CPU照原队列记录。core列表/资源审计由父线程按原host合同核对：

```powershell
python "$r14Root/src/queue_independent.py" --root $supplementRoot --legacy-repo $legacyRoot --manifest "$supplementRoot/configs/observable19_exposed24_v1.json"
```

单任务原C2入口（只有父线程决定执行）为同一个 wrapper 不带 `--check-imports`；运行前已有输出目录会拒绝覆盖。所有新任务输出在 supplement/RESULTS，不改原预测或原评分。尚未评分，也没有新科学/新颖性结论。
