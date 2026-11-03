# R6：真实 FFT-TSOM 计算接线

状态：**完整两例四任务已执行并封存评分**。新入口是 `src/cie_fft.py::run_cie_fft`；旧 baseline、R2/R4/R5 solver 与 R6 primary dense benchmark 均未改动。完整长轨迹与实际费用见本轮RESEARCH_REPORT_ZH第5节、RESULTS/scoring_fft及RESULTS/analysis/FFT_ENGINEERING_MATCHED.json；下列有限接线检查单独保留，不冒充完整结果。

## 接口与来源

```python
from research.research_max_struct_r6.src.cie_fft import run_cie_fft

cfg = dict(existing_frozen_config)
cfg.update(arm='a_material',  # 或 cie_dual
           fft_geometry={k: public[k] for k in ('coords', 'radius', 'k')})
result = run_cie_fft(model, public['Y'], m0, cfg)
```

`result` 保留 `m/candidates/history/stage_states/cost/last_joint_m/status/receiver_rank/source_metrics`，原 `iteration_callback/stage_callback/receipt_path/cost_callback` 可用。调用者继续使用已有原始物理 corrector。`cie_dual` 默认 sourcewise PR；`a_material` 默认 global PR。若显式传入 `sourcewise_cg` 则遵从调用者，sourcewise 要求对角 source metrics。只支持这两条 CIE arm 和 PR；不将 bounded exact-current SVD 诊断或 CSI 冒充 FFT arm。

必须传入真实 uniform Cartesian geometry，不能从任意 dense G 猜测 Green 核。Geometry 与原 `model.G` 须来自同一原 disk-pulse assembler。辅助传播用新 FFT action；`model.state`、接收 `S`、候选原始 Phi、原材料 box/prior、readout 与物理 corrector 保持原模型。material prior 使用原 callable，未加 coefficient ridge，未换材料参数化。

算法循环取自 R5 `cie_factorial.py`（继承 R2 global loop），CIE sourcewise readout 调用原 FBG-R3 `_readout(J,Q,...,ws)`；A readout 直接调用 `onehop_material_readout`。PR+、非下降方向重启、解析二次步长、全局/逐源 coefficient stop、DPCS affine offset、beta current/material 分离与 Fourier 零填充 transport 都保留。矩阵乘法顺序改变带来浮点扰动，因此不承诺长非线性轨迹逐位相同或相同 endpoint。

## 真正的 action 与 adjoint

`src/fft_tsom_ops.py::NativeFourierAPCS` 实现

\[
U A=Q_d\,\mathrm{IFFT}_{\rm ortho}(M A),\qquad
U^*z=M^*\,\mathrm{FFT}_{\rm ortho}(Q_dz),\quad Q_d=I-V_dV_d^*.
\]

`M` 按历史 four-corner mask 同顺序 scatter/gather。native Fourier coefficient 的 Euclidean CG 度量未被 orthonormalization 改写。DPCS 投影仍为低秩 dense 操作，接收与接收伴随仍为 dense S/S*；不能把整个 solver 的每项都标成 O(N log N)。

`DiskPulseToeplitzGreen` 用 native grid lag 构造实际 disk-pulse 系数：非零 lag 是 `.5j*pi*k*radius*J1(k*radius)*H0(k*distance)`，零 lag 是 `.5j*pi*k*radius*H1(k*radius)-1`。与旧 `research_z_alg_r4_diversity/src/physics.py::build_operators` 相同，包含有限 disk-pulse 自项。对每轴补至至少 `2*n-1`（24 网格选48），负 lag 映射到补零 embedding 尾部；输入先补零，卷积后只裁原窗口。这是有限 Toeplitz 作用，没有用原24网格周期核替代原物理。adjoint 用 embedding spectrum 的共轭，配同样的补零/裁剪。

固定 m 下，CIE `p=beta/(1+beta*m)`，新的 current block 使用

\[
H a=p\odot(Ua-m\odot G(Ua)),\qquad
H^*z=U^*\{\bar p\odot z-G^*(m\odot\bar p\odot z)\}.
\]

因此 `U/GU/H/SU` 全部换成 action/adjoint，不仅把最后一次 `U@A` 换成 FFT。独立 `FixedCurrentBlock` 提供一般固定材料 `P/P*` residual composition 和 `readout(J,Q)` 回调；完整两arm入口用同一已知算子组合，并直接调用原 material readout。A 在这些实验中位于材料读出位置；FFT 不是新 A 科学公式。

## 数值零作用的明确差别

历史 `fourier_apcs` 显式生成 `Z=Qd F`，再用 SVD 仅截掉 roundoff-sized singular actions，最后仍重建到原系数坐标。新 action 不构造 Z、不做 APCS SVD，而直接实现 `Qd IFFT(mask A)`。因此二者在精确算术下是同一 native action，在浮点算术下并非同一矩阵。

`audit/FFT_WIRING_RECEIPT.json` 的同一真实24×24 public geometry，DPCS rank17：

| mask | native系数 | 历史APCS rank | 历史去除数值零作用 | 随机action相对差 |
|---|---:|---:|---:|---:|
| mf10 |400|400|0|7.82e-15|
| mf12 |576|559|17|7.87e-15|

mf12 的单位 DPCS null coefficient 测得：native current 绝对范数1.24e-15，历史SVD-cleaned current5.59e-15，两者绝对差5.74e-15。近零输出只看绝对误差；相对近零输出误差没有稳定解释。主1500更新计划止于mf10，mf12只属于显式追加 continuation。该记录没有声称 native numerical null action 完全为零，也没有重新把系数换成正交坐标。

## 有限实接线回执与成本

复现新接线检查：

```sh
python research/research_max_struct_r6/experiments/fft_wiring_check.py
```

检查仅打开已暴露开发弯条的公共 `G/E0/S/Y/grid_shape/coords/radius/k/host`，没有读取 GT。一次随机 Green action 与原 dense G 的相对差2.46e-15；新 Green adjoint pairing 差7.31e-15。随后每arm各做**一次**真实mf10 current更新、原readout和原物理候选验收（不是整库重复恒等式检查）。

| arm | current相对差 | 材料相对差 | dense wall/s | FFT wall/s |
|---|---:|---:|---:|---:|
| cie_dual |2.02e-15|5.18e-15|0.05200|0.02214|
| a_material |2.01e-15|4.21e-11|0.06488|0.03014|

每个执行均付2次原物理LU/24 complex RHS。此处 wall 包含 solver setup/readout/原物理候选验证，外部模型创建/公共算子导入不在每solver窗内；整个有限检查 wall0.276s。结果是新接线资格，不是长轨迹、完整反演速度或 GT recovery 结论。初始小回执的 `S_columns` 未计 post-update history 的一批12列接收乘法；源码已补计；当时sourcewise scalar readout也沿用旧计数模板多记2次whitening products，源码现按实调用修正。FFT transform及LU/RHS计数不受影响。

以mf10为例，历史 U/GU/H/SU 四组数组合计11,366,400 bytes；mf12为16,367,616 bytes。新入口不分配这些数组；保留48×48 complex spectrum36,864 bytes，mask系数分别76,800/110,592 bytes。还存在 FFT padded 临时workspace、Vd、stage/iteration snapshots及原物理 dense G/LU缓存；这些数组数字**不是**进程峰值内存或整个反演的内存收益。每次实跑 receipt 记录 batched FFT调用与独立2D transforms、dense接收列数、读出成本、LU/RHS、失败及wall/cpu。

本机无新依赖、无CUDA/远端作业、无GT、无commit/push。已有CUDA catalog/REPRODUCE 已读；本实现是与原 scalar-TM uniform disk-pulse 物理匹配的 CPU scipy.fft 接线，未把现成3D vector FFT当作本问题的替代算子。

## 完整独立工程运行器（已由父线程统一调度完成）

`experiments/run_fft_case.py` 复制主 `run_case.py` 的public加载、observer、原corrector、快照、封存和输出合同，只替换辅助solver调用，并从合并public取coords/centers、radius、k。差异存于 `audit/FFT_RUNNER_DIFF.txt`。主dense文件没有修改。

`configs/fft.json` 引用primary配置来源，完整继承3000辅助更新（3×500 + mf12×1500）、3000原corrector及100/500/1000/2000/3000 accepted checkpoints；仅增加独立backend/provenance。`configs/FFT_manifest.json` 固定旧b20d和历史obs_b2449bfb61f845bf、各cie_dual/a_material，共四任务；taskid加`__fft`，输出到`RESULTS/fft`。例如：

```sh
python research/research_max_struct_r6/experiments/run_fft_case.py \
  --manifest research/research_max_struct_r6/configs/FFT_manifest.json \
  --task-id obs_b20d5a55677d4eca__cie_dual__fft
```

实现worker仅做语法与CLI检查，没有重复启动完整作业；随后父线程实际完成四任务并独立封存评分。峰值RSS记整个child process high-water（macOS ru_maxrss按bytes；其他Unix按KiB换算）；它包括import、原dense物理model与snapshots，不能当作FFT独有workspace或直接和不同记账窗口的dense RSS相比。两例A的总CPU分别59.02→49.76秒、70.58→54.44秒，均保持恢复；同FFT下Dual更便宜。完整wall/FFT/原Phi/终点差异在正式报告和matched表，不能用早期接线PASS声称长轨迹完全一致。
