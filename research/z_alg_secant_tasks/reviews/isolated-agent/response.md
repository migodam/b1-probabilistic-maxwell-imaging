原回复。文中的本机文件链接已换为公开研究资料或附件来源说明。

# Z-ALG 独立反方证据审稿报告

> 审稿约束确认：本审查全程只读，未向工作区写入或修改任何文件。以下所有数值均由只读 Python 直接从工作区 `.npz` 二进制数组、`.json` 区间记录与源码公式独立复算得出。

---

## 1. 工具运行状态与审计口径说明

- **可直接运行**：只读终端命令与系统 Python 3.13.13（标准库 `zipfile` / `struct` / `json` / `decimal` / `math` / `cmath`），以及文件读取工具。
- **无法直接运行的依赖项（逐项记录）**：
  1. `numpy`、`scipy`、`python-flint`、`mpmath` 在隔离环境的所有系统 Python（`/opt/homebrew/opt/python@3.13/libexec/bin/python3`、`/opt/homebrew/bin/python3.12`、`/usr/bin/python3`）中均**未安装**（`ModuleNotFoundError`）。
  2. 因此，无法直接运行原工程的 192 位 Arb 认证脚本（如 [full_port_radius.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/full_ports/full_port_radius.py)、[certify_separated.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/certify_separated.py)、[interval_physics.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/cert/interval_physics.py)）以及依赖 SciPy 的 [cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py)、[replay.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/replay.py)。
- **替代独立复算手段**：
  1. 用纯 Python 标准库直接解析 `.npy`/`.npz` 二进制浮点与复数数组；
  2. 用 60 位精度 `decimal` 配合显式 `ROUND_FLOOR` / `ROUND_CEILING` 有向舍入区间算术，独立复算全部完整跨端口区间半径；
  3. 用纯 Python 实现 Bessel/Hankel 级数与渐近展开 + 复数高斯消元，独立重建离散脉冲正演 `FullModel`（[separated.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/separated.py#L83-L157)）与连续圆柱 Fourier–Bessel 正演 `Cylinders`（[cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py#L36-L111)），并与存档 `.npz` 交叉校验：
     - `Pulse` vs [conditioning/materials_and_ports.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/materials_and_ports.npz) (`L2_R30__AA`): Frobenius 误差 **$2.770\times 10^{-16}$**；
     - `Cylinders(M=4)` vs [continuum/materials_and_ports.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/materials_and_ports.npz) (`M16_L2_AA`): Frobenius 误差 **$2.051\times 10^{-15}$**；
     - 36 盘历史宿主 `Pulse` vs [historical_inputs.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/historical_inputs.json): $\|F_p(A)-y_-\| = 2.010\times 10^{-15}$，$\sigma_{\min}(J_p(A)) = 6.707108\times 10^{-6}$、$\sigma_{\min}(J_p(B)) = 7.540188\times 10^{-6}$（与 [conditioning/README.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/README.md#L9) 报告的 `0.00000670711` / `0.00000754019` 逐位一致）。

---

## 2. 从原始区间与复数数组直接复算的主张与反方审计

### 2.1 主复算对象：完整跨端口复数矩阵与 10% 有限噪声球主张

**待审主张**：[full_ports/README.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/full_ports/README.md) 与 [full_ports/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/full_ports/results.json) 基于 [krawczyk.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/results/krawczyk.json) 与 [conditioning/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/results.json) 中的复区间外包，按共同中心矩阵 $Y_*$（对角块为重复本地目标 $y_*$，跨端口块为各根区间中点均值 $\bar{Y}_{\text{cross}}$）计算：
$$
R_{\text{lo}} = \frac{1}{2}\max_{\alpha\ne\beta}\|F_{\text{all}}(m_\alpha)-F_{\text{all}}(m_\beta)\|_F,\qquad
R_{\text{hi}} = \max_{\alpha}\|F_{\text{all}}(m_\alpha)-Y_*\|_F,\qquad
\rho = \frac{R_{\text{cov}}}{\|Y_*\|_F}.
$$
据此声称：在 10% 联合 Frobenius 噪声球下，历史 6 槽两/三组件（$40{,}000\lambda / 200{,}000\lambda$）分别可共球容纳全部 4 根 / 8 根，而偏置 2 槽（$30\lambda$）全 4 根 / 全 8 根被排除（三组件甚至逐对分离）。

**独立 60 位有向舍入区间复算结果**（另补充未纳入 `full_ports` 的**跨口自身归一化**、**直达入射场 $E_{\text{inc}}=\frac{i}{4}H_0^{(1)}(kd)$ 归一化**，以及连续圆柱 `M=16` 复算）：

| 组别与配置 | 复算覆盖半径 $[R_{\text{lo}}, R_{\text{hi}}]$（与 Arb 差） | $\|Y_*\|_F$ / $\|Y_{\text{cross}}\|_F$（跨口占比） | 占全散射信号 $R_{\text{cov}}/\|Y_*\|_F$ | **占跨端口信号 $R_{\text{cov}}/\|Y_{\text{cross}}\|_F$** | **占直达入射场 $R_{\text{hi}}/\|E_{\text{inc}}\|_F$** |
|---|---|---|---:|---:|---:|
| 历史 6 槽 `L2_R40000` (4 根) | $[7.836535, 7.861144]\times 10^{-6}$ ($\le 1.7\times 10^{-21}$) | $0.165823$ / $4.5466\times 10^{-4}$ (**0.274%**) | $0.004726\%$–$0.004741\%$ | **$1.724\%$–$1.729\%$** | $0.00462\%$ ($\|E_{\text{inc}}\|_F=0.17004$) |
| 历史 6 槽 `L3_R200000` (8 根) | $[6.326714, 6.329245]\times 10^{-6}$ ($\le 8.5\times 10^{-22}$) | $0.203091$ / $3.2152\times 10^{-4}$ (**0.158%**) | $0.003115\%$–$0.003116\%$ | **$1.968\%$–$1.969\%$** | $0.00304\%$ ($\|E_{\text{inc}}\|_F=0.20825$) |
| 偏置 2 槽脉冲 `L2_R30` (4 根) | $[0.00699331, 0.00753255]$ ($\le 8.7\times 10^{-19}$) | $0.025836$ / $0.003813$ (**14.76%**) | $27.068\%$–$29.155\%$ | **$183.4\%$–$197.5\%$** | **$8.727\%$** ($\|E_{\text{inc}}\|_F=0.08632$) |
| 偏置 2 槽脉冲 `L3_R30` (8 根) | $[0.00998759, 0.01009783]$ ($\le 1.7\times 10^{-18}$) | $0.031870$ / $0.006017$ (**18.88%**) | $31.339\%$–$31.685\%$ | **$166.0\%$–$167.8\%$** | **$9.377\%$** ($\|E_{\text{inc}}\|_F=0.10769$) |
| 连续圆柱 `M=16, L2_R30` (4 根) | $[0.0063363, 0.0067707]$（`.npz` 复算） | $0.024400$ / $0.003547$ (**14.54%**) | $25.97\%$–$27.75\%$ | **$178.6\%$–$190.9\%$** | **$7.844\%$** |
| 连续圆柱 `M=16, L3_R30` (8 根) | $[0.0090904, 0.0091620]$（`.npz` 复算） | $0.030094$ / $0.005612$ (**18.65%**) | $30.21\%$–$30.44\%$ | **$162.0\%$–$163.3\%$** | **$8.508\%$** |

**算术结论 vs 物理审稿结论**：
1. **代码算术完全准确**：所有四组 Arb 区间端点与我们的 60 位有向舍入复算吻合至 $10^{-18}$–$10^{-21}$。
2. **反方核心发现一（未白化全局范数制造的尺度假象）**：在历史 6 槽 `L2_R40000` 与 `L3_R200000` 中，局部散射信号为 $\Theta(1)$ 而跨口信号按二维柱面波衰减为 $\Theta(R^{-1/2})$，导致整个跨口块仅占 $\|Y_*\|_F$ 的 **0.274%** 与 **0.158%**。所谓“10% 全矩阵噪声球容纳全部 4 根/8 根”，其噪声半径（$0.01658$ / $0.02031$）是**整个跨口信号本身的 36.5 倍与 63.2 倍**——任何 $R\to\infty$ 的系统只要使用未白化的全局 Frobenius 球，都会把原本与 $R$ 无关的 $\Theta(1)$ 跨口相对分离度（`1.73%` / `1.97%`）压成 $O(R^{-1/2})$ 并必然落入固定百分比球内。
3. **反方核心发现二（直达入射场扣除使 2 槽 30$\lambda$ 的 10% 判决直接反转）**：在 2 槽 $30\lambda$ 阵列中，散射场仅占直达入射场的 `29.9%`（跨口散射仅占跨口直达场的 `18.5%`）。若噪声或背景扣除标定误差按总场/直达场 $\|E_{\text{inc}}\|_F$ 的 10% 计，则脉冲与连续圆柱的覆盖半径分别跌至 **`8.73% / 9.38%`**（脉冲）与 **`7.84% / 8.51%`**（连续圆柱 `M=16`），**全部落入 10% 球内**，“排除全部根”的结论直接翻转。

---

### 2.2 补充复算：高维复数数组与任务约束 Gauss–Newton（`study/highdim/`）

直接从 [r10_twod_0001.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/r10_twod_0001.npz)、[austria_eps4_400MHz.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/austria_eps4_400MHz.npz)、[OFFLINE_SCORE.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/OFFLINE_SCORE.json) 与 [GN_BUDGET_COMPARISON.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/GN_BUDGET_COMPARISON.json) 复算关键数字：

- **24² 环 (`r10_twod_0001`, 同网格 24²→24²)**：$\|Y_{\text{arch}}\|=1.00249709$，$r=0.01051718$（$1.0491\%$）。对 `r12_g19__projected_task__025`：任务差 `0.376464`（`3.611687` vs `3.235223`），端间场距 `0.000375`（`0.03563 r`），两端距旧 $Y$ 为 `0.99094 r` / `0.98641 r`，偶次漂移 `0.37286 r`。但在材料能量分布上，**`49.54%` 的 $\|\Delta m\|_2^2$ 落在 ROI 之外**（ROI 内 $L^2=8.4572$，ROI 外 $L^2=8.3790$）。同预算对照下，弱投影与普通任务约束 GN 跨度完全打平（`0.37500` vs `0.37500`）。
- **32² Austria (`austria_eps4_400MHz`, 跨网格 150²→32²)**：$\|Y_{\text{arch}}\|=10.60397148$，R12 容差 $12\%$（$r=1.272477$），原 R8 噪声上界 $9.9292\%$（$r=1.052892$）。复算 `cie_dual__projected_task__025`：任务差 `0.089977`，端间场距 `0.003937`（`0.00309 r`），两端距旧 $Y$ 为 `1.040413 / 1.040332`（均在原 R8 球内），但 **`48.81%` 的扰动能量在 ROI 外**且端点触碰非物理盒顶 `20.0000`（原 Austria 对比度上限仅 3）。更重要的是：
  1. **粗网格真值本身严重出球**：[OFFLINE_SCORE.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/OFFLINE_SCORE.json#L22) 显示 Austria 粗网格真值相对残差高达 **`22.526%`**（远超 `9.93%` 与 `12%` 球），且三份 SRAP 观测的粗网格真值残差为 **`13.714%`、`14.946%`、`11.816%`** 全部出球——这意味着在跨网格高维场景中，**离散建模误差已经主导了噪声球语义**，在粗网格球内构造的两点下界并不能约束真实细网格/连续物理目标的恢复误差；
  2. **算法优势反转**：在 Austria 原 R8 精确噪声球内，同物理调用数（7 LU + 1 FullJ）下标准任务约束 GN 跨度为 **`0.13125`**，反而比弱投影（`0.08750`）大 **+50%**。

---

## 3. 核心审查：连续圆柱与离散种子能否区分「物理效应」「模型误差」与「测量几何」？

**结论：不能。当前连续圆柱对照（[cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py)）与两套离散种子（历史 6 槽、[conditioning.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/conditioning.py) 偏置 2 槽）在设计上把三者严重耦合，且我们的独立复算已直接抓出两处决定性反证。**

### 3.1 物理效应 vs 模型误差：历史 36 盘宿主在连续圆柱下 Fold 局部消失，且坐落在单盘 Mie 共振上

离散脉冲模型（[separated.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/separated.py#L115-L118)）假设每个半径 $a=0.04\lambda$ 圆盘内部电场恒定，给出单盘单极散射系数
$$
c_0^{\text{pul}}(\chi)=\frac{u\chi}{1-\chi g_{00}},\qquad u=\frac{i\pi ka}{2}J_1(ka),\quad g_{00}=\frac{i\pi ka}{2}H_1^{(1)}(ka)-1,
$$
而精确二维 TM 连续介质圆柱的 $n=0$ 阶 Mie 散射系数为（$k_i = k\sqrt{1+\chi}$）：
$$
c_0^{\text{cyl}}(\chi)=\frac{k_i J_0'(k_i a)J_0(ka)-k J_0(k_i a)J_0'(ka)}{k J_0(k_i a)H_0^{(1)\prime}(ka)-k_i J_0'(k_i a)H_0^{(1)}(ka)}.
$$

我们直接计算了不同 $\chi$ 下的单盘误差与**历史 36 盘宿主（[historical_inputs.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/historical_inputs.json)）的连续圆柱重放**（该项目自身从未在 36 盘宿主上运行 `cylinders.py`，并在 [06_CLAIM_LEDGER_AND_NEXT_STEPS.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md) 列为未解决缺口）：

1. **单盘单极 Mie 共振与 Richmond 脉冲失效**：
   - 在 $a=0.04\lambda$ 处，精确连续圆柱的**第一单极 Mie 共振峰**（达到无损光学定理幺正上限 $|c_0^{\text{cyl}}|=1$）位于 **$\chi = 18.36$**（$k_i a = 1.1058$，实算峰值 $|c_0^{\text{cyl}}| = 0.99999996$）。
   - 历史 6 槽种子 `seed_A` / `seed_B` 的 6 个活跃对比度中有 3 个直接压在共振峰上或跨过共振峰：**$\chi_0 = 18.88/19.12$、$\chi_4 = 19.26/19.33$、$\chi_5 = 26.17/25.77$**（盘内 $J_0(k_i a)$ 跌至 `0.708`–`0.615`，盘内场根本不是常数！）。
   - 对应单盘单极相对误差：$\chi=8$ 时 **`10.93%`（$-3.56^\circ$）**，$\chi=19$ 时 **`15.69%`（$-9.02^\circ$）**，$\chi=26.17$ 时 **`13.82%`（$-7.52^\circ$）**。

2. **历史 36 盘宿主在连续圆柱 PDE 下的决定性反证**（`Cylinders` 在 $M=1\to 2$ 已收敛至 $<3\times 10^{-6}$）：
   - **同材料模型误差高达 ~10%**：$\|F_{\text{cyl}}^{(M=2)}(A)-F_{\text{pul}}(A)\| = \mathbf{0.012150}$ (**`9.87%`**)，$\|F_{\text{cyl}}^{(M=2)}(B)-F_{\text{pul}}(B)\| = \mathbf{0.012171}$ (**`9.90%`**)。这一离散化误差比 $R=40{,}000\lambda$ 跨口覆盖半径（$7.86\times 10^{-6}$）大 **1,546 倍**。
   - **两根在连续模型中不再等数据**：$\|F_{\text{cyl}}^{(M=2)}(A)-F_{\text{cyl}}^{(M=2)}(B)\| = \mathbf{0.000183}$（`0.15%`），是 $L=2, R=40{,}000\lambda$ 最小跨口分离（$1.01\times 10^{-5}$）的 **18 倍**、$L=3, R=200{,}000\lambda$（$5.19\times 10^{-6}$）的 **35 倍**。
   - **关键证据——Jacobian 行列式变号消失与伙伴根坍缩**：

| 评估点（36 盘宿主） | 离散脉冲 `Pulse`：$\det(J)$ / $\sigma_{\min}(J)$ | 连续圆柱 `Cyl(M=1)`：$\det(J)$ / $\sigma_{\min}(J)$ |
|---|---|---|
| `seed_A` | **$+4.0289\times 10^{-19}$** / $6.7071\times 10^{-6}$ | **$-1.0796\times 10^{-17}$** / $1.7075\times 10^{-4}$ |
| 中点 `(A+B)/2` | **$+1.0908\times 10^{-20}$** / $1.8712\times 10^{-7}$ | **$-1.0916\times 10^{-17}$** / $1.7680\times 10^{-4}$ |
| `seed_B` | **$-4.2554\times 10^{-19}$** / $7.5402\times 10^{-6}$ | **$-1.1083\times 10^{-17}$** / $1.8415\times 10^{-4}$ |

   在离散脉冲中，`seed_A` 与 `seed_B` 跨过 $\det(J)=0$ 的折叠面；而在连续圆柱 `Cyl(M=1)` 中，**整个线段上 $\det(J)$ 恒为负且几乎不变**，中点 $\sigma_{\min}$ 暴涨 **945 倍**！以 `seed_B` 为初值、以 $F_{\text{cyl}}(\text{seed}_A)$ 为目标执行 Newton（`M=0` 与 `M=1`），均在 3 步内**直接坍缩到 `seed_A`**（$\|z_{\text{stop}}-\text{seed}_A\|=0.000000$，残差 $9.08\times 10^{-15}$）；若以脉冲观测 $y_-$ 为目标，两初值在 `Cyl(M=1)` 中均不收敛（残差卡在 `8.10e-3` / `6.85e-3`）。

3. **为什么 2 盘连续圆柱对照 [cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py) 看起来“验证了”离散机制？**
   - [cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py#L151-L154) 对每个角阶 $M$ 均**重新用自身正演生成目标 `one.target = y_a[one.local]`**，完全规避了跨模型误差；
   - 我们对同一份 `L2_R30` 脉冲根 $m$ 做同材料双模型前向差实算：$\|F_{\text{cyl}}^{(M=4)}(m)-F_{\text{pul}}(m)\|_F = \mathbf{0.002136\text{–}0.002426}$（**`8.70%–9.59%`**），其中绝大部分来自 `Pulse → M=0` 的**径向内场近似误差**（`full_norm` 从 `0.02584` 跌至 `0.02411`，即 `-6.7%`），而 `M=0 → M=16` 仅变化 `+1.2%`；
   - 该 `8.7%–9.6%` 模型误差是 `0.013λ` 对称破缺信号（`0.000451`）的 **4.7–5.4 倍**，是 $30\lambda$ 局域跨组件耦合残差（`0.000898`）的 **2.4–2.7 倍**！2 盘连续模型之所以在 $M=0..16$ 均保住 4/4 与 8/8 根，只是因为 `rx_x=0` 时两盘具有精确镜像对称（两个模型共享同一拓扑对称性），并不能证明非对称强对比宿主（如 36 盘历史 Fold）能跨过 10% 径向模型误差而存活。

### 3.2 物理效应 vs 测量几何：旋转阵列平移轴使全端口覆盖半径坍缩 133 倍

在 [conditioning.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/conditioning.py) 与 [cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py) 中，单组件两圆盘位于 $(\pm 0.18\lambda, 0)$（沿 $x$-轴排列），本地源/收位于 $(0,-0.85\lambda)$ 与 $(0.013\lambda, 0.95\lambda)$（沿 $y$-轴前向照明），而多组件恰被安排**沿 $x$-轴平移** $(R, 0)$。

这意味着跨组件波束恰好以 $90^\circ$ 侧视入射两盘连线，双盘间产生最大光程差 $0.36\lambda$（相位差 $130^\circ$），使 `AB` 与 `BA` 跨口实部直接反号（`+0.00453` vs `-0.00403`）。这究竟是“$30\lambda$ 下 2 槽种子的物理分辨性”，还是**人为选定的平移方向几何红利**？

我们保持**完全相同的组件、对比度盒、局部收发点（含 `rx_x=0.013`）与间距 $R=30\lambda$**，仅将第二组件的平移矢量从 $(30\lambda, 0)$ 改为 $(0, 30\lambda)$（沿 $y$-轴）与 $(21.21\lambda, 21.21\lambda)$（$45^\circ$），独立求解完整耦合 4 根并计算全端口覆盖半径：

| 模型与配置（$L=2, R=30\lambda$） | 收敛根数 / $\min\|\Delta m\|_2$ | $\|Y_*\|_F$ / $\|Y_{\text{cross}}\|_F$ | 全端口覆盖半径 $[R_{\text{lo}}, R_{\text{hi}}]$ | 占全散射信号 | 占跨口信号 | 10% 噪声球判决 |
|---|---|---|---:|---:|---:|---|
| 脉冲 `rx_x=0.013`, **$x$-轴平移（原文）** | 4/4 / $8.3430$ | $0.025836$ / $0.003813$ | $[0.006993, 0.007533]$ | **27.07%–29.15%** | 183.4%–197.5% | 排除全 4 根 |
| 脉冲 `rx_x=0.013`, **$45^\circ$ 平移** | 4/4 / $7.6134$ | $0.026164$ / $0.005619$ | $[0.004348, 0.004426]$ | 16.62%–16.92% | 77.4%–78.8% | 排除全 4 根 |
| 脉冲 `rx_x=0.013`, **$y$-轴平移** | 4/4 / $7.1576$ | $0.026953$ / **$0.008573$** | **$[0.000056, 0.000058]$** | **0.209%–0.213%** | **0.66%–0.67%** | **0.22% 球即容纳全 4 根** |
| **连续圆柱 `M=4`, `rx_x=0.013`, $y$-轴平移** | 4/4 / $7.0719$ | $0.025431$ / **$0.007999$** | **$[0.000050, 0.000051]$** | **0.197%–0.202%** | **0.63%–0.64%** | **0.21% 球即容纳全 4 根** |
| 脉冲 `rx_x=0.0`, **$y$-轴平移** | 4/4 / $7.4850$ | $0.026860$ / $0.008536$ | **$[0.000004, 0.000004]$** | **0.014%** | **0.04%** | `AA`–`BB` 跨口差严格为 **0** |

**铁证**：仅把阵列平移方向从 $x$ 旋转 $90^\circ$ 到 $y$，跨口总信号反而**增强 2.25 倍**（`0.00355 → 0.00800`），而连续圆柱 `M=4` 的全端口覆盖半径却**暴跌 133 倍**（从 `27.75%` 坍缩至 **`0.202%`**），`min_pair` 暴跌 54 倍（从 `7.49%` 至 **`0.134%`**）。这证明 [full_ports/README.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/full_ports/README.md) 中“偏置两槽 30$\lambda$ 在 10% 噪声下排除全四根、历史六槽 40,000$\lambda$ 落入 10% 球”的对比，**实质上混杂了（i）阵列平移轴相对局部奇异方向的方位角、（ii）跨口/本地幅度比未白化、（iii）近折叠种子条件数三者，根本无法据此归因于单一物理机制**。

---

## 4. 最省成本的证伪实验与更有价值的原创命题

### 4.1 最省成本、可能直接推翻当前解释的实验（改动 <20 行代码，单核 <2 秒 CPU）

利用现有 [cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py) 即可在 2 秒内完成两项低成本决断检验，直接检验当前两条主线是否成立：

1. **检验 1（击穿历史 6 槽 Fold 的连续 PDE 合法性）**：将 [historical_inputs.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/historical_inputs.json) 的 36 盘坐标直接输入 `Cylinders(order=0,1,2,4)`，从 `seed_A/B` 出发解增广转折点方程 $\Gamma_F(z, a, \lambda)=0$（[06_CLAIM_LEDGER_AND_NEXT_STEPS.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md)）并延拓圆盘半径 $a \to 0.04\lambda$。
   - **推翻判据**：若在连续圆柱 PDE 中，该历史 Fold 要么移出 $[0.5,30]^6$ 合法盒，要么仅通过依赖单盘 Mie 共振（$\chi\approx 18.36$ 附近）重建，则此前 `z_alg_host_fold` + `multiplicity` 对“非共振真实连续介质多根”的叙事须彻底降格为“高对比离散脉冲近似的反例”。
2. **检验 2（同构连续圆柱下的方位角扫描与白化噪声合同）**：在 `Cylinders(M=8, L=2,3, R=30)` 中，将组件平移向量参数化为 $R(\cos\phi,\sin\phi)$ 扫描 $\phi\in[0,\pi/2]$，同时报告三种口径的覆盖半径：（a）未白化散射阵 $\|Y_*\|_F$、（b）**逐通道白化跨口比** $R_{\text{cov}}/\|Y_{\text{cross}}\|_F$、（c）含直达场标定误差 $\|E_{\text{inc}}\|_F$；并加入**跨模型盲反演**（用 `M=16` 生成数据、用 `Pulse` 反演）。
   - **推翻判据**：若在 $\phi=\pi/2$ 处 $L=3, R=30\lambda$ 仍出现全部 8 根落入 $<1\%$ 白化跨口噪声球，同时跨模型盲反演误差（~9%）超过 $\phi=\pi/2$ 的跨口分离 100 倍以上，则说明“跨端口自动消歧”与“脉冲-连续等价”两条当前解释均被证伪。

### 4.2 更有价值且可证伪的原创命题（及与最接近先例的界限）

把目前割裂的两条线——Pro 原文 [01_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md) 的连续小区域 $O(\varepsilon^2)$ 空间角 vs $O(\varepsilon^4)$ 频率阶数，与 [MATHEMATICAL_THEORY.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/MATHEMATICAL_THEORY.md#L86-L143) 的多组件全端口多根——合流为一条**统一、物理口径干净（按逐通道信噪比白化、避开单盘 Mie 共振）、且直接可证伪的定量命题**：

> **候选原创命题（有限孔径多组件连续 TM 散射中的「跨口几何相变」与白化不可辨界）**
> 在二维无损连续标量 TM（精确出射 Hankel 核、均匀连续圆柱或有界 Lipschitz 夹杂，且满足非共振电尺寸上界 $k\sqrt{1+\chi_{\max}}\,a \le c_0 < 1$ 以排除单盘 Mie 共振）中，将 $L$ 个具有两个内部正则等局部数据分支的组件按方向 $\hat{e}_\phi=(\cos\phi,\sin\phi)$、间距 $R$ 排布，并记录**全部 $L\times L$ 复散射矩阵**，但以**逐通道白化相对跨口信噪比** $\eta_{\text{cross}}:=\frac{\|F_{\text{cross}}(m_\alpha)-F_{\text{cross}}(m_\beta)\|_F}{\|Y_{\text{cross}}\|_F}$ 度量可辨性（从而消除 $R^{-1/2}$ 几何衰减带来的未白化假象）：
> 1. **非共轴相（通用角 $\hat{e}_\phi \cdot (r_2-r_1) \notin \frac{\pi}{k}\mathbb{Z}$）**：随着 $R\to\infty$，任意一对异分支的**白化跨口分离 $\eta_{\text{cross}}$ 收敛于严格正的 $\Theta(1)$ 常数**（不随 $R$ 衰减！未白化覆盖半径为 $\Theta(R^{-1/2})$）。此时“增大 $R$ 使全端口落入固定白化噪声球”**不可能成立**；
> 2. **共轴退化相（$\hat{e}_\phi \perp (r_2-r_1)$，即阵列沿单组件一阶盲方向排列）**：跨口单次散射首阶路径差为零（或受破缺偏移 $\delta_{\text{rx}}$ 控制），导致**白化跨口分离 $\eta_{\text{cross}}$ 从 $\Theta(1)$ 相变为 $O(\delta_{\text{rx}}/R + R^{-1/2})$**（未白化从 $\Theta(R^{-1/2})$ 相变为 $O(R^{-1})$）。由此**无需超天文间距 $R\sim 10^4\lambda$，在中等物理间距（如 $R=30\lambda$）的连续 PDE 完整跨端口测量中，即可在 $<1\%$ 白化噪声与有限发射功率下稳定容纳全部 $2^L$ 个相异正材料解释**；
> 3. **可证伪预测**：在 $R\in[10\lambda, 300\lambda]$ 的连续圆柱扫频/扫角实验中，将对数斜率 $\frac{d\log \eta_{\text{cross}}}{d\log R}$ 对 $\phi$ 作图，必须在 $\phi=\pi/2$ 处测到斜率从 $0$ 跳变至 $-1$（当 $\delta_{\text{rx}}>0$ 且主导时）或 $-1/2$（当 $\delta_{\text{rx}}=0$ 的互散射主导时），且白化跨口覆盖半径比 $R_{\text{cov}}(\phi=0)/R_{\text{cov}}(\phi=\pi/2)$ 随 $R$ 线性或根号发散。

**与最接近先例的精确界限**：

| 最接近先例 | 先例已覆盖的范围 | 本命题的严格边界与新增可证伪贡献 |
|---|---|---|
| **Griesmaier–Schätzle (2024)**（[PRIOR_ART.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/PRIOR_ART.md#L18)） | 给出二维 Helmholtz 两组散射路径展开与 $O(R^{-1})$ 远场算子分离稳定性常数，前提是**全入射×全出射远场算子**与全局 Born 条件。 | 不含“有限单源/单收局部多根拼接”，更未给出全跨端口白化信噪比随平移方位角 $\phi$ 的 $\Theta(1)\to O(R^{-1})$ 几何相变与 $2^L$ 实正材料共球界。 |
| **Everett (1996, JCP; 2026, GJI)** 与 **Vasco (2000, GJI)**（[05_LITERATURE_AND_PRIOR_ART.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md)） | 离散有限差分 Helmholtz/扩散系统的多项式化、Gröbner 商环代数根数上界（含复/负根）、同伦寻多根。 | 限于小规模离散代数系统，无开域连续 PDE 出射 Hankel 跨组件耦合、无实正材料盒内 $2^L$ 正则根下界与白化全端口噪声几何相变。 |
| **Kuo–Novikov (2025, Inverse Problems)** 与 **Bonnet-Ben Dhia–Chesnel–Nazarov (2015)** | 固定入射方向点势节点隐身构造（改变或增加支撑/退化系数方向）及有限方向小扰动不可见性。 | 依赖总场节点隐身或微扰隐身；本命题限定**固定支撑、无节点、非单盘 Mie 共振的实正介电圆柱盒**，并同时刻画“何时全跨端口仍不可分”与“最小角度破缺恢复率”。 |
| **Beinert–Plonka (2015)** 与 **Novikov–Xu (2026)** | 强分离支撑的 $2^{N-2}$ 指数非唯一性（二选一组合）。 | 数据为**无相位 Fourier 强度**、线性无耦合正演；本命题针对**完整复数相干全波跨端口矩阵**，必须处理全阶非零 Hankel 多次散射耦合。 |

---

## 5. 难以核查的缺口清单

1. **192 位 Arb 认证源码的现场复跑缺口**：因环境缺少 `python-flint`，我们复算了 [krawczyk.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/results/krawczyk.json) 与 [conditioning/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/conditioning/results.json) 输出的端口区间到 [full_ports/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/full_ports/results.json) 的全链路（一致至 $10^{-18}$–$10^{-21}$），并用独立双精度正演核对了脉冲与连续圆柱解，但**未能现场重跑** [certify_separated.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/certify_separated.py) 与 [interval_physics.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/cert/interval_physics.py) 内部的 192 位 Krawczyk 收缩步本身。
2. **连续圆柱的无限角阶截断尾界缺口**：[cylinders.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/continuum/cylinders.py#L215) 明确仅做谱截断（$M\le 16$）与 48 采样点边界跳差数值核对，缺乏解析尾界与区间认证，因而“连续圆柱精确 PDE 存在 $2^L$ 根”目前仅在极小对比度极限有存在性草稿（[MATHEMATICAL_THEORY.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/MATHEMATICAL_THEORY.md#L127-L143)），在有限大对比度 $[0.5,30]$ 尚无严格证明。
3. **高维真值文件与上游历史工作树缺失**：[OFFLINE_SCORE.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/OFFLINE_SCORE.json) 与 [separated.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/separated.py#L21) 引用的外部路径（`research/delegated/...`、`sealed_truth/...`、`srap_cie_array_study_20261121/...`）不在当前隔离工作区内，无法独立核验细网格（96²/150²）真值生成、R8 噪声方向向量 $u$ 的原始采样与历史 36 盘宿主原证书 `point_fibers.json`。
4. **Pro 原文显式常数与文献全文缺口**：[01_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md) 第 4 节的 $1/8, 1/16$ 常数及第 6 节统一频率余项 $O_{C_x^2,I}(\varepsilon^6(1+|\log\varepsilon|)^2)$ 仅有解析草稿、无对应代码或机器计算核对；且 [PRIOR_ART.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/multiplicity/PRIOR_ART.md#L38) 注明 Everett (1996)、Dmitriev–Rumyantseva (2021)、DeFacio–Kim (1992) 与 Everett (2026) 部分正文未能完整获取全文，文献排重仍留盲区。
