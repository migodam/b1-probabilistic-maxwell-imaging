原回复。文中的本机文件链接已换为公开研究资料或附件来源说明。

下面基于隔离工作区中 [RESEARCH_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/RESEARCH_REPORT_ZH.md)、[MATHEMATICAL_THEORY.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/MATHEMATICAL_THEORY.md)、[box_modulus.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/bounds/box_modulus.py)、[task_caps.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/bounds/task_caps.py)、[array_secant.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/array/array_secant.py)、[replay.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/replay.py) 及两份原研究 [01_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md) / [01_EXECUTIVE_VERDICT.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md)，并在不改动工程的前提下，直接对冻结数据 [physical_arrays.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/array/results/physical_arrays.npz)、[public_vie_inputs.npz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/public_vie_inputs.npz) 和 [results/](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/OFFLINE_SCORE.json) 做了一组独立可复算的验证计算（正演重放相对误差均在 $3.1\times10^{-15}\text{–}8.5\times10^{-15}$），给出方向纠正与两个实质不同的物理候选。

---

## 一、用自检证据纠正三个可能走偏的方向

### 1. 现稿主线滑向了“扩张孔径的负向计数”，而正向恢复被困死在材料空间的 Neumann/细分墙里
- **已核对事实**：核对 [bounds/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/bounds/results.json) 全部 72 个条件：**72/72 的整盒界给出了 0% 收紧**（`best_diameter_upper == prior_task_diameter`，其中高对比 36 组里有 4 组连绝对 Green 优控级数 $|H|D_r$ 都不收缩）。为了在 8 槽低对比（$\pm 0.75$）挤出 **1%–5%** 的收紧，[refine_caps.py](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/bounds/refine_caps.py) 不得不对 6 个外围坐标做指数级盒子细分（较厚层扫满 4095 个盒仍有 866 个未决叶），且 [gain_pressure.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/bounds/gain_pressure.json) 显示 0.77% 共享增益漂移就能抹平它。
- **为什么走偏**：两份原研究里最有物理穿透力的抓手是**光学定理 / 能量守恒**（[01_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md)）。但进入 `study/` 的中尺度与高维实验后，方法退回到了把散射矩阵 $Y\in\mathbb C^{N_r\times N_s}$ 拍扁成无结构实向量 `pack(y)`、在高度非凸的材料坐标 $m$ 上做单点 Jacobian / Neumann 级数展开和空间二分。这必然在高对比、高维（576/1024 像素）撞上指数墙；于是稿件只能转去写把组件拉远到 $30\lambda\text{–}200{,}000\lambda$ 的 $2^L$ 根计数（`multiplicity/`），离“固定孔径下高维恢复/选接收”越来越远。

### 2. 一个必须提前踩刹车的高维直觉误区：把介观光学的“相干波前聚焦（GWS）”直接搬到高维材料反演，发射端收益几乎被封死，单模聚焦反而会变差 2–3 倍
- **已核对事实**：在相同发射孔径、相同总发射能量 $\operatorname{Tr}(C_w)=P$ 和相同每通道噪声 $\sigma$ 下，我们对全波完整 Fisher 矩阵 $\mathcal I(C_w)=\sigma^{-2}\operatorname{Re}((P^\dagger P)\odot(EC_wE^\dagger)^T)$ 用 Frank–Wolfe 求解了全局最优发射协方差 $C_w^*\succeq 0$（以及忽略外围像素的单目标聚焦波前 `naiveGWS`），对比逐源均匀发射 $C_w^{\rm unif}=(P/N_s)I$：

| 场景（相同孔径、相同总发射能量 $P$） | 均匀逐源 $\mathrm{std}$ | 单目标聚焦 `naiveGWS` | 全局最优 $C_w^*$ | 最优/均匀增益 |
|---|---:|---:|---:|---:|
| `disks08_low`（$N=8 < 2N_r=12$，`rich`） | 0.1516 | —（rank-1: 0.1012） | **0.0932** | **1.63×**（95% 能量进 1 模） |
| `r10_twod_0001`（$N=576$，12×48，盒先验） | 0.2174 | **0.4658（差 2.14×）** | 0.2145 | **1.01×** |
| `austria_eps4_400MHz`（$N=1024$，20×40，盒先验） | 0.0537 | **0.1578（差 2.94×）** | 0.0509 | **1.06×** |
| `s2d_007__pw180_10x20`（$N=1024$，10×20） | 0.0753 | **0.1531（差 2.03×）** | 0.0724 | **1.04×** |
| `s2d_007__pw360_10x20`（$N=1024$，10×20） | 0.0483 | **0.0997（差 2.06×）** | 0.0468 | **1.03×** |

- **物理解释与纠正**：当自由度 $N \gg 2N_r$ 时，把全部发射能量集中到 1 个聚焦波前会令照射矩阵 $E C_w E^\dagger$ 掉秩到 1，外围 500–900 个未知干扰像素立刻能完美抵消 ROI 信号；即使全局最优 $C_w^*$ 也不得不把能量铺开到几乎所有发射角（如 Austria 前 6 个特征值为 `[3.18, 2.65, 2.13, 1.84, 1.65, 1.46]`，总和 20），因此**在高维未知介质里，固定孔径下仅靠发射端波前整形几乎没有大收益（1.01×–1.08×）**。真正有大收益的两个抓手是：**（1）切换表示——在电流空间用无源守恒律凸对偶做全局非线性恢复/认证；（2）在发射能量不变前提下，优化接收端空间采样布局**。

### 3. 把 $Y$ 拍扁丢掉了高维强散射下极度稳健的算子守恒量
- **已核对事实**：在 [austria_eps4_400MHz](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/highdim/results/austria_eps4_400MHz.json)（1024 像素、10% 复噪声）中，4 个都落在原数据球内的不同历史求解器终点（`cie_dual`, `cie_single`, `global_dual`, `paper_nie`）的锐利 ROI 均值 `sharp_ROI` 剧烈波动 **44%**（`1.8554` 到 `2.6693`，跨度 `0.814`），甚至线性化的 Backus–Gilbert 核均值 `q_BG^T m` 在强非线性分支 `global_dual` 上仍漂移 **18%**（`1.9557` vs `2.2998`）；但**精确离散光学定理量** $\operatorname{Im}\operatorname{Tr}(E_0^\dagger J)=\operatorname{Tr}(J^\dagger(\operatorname{Im}G)J)$ 在这 4 个终点上被死死锁在 **`2008.02–2010.50`（仅 0.12% 展宽！）**（`r10_twod_0001` 的 `r12_dr` 与 `r12_g19` 上更是 `21.740` vs `21.741`，仅 **0.005%**，且其 `q_BG^T m` 跨度比 `sharp_ROI` 缩小了 **12.6×**，从 `0.1695` 降到 `0.0134`）。同时，Austria 的 $\operatorname{Tr}\operatorname{Re}(J^\dagger GJ)\approx 14{,}968$ 是 Born 项 $\operatorname{Tr}\operatorname{Re}(E_0^\dagger J)\approx -7.38$ 的 **2000 倍**，说明任何依赖弱散射/Neumann 级数的界都注定失灵，而建立在 $J$ 的二次守恒律上的方法才匹配真实物理。

---

## 二、两个实质不同的相邻学科物理候选

### 候选一（更值得做，主推）：纳米光子学 / 电磁基本极限中的“电流空间无源与局域守恒律凸对偶（Passivity / Local Conservation-Law Dual）”

**物理思想来源**：纳米光子学与计算电磁学中的基本散射极限理论——把非凸的材料参数 $m\in\prod_i[a_i,b_i]$ 消去，将全波 Lippmann–Schwinger 方程 $J_s = D_m(E_{0,s}+GJ_s)$ 转化为**诱导极化电流 $J$ 空间中的精确二次能量/无功功率守恒律**（局域与全局光学定理 + 对比度盒无源不等式），用拉格朗日对偶（单个 Hermitian 正定性 / Cholesky 检验）绕开 $2^N$ 盒子细分与 Neumann 级数收缩条件，直接给出全多次散射、高维、全局非线性有效的证书。

#### 1. 可算的核心量
记 $E_s = E_{0,s}+GJ_s$（关于 $J_s$ 仿射），盒中心 $\alpha_i=(a_i+b_i)/2$，半宽 $\beta_i=(b_i-a_i)/2$。对任意实正材料 $m_i\in[a_i,b_i]$，在每个像素 $i$ 上、跨全部 $N_s$ 次照射（记 $J_i,E_i\in\mathbb C^{N_s}$）严格成立两条无近似二次守恒律：
1. **跨源实材料/局域光学定理等式**（强制全部 $N_s$ 次照射看到**同一个实数** $m_i$）：对任意 Hermitian 乘子 $N_i=N_i^\dagger\in\mathbb C^{N_s\times N_s}$，
   $$\operatorname{Im}\!\big(J_i^\dagger N_i E_i\big)=0.$$
2. **跨源对比度盒无源不等式**：对任意半正定乘子 $M_i\succeq 0\in\mathbb C^{N_s\times N_s}$，
   $$(J_i-\alpha_i E_i)^\dagger M_i (J_i-\alpha_i E_i) - \beta_i^2 E_i^\dagger M_i E_i \le 0.$$

将它们与数据拟合 $\|SJ-Y\|_F^2$ 合成关于 $J$ 的精确二次拉格朗日函数并求解析下确界，得到**光滑凹对偶下界**：
$$\min_{m\in\prod_i[a_i,b_i]} \|S J(m)-Y\|_F^2 \;\ge\; \underline d^2(a,b;S,Y) := \sup_{\substack{M_i\succeq 0,\; N_i=N_i^\dagger\\ H(M,N)\succ 0}} \Big\{\|Y\|_F^2 + \langle E_0, D_{ab M}E_0\rangle - r(M,N)^\dagger H(M,N)^{-1} r(M,N)\Big\}.$$
- 只要找到一组可行乘子满足 $H(M,N)\succ 0$ 且 $\underline d(a,b;S,Y)>\delta$，即**一次 Cholesky 分解完成对整个 $N$ 维非线性假设盒 $\prod[a_i,b_i]$ 的全局排除**——无需初值、无需线性化、无需 $\rho(|H|D_r)<1$。
- **可直接恢复的高维材料量**：
  - 全局峰值对比度严格下界 $\max_{i\in\Omega} m_i > b_{\rm lb}$（令所有像素 $m_i\in[0,b]$ 二分求解）；
  - 区域存在性/对比度下界 $\max_{i\in\mathrm{ROI}} m_i > b_{\rm ROI}$（在给定外围背景杂波上界 $m_{\rm out}\in[0,b_{\rm out}]$ 下排除空 ROI 或低对比 ROI）；
  - 区域极化储能 $\sum_{i\in\mathrm{ROI}} m_i|E_i|^2 = \sum_{i\in\mathrm{ROI}}\operatorname{Re}(J_i^\dagger E_i)$（本身就是 $J$ 的二次型，可直接进目标函数）；
  - **定量预报加哪种接收**：$H(M,N)$ 对 $S^\dagger S$ 是显式仿射的，改变/新增接收集合直接给出非线性全局下界 $\underline d$ 的定量跃升。

#### 2. 依赖的数据
仅需冻结的测量矩阵 $Y$、确定性噪声半径 $\delta$、已知背景线性算子 $(G,E_0,S)$（或已知宿主修饰后的 $(G_{\rm host},E_{\rm host},P_{\rm host})$）和待检材料盒 $[a_i,b_i]$。**不依赖已有重建中心 $m_0$、不依赖真值、不依赖局部 Jacobian。**

#### 3. 已核对的正例（直接复算自工作区数据）
- **576 维高维全波（`r10_twod_0001`，12 tx × 48 rx）**：
  - 仅用 **2 个全局标量乘子** $(\mu_0,\nu_0)$（全 12 个源共用一次 $576\times 576$ Cholesky，约 15 ms，**零线性化、零细分**），直接全局证出 **$\max_{i\in\{1..576\}} m_i > 2.7095$**（在 $b=2.7095$ 处 $\underline d=0.0106>\delta=0.0105$；离线真值 ROI 均值为 `3.447`）。
  - 仅用 **3 个标量乘子** $(\mu_{\rm ROI},\mu_{\rm out},\nu_0)$ 对中央 64 格 ROI 做全局排除，允许全部 **512 个外围像素在 $[0, b_{\rm out}]$ 任意独立变化**：
    - `outer <= 3.0` 时：`ROI <= 0.0` 得 `LB/rad = 7.75`，`ROI <= 0.5` 得 `LB/rad = 3.72`（**全局排除空 ROI 并认证 $\max_{\rm ROI} m > 0.5$**）；
    - `outer <= 2.5` 时：`ROI <= 2.0` 得 `LB/rad = 9.58`（**全局认证 $\max_{\rm ROI} m > 2.0$**）。
- **固定发射孔径与发射预算下，接收端布局的巨大定量收益**：
  - `r10_twod_0001`（固定全部 12 个发射源和能量），从 1/4 环（12 rx）出发：
    - 1/4 环（12 rx）：全球非线性下界 $\max m > 1.333$；
    - +4 个相邻接收（16 rx）：$\max m > 1.444$（仅增 `+0.111`）；
    - **+4 个 Schur 优选接收（16 rx，零额外发射能量）**：**$\max m > 1.914$（增 `+0.581`，是相邻增设的 5.2 倍！）**；
    - 全环（48 rx）：$\max m > 2.708$。
  - `austria_eps4_400MHz`（1024 像素，固定 20 tx）：1/4 环（10 rx）`0.162` → +4 相邻 `0.191`（`+0.029`）→ **+4 Schur 优选 `0.214`（`+0.052`，1.79 倍）** → 全环 `0.375`。
  - `s2d_007`（1024 像素，同为 10 tx × 20 rx）：半环 `pw180` 为 `0.3465`，全环 `pw360` 提升至 **`0.4633`（+34%）**。
- **中尺度已知宿主（`disks08/12/16`，低/高对比全部 6 组）**：
  - 即使单源解耦对偶，也在 **6/6 组**（含 `box_modulus.py` 失效的全部高对比组）上零细分全局排除 `ROI = 0.5`（外围 $\pm 0.75$），`LB/delta = 1.01–3.74`。
  - 在 `disks08_low` 最难的 `pos_face_0.75`（[faces/results.json](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/faces/results.json) 实际最小全波残差仅 `1.0636 δ`）上，加入跨源耦合 $(M_i,N_i)$（$32\times 32$）后，**零盒子细分**就把 `LB/delta` 从 `0.4764` 推到 **`0.7333`**。

#### 4. 会失败的情形（已核对的物理边界 + 猜想）
1. **严格无损介质的准静态/暗共振天花板（已核对）**：我们直接算出了两个网格的实部 Green 最大特征值：
   - 24² 网格 $G_0$：$\lambda_{\max}(\operatorname{Re} G_0)=0.6282 \Rightarrow b_{\rm crit}=1/\lambda_{\max}=1.5919$（$\operatorname{Im}G_0$ 仅 17 个特征值 $>10^{-6}$）；
   - 32² 网格 $G_1$（电尺寸更大，2 m @ 400 MHz）：$\lambda_{\max}(\operatorname{Re} G_1)=2.9264 \Rightarrow b_{\rm crit}=0.3417$（$\operatorname{Im}G_1$ 仅 35 个特征值 $>10^{-6}$）。
   - 在严格无损（$\operatorname{Im}\chi=0$）时，一旦 $b_{\max}$ 明显超过 $b_{\rm crit}$，高空间频率的**非辐射（倏逝）暗电流**同时满足 $SJ\approx 0$ 与 $(\operatorname{Im}G)J\approx 0$，导致 $I-b_{\max}\operatorname{Re}G$ 不定，简单 2 标量对偶随之归零（`r10` 在 $2.71$ 以上，`austria` 在 $0.38$ 以上）。**【猜想，待检验】**：引入实际物理中微小的无源损耗 $\operatorname{Im}\chi\ge\epsilon''>0$、启用逐像素/跨源 $(M_i,N_i)$ 耦合、或叠加平滑度二次约束，应能大幅抬高该天花板，但在严格无损大电尺寸高对比下能抬高多少尚未验证。
2. **线性区域均值 $q^Tm$ 不是 $J$ 的原生二次型**：此方法原生最紧地界定“整盒排除 / 峰值对比度 / 极化储能”，若要直接紧界线性均值 $q^Tm$，还需结合少量分层或辅助变量（不如整盒排除干净）。
3. **跨网格模型误差与未知增益**：若粗网格参考本身超出噪声球（如 `austria` 150²→32² 离散误差 `22.5%` > 原 `10%` 噪声），必须在 $\delta$ 中显式计入建模误差或将增益并入对偶，否则证书仅对当前离散算子成立。

#### 5. 最接近的原研究
- **外部文献**：Kuang & Miller (*PRL* 2020, "Computational bounds to light–matter interactions via local conservation laws")；Angeris, Vučković & Boyd (*ACS Photonics* 2019 / *PRL* 2021)；Molesky, Chao & Rodriguez (*PRX* 2020)；Gustafsson et al. (*NJP* 2020)。**本质区别**：这些工作全部面向**正向光子器件设计**（在无测量数据下最大化某响应的物理上限），几乎未用于**逆散射数据一致性** $\|SJ-Y\|_F\le\delta$ 下的高维材料下界认证，也未引入多源逆散射特有的**跨源 Hermitian 同材料守恒律** $J_i^\dagger E_i = E_i^\dagger J_i$。
- **本项目原研究**：最直接承接 [01_FULL_REPORT_ZH.md](https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/z-alg-secants-20261122/research/z_alg_secant_tasks/SOURCES_AND_SCOPE.md) 的光学定理主线，但将其从“自由空间、两个极小 Born 圆盘、单前向通道”一次性推进到了“任意高维、强多次散射、全矩阵有限端口”。

#### 6. 最小的配对实验
在冻结的 `r10_twod_0001`（576 px）、`austria_eps4_400MHz`（1024 px）和 `disks08_low/high`
