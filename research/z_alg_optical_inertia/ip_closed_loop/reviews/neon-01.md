以下意见基于对隔离工作区 `closed-loop-26/first/` 及其关联背景文件的逐行代码审读与独立数值复算（未修改工作区任何文件）。下文严格区分**已核对结果**与**猜想/方向建议**。

---

## 一、核心诊断：当前闭环研究卡在哪里？

先给结论：**目前这套 `closed-loop-26/first/` 的数学骨架大体正确，但实验设计把真正有物理意义的新现象“稀释”成了已知经典算法的重演。**

已核对的三个关键事实：
1. **配对光束优化被强基线几乎完全覆盖（`pair/main_32` 与 `beams/`）**：
   - 在 [proposal_arrays.npz](../pair/main_32/proposal_arrays.npz) 和 [proposal.json](../pair/main_32/proposal.json) 的全部 4 个算例（Case 0–3）上，我用纯 Python 复算了 `known_pair_coherent`（两中心差矩阵 $F_1-F_0$ 的最大右奇异向量）与 `positive_quadratic`（双二次型对偶最优光束）的夹角及精确对偶下界 $\sqrt{a^* H a}/r_{\text{new}}$：

| Case | 夹角 $\angle(a_{\text{coh}}, a_{\text{pos}})$ | 中心上界 $/r$ | `known_pair_coherent` 在正定型下的 $d/r$ | `positive_quadratic` 的 $d/r$ | `indefinite_quadratic` 的 $d/r$ | 正定型椭球半宽 $/r$ |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | $0.346^\circ$ ($|\cos|=0.999982$) | $2.4000$ | `2.382495` | `2.382512` (**+0.0007%**) | `2.328633`（更差） | $(0.0082, 0.0093)$ |
| 1 | $0.893^\circ$ ($|\cos|=0.999879$) | $2.4000$ | `2.375888` | `2.375987` (**+0.0042%**) | `2.253394`（更差） | $(0.0124, 0.0118)$ |
| 2 | $0.382^\circ$ ($|\cos|=0.999978$) | $2.4000$ | `2.379451` | `2.379461` (**+0.0004%**) | `2.304855`（更差） | $(0.0107, 0.0098)$ |
| 3 | $0.050^\circ$ ($|\cos|=0.9999996$) | $2.4000$ | `2.386722` | `2.386726` (**+0.0002%**) | `2.231641`（更差） | $(0.0062, 0.0071)$ |

   - 在旧结果控制的 20 个 `main_32` 算例（[beams/main_32/results.json](../beams/main_32/results.json)）中，两类盒子的收缩参数 $\eta_{\text{upper}}\in[0.00285, 0.0116]$（中位数 $0.00402$），盒子宽度小到材料几乎退化为两个已知点；`source_shaped` 与 `direct_coherent` 的前瞻全盒分离间隙在 **20/20 算例上差值全为 `0.00000`**（精确一致），光学排除 `source_optical_exclusion` 为 **0/20**。
   - 在 12 个 `confirm_32` 算例（[beams/confirm_32/results.json](../beams/confirm_32/results.json)）中，$\eta_{\text{upper}}\in[0.0127, 0.0611]$，`direct_coherent` 与 `source_shaped` 事后排除数完全相同（均为 6/12），前瞻分离间隙全为负（0/12），二者最大差值仅 $+0.020r$（行 300：`-0.645` vs `-0.665`），光学排除仍为 **0/12**。
2. **填充 FFT 扇区地图有效，但尚未触及不定状态物理（`maps/`）**：
   - 在 [sector_32_solid_primary.json](../maps/results/sector_32_solid_primary.json) 中，当弱背景 $b=0.03125$ 时，旧的 Neumann/强制性方法（[map_32_solid_primary.json](../maps/results/map_32_solid_primary.json)）是 **0/1024** 排除，而 [sector_map.py](../maps/sector_map.py) 排除了 **1017–1021 / 1024** 个中心——这是本批文件里**最实在的算法进展**。
   - 但所有 8 组场景中 `negative_state_excluded` 全部为 **`0`**（~8160 个被排除中心里没有一个用到负惯性计数 $n>0$），原因在于 `[3.75, 4.25]` 区间太窄（相对半宽仅 6.25%），局部移位算子几乎正定。
3. **实测实验室表（`lab/`）既未测本方法，反倒暴露了经典 SVD 相干波束在跨孔径外推时的脆弱性**：
   - [lab/measured.py](../lab/measured.py) 只对两个**已知实测目标**的半圈接收差矩阵 `delta[train]` 做 SVD 取第一右奇异向量，与材料类、二次型外包络、噪声球和物理模型完全无关。
   - 在 [lab/results.json](../lab/results.json) 和 [lab/held_ports.csv](../lab/held_ports.csv) 中，Fresnel 仅用 2 个对向发射口（`0°, 180°`），中位增益仅 `1.088`；而 Manitoba（12 发、6 训 6 测）在 3.5 in 校准下仅赢 `45/84`（53.6%，`empty_last` 为 `42/84=50.0%`，5.0 GHz 更是 **0/12 全败**，中位比 `0.715`），换到留出 2 in 钢柱校准后跌至 **`32/84`（38.1%，中位比 `0.865–0.877`，劣于单口发射）**。

---

## 二、逐项数学与实现核查（已核对结果）

### 1. 双二次型距离对偶的精确条件（[paired.md](../theory/paired.md) & [distance.py](../pair/distance.py)）
- **弱对偶与复 SDP 秩一精确性**：核对成立。固定 $a$ 时，最小化 $\|y_1-y_2\|^2$ 满足两个复 Hermitian 二次约束 $[a;y_i]^*Q_i[a;y_i]\le 0$，提升到 $X=[y_1;y_2;1][y_1;y_2;1]^*\succeq 0$ 后仅有至多 3 条实仿射约束。在**复** Hermitian 空间中极点秩 $r$ 满足 $r^2\le 3\Rightarrow r=1$（注：实对称空间中 $r(r+1)/2\le 3$ 允许 $r=2$，故此精确性依赖复数域结构）。存在 $\lambda_1,\lambda_2>0$ 使 $D=\begin{bmatrix}I+\lambda_1 R_1&-I\\-I&I+\lambda_2 R_2\end{bmatrix}\succ 0$ 确实保证提升后目标子水平集有界，从而由 Slater 条件 + 强对偶 + 极点存在性得到精确性。由于双上确界 $\sup_a \sup_\lambda = \sup_\lambda \sup_a$ 无条件成立，换序取极值也成立。
- **容易被忽略的四个隐含前提与实现漏洞**：
  1. **严格可行性（Slater）与空纤维陷阱**：定理要求对固定 $a$，两类纤维 $\{y_i:[a;y_i]^*Q_i[a;y_i]\le 0\}$ 均严格非空。若旧数据噪声半径估计过小、或某一类已被旧数据排除导致某方向上纤维为空，对偶值 $a^*Ha$ 会随 $\lambda_i\to\infty$ 发散，[distance.py](../pair/distance.py) 的 `optimize()` 从未验证两类纤维非空。
  2. **$\lambda$ 参数化截断与非凸局部
