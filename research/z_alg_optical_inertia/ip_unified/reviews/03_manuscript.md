> 独立顾问原意见，内部文件引用已转为文件名。文中的候选结论属于意见发表时的判断；最终取舍及区间复核见 [答复](RESPONSE_ZH.md)。

### 审稿与创意顾问总评（*Inverse Problems* 审稿视角）

**总体决定建议：大修（Major Revision）**。

本文围绕“有限孔径电磁散射复数据能否在重建之前排除整个低对比度介质类”这一高价值反问题，构建了基于二次型惯性交错（rank interlacing）、局部本构行补偿、反向障碍（Proposition 4.1）以及有损接收端椭球对偶（Proposition 3.2）的严谨数学框架。经只读核查 manuscript_ip.tex 及其全部独立分节（model.tex、contrast.tex、physics.tex、experiments.tex、appendix.tex、references.tex），并交叉复核三组自选原始数值文件（u100_best_strict.json、u100_best_strict_source.json、finite_0_directional.json 及其 22 行配套验证文件 finite_0_pure_shift_verified.json），确认 Table 2 的 18 个核心比值及主要向外舍入区间端点数理自洽、治学严谨。

但在以下**四个关键科学问题**上，稿件存在**表述过冲、既有证据未用足、机制诊断未触底以及两大算法主线人为割裂**的决定性缺陷；修正它们不仅决定 *Inverse Problems* 审稿成败，还可直接把本文从“诊断旧松弛的失败”推进到“更强且统一的物理下界理论”。

---

### 四项决定性科学审稿意见（区分【已核对证据】与【推论/猜想】）

#### 1. 中心项（$J_0, E_0, M_Y, H_\eta$）：正文定位与真实 Helmholtz 表现严重倒挂
* **【已核对证据】**
  1. appendix.tex 的二元无源反例经逐项验算完全正确（$K=\diag(-1/2,3/2)$，$H=199/48$，$M_Y=+3/2 \implies \delta^2\ge 127/48$），但它是**非 Helmholtz 代数构造**。在所有真实 Helmholtz 算例上，中心项从未改变任何判定：52 行中心感知 SDP 在 Crescent and plate 上仅得 $0.30104r$，反低于纯秩界 $0.30105r$（experiments.tex）；全区域缺陷界（physics.tex）耗时 159 s 仅将 Table 2 六个目标的 $d/r$ 提升 $0.00006\text{--}0.00016$（相对增幅 $<0.08\%$）。
  2. **可复算根因（直取原始数值）**：
     * **悬崖正谱隙锁死权重**：u100_best_strict.json 显示所复用的权重 `u100_best` 正谱隙仅 $g_+/\mu = 0.001953125 = 2^{-9}$（比负谱隙 $0.03017$ 小 15 倍）。在该算子上扣除均匀上界 $\omega(2\|I-\frac{u}{2}G\|^2 + \frac{u^2}{2}\|G\|^2)$（其中 `green_norm_upper` $=10.5585$），在 $\omega/\mu = 10^{-3}$ 时正谱隙已坍缩至 $2^{-20}=9.54\times 10^{-7}$（u100_best_strict_source.json），彻底锁死了 $\omega$ 的提升空间。
     * **更深层的模式尺度失配**：核查 u100_best_strict_source.json 六个目标的 `full_template_defect_lower`（$m_j$）：巨额本构缺陷集中在**被 $\Tail_{13}$ 丢弃的前 13 个源模态**（$m_1 \in [1.49\times 10^4,\,4.05\times 10^4]$）和**底部深度负模态**（$j=29\text{--}32$，$m_j \in [214,\,1300]$），而在真正求和的**保留正模态**（$j=14\text{--}20$）上仅 $m_j \in [2.9,\,28.1]$。以 Slotted block 为例，即便将 $\omega/\mu$ 强行放大 $1000\times$ 至 $1$ 且令 Schur 折损 $\rho=0$，$\mu\sum_{j=14}^{17} m_j = 1.354\times 10^{-5}\times 23.06 = 3.12\times 10^{-4}$，仍比 $r^2 = 7.83\times 10^{-3}$ 小 25 倍；反过来，前 13 个模态高达 $10^5$ 量级的 `endpoint_energy_upper` 通过式 (4.6) 的 $\rho$ 把保留模态本就微小的增量再腰斩（Slotted block $j=14$：$\omega m_{14}=1.498\times 10^{-7}$，实际增量仅 $7.62\times 10^{-8}$）。
     * **试探中心 $\widehat J$ 的精确代数退化**：将 $\widehat J = J_0 + \Delta J$ 代入 physics.tex 式 (4.4)，可得精确恒等式 $H_\eta = H + \Delta J^*(K - \eta^{-1}K^2)\Delta J$。在 $\eta = 2^{-20}\mu$ 下，任何沿 $|\lambda|>\eta$ 方向的偏移都严格拉低 $H_\eta$；且数值日志中 `stationarity_deduction` 仅 $6.4\text{--}11.4\times 10^{-14}$，表明式 (4.4) 实际上只是数值求解 $K\widehat J \approx L$ 的误差包络工具，放在正文易被误读为可自由优化的算法变量。
* **【推论/猜想】**
  在标量 $K$ 固定、仍丢弃前 $\nu$ 个源模态的前提下，全细胞中心缺陷对真实大域无损排除基本无望；中心项若想在 Helmholtz 问题中发挥实质作用，只能通过**(a) 反向障碍定理中的负特征值吸收（见创意 1）**，或**(b) 放弃悬崖权重并改用非均匀对角权重以保护近零正谱隙**来验证，后者效果目前尚属待验证猜想。

#### 2. “52 像素”表述：摘要过冲、术语冲突、网格依赖，且低估了自身已有数据
* **【已核对证据】**
  1. **摘要过冲与术语混用**：manuscript_ip.tex 摘要称 *"no selection of up to 52 pixel balances can exclude the cap for the tested weights"*，遗漏了 Proposition 4.1 仅针对**定理 3.1 的一致对角非负框架**（contrast.tex 式 (3.3)）这一核心限定（Table 2 标题准确写了 *diagonal local compensation bound*，但摘要删掉了 *diagonal* 与 *uniform frame*）。耦合正定度量仅在 3 个特定嵌套子集（8、26、52 行）上验证，精确 $x$-依赖消元（contrast.tex）也不受此限。更严重的是，model.tex 与 experiments.tex 把 *balances* 明确定义为**改变 $K$ 本身的物理功率/扇区乘子 $D_\mu, D_\lambda$**，而摘要的 *pixel balances* 却指代固定 $K$ 下的本构行框架，极易误导审稿人。
  2. **网格依赖性**：52 像素在 $96^2$ 网格上对应**面积占比 $52/9216 = 0.564\%$**；在网格加密 $h\to 0$ 时，单像素迹按 $O(h^2)$ 衰减，固定像素数 $k$ 并非连续不变量（真正网格无关的是 experiments.tex 的“全 9,216 像素仍 $c\le 8$”，即 100% 面积上限）。
  3. **已有计算被低估（可零成本复算立即升级）**：experiments.tex 仅凭 11 维负子空间上 52 行迹 $<0.95402$ 得出 $c\le 2$，使 Double chevron 的 52 像素上限逼近阈值（$0.97812r$）。**然而 u100_best_strict.json 已直接给出完整 13 维负空间（$d=\nu=13$）的 52 行迹上界 $\beta = 1.75036246 < 2$（`negative_energy_lower` $= 1 - 3.4\times 10^{-14}$）！**直接套用正文自己的 Proposition 4.1 公式 $c \le \nu - d + \lfloor\beta\rfloor$，即得
     $$c \le 13 - 13 + \lfloor 1.75036246 \rfloor = \mathbf{1} \quad (\text{而非 } c \le 2)！$$
     进一步由降序项平均值上界 $b_{(52)} \le S_{52}/52 = 0.0336608$，立得 $S_{59} \le 1.75036246 + 7\times 0.0336608 = 1.98599 < 2\times\text{negative\_energy\_lower}$，即**仅凭现有 JSON 数据，已严格证明任意至多 59 个像素均满足 $c \le 1$**。仅 $m=11\to 12$ 一项就将 Double chevron 的 $\Tail_m(H)$ 从 $0.0020737$ 降至 $0.0011942$（下降 42%），Crescent and plate 从 $0.0027183$ 降至 $0.0015876$（下降 42%）。
  4. **Proposition 4.1 符号中心上限的结构性松弛**：physics.tex 用 $\tr((A+P)_+)\le \tr(A_+)+\tr(P)$ 彻底扔掉了 $Q^*HQ$ 的负特征值（Double chevron 上 $\sum_{j=21}^{32}(-h_j)\approx 0.0119 > r^2=0.0109$），而 u100_best_strict_source.json 表明中心缺陷的绝大部分恰集中在这些深度负特征向量（$j=29\text{--}32$）上。
* **【推论/猜想】**
  将 $c\le 1$ 与下文创意 1 的负谱 Schur 吸收引理结合后，符号中心不可能结论预计可远超 52 像素（推测可达百像素量级），具体新阈值需在现有中心矩阵上复核后方能确认。

#### 3. 损耗共标量严格分离：数学严谨但实证裕度过脆，且割裂了全文两条主线
* **【已核对证据】**
  1. **数学合法性确认**：经 Woodbury 恒等式复核，对纯耗散约束 $(J-J_c)^*M(J-J_c)\preceq C_B$ 施加全电流单标量惩罚 $\tau$ 所得源矩阵严格等于 $H(\tau)=\tau[Z_Y^*(R_S+\tau I)^{-1}Z_Y-C_B]$；因此 Proposition 3.2 的 $\sup_{\tau\ge 0}\tr[H(\tau)]_+$ 确实严格涵盖纯耗散下所有共标量惩罚与任意源子空间投影器。
  2. **主推算例过于脆弱，而同一份数据中藏着 30 倍更稳健的分离**：
     * 摘要与 experiments.tex 主推的严格分离出现在小网格 $12\times 12$（$ka=3.55$）、高损耗角 $t=0.2$、$u=1.0$，且**仅在两个目标中的第一个**成立（目标 2 的共标量本身已达 $0.07308/0.06748 = 1.083r > 1$，finite_0_directional.json）。目标 1 的全 $\tau$ 区间上界为 $0.99873r$（finite_0_pure_shift_verified.json，点最优 $0.99344r$），距阈值仅 **$0.127\%$**，略换噪声实现即可能翻转。
     * 反观同一文件 finite_0_directional.json：在损耗低 4 倍的 **$t=0.05,\,u=0.75$**（目标 1）上，共标量点最优仅 **$0.8025r$**，而 `best_common_shift`（仅 0.85 ms、已区间验证）已达 **$1.01324r > 1$**（接收端椭球 $1.02003r$），分离裕度宽约 30 倍，但作者并未对该算例做 101 区间全 $\tau$ 封口证明。
  3. **一处已核实的计时笔误与偏弱的基线对比**：experiments.tex 称首个正例的 `right_singular` 耗时约 **0.5 ms**、直接迭代投影 **0.319 s**。核查 finite_0_directional.json（$t=0.2,u=1.0$ 目标 1），`receiver_reference` 确为 `0.31854 s`（且因触及 `iteration_limit: 2500` 中止，非收敛最优），但 `right_singular` 实为 **`0.88 ms`**（`0.000878 s`，另加验证 `20.8 ms`）；且将启发式基与未加速一阶迭代投影比耗时，说服力不如对比 $10\times 10$ 对偶矩阵 $\Lambda$ 的标准内点/牛顿法。
  4. **核心方法割裂**：experiments.tex 自己展示，在 $t=0.2,u=1.25$ 时“全电流耗散平移 + 扇区”得 $1.28439r$，几乎是纯接收端椭球（$0.663r$）的 2 倍；而在大域（$1.25\text{ m}$ 与 $96^2$）上纯接收端椭球因丢弃扇区/相位律完全崩塌（$\le 0.151r,\,0.395r$）。然而论文从未将“相干方向独立定价”与“扇区/相位二次型”结合。
* **【推论/猜想】**
  纯耗散投影在大域失效并非“相干独立定价”的固有极限，而是丢弃实部扇区与固定相角等式所致；按创意 2 统一两者后，有望在 $1.25\text{ m}$ 与 $96^2$ 有损域上显著抬高下界并产生大域分离（属待验证算法猜想）。

#### 4. 连续谱与连续算子（$L^2(\Omega)$）适用范围：Weyl 律边界、面积-周长相变与 2D TM 连续迁移
* **【已核对证据】**
  1. **谱类型与 Weyl 常数核验**：在有界域 $\Omega\subset\R^2$ 上，$H_G, R, S^*S$ 均为紧算子，$K/\mu = I + \text{compact}$ **并无连续谱**（本质谱仅单点 $\{\mu\}$），负惯性 $\nu$ 在连续 $L^2(\Omega)$ 上仍有限；physics.tex 式 (4.7) 的主常数 $\frac{k^2|\Omega|}{4\pi}u = \pi(W/\lambda)^2 u$ 推导准确。但 appendix.tex 证明有两处可立即澄清与简化：应指明 $A=k^2(1-\Delta_{\R^2})^{-1}$ 定义在全空间（其与 $H_G$ 的核差主奇性为 $\frac{k^2+1}{8\pi}r^2\log r$，属 $-4$ 阶迹类算子，余项 $N(\varepsilon)=O(\varepsilon^{-1/2})=o(\varepsilon^{-1})$ 无需重型工具），且“平滑域 $\Omega$”条件过强——夹逼证明仅需 $|\partial\Omega|=0$，刚好涵盖正文所有正方形实验域。
  2. **Weyl 律在实测物理区间的适用边界**：
     * 式 (4.7) 仅在**固定 $k$、固定紧扰动 $C$、$u\to\infty$** 下成立。而在全空间傅里叶符号下，$\widehat{H_G}(\xi)=\operatorname{P.V.}\frac{k^2}{|\xi|^2-k^2}$ 在传播盘 $|\xi|<k$ 上为**负**，仅在**倏逝环带 $k < |\xi| < k\sqrt{1+u}$** 上导致 $1-u\widehat{H_G}<0$，且在 $|\xi|=k$ 处存在极点。
     * 在实际大电尺寸反问题（**固定 $u\sim O(1)$、大 $ka$**）中，$k^2(1-\Delta)^{-1}$ 替换不再控制 $|\xi|=k$ 处的极点，且优化后的 $C = \mu^{-1}S^*S + (\lambda/\mu)R$ 随 $u$ 动态变化并稳定至多 $p=\rank(S^*S)$ 个辐射模态。核对本文自身数据：在 $96^2$（$W/\lambda=2.67$）、$u=1$ 处，$\frac{k^2|\Omega|}{4\pi}u \approx \mathbf{22.35}$，而实测 $\nu = \mathbf{13}$（高估 72%）；$u=1\to 1.25$ 实测斜率 $\Delta\nu/\Delta u = 8$ 远低于 Weyl 斜率 $22.35$。
  3. **定量解释全文经验相图的“面积-周长相变”**：内部负模按面积 $\propto (ka)^2 u$ 增长，而可观测相干源自由度按周长 $\sim 2ka$（Bucci–Franceschetti）饱和，故纯秩界存在结构性失效边界 $u_{\rm crit}\sim O((ka)^{-1})$——这与 Table 1 中 $24^2$（$ka=4.27$，可排 $u=5$）、$96^2$（$ka=11.85$，可排 $u=0.84$ 但 $u=1$ 失效）、Fresnel 8 GHz（$ka=23.7$，$u=0.5$ 失效）完全吻合。
  4. **2D TM 连续箱 $L^\infty(\Omega;[0,u])$ 的关键理论补充（已验证推导）**：manuscript_ip.tex 仅将连续介质前向离散误差列为开放需求。事实上，对于本文采用的 **2D TM 标量 Helmholtz 方程**，因 $G:L^2(\Omega)\to L^2(\Omega)$ 紧且不存在均匀化（homogenization）张量偏移，结合 $\mathcal K_u=L^\infty(\Omega;[0,u])$ 的弱-$*$ 紧性与 Rellich 唯一延拓，可直接反证得到一致预解界 $\sup_{\chi\in\mathcal K_u}\|(I-GM_\chi)^{-1}\|_{L^2\to L^2}<\infty$，从而保证**粗糙有界类 $L^\infty(\Omega;[0,u])$ 上存在一致 $O(h)$ 的连续-离散误差界**。但必须明确警示：这一性质**不能**推广至 2D TE 或 3D 矢量 Maxwell（后者存在强均匀化效应，弱-$*$ 极限不保持本构关系）。

---

### 三个具体可执行改动与算法创意（标出精确文件与段落）

#### 改动/创意 1：修正“52 像素”表述、零重算升级至 $c\le 1$（至多 59 像素），并引入 3 行“负谱 Schur 吸收引理”突破障碍天花板
* **涉及文件与段落**：
  * manuscript_ip.tex（Abstract 第 6 句）；
  * physics.tex（§4.2 Proposition 4.1 及其证明）；
  * experiments.tex（§5.2 第 3–4 段及 Table 2）。
* **具体可执行内容**：
  1. **精确措辞与网格无关化**：将摘要中的 *"up to 52 pixel balances"* 改为 *"any selection of up to 52 cells ($0.56\%$ of the search area) under Theorem 3.1's uniform diagonal frame at the tested weights"*，彻底消除与 model.tex 物理功率 *balances* 的术语冲突；并在 §5.2 明确指出“52 像素”仅是 $96^2$ 网格的面积分数，而全区域 $c\le 8$ 才是连续极限不变量。
  2. **零重算立即加强主结论（已核对）**：直接引用 u100_best_strict.json 中已算出的 $d=13$ 迹上界 $\beta=1.75036246<2$，按 Proposition 4.1 将 52 像素界从 $c\le 2$ 立即升级为 **$c\le 1$**（并由排序单调性自然推至 **59 像素 $c\le 1$**），同步压低 Table 2 的 `52-cell limit`。
  3. **核心数学创意（3 行负谱 Schur/Young 吸收引理，已验证证明）**：在 Proposition 4.1 式 (4.3) 中，将 $Q=[Q_+, Q_-]$ 拆分为非负/浅负块 $Q_+$ 与深度负特征值块 $Q_-$（记 $Q^*HQ=\diag(A_+,-A_-)$ 且 $A_-\succ 0$）。对任意 $P = Q^* D_{x^0}^* D_{x^0} Q \succeq 0$，只要在深度负块上满足 $P_{--}\preceq \theta A_-$（$\theta\in0,1)$，可逐像素由 $\sum_{\text{largest }k}\bar\omega_i\|A_-^{-1/2}(d_i(x_i^0)Q_-)^*\|_2^2 \le \theta$ 严格验证），对任意 $z=(z_+,z_-)$ 利用 $P\succeq 0$ 蕴含 $2\Re(z_+^*P_{+-}z_-)\le \frac{\theta}{1-\theta}z_+^*P_{++}z_+ + \frac{1-\theta}{\theta}z_-^*P_{--}z_-$，立得
     $$z^*(A+P)z \le z_+^*\!\left(A_+ + \frac{1}{1-\theta}P_{++}\right)\!z_+ \implies d_{\rm mix}^2 \le \Tail_m(H) + \frac{1}{1-\theta}\sum_{\text{largest }k}\bar\omega_i\|d_i(x_i^0)Q_+\|_2^2.$$
     **【猜想/预期收益】**鉴于 [u100_best_strict_source.json 已证实 $j=29\text{--}32$ 深度负模态的缺陷远大于保留模态，该引理将把深度负模态上的虚假中心增益完全吸收进 $-A_-$，预计可在不重算前向算子的前提下大幅压低 Table 2 第 3 列并将符号中心障碍显著推离 52 像素。

#### 改动/创意 2：提出“相干分裂-零惯性混合定价”统一定理（打通全文两半），并加固损耗分离证据链
* **涉及文件与段落**：
  * contrast.tex（§3.2 第 4 段及 Proposition 3.2）；
  * experiments.tex（§5.3 第 1–3 段）；
  * manuscript_ip.tex（Abstract 第 8–9 句）与 manuscript_ip.tex（§6 第 3 段）。
* **具体可执行内容**：
  1. **改变研究方向的统一算法定理（已验证 4 行证明）**：利用 Frobenius 范数的酉不变性及固定损耗角扇区 $Q_r\preceq 0$、相角等式 $Q_{\rm phase}=0$ 与耗散不等式 $(J-J_c)^*M(J-J_c)\preceq C_B$ 均在 $s\times s$ 源矩阵序下成立的事实，对任意正交相干源基 $U=[A_1, Q_2]$（其中 $A_1=[a_j]_{j\in\mathcal S_1}$），将 $K_j(\mu_j,\lambda_j,\tau_j)$ 保持严格正定（$\nu_j=0$，对任何 $t>0$ 只要 $\tau_j$ 足够大必然存在）的方向与共享负惯性的子空间 $Q_2$ 联合，立得
     $$\delta(u;Y)^2 \ \ge\ \sum_{j\in\mathcal S_1}\max_{\substack{\mu_j,\lambda_j,\tau_j\\ K_j\succ 0}} a_j^* H(\mu_j,\lambda_j,\tau_j) a_j\ +\ \Tail_{\nu(\mu,\lambda,\tau)}\!\big(Q_2^* H(\mu,\lambda,\tau) Q_2\big).$$
     该定理在 $\mathcal S_2=\emptyset,\mu_j=\lambda_j=0$ 时精确退化为 Proposition 3.2，在 $\mathcal S_1=\emptyset$ 时精确退化为整电流耗散-扇区秩界（$1.284r$），从而**严格支配全文现有的两种独立方法**，且仍只需低维标量/小锥优化，无需大尺度 SDP。**【猜想/研究方向价值】**在 $1.25\text{ m}$ 和 $96^2$ 等大域有损问题中，由于不再抛弃实部扇区与相角约束，该框架有望突破纯接收端椭球 $0.15r\text{--}0.395r$ 的极限，使“相干独立定价”真正跨越 $12\times 12$ 玩具尺度。
  2. **加固 §5.3 实证并修正笔误**：将全 $\tau$ 严格上界证明补充/主推至 finite_0_directional.json 中已呈现 $0.8025r$ vs $1.01324r$（裕度宽 30 倍且损耗角更贴近真实介质）的 $t=0.05,\,u=0.75$ 算例，避免全文核心分离结论孤悬于 $t=0.2,\,u=1.0$ 仅 $0.127\%$ 的边界裕度；同时将 experiments.tex 的 `0.5 ms` 更正为 `0.88 ms`，并说明 `0.319 s` 来自撞及 `iteration_limit: 2500` 的一阶投影。

#### 改动/创意 3：重塑中心缺陷章节定位，并严格界定连续算子（Weyl 律与 2D TM 连续迁移）适用范围
* **涉及文件与段落**：
  * physics.tex（§4.3 与 §4.4）；
  * appendix.tex（Weyl 律证明段）；
  * manuscript_ip.tex（§6 Discussion 末段）。
* **具体可执行内容**：
  1. **瘦身并诚实重定位 §4.3**：将式 (4.4) 中的 $\widehat J$ 与 $-\eta^{-1}(K\widehat J-L)^*(K\widehat J-L)$ 移至 appendix.tex 归位为“非精确驻点数值验证”，并在正文 §4.3 明确给出为什么 `u100_best` 上全细胞标量缺陷仅增 $10^{-4}r$ 的两条机制诊断（悬崖正谱隙 $g_+/\mu=1.95\times 10^{-3}$ + 缺陷集中在被 $\Tail_{13}$ 丢弃的前 13 模态）；若作为后续研究建议，可指出应改在 $g_+/\mu=O(1)$ 的非悬崖算子上测试**空间加权对角缺陷** $D_X^*\diag(w)D_X$（仍保留逐像素闭式投影）以避开近零正特征空间（属待验证猜想）。
  2. **严格界定 §4.4 与 §6 的连续算子适用域（已验证数学事实）**：
     * 在 physics.tex 将“smooth domain”放宽为 $|\partial\Omega|=0$（覆盖正文正方形网格），并在 appendix.tex 明写 $A=k^2(1-\Delta_{\R^2})^{-1}$ 及其差算子的 $-4$ 阶迹类性质，同时点明有界域上 $K/\mu$ 仅有离散谱累积于 $\mu$、无连续谱。
     * 在 §4.4 明确区分“固定 $k,\,u\to\infty$ Weyl 律”与“固定 $u$、大 $ka$ 下的倏逝环带 $k<|\xi|<k\sqrt{1+u}$”，指出 $96^2,u=1$ 时实测 $\nu=13$ vs Weyl 预测 $22.35$ 的差距来自随 $u$ 变动的 $S^*S/\mu + (\lambda/\mu)R$，并写明由面积 $(ka)^2u$ 对抗周长 $2ka$ 导出的结构性相变律 $u_{\rm crit}\sim O((ka)^{-1})$。
     * 在 manuscript_ip.tex 补充 **2D TM 极化的独特数学红利与物理边界**：在 2D TM 标量 Helmholtz 下，因弱-$*$ 紧性与无均匀化，粗糙类 $L^\infty(\Omega;[0,u])$ 上存在一致预解界与一致 $O(h)$ 连续-离散前向误差界（使有限网格排除向连续介质迁移在理论上可闭环）；但该结论**不适用于** 2D TE 或 3D 矢量 Maxwell 散射。

