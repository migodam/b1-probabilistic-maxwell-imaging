# 论文与端点判读

原意见由独立 CLI 研究会话形成。下文保留其数学提案、当时可见数据与未复核判断。作者接受、补算和纠正的结论见[作者回应](RESPONSE_ZH.md)。

以下意见基于对 `guided_20261125/final/` 内最新一轮稿件（含 14:40 更新版 `theory.tex`、`results.tex`、`technical.tex` 及 `same_site.json`、`cie_with_corrector_*.json`）与全部原始数据、对照脚本的逐项复算核查。全文严格区分**【已核对】**（直接由源码/数据/代数验证确认）与**【猜想/建议】**（尚未跑实验的新算法或实验构想）。

---

## 一、六大核心审查项：已核对的结论与数据/公式纠错

### 1. 奇异 $K$ 时带边矩阵正负计数：理论 $b(t)$ 与实用 $\np$ 检验存在本质不对称

**【已核对】代数机制与反例**
设 $K\in\C^{n\times n}$ 奇异，$\dim\ker K = k_0$，记不一致零空间秩 $r_1 = \operatorname{rank}(L^*|_{\ker K})$，**一致（持久）零空间维数** $d_0 = \dim(\ker K \cap \ker L^*) = k_0 - r_1$，记缩减广义 Schur 补 $H_{\rm red} = (C - L^* K^\dagger L - t I_s)|_{\ker(L|_{\ker K}^*)}$（维数 $s-r_1$）。由块合同变换，$\mathcal A - t P_s$ 的精确惯性为：
$$\nn(\mathcal A - t P_s) = \nn(K) + r_1 + \nn(H_{\rm red}),\qquad \np(\mathcal A - t P_s) = \np(K) + r_1 + \np(H_{\rm red}),\qquad n_0(\mathcal A - t P_s) = d_0 + n_0(H_{\rm red}).$$
由此立刻暴露出 Theorem 2.1（theory.tex）中两个准则在奇异 $K$ 处的**分裂**：
1. **积分计数 $b(t) = [s - \nn(\mathcal A - t P_s)]_+$**：化简为 $[\np(H_{\rm red}) + n_0(H_{\rm red}) - \nn(K)]_+$，**不受 $d_0$ 损失**，数学上完全严格。
2. **单阈值正特征值计数 $\np(\mathcal A - t P_s) \ge n + j$**（亦即 run.py 实际执行的排除判据）：化简为
   $$\np(\mathcal A - t P_s) - n = \np(H_{\rm red}) - \nn(K) - d_0.$$
   - 在**一般非一致奇异穿越**（$d_0=0, r_1>0$）时，$\mathcal A - t P_s$ 保持非奇异（增加 $r_1$ 对 $\pm 1$ 双曲特征值），$\np\ge n+j$ 确实有效且避免了求 $K^{-1}$。
   - 但在**存在持久零方向**（$d_0 > 0$）时，$\np$ 检验会遭受与负特征值同等级的**虚幻 $-d_0$ 惩罚**；更关键的是，**外向舍入区间算术（Arb/binary64 误差界）无法直接证明 $n_0$**（精确零特征值只能包络为 $[-\epsilon,+\epsilon]$，可能判为负），因此任何严格数值验证都只能退化到 $\np$ 检验，而不能直接实现理论的 $s-\nn$ 检验。
3. **直接在稿件自带的有理奇异例上验证失效**：在 technical.tex 的例子（$n=1, s=2, h=W=1, V=0, T=5/4, E_c=(3/4,0), P=(3/4,0)^\top, y=\diag(-5/4,2)$）中，经核算确有 $K=0, L=(0,0), C=\diag(1,4)$（$d_0=1, r_1=0$），且 $\Delta=-1$ 达到残差平方 $5$。此时 $\mathcal A - t P_s = \diag(0, 1-t, 4-t)$：
   - 理论积分 $\int_0^\infty b(t)\,dt = 2\times 1 + 1\times 3 = 5$（准确）。
   - **但单阈值检验 $\np(\mathcal A - t P_s) \ge 1 + j$**：当 $j=1$ 时仅对 $t<1$ 成立（最多排除 $r^2<1$，损失 80% 的真实下界 $5$，甚至连单特征值 $4$ 都够不到）；当 $j=2$ 时因永久零特征值，$\np\le 2 < 3$ 对所有 $t\ge 0$ 恒失败！若只取第 2 个源（$s=1, C=4$），$b(t)$ 给出下界 $4$，而 $\np$ 检验给出的下界为 **$0$**。
4. **实证覆盖缺口**：restricted.json 虽测试了近奇异 $K_U$（`state_float_min_abs` $\approx 1.16\times 10^{-8}$, $\nu=13$），但 4 个受限空间全部未排除（`positive_count_lower = 93..100` vs `required = 101..132`）；confirmation.json 中 98 次排除全部在 `j=1, trials=1` 通过；native_l4.json 中 2 次排除均为 $\nu=0$。也就是说，**全文尚无一例实证排除是在真正奇异或近奇异 $K$ 上触发的**。

---

### 2. 区间套入的权重传递：证明措辞瑕疵、向上舍入违规、行序错配与未开发自由度

**【已核对】**
1. **Proposition 2.2 证明笔误**（theory.tex）：
   - 证明写 *"For $\alpha=1$ inclusion implies equality of the intervals."*，表述有误：一般 $\alpha=1$ 并不蕴含区间相等（例如同心子区间 $m=m', h<h'$ 时 $\alpha=1$ 同样满足 $Q_I - Q_{I'} = (h'^2-h^2)e^*e\succeq 0$）。真正需要说明的是：在取定 $\alpha=h'/h$ 的前提下，$\alpha=1$ 等价于 $h=h'$，此时由 $I\subseteq I'$ 才推出 $I=I'$（从而 $M(1)=0\succeq 0$，补足 Sylvester 主子式 $\alpha-1>0$ 的边界情形）。
2. **精确可行 $\alpha$ 区间与端点相切退化**：
   - 二次判别式 $\det M(\alpha) = -h^2\alpha^2 + (h^2+h'^2-(m-m')^2)\alpha - h'^2 \ge 0$ 在整段闭区间 $[\alpha_-, \alpha_+]$ 上成立，且由韦达定理 $\sqrt{\alpha_-\alpha_+} = h'/h$（文中取的是几何中点）。
   - 当子区间与父区间**共享端点**时（如 $[2,6]\subset [0,6]$，$|m-m'|=h'-h=1$），$\det M(\alpha) = -(\alpha h - h')^2$，此时 $\alpha=h'/h=1.5$ 是**唯一**满足逐行优势的标量。
   - 当子区间**严格内嵌**时（如 $[3,5]\subset [0,6]$ 得 $\alpha\in[1.1459, 7.8541]$；$[3.5,4.5]\subset [0,6]$ 得 $\alpha\in[1.1296, 31.8704]$），存在很大的 $\alpha$ 优化空间。
3. **秩-$2k_p$ 更新公式的行序错配**（theory.tex）：
   - 式 \eqref{eq:update} 的代数本身经数值核对**正确**（$\mathcal H' = \mathcal H_b + A_p^* Z_p^{-1} A_p$ 与 $\nn(K') = \np(Z_p) - \np(D_p^{-1})$ 均精确成立），但文字写 *"stack the current and Green rows on the changed cells as $T_p$, the incident rows as $e_p=[0;B_{\rm patch}]$, and put $D_p=\mathrm{diag}_{\rm blocks}(D_i)$ in matching row order"* 自相矛盾：若按块堆叠使 $e_p=[0; B_{\rm patch}]$，则 $D_p$ 是 $2\times 2$ 块对角矩阵而不是 $\mathrm{diag}_{\rm blocks}(D_i)$；若按像素交错排列 $(j_i, e_i)$ 使 $D_p = \mathrm{diag}_{\rm blocks}(D_i)$，则 $e_p$ 必须交错为 $(0, B_1, 0, B_2,\dots)^\top$。
4. **单调性核对通过，但 results.tex 存在 2 处违反“向下舍入”原则的数值错误**：
   - 将 nested.json 与 mechanism/results.json 逐项比对：在全部 8 组 `(site, b, mu)` $\times$ 3 个子宽度 $\times$ 8 个读数上，权重传递后的子盒下界均严格 $\ge$ 父盒 $[0,6]$ 下界，且在 $b=1/8$ 时负惯性确实从父盒 $\nu=6$ 降至子盒 $\nu=5,4,4$（而未传递权重的匹配 $w$ 反而升至 $\nu=8$，导致排除数从 6/6 跌至 1–3/6）。
   - **违规 1（严格下界向上舍入）**：results.tex 写 *"with lower endpoints from $2.36r$ to $5.43r$"*。核查 nested.json，`site=2, b=0.125, mu=0.0009765625, w=2` 上 6 个主目标的最大严格下界为 **`5.42283978`**（最小为 `2.36352302`）。严格下界向下舍入应为 **`5.42r`**，写 `5.43r` 构成越界（不安全下界）！
   - **违规 2（同位点低幅值对照也向上舍入，且未交代跨 $\mu$ 取优）**：
     - same_site.json 和 same_site.json 的精确值为 `0.2427235` 与 `0.2772215`，results.tex 写成 `0.243r` 和 `0.277r`——前者向上舍入，应改为 **`0.242r`**（顺便 noted：`same_site.json` 中 4 条未排除项出现未截断的负非规格化浮点 `-4.9e-322`，公开数据前建议截断至 `0.0`）。
     - 同一段中两个控制目标的错误位点下界 `2.518r`（`site=1, mu=0.0009765625`）与 `0.496r`（`site=2, mu=0.000244140625`，在 `0.0009765625` 下仅为 `0.2967r`）分别取自两个不同标量 $\mu$，虽对下界合法（$\max_\mu$），文中应明言“在两个预设标量基权重中取较大者”。

---

### 3. 原生非收缩例到底排除了什么：5 行表里只有 1 个真正分离，且高度依赖精确外部真空

**【已核对】逐行拆解 Table 2（paired_summary.csv 与 results.tex）**

1. **前 4 行均不足以证明“需要 $\nu>0$ 的同支持高对比度分支排除”**：
   - 第 1、3、4 行（`r10_twod_0001`, `optfresh_01`, `optfresh_02`）：把盒下界放宽到 $0$（同支持 `[0, c_i+w]`）后**依然被强排除**（`56.14r`、`3.34r`、`6.27r`），说明排除主要来自**支持集与外部固定真空的失配**，而非排除了高对比度分支本身。
   - 第 2、3、4 行（`shiftbox32_01`, `optfresh_01`, `optfresh_02`）：同一盒上**正定二次型（$\nu=0$）的下界反而全面高于 $\nu=9$ 的不定二次型**（`12.23 > 10.70`, `25.05 > 18.94`, `23.20 > 21.98`）。核查 mechanism/verified.json 可知，原因是这三例的 9 个负扇区方向均对接收阵可辐射（不在 $\ker P$ 内），把标量 $\mu$ 缩小 $256\times$ 后 $P^*P$ 就能压正全部负特征值。
2. **第 5 行（`control_disk075`，1153 格，$w=2.5$，$\nu=7/32$，$d/r=1.7284581$）是全文唯一的 $N=1$ 真正非收缩同支持分离例**，但它精确排除了什么、**没排除什么**必须说透：
   - **确已证明**：在**全部 $9216-1153=8063$ 个外部像素固定为精确真空**的前提下，以 Born 伴随峰值定位的 1153 格圆盘支持上，族 $\mathcal B = \{\max(0, c_i - 2.5) \le \chi_i \le c_i + 2.5\}$（名义锥形峰值为 4，顶点处为 $[1.5, 6.5]$）中**不存在任何全波解**落在 2% 噪声球内（下界 $1.7284581r$）；与此同时，真实的 416 格、峰值 $0.625$ 圆盘位于**同一 1153 格支持内**并满足全波残差 $\le 0.999950r$。此外，均匀上界扰动 $+2.5 I$ 确实落在盒内且经验证 $\rho(2.5T)>1.0864$（Neumann 级数真发散，不只是 2-范数充分条件 $w_{\rm Neumann}\approx 0.48605$ 的 $5.1435\times$ 超出）；且因存在接收零空间的“暗负模”（扇区值 $\le -0.2220$，辐射 $\le 1.474\times 10^{-17}$），任何**标量权重 $W=\mu I$**（含 $|\tau/\mu|<1.506\times 10^{16}$ 的全局光学项）都不可能使 $K\succ 0$。
   - **尚未证明 / 摘要与正文需收紧的四点**：
     1. **截断非对称盒被简写为“半宽 2.50”**：disk_radius.json 明确盒定义为 `max(0, c_i - w) <= m_i <= c_i + w`。在边缘格 $c_i < 2.5$ 处下界被截断在 $0$，真实半宽是 $(c_i+2.5)/2 < 2.5$、中点发生偏移。摘要（manuscript_ip.tex）仅写 *"centred at peak four with half-width 2.50"* 会让读者误以为所有活跃格均未截断（且限制像素 `(33,49)` 处 `c_i=3.864` 而真值圆盘为 `0.0`，恰恰体现高分支在该核心区域强制了严格正下界 $\ge 1.364$）。
     2. **外部真空假设是决定性的**：一旦在 same_site.json 中允许 9152 个外部背景像素在 $[0, 1/16]$ 或 $[0, 1/8]$ 自由变化，同支持上的峰值 0.625 与 $[2,6]$ 分支即**完全无法分离**（下界仅 `0.242r`–`0.384r`）。
     3. **Table 2 中 `Pos.` 列写 `---` 易被误读为“所有正定二次型均不可能”**：disk_radius.json 仅否定了**标量 $W=\mu I$**，并未否定一般逐格对角权重 $W\succeq 0$。
     4. **固定非零宿主实验措辞有误**（results.tex vs paired_summary.csv）：正文写 *"retaining the full background gives lower bounds $29.006r$ and $27.730r$"*，极易误导为“保留了 9216 格未知背景”。实际 `paired_summary.csv` 显示 `active_cells=515, unknown_cells=437, fixed_nonzero_cells=78`——其余 8701 格仍是**固定真空**，只是保留了 78 格已知对比度 `3.25` 的宿主散射相互作用（且表中另有 `inverse_center_illumination` 权重给出 `26.355r` 和 `26.098r`，以及 `peak_lower=3.25, peak_upper=6.5, box_half_width=2.5` 的下界截断信息未在正文说明）。

---

### 4. 证明和单图精修比较的计时：Table 4 漏报 2.882 s Arb 重验成本，且两组比较的停止规则不对称

**【已核对】**
1. **Table 3（前瞻 8 目标，confirmation.json 与 confirmation_online.json）数值精确无误，但有三大隐含结构必须挑明**：
   - **在线耗时与 Arb 严格证明耗时相差一个量级**：`Test surviving boxes`（`fast_boxes`）中位耗时 `54.6 ms`（`179.5 LU`）**仅含 binary64 误差界**（其注释在 run.py 明确依赖 BLAS 精度假设）；一旦计入独立的 Arb 复球重验（中位 `362.7 ms`），每目标中位总耗时升至 **`419.1 ms`（~0.42 s）**，是 `Points until first fit`（`top_early`，`37.1 ms`）的 **11.3 倍**、甚至比穷举全部 16 盒（`all_cluster`，`82.5 ms`）慢 **5.1 倍**。
   - **停止规则不对称**：逐案核查 `confirmation_online.json` 发现，`fast_boxes` 在 **8/8 个目标上都精修了全部存活盒**（`surv == fb_visited`，共 2–7 盒，中位 3.5 盒），而 `top_early` 在 7/8 目标上**精修 1–2 盒即早停**。因此 `54.6 ms` vs `37.1 ms` 并不只是筛盒开销（`fast_boxes` 筛盒仅 `10.93 ms` vs `top_early` 64 样本正演 `4.25 ms`），更反映了“精修全部存活盒”与“早停找首解”的任务差异！
   - **真正胜出场景与宽盒主导现象**：
     - 在唯一未于 42 步内拟合噪声球的困难目标 `confirm49_20261125106` 上，`top_early` 退化为遍历全部 16 盒（`74.6 ms / 577 LU`），而 `fast_boxes` 先排除 14/16 盒，耗时仅 **`56.8 ms / 119 LU`**（**反超**，LU 省 4.8 倍）。
     - 在全部 8 个前瞻目标上，覆盖真值的盒索引仅为 `{0: 8, 15: 2}`，最终胜出盒在 8/8 目标上**全部是宽先验盒 `box 0`**——这就是三组方法图像/留出误差中位数完全一致（`10.02%`, `1.16%`, `1.67%`）的直接原因，正文应坦承“前瞻集内窄盒均未覆盖真值，三臂均收敛到同一宽先验盒 0”。
2. **Table 4（新增 12 目标 CIE+校正器，results.tex）漏报了 `2.882 s` 的 Arb 重验耗时！**
   - 核查 cie_with_corrector_online.json 与 cie_with_corrector_offline.json：
     - Table 4 报告 `All boxes = 4.963 s`、`Point rank until fit = 0.890 s`、`Neumann survivors = 6.456 s`、`Bordered-test survivors = 1.892 s`。
     - 但 `1.892 s` 中只含 `0.103 s` 的 binary64 在线筛盒，**未计入 `cie_with_corrector_offline.json` 中 `arb_checked_shifted_exclusions` 的独立 Arb 复球重验总耗时 `2.882 s`**！
     - 若按与 Table 3 同等口径加上 Arb 重验，Table 4 中 `Bordered-test survivors` 总耗时变为 **`4.774 s`**——几乎与穷举 `All boxes`（`4.963 s`）打平，且比 `Point rank until fit`（`0.890 s`）慢 **5.36 倍**！不补注这一数字，审稿人对照附件 JSON 会立刻质疑计时口径前后不一。

---

### 5. 真实 Flow 实验与“72 个物体”口径

**【已核对】**
1. **真实 Flow 实验（native_l4.json 与 results.tex）**：
   - **两处被排除盒均为“正定且收缩”**：36 个盒中仅排除 2 个（`140000` box 1: `3.2943r`; `140111` box 9: `8.4557r`），核查 `native_l4.json` 显示二者均为 `negative_inertia = 0` 且 **`strict_kappa_upper = 0.5101, 0.5155 < 1`**！Neumann 未排除它们不是因为非收缩，而是三角不等式范数界过松；因此 Flow 实验**没有进入非收缩/负惯性区间**，其真正科学价值是**审计生成先验的覆盖失效与严重非唯一性**（0/36 盒含真值，却有 15/36 盒和 25 对不交盒对在 1% 误差球内找到可行解，材料误差高达 26.35%–79.36% 而留出接收误差 <1%）。
   - **34 s vs 21.36 s 的账目拆分缺失**：results.tex 写 *"The sector search and verification here take about 34 s for the three readings"*，而 `native_l4.json` 的 `total.wall_seconds` 仅为 **`21.358 s`**（`setup 2.801 s` + 仅 2 个排除盒的严格验证 `8.886 s + 9.619 s`），其余 34 个 `RETAIN_UNRESOLVED` 盒的浮点扫描时间（约 `12.6 s`）未记录在 `native_l4.json` 中。应在正文显式拆分为“约 13 s 浮点预筛 + 21.36 s 对 2 个排除盒的严格区间验证”。另外，box-first 总耗时（`29.87 + 21.50 + 27.72 = 79.09 s`）慢于穷举精修（`53.88 s`）且是点排序早停（`23.67 s`）的 3.3 倍，稿件已坦诚给出数字，值得肯定，但应在结论中明确指出“在 32² 全域流先验下，盒筛选不是加速手段，而是多解审计工具”。
2. **“72 个物体 / 121 个视角”口径审查（objects.csv 与 Table 1）**：
   - **视图计数矛盾（120 vs 121）**：核算 `objects.csv` 共 73 行，其中 `independently_counted=True` 的 72 个物体合计 **`120` 个 `observation_views`**，加上被标记为 `False` 的单独 6 参数多解例 `GENERAL6_00` 的 `1` 个视图才等于 `121`。而 Table 1 标题（results.tex）说 *"The 121 recorded acquisition views are paired conditions. A separate six-parameter observation tests ambiguity."*，把 `GENERAL6_00` 既说成“单独的”又算进了 121。
   - **新增 `same_site.json` 对照未入表**：results.tex 新加入两个峰值 `0.625` 的同位点低幅值读数，并称 *"in addition to the registered acquisition views"*。由于对比度已从 `3.884 / 4.470` 变为 `0.625`，物理上属于新的介质分布 $\chi(x)$，建议明确标注为“基于已有两个位点支持集的 2 组额外低幅值配对控制读数（不计入 72 个独立样本总数）”。
   - **警惕“凑 72 样本”的审稿观感**：72 个样本横跨 6 种网格（$12^2$ 到 $96^2$）、3 种本构模型、5 种噪声合同，其中**真正前瞻盲测仅 8 个（11%）**、**实测暗室物体 7 个（10%，且无加载场景正演误差界）**、**核心非收缩暗模同支持分离例仅 1 个（`control_disk075`）**、**真实 Flow 仅 3 个**。切忌把“72”当作统一基准规模来强调，建议在 Section 4.1 开头直接按**核心论断的有效样本量**分层列出（$N=28$ 大域 Born/全局帽；$N=1$ 非收缩暗模同支持分离；$N=6+2$ 未知背景双位点传递；$N=8$ 前瞻候选池；$N=3$ Flow 审计；$N=7$ 条件暗室）。

---

## 二、需直接修改的稿件字句清单（按文件定位）

| 文件与行号 | 原文问题 | 建议修改字句（可直接替换） |
|---|---|---|
| manuscript_ip.tex | *"At a specified error tolerance, an inertia count on a bordered matrix tests the entire interval, including when the current block is singular."* 未区分 $n_-$ 与 $n_+$ 在持久零空间下的差异；*"centred at peak four with half-width 2.50"* 省略了非负截断与外部真空前提。 | 改为：*"At a specified error tolerance, an inertia count on a bordered matrix avoids inverting a singular current block; its negative-index count handles persistent zero modes directly, whereas strict positive-count certification retains full power at generic rank drops and loses any persistent nullity."* 并将圆盘句改为：*"Under fixed exterior vacuum on the same 1,153-cell support, every member of the nonnegative clipped class $[\max(0,c_i-2.5),\,c_i+2.5]$ about a peak-four taper misses the error ball ($>1.72r$, seven negative current modes)..."* |
| theory.tex | Theorem 2.1 未说明 $\np(\mathcal A-tP_s)\ge n+j$ 对 $d_0=\dim(\ker K\cap\ker L^*)$ 的敏感性。 | 在证明末尾补一句：*"When $K$ is singular with $d_0=\dim(\ker K\cap\ker L^*)$ and $r_1=\operatorname{rank}(L^*|_{\ker K})$, congruence shows $\nn(\mathcal A-tP_s)$ increases only by $r_1$, so $b(t)$ is unaffected by $d_0$, whereas $\np(\mathcal A-tP_s)-n$ is reduced by $d_0$ because $\mathcal A-tP_s$ inherits nullity at least $d_0$."* |
| theory.tex | *"For $\alpha=1$ inclusion implies equality of the intervals."*（逻辑不精确） | 改为：*"For the transported choice $\alpha=h'/h$, the boundary value $\alpha=1$ means $h=h'$, in which case inclusion implies $I=I'$ and the difference matrix vanishes. More generally, $(h'/h)$ is the geometric mean of the feasible interval $[\alpha_-,\alpha_+]$ where the determinant is nonnegative, collapsing to the unique valid scalar when $I$ shares an endpoint with $I'$."* |
| theory.tex | *"the incident rows as $e_p=[0;B_{\rm patch}]$, and put $D_p=\mathrm{diag}_{\rm blocks}(D_i)$ in matching row order."*（行序自相矛盾） | 改为：*"stack the interleaved pairs $(j_i,e_i)$ on the changed cells as $T_p J + e_p$ (with $e_p$ zero in each current row and $B_i$ in each field row), and put $D_p=\mathrm{diag}_{\rm blocks}(D_i)$."*（或改写为按块堆叠对应的 $2\times 2$ 对角块矩阵 $D_p$）。 |
| results.tex | *"The 121 recorded acquisition views are paired conditions. A separate six-parameter observation tests ambiguity."*（120 vs 121 冲突） | 改为：*"The 72 specimens contribute 120 recorded acquisition views across paired conditions; a separate six-parameter observation adds one ambiguity view (121 registered views in total), and two extra low-amplitude paired site controls are reported separately."* |
| results.tex | *"retaining the full background gives lower bounds $29.006r$ and $27.730r$"*（误导为未知全背景） | 改为：*"Two further experiments placed a separate, known 78-cell dielectric host of contrast 3.25 outside the 437 uncertain cells (with the remaining 8,701 exterior cells fixed at vacuum); eliminating the known host through the active-cell Green operator gives lower bounds $29.006r$ and $27.730r$ under identity weighting (and $26.355r$, $26.098r$ under inverse-illumination weighting), both with $\nu=7$ and 32 illuminations."* |
| results.tex | Table 2 中 `Pos.` 列为 `---` | 改为 `\text{obstr.}^*` 并在表注说明 `*Obstructed for every scalar W=\mu I (and |\tau/\mu|<1.506\times 10^{16}); general diagonal W not ruled out.` |
| results.tex | *"from $2.36r$ to $5.43r$"*、*"$0.243r$ and $0.277r$"*（两处严格下界向上舍入） | 将 `5.43r` 改为 **`5.42r`**，将 `0.243r` 改为 **`0.242r`**；并注明控制例的 `2.518r` 与 `0.496r` 取自两个预设基权重 $\mu\in\{2^{-12}, 2^{-10}\}$ 中的较大下界。 |
| results.tex | Table 3 正文未明言 `fast_boxes` 精修所有存活盒、且 8/8 均由宽盒 `0` 胜出 | 在第 69 行补一句：*"The point-ranked arm stops at the first feasible box (1--2 boxes on seven targets; all 16 on the hard target, where the bordered test is faster at 56.8 vs 74.6 ms and 119 vs 577 LU), whereas the bordered arm refines all 2--7 surviving boxes; on all eight targets the broad prior box 0 wins."* |
| results.tex | Table 4 漏报 `2.882 s` 的 Arb 重验耗时 | 在第 78 行和 Table 4 表注补明：*"The reported 1.892 s uses online binary64 error bounds; independent Arb ball replay of all shifted exclusions adds 2.882 s across the 12 targets (4.774 s total, comparable to exhaustive 4.963 s)."* |
| results.tex | *"both with positive state curvature"*、*"about 34 s"* | 补注 *"both with positive state curvature and contractive operator norm ($\kappa \le 0.52$)"*，并将 *"about 34 s"* 拆明为 *"about 34 s in total (roughly 13 s for the 36-box float scan plus 21.36 s for strict Arb verification of the two excluded boxes and shared setup)"*。 |
| technical.tex | 有理奇异例只报告了 $b(t)$ 积分为 5，掩盖了 $\np$ 单阈值检验仅得 $<1$ | 在该段末尾补一句：*"At the same time, because $\mathcal A-tP_s=\diag(0,1-t,4-t)$ retains a zero eigenvalue, the strict positive-count test $\np(\mathcal A-tP_s)\ge 1+j$ certifies only $t<1$ for $j=1$ and never triggers for $j=2$, illustrating the $-d_0$ gap between $s-\nn$ and $\np$ at persistent singularities."* |

---

## 三、一项不夸大的 *Inverse Problems*（IP）定位

基于全部经核查的正面与负面证据，最稳健、抗审稿攻击的 IP 定位应是：

> **“面向有限孔径非线性电磁逆散射的可认证材料假设否决（Certified Material-Class Falsification），作为重建求解器的正交验证层而非加速替代品。”**

具体界定三条清晰边界（建议直接在摘要末句、Introduction 末段和 Discussion 呼应）：
1. **核心数学贡献（强项）**：从共材料逐格实区间约束导出统一的**不定 Hermitian 二次型惯性下界框架**（Theorem 2.1 的带边矩阵计数 + 非奇异 Schur 尾和），配套**非同心区间嵌套的显式权重传递**（Proposition 2.2）与**未知弱背景下的秩-$2k_p$ 局部补丁更新**（\eqref{eq:update}）。它能在全波、多照明、超越 Neumann 收缩半径乃至存在暗负模的强散射体制下，给出数据空间严格距离下界。
2. **明确不宣称什么（防过度包装）**：不宣称比“单点全波排序+早停精修”更快得到**单张**可行图像；不宣称提供连续介质或暗室加载场景的先验建模误差界；不宣称能用线性 Born 否决或负模范数下界指导无约束非线性优化器的初始化或步长选择。
3. **实际不可替代的应用生态位**：(a) **物理断言的整类否决**（排除一整片高对比度或错误位点假设，而不仅是反驳单个点解）；(b) **生成式/学习先验的覆盖与非唯一性审计**（如 Flow 实验中揭示 0/36 覆盖与 25 对不交可行分支）；(c) **困难/无早停命中实例的死分支剪枝**（如前瞻例 `confirm49_20261125106` 中 LU 降为点排序的 $1/4.8$）。

---

## 四、能改变研究方向的算法思想与低成本决定性实验（区分已核对事实与猜想）

以下 4 个方向均直击上述已核对瓶颈，可在后续研究或大修中以极低计算成本显著改变方法的能力上限：

### 思想 1：用“纯双精度/Rump 区间乘法（Ozaki 分裂）+ 认证 $LDL^*$ 惯性”替代 104-bit Arb 标量循环，消灭 10 倍验证墙【猜想】
- **【已核对依据】**：在 Table 3 和 Table 4 中，真正的计时杀手不是带边矩阵筛盒本身（binary64 仅 `10.9 ms`/目标、`0.103 s`/12 目标），而是 Arb 逐元素高精度复球重验（`362.7 ms`/目标、`2.882 s`/12 目标），使严格证明比穷举精修还慢。
- **【算法猜想】**：对 $53\times 53$ 乃至 $200\times 200$ 矩阵，完全不需要调用 104-bit Arb 做稠密 $Q^* A Q$。利用 **IEEE-754 定向舍入（`fesetround`）的 Rump 中点-半径区间矩阵乘法**（仅需 2–3 次 BLAS `zgemm` 即可得到数学上 100% 严格的含舍入外包络，无需假设 BLAS 无误差！），再对小矩阵做带外向舍入补偿的 Cholesky/$LDL^*$ 惯性检验。
- **【预期改变与验证实验】**：预计可将 Table 3 的严格复球验证从 `362.7 ms` 压到 **`< 5 ms`**，将 Table 4 的总严格耗时从 `4.774 s` 压回 **`~1.95 s`**（首次在**包含完全严格证明**的前提下显著快于穷举 `4.963 s`！）。验证实验只需在现有 `confirmation.json` 的 98 个盒上对比定向舍入 BLAS 认证与 Arb 认证的余量与耗时（数秒内完成）。

### 思想 2：“先早停命中、后惰性带边否决”混合流水线（Hit-First, Certify-Rest）【猜想】
- **【已核对依据】**：在 `confirmation_online.json` 中，`top_early` 仅花 `4.25 ms` 扫 64 点并在 7/8 目标上精修 1–2 盒即命中（`37.1 ms`）；但在无早停命中的目标 `106` 上浪费了 `577 LU`（`74.6 ms`），而 `fast_boxes` 只用了 `119 LU`（`56.8 ms`）。同时 `fast_boxes` 之所以比 `top_early` 慢，一半是因为它强行精修了全部 2–7 个存活盒。
- **【算法猜想】**：改为两阶段自适应策略：(i) 先花 `4.25 ms` 做点排序，仅精修排第 1（最多前 2）个盒；(ii) **若命中噪声球且用户只需一解，立即结束（`~37 ms`）；若前 2 盒未命中（如目标 `106`）或用户要求对未访问分支出具“非可行证明”，再触发带边矩阵惰性筛除其余盒**。
- **【预期改变与验证实验】**：在 8 个前瞻目标上，寻找单解的中位时间应与 `top_early` 持平（`~37 ms`），同时在困难目标 `106` 上保持 `fast_boxes` 的 `~57 ms` 优势；若需证伪其余分支，也无需重复精修已筛除盒。该策略可直接复用 `confirmation_online.json` 中的逐盒数据零成本仿真评估。

### 思想 3：释放 Proposition 2.2 的自由度（逐格优化 $\alpha_i\in[\alpha_{i,-},\alpha_{i,+}]$ + 局部光学乘子 $V_p$）与对角 $W$ 暗模消除实验【猜想】
- **【已核对依据】**：
  1. 在 `control_disk075` 上，标量 $W=\mu I$ 被暗模阻断无法正定，但理论允许逐格对角 $W\succeq 0$ 与 $V$；
  2. 在 `same_site.json` 中，未知背景 $b=1/8$ 下同支持峰值 `0.625` 与 $[2,6]$ 仅得 `0.242r`–`0.277r`；
  3. 对严格内嵌子区间（如 $[3,5]\subset [0,6]$），$\alpha_i$ 并非只能取 $h'/h$，而是可在宽区间 $[\alpha_{i,-},\alpha_{i,+}]$ 内逐格选择且仍保证优于父盒。
- **【算法猜想与决定性小实验】**：
  1. **实验 A（暗模是否可被非标量对角 $W$ 消除）**：在 `control_disk075`（1153 格，仅做浮点探索即可，数秒级）上，针对已找到的暗负模 $u_{\rm dark}$ 构造非均匀对角权重 $W_i \propto (|u_{{\rm dark},i}|^2 + \epsilon)^{-1}$ 或以一阶次梯度最大化 $\lambda_{\min}(K(W,V))$，检验是否存在对角 $(W,V)$ 恢复 $\nu=0$。若**存在**，说明暗模阻断只是标量权重的局限，对角加权可彻底提升下界；若**连小型 SDP 也证实不存在**，则把原稿的“仅限标量 $W$”升级为真正物理级的“全对角二次松弛暗模障碍定理”！
  2. **实验 B（$2k_p$ 缩减空间上的局部 $(W_p, V_p)$ 寻优）**：在 `same_site.json` 的 $8\times 8$ 补丁上，矩阵维度仅 $2k_p = 128$，保持背景权重不动、仅优化补丁 64 格的局部 $(\alpha_i, v_i)$（仍在 \eqref{eq:update} 秩-$2k_p$ 框架内），直接判定 `0.242r`/`0.277r` 究竟是“固定均匀权重不够紧”还是“$b=1/8$ 未知背景下的根本松弛极限”。

### 思想 4：奇异 $K$ 的“持久零空间约化（Deflation）”与面向 Flow 未决盒的分支定界传递【猜想】
- **【已核对依据】**：(a) 我们证明了 $\np(\mathcal A-tP_s)$ 在 $d_0=\dim(\ker K\cap\ker L^*)>0$ 时损失 $d_0$ 个源模式，即使理论 $b(t)$ 不受影响；(b) 真实 Flow 实验中有 19/36 个盒处于 `RETAIN_UNRESOLVED`（且 0/36 含真值）。
- **【算法猜想】**：
  1. 对已知结构性零空间（如局部 $W_i=0$ 导致的暗零方向或接收零空间内精确零模），在调用 $\np$ 检验前先用正交基约去 $U_0\subseteq \ker K\cap\ker L^*$（在 $U_0^\perp\oplus \C^s$ 上检验 $\np \ge (n-\dim U_0) + j$），即可在**保持外向舍入严格性**的同时彻底消除 $-d_0$ 虚幻惩罚。
  2. 对 Flow 实验中 19 个未决盒，与其在粗簇盒上反复调参，不如沿各簇粒子主成分（PCA 第一方向）二分 1–2 层，并直接利用 **Proposition 2.2 的显式权重传递**复用父盒算子（避免对子盒重算全部转移矩阵），检验能否把 19 个未决盒拆分成“可严格排除子区”与“真多解子区”，从而把 Flow 实验打造成真正的**生成先验可行域剖分（Set Inversion）**示范。
