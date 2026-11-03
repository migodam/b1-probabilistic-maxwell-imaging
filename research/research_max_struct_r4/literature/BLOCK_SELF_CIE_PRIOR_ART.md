# Block/self CIE：最接近的原文公式与证据边界

检索日期：2026-11-03。任务仅为 prior-art extraction；不访问本轮数值结果/GT，不判定核心等价定理或原创性。A=已读可定位全文；C=仅身份核验。**已核实：局部 Green 项抽取、空间变参数 NIE/CIE、自作用极化率重整化、一般矩阵 T 重整化，以及 Y₀-CSI/NIE 混合均有直接先例。尚未核实同一四单元块、同一有限残差目标与更新算法。**

## 1. CS-EB 与任意空间 β 的 NIE：直接分母先例

Bevacqua–Isernia，2020 预印本/2021 OJAP 综述，A；[作者全文](https://arxiv.org/pdf/2012.15252)，[正式 DOI](https://doi.org/10.1109/OJAP.2021.3063248)。页码以下按 12 页作者版本。原文使用 \(W=\chi E_t\)、内部算子 \(A_i\)。

**CS-EB，§III.A，p3，Eq5–6：**

\[
 W=pE_i+pA_i^{CS-EB}[W],\quad
 A_i^{CS-EB}=A_i-f_D I,\quad
 f_D(r)=\int_D G_b(r,r')\,dr',\quad
 p(r)=\frac{\chi(r)}{1-\chi(r)f_D(r)}.
\]

其 \(f_D\) 是**域积分形成的逐点乘法函数**，不能直接写成离散 Green 矩阵的对角元。综述追溯 Isernia–Crocco–D’Urso 2004 与 D’Urso–Isernia–Morabito 2010；本轮未获取这两篇原文，因此上述公式依据是已读综述。原文也讨论将参数化的域尺寸用于构造等价方程族。

**NIE，§III.B，p4，Eq9–10：**

\[
 \beta(r)W=R(r)\beta(r)W+
 R(r)\{E_i+A_i[W]\},\qquad
 R(r)=\frac{\beta(r)\chi(r)}{1+\beta(r)\chi(r)}.
\]

原文明确允许不同空间函数 \(\beta(r)\)，要求分母非零；并指出增大 β 后 R 接近支撑指示函数，会恶化从 R 提取材料的步骤。Eq11–12 的算子简写涉及 β 与积分算子的位置；**空间变 β 的映射宜从 Eq9 核对，不可把乘法与积分算子自动交换。**

**与当前 self-only 的关系（对照推断，交由主线程核准）：** 若 \(G_{\rm in}=\operatorname{diag}G\)，当前 \(B=\beta I-G_{\rm in}\) 是已知逐点乘法。其分母 \(I+XB\) 与空间参数 NIE 极为接近；仍须逐项比较外部 β 的残差分子、目标归一化、实际优化坐标和 material/current 更新。分母同族不能自动推出有限罚函数或迭代轨迹相同。原文范围是逐点 β，不是一般非对角块。

## 2. 已落地 CIE-I / SOM：不是只有代数重写

Zhong–Xu，*Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*，2019，A；[出版全文](https://doi.org/10.3390/jimaging5020027)，[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/)。

本地已存出版 XML：`research/research_fbg_r3_experiments/audit/literature_cie2019.xml`。§2.1 Eq3/5 使用与上式相同的
\[
 \beta I_l=R\beta I_l+R(E_l^{inc}+G_D I_l),\qquad
 R=\beta\chi/(1+\beta\chi);
\]
正文明确 β 可为常数或随位置变化。§2.2 Eq13–15 对数据残差与 CIE 电流残差平方作归一化求和，§2.3 在 SOM/FFT 子空间中交替更新辅助电流系数与 R，并采用不同 β/子空间尺寸的分阶段策略。因而“允许不满足本构关系的辅助电流，并改变本构残差尺度”本身已有实际反演实现。这里没有读到把 β 替换为真实四单元 Green 块的同一算法。

NIE 最初论文 Zhong–Lambert–Lesselier–Chen，2016，DOI [10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492)，身份已由 Crossref 核验；**原文仍未取得**，本轮不能声称逐式检查了这篇原始论文。

## 3. DCTMC II：离散自项与一般矩阵重整化的强先例

Levinson–Markel，*Solution of the nonlinear inverse scattering problem by T-matrix completion. II. Simulations*，2016，A；[作者全文](https://arxiv.org/pdf/1505.06777)，[正式 DOI](https://doi.org/10.1103/PhysRevE.94.043318)。§II，p3，Eq17–24：

\[
 d_n=\alpha_n\left(e_n+\sum_m\Gamma_{nm}d_m\right),\quad
 \alpha_n=\frac{h^3\chi_n}{1-(kh)^2(\xi+ikh)\chi_n},
 \quad \Gamma_{nm}=(1-\delta_{nm})G_0(r_n,r_m),
\]
\[
 T[V]=(I-V\Gamma)^{-1}V,\qquad AT[V]B=\Phi.
\]

这里确实把离散自作用吸收到逐体素极化率，余下 \(\Gamma\) 对角为零；§II 原文说明辐射修正和三维立方体离散假设，不能把系数搬到当前二维 TM 网格。

**p3，Eq24 之后的未编号式**还明确给出：
\[
 (I-V\Gamma)^{-1}V=(I-V'\Gamma')^{-1}V',
 \qquad V'=PV,\quad
 \Gamma'=\Gamma-V^{-1}(I-P^{-1}),
\]
其中 P 可为任意可逆矩阵；原文另注明 V 的可逆性并非实际必要条件。若 D 是 \(\Gamma\) 的对角部分，原文选择 \(P=(I-VD)^{-1}\) 消去该对角。

这是**逆散射方法论文中的明确一般恒等式和自项实现**，比泛指 Born/T-matrix 传统更接近当前候选。一般矩阵变换可作为主线程检查块重求和的先例起点；本轮没有把其所有矩阵选择都认作已实现的 block-CIE。DCTMC 的反演变量/数据完成及局域性投影与当前辅助电流罚目标不同，恒等式也不证明二者优化效果等价。

## 4. Y₀：需与 self/邻接块抽取分开

Bevacqua–Isernia，2020 预印本/2021 TAP，A；[作者全文](https://arxiv.org/pdf/2006.13300)，[正式 DOI](https://doi.org/10.1109/TAP.2021.3060147)。§III.A–B，pp3–4，Eq5–10：

\[
 G_b=-j\frac{k_b^2}{4}J_0-\frac{k_b^2}{4}Y_0,\qquad
 W=\chi\widehat E_i+\chi A_i^{Y_0}W,\qquad
 \widehat E_i=E_i-j\frac{k_b^2}{4}F_{J_0}.
\]

\(F_{J_0}\) 从包围目标的测量曲线数据计算；剩余算子保留 Y₀ 核。§V，p6，Eq13 将该状态残差与数据残差放入 CSI 目标并优化。**这是利用观测数据计算辐射贡献的 Green 核分解，不是单个网格自项或四单元邻接块。**上述 OJAP 综述 §V.A，p8，Eq27–28 进一步实现 Y₀-NIE-CSI 混合，所以“将两个重写家族结合”也有实际先例。

## 5. 留给主线程的精确比较边界

| 项目 | 已核实先例 | 仍需核对的当前差异 |
|---|---|---|
| 对角 self-only | CS-EB 局部乘法抽取；空间 β 的 NIE/CIE；DCTMC 明确离散自项重整化 | 外部 β 分子及固定行权重、变量约束、目标归一化与优化坐标 |
| 真正多单元块 | DCTMC 给出一般可逆 P 的 T 不变恒等式 | 尚未核实真实 disjoint4-cell Green 块进入同一 CIE 罚目标的既有实现 |
| Y₀/CS-EB/NIE | 已有核分解、局部项抽取及 CSI 混合实现 | 当前块并无 Y₀ 所用的数据可计算辐射分量，不能仅凭“Green 分解”混同 |
| 实验解释 | 本任务没有读取/重跑当前试验 | self-only 对照的效果、块内回路是否必要、收益是否超过成本，由主线程判断 |

有限检索没有找到同式实现，不等于不存在；这里不给“first/new”的结论。所有旧 prior-art 产物保留。

## 检索与回执

本轮 5 条定向 web 查询；其中不完整数字 DOI 查询返回大量无关结果，已记作失败，没有据此取证。重点来自旧种子原文与其引用链。ScholarQA 单并发、无 snippets 的 collect 再次 HTTP429；verify 亦单次执行，具体状态见回执。按主线程已授权的 fallback，使用作者预印本、出版 XML 与 Crossref 身份核验；不隐称 Semantic Scholar 成功。

- `receipts/BLOCK_SELF_CIE_SEARCH_LOG.json`：检索式、证据层级、原文位置和未决项。
- `receipts/scholarqa_block_collect.json`、`scholarqa_block_verify.json`：实际返回。
- `receipts/block_downloads.json`、`crossref_block_*.json`：访问/身份回执。
- 完整第三方 PDF、提取文本和网页原始返回保存在忽略的 `local_only/`，不作为拟公开仓库正文复制。

