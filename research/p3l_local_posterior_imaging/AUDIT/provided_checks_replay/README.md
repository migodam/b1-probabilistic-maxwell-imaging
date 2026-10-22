# P3-L 独立数学核验：0 Maxwell

本目录是报告的可复现数学控制实验。它核验概率定义、积分公式和接受判据，**没有运行 Maxwell/PDE forward、adjoint、材料反演或付费物理实验**。这里的正结果不能作为 Maxwell 后验已校准、材料风险已改善或求解器已加速的证据。

## 文件与复现

- `run_math_checks.py`：完整计算与作图脚本；执行 `python run_math_checks.py`。
- `results.json`：四类核验的完整数值、积分误差估计、浮点容差、软件版本及随机种子。
- `summary.txt`：机器生成的简短结果。
- `01_profile_boundary.png`、`.pdf`：强方向体积和活跃边界的反例。
- `02_certification_prior.png`、`.pdf`：MH 区间核验和 P1 decoder 高斯矩。

本次环境：Python 3.12.14，NumPy 2.3.5，SciPy 1.17.0，Matplotlib 3.10.8。随机种子固定为 20261019；图中采用英文标签。脚本不访问网络、不导入 Maxwell 求解器。`quad` 的误差是数值积分器报告的估计，Taylor 余项是解析上界；本实验使用普通浮点运算，未声称完成严格区间算术证明。

## 1. 条件优化 ridge 不能替代强方向边缘化

在整个实平面上取

\[
U(a,b)=\frac12 a^2+\frac12 e^{2ka}(b-a^2)^2,\qquad k=0.8.
\]

强方向条件极小值为 \(b^*(a)=a^2\)，所以把 \(e^{-U(a,b^*(a))}\) 当作关于 \(da\) 的密度再归一化，会得到 \(N(0,1)\)。但真正积分结果为

\[
\int_{\mathbb R}e^{-U(a,b)}\,db
=\sqrt{2\pi}\,e^{-a^2/2-ka},
\]

因此 \(a\sim N(-k,1)\)，且

\[
Z=2\pi e^{k^2/2},\qquad
\mathbb E[a]=-0.8,\qquad
\mathbb E[b]=\mathbb E[a^2]=1+k^2=1.64.
\]

用 profile 分布平均 ridge 则给出 \(\mathbb E[b]=1\)。这是对强方向条件体积的遗漏，不是候选发现数量不足。由于本例强条件正好是 Gaussian，Laplace 的 \(H_{bb}^{-1/2}=e^{-ka}\) 因子在这里完全精确。

数值积分结果：\(Z=8.652748934984757\)，\(\mathbb E[a]=-0.8000000000000002\)，\(\operatorname{Var}(a)=0.9999999999999999\)，\(\mathbb E[b]=1.6400000000000001\)。另外在 \(a=-3,-1,0,1,3\) 五个位置直接对原始 \(b\) 坐标积分，均与闭式表达一致。

## 2. 活跃边界处必须保留线性项和截断质量

在 \(b\in[\ell,u]\) 上考虑 \(\exp(-gb-Hb^2/2)\)，\(H>0\)。令

\[
\mu=-g/H,\quad s=H^{-1/2},\quad
\alpha=(\ell-\mu)/s,\quad\beta=(u-\mu)/s,\quad
M=\Phi(\beta)-\Phi(\alpha).
\]

完成平方给出

\[
Z=\sqrt{2\pi/H}\,e^{g^2/(2H)}M,
\]

\[
\mathbb E[b]=\mu+s\frac{\varphi(\alpha)-\varphi(\beta)}{M},
\]

\[
\operatorname{Var}(b)=s^2\left[
1+\frac{\alpha\varphi(\alpha)-\beta\varphi(\beta)}{M}
-\left(\frac{\varphi(\alpha)-\varphi(\beta)}{M}\right)^2
\right].
\]

无限端点处按极限解释 \(\alpha\varphi(\alpha)\) 等项。脚本核验下边界、上边界和两个有限区间共四例。第一例 \(g=1.25,H=3,b\ge0\) 的 constrained mode 是 0，梯度 \(g\) 非零；正确保留 Gaussian mass 为 0.2352432110，积分为 0.4417168099，条件均值为 0.3379646399。若直接丢弃线性项，积分被高估约 63.82%。这说明在 native uniform prior 的活跃边界附近，不能无条件套用以驻点为中心、无限支撑的 Gaussian Laplace 公式。

## 3. 同随机数的 MH 区间判据

这是一个明确标记为 `ANALYTIC_TOY_CONTROL_NOT_MAXWELL` 的一维控制。设

\[
f(x)=0.15x^2+\sin x,\quad y=0.25,\quad\sigma=0.6,
\quad p_0=N(0,1),\quad x'=x+0.8\xi.
\]

对 \(\sin x\) 采用奇数阶 \(d\) 的 Taylor 多项式，解析上界为

\[
|\sin x-P_d(x)|\le |x|^{d+2}/(d+2)!.
\]

因为随后那个偶数阶系数为零，这是 Taylor 公式给出的全实轴有效余项界。设白化预测误差界为 \(\epsilon\)，近似 residual 为 \(\widetilde r\)，则真实 Gaussian misfit 位于

\[
L=\tfrac12\max(0,|\widetilde r|-\epsilon)^2,
\qquad U=\tfrac12(|\widetilde r|+\epsilon)^2.
\]

若 prior 与提案比的精确对数项为 \(c\)，则真实 MH log ratio 被夹在

\[
\ell_-=c+L(x)-U(x'),\qquad
\ell_+=c+U(x)-L(x').
\]

沿用同一个均匀随机数，\(\log u<\ell_-\) 接受，\(\log u>\ell_+\) 拒绝，否则细化；各层区间取交集。当上界有效且提案相同，这一判据给出和精确 MH 相同的接受决定。对给定 \((x,x')\)，初次必须细化的条件概率为

\[
e^{\min(\ell_+,0)}-e^{\min(\ell_-,0)}
\le\min(1,\ell_+-\ell_-).
\]

脚本先运行精确参考链，再在其 50,000 组共享提案与均匀随机数上验证上述判据。**参考真值在本次控制中已经计算过**；报告中的 fallback 次数是区间算法的逻辑计数，不能理解为实测物理计算成本。

| 初始 Taylor 阶数 | 实测初始细化比例 | 平均条件细化概率 | 平均区间宽度概率上界 | 接受决定差异 |
|---:|---:|---:|---:|---:|
| 1 | 0.217800 | 0.219847 | 0.545443 | 0 |
| 3 | 0.055200 | 0.055049 | 0.195428 | 0 |
| 5 | 0.007300 | 0.007314 | 0.050444 | 0 |
| 7 | 0.000640 | 0.000687 | 0.015035 | 0 |

使用阶数 1、3、5、7、9、13 逐步细化，本次所有判据在最终层已决定，无逻辑 fallback。最大浮点区间包含偏差为 \(8.88\times10^{-16}\)，低于脚本采用的 \(10^{-10}\) 容差。控制不认证采样收敛、ESS、Maxwell surrogate 界、后验校准或加速比。

## 4. P1 smootherstep 高斯矩与原 prior 精确均值

采用 P1 真正的 decoder：

\[
\phi=\mathbf1_{\rm ROI}h(f/0.5),\quad
h(s)=6t^5-15t^4+10t^3,\quad t=\operatorname{clip}((s+1)/2,0,1).
\]

对 \(f\sim N(\mu,s_f^2)\)，内部区间上的 \(h^p\) 是次数 \(5p\) 的多项式。脚本用截断 Gaussian 矩递推计算 \(p=1,2\) 的解析表达，并用独立一维 `quad` 检验：

\[
I_0=\Phi(\beta)-\Phi(\alpha),\quad
I_1=\varphi(\alpha)-\varphi(\beta),
\]

\[
I_j=\alpha^{j-1}\varphi(\alpha)-\beta^{j-1}\varphi(\beta)+(j-1)I_{j-2}.
\]

然后将 \(f=\mu+s_fz\) 展开，并加入上饱和区的概率质量。五组参数是数学控制中的选定 Gaussian marginals，不是新的 Maxwell 对象。

| \(\mu\) | \(s_f\) | \(\mathbb E\phi\) | \(\mathbb E\phi^2\)（quadrature） |
|---:|---:|---:|---:|
| 0 | 0.15 | 0.5 | 0.307318946207440 |
| 0 | 0.50 | 0.5 | 0.417808135003471 |
| 0 | 1.00 | 0.5 | 0.457369171224276 |
| 0.25 | 0.35 | 0.734315743069388 | 0.643508047681417 |
| −0.40 | 0.60 | 0.262557724905460 | 0.206060478139833 |

最大解析式与数值积分绝对差约 \(1.11\times10^{-12}\)。对 P1 原 prior，各 cell 的 Gaussian field 都是零均值；由 \(h(-s)=1-h(s)\) 得 ROI 内 \(\mathbb E\phi=1/2\)，与各 cell 的方差无关。已知

\[
\chi=0.1+0.04i+\phi(a+ib),\quad
a\sim U[0.8,1.8],\quad b\sim U[0.03,0.18],
\]

且材料幅度独立于 Gaussian field，因此

\[
\boxed{\mathbb E_{p_0}\chi(x)=0.75+0.0925i\quad(x\in\mathrm{ROI})},
\]

ROI 外均值为背景 \(0.1+0.04i\)。对应 ROI 内 \(\mathbb E\epsilon_r=1.75+0.0925i\)。这是可精确实现的 **zero-physics prior-mean 基线**，没有条件于观测，不能称为 posterior mean。若按 P1 的 ROI 256 cells、总网格 576 cells 取空间平均，其值为 \(0.3888888889+0.0633333333i\)。P1 参数来源是本次已读项目证据 `evidence/p1/MODEL_AND_PRIOR.md`；脚本把这些常数显式写出，不依赖该证据文件才能运行。

## 报告图注建议

**图 1：** 两个可解析反例展示概率质量不能由 ridge 似然代替。A：\(U(a,b)=a^2/2+e^{1.6a}(b-a^2)^2/2\) 的 normalized profile 是 \(N(0,1)\)，正确弱边缘是 \(N(-0.8,1)\)。B：\(e^{-1.25b-1.5b^2}\mathbf1_{b\ge0}\) 的边界 mode 为 0，但均值为 0.338；完成平方后的 Gaussian 仅有 0.2352 概率质量保留在支撑内。

**图 2：** 纯数学控制，不含 Maxwell 评估。C：一维 sine forward 的 Taylor 误差界产生 likelihood 区间；50,000 个共享提案与随机数上的接受判据与精确参考一致，初始细化率随余项缩小而降低。D：P1 饱和 decoder 的 Gaussian 一阶、二阶矩与独立数值积分一致；零均值 Gaussian field 下的 \(\mathbb E\phi=1/2\) 是严格对称性结论。
