# 共享白化、材料链式导数与原噪声标定的独立只读审查

结论：**PASS（已读版本的代数与排列一致性）**。本审查未修改共享文件、未调用正演、未读取真值数组。接口核验中已付费的单个真实支持出生状态另计 80 RHS；它确认两种精确方程写法的接口一致，不是科学成败判据。

## 源优先排列

[`shared/core.py`](../shared/core.py) 第 37—41 行中，`pack(F)` 将 `F.real.T` 与 `F.imag.T` 在接收维连接，所以每个源依次输出 32 个实部与 32 个虚部，再进入下一个源。`white_jac(J)` 把复 Jacobian 的 256 行解释成 8×32，并同样逐源连接实部与虚部。二者一致。

`ExactSolver` 第 67 行的 `einsum('ri,is->sri', ...)` 明确令源在最外层，故其 `(256,256)` 占据率 Jacobian 与上述排列一致。第 68 行 `J_delta=(...).T.reshape(-1)` 也按源优先排列。

## 白化与原噪声

原 [`TARGET.json`](../../p3l_r4_conditional_transport/configs/TARGET.json) 第 47—55 行冻结圆对称复高斯噪声、固定绝对 \(\sigma=0.00119211798514519\)、实虚各 \(\sigma^2/2\)，并指定逐源 `sqrt(2)/sigma Re32/Im32` 白化。`core.white` 与之相同。

[`shared/prepare_digits.py`](../shared/prepare_digits.py) 第 43 行用 \(\eta=(\sigma/\sqrt2)(\xi+i\zeta)\)，\(\xi,\zeta\) 为独立实标准正态，因此每个实分量的方差恰为 \(\sigma^2/2\)。这是沿用原绝对仪器噪声，不是按数字重新归一 SNR。

`core.phi` 第 43 行采用 \(\Phi=\|F-y\|_F^2/\sigma^2\)。于是
\[
\Phi=\tfrac12\|W_\sigma(F-y)\|_2^2.
\]
`invert_digits.py` 的 `.5*r@r` 数据目标与此一致，不存在遗漏或重复的 2。256 个复观测对应 512 个白化实观测；噪声真均值下 \(2\Phi\) 为 512 自由度卡方，而 \(\mathbb E\Phi=256\)。这只说明标定，不保证每个数字观测均落在给定分位数门内。

## 精确物理 Jacobian

`core.ExactSolver` 第 63—68 行使用电流方程
\[
M=I-\delta AC,\quad h=M^{-1}\delta AE,\quad e=E+Ch,\quad T=OM^{-1}.
\]
其全微分为 \(dF=T\,dD\,e\)，所以
\[
\partial_{a_i}F_{rs}=\delta T_{ri}e_{is},\qquad
F_\delta=T(Ae).
\]
这正是两个返回量的实现；接收端求解使用转置而非共轭转置以构造解析复导数，随后才转换为实参数白化 Jacobian。

## N / P 材料链式导数

[`shared/invert_digits.py`](../shared/invert_digits.py) 第 19—30 行中，属性映射为
\[
\delta=.8+\mathcal N(x_\alpha)+i[.03+.15\mathcal N(x_\beta)].
\]
因此两列导数是 \(F_\delta\varphi(x_\alpha)\) 和 \(i\,.15 F_\delta\varphi(x_\beta)\)，实现中正是这两项，没有遗漏实部区间长度 1 或虚部区间长度 .15。

N 模型直接使用占据率列；P 模型的
\[
a=s\!\left(\frac{\mu+A_{ROI}z}{\tau}\right),\qquad
\frac{\partial a}{\partial z}=s'\!\left(\frac{\mu+A_{ROI}z}{\tau}\right)A_{ROI}/\tau
\]
与原 `src/p1_material/prior.py` 的 quintic decoder 和温度缩放一致。`Problem.jac` 先乘该实链式矩阵，再附加属性两列，最后仅白化一次。P 返回 `(512,66)`，N 返回 `(512,258)`。

## 本审查的范围

本结论不验证优化是否达到全局最优，不解释数字重建质量，也不把近似输出残差界升级为 Jacobian 误差界。当前受控求解器在请求 Jacobian 时使用精确场和精确导数，并实际支付 40 RHS256；该费用必须保留在反演总账。
