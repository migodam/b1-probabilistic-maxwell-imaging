# Track A：响应概率坐标、正质量饱和层与有限预算结果

本阶段完成了三个解析数学单元、3个冻结二维子空间、6个独立 prior-predictive parents、3档配对噪声的18份实际 TM24 观测及原坐标积分尝试。结论是 **局部响应积分在正常区域可以准确，但完整概率覆盖与同精度成本门没有通过，扩大保持 HOLD**。没有进入 full66，也没有将局部 inverse 成功当成全局生成密度。

最明确的失败机制来自实际 decoder：饱和平台包含正先验质量，其响应秩与过渡区不同。连续 rank1 chart 漏掉这些层，不能靠更小 root residual、失败节点重归一化或 identity POU 补回来。该机制属于 decoder/潜变量表述，不是新发现的 Maxwell 传播零空间。

## 目标与独立性

冻结活动坐标分别为 S0=(0,1)、S1=(0,64)、S2=(12,27)，其余64个 Gaussianized 坐标为公开零值。活动 prior 为 **独立标准 Gaussian 条件于 [-4,4]²**，归一化为 [Φ(4)−Φ(−4)]²。它们是新二维小模型，绝不是 full66 posterior 的边缘。

原 P1 `MaterialPrior`、24×24 VIE、8 sources×32 complex receivers、source-major Re32/Im32 和 σ=0.00119211798514519 保持不变。数据按照 dσ=f(x*)+σε 成对生成，σ比例为1、1/2、1/4；同 parent 三档共用一份 ε。18行对应6个独立 parent，不做把18当独立样本的显著性检验。没有读取 B 的锁定 truths。实际接口和成本均经过共享 `physics_adapter.py`。

## 解析数学单元：PASS

`math/A1_MATH_RESULTS.json` 含19组扫描、7个检查，另有3项积分与稳定尾部单元测试。

- 条件体积 f(a,b)=k(a)b：原坐标与响应坐标归一化和均值一致。σ=.125时 profile 的弱边缘 TV=.213906，删 determinant 的 TV=.216927，偏差并不随噪声消失。
- f(x)=x²：原坐标与包含奇异体积的双分支响应积分一致。y=0的厚度随 sqrt(σ)；仅保留一个对称分支的 TV=.5。
- 高条件数线性 Gaussian：σ=.001时，被错误补为 prior 的弱方向 TV=.921656，双向 KL 为2.9023和446.5982。线性 GN/RTO 精确等式仅是解析单元，不冒充非线性 Maxwell baseline。

## 实际 Maxwell：正常区与混合秩必须分开

### 正常的属性–几何 coupling

S1 的 regular parent P1 有可信的数值原坐标收敛趋势。原坐标使用全局 GL24/48 与 mode-adapted latent GH15/21，通过已知正 Gaussian partition 共同积分。响应积分则取属性方向为强坐标、场 offset 为弱坐标，保留完整 prior、全部512维残差、σ/|B| 和外层弱变量体积。

|σ比例|原坐标 Z 改变量|原坐标材料均值变化/SD|响应与原坐标材料均值差/SD|GN TV 数值估计|Laplace TV 数值估计|
|---|---:|---:|---:|---:|---:|
|1|.001051|.002167|.024948|.217107|.216357|
|1/2|.00002756|.00003641|.002981|.117843|.117717|
|1/4|.0000001716|.0000005341|.0003388|.059491|.059487|

在 σ/4，GN/Laplace 材料均值误差约.1075 SD，响应积分约.000339 SD。这里可以说响应积分解决了局部非 Gaussian 的材料矩误差；不能说它已经满足全局 posterior TV/branch/覆盖门，也没有证明总成本节约20%。Gaussian 密度指标是相对于收敛数值参考的积分估计；σ/4的 q 积分检查约.999973，仍不是带严格舍入包络的证书。

局部 QR-RTO 已实际实现：固定 augmented map h=[x;(f−d)/σ] 的 QR 基，真实 forward 求 inverse，付费两方向 JVP 计算 |det(QᵀDh)|。regular coupled 三档均8/8 roots 成功。成功集与分支选择的全局归一化仍未建立，故只报告 **LOCAL_NUMERICAL_RTO_ATTEMPT**，不把 ideal logq 当合格 unconditional density。

### 正质量 host 平台：rank0、rank1、rank2 三层

S1 中实际场为 .65 x0，温度为.5，故三个原坐标区域为：

|区域|decoder 结构|一般响应秩|处理|
|---|---|---:|---|
|x0≤−10/13|所有 ROI 为 host，属性也不可见|0|精确 Gaussian prior 体积乘一次常材料 forward|
|−10/13<x0<10/13|场与属性共同变化|通常2|二维原坐标参考及1强+1弱响应积分|
|x0≥10/13|纯 inclusion，场坐标为 gauge，属性可见|1|精确场 prior 体积乘属性1D积分|

下层 prior 质量为.22086048232212066，上层质量相同。host 材料方差为零并不表示其 latent 条件方差为零。V3 已在原坐标 summary 中保留层内 latent 方差；首个修复前 response covariance 单独保持不合格，未用于主结论。

旧 mode-partition 参考在 host parent 上仍未收敛，不能引用其 .70–.78 host posterior mass。零新增物理的直接 GL24/48 readback 保留精确 host 层，得到更稳定结果：

|σ比例|host posterior mass 数值估计|Z 改变量|材料均值变化/SD|
|---|---:|---:|---:|
|1|.926520|.0001706|.01706|
|1/2|.941277|.0012468|.02142|
|1/4|.953144|.0020950|.00521|

前两档仍未达到1% SD材料均值门。该 readback 的 pure 层在现有 float64 求积中 underflow；这里没有将其升级为严格漏质量上界。数值趋势支持“噪声缩小时 host 原子质量增加”，不能声称精确 branch certificate。

在接缝附近 H(−1+ε)=1.25 ε³+O(ε⁴)，所以一般正则物理方向的 seam 厚度为 σ^(1/3)，响应体积有 s^(−2/3) 型奇性。主线程另做了已计费 Maxwell action probe：实际小 h response slope=2.996758、材料方向 Jχ norm≈1.777、latent J 在 host seam 为零。该证据把 decoder 退化与非零物理材料作用分开；参见上级 `audit/seam_maxwell/RESULTS.json`。

### 高阶场方向与正常 rank2 chart

S2=(12,27) 的6份观测均完成原坐标9/13 cubature、GN/LIS、observed-Hessian Laplace 及每观测4个实际 QR-RTO local root 尝试。mode ROI 中确有 host/pure 平台单元，例如末个 case 两者各.06640625，过渡区.8671875；这只叫 saturation exposure，没有宣称发现物理 fold。

其中 parent P0 的三档噪声完成独立的 rank2 response chart，3×3→5×5求积。σ/4时原/响应 Z 相差4.82e−6，材料均值相差1.59e−6 SD；该 chart 没有弱变量。P1 的额外 rank2 求积触发自设32680状态软停止，partial 文件保留；其原坐标与 baselines 已完成。S2 最后一档的 GN/Laplace TV 数值估计约.00434/.00431，说明简单 Gaussian 在该区域已经很强。

S0 的 saturation parent P0 三档原坐标积分均未收敛；P1 的三档数值矩较稳定。全部状态逐行列于 `OBSERVATION_SUMMARY.csv`，不会用正常区结果覆盖失败区域。

## 概率门与 baseline 合法性

精确公式不等于积分证书；离散求积测度与连续 posterior 直接比较的 TV 本身为1。本阶段只将 TV/KL 用在可归一化连续 Gaussian 近似密度，且限定到参考收敛与 q 积分检查合格的 case。旧参考下出现负 KL、TV>1 或 q 积分偏离的数值全部视为 **INVALID_REFERENCE_METRIC**，不作方法胜负依据。

Chart 的归一化、材料矩和分支质量是数值积分比较。有限 inverse 残差或有限 B 样本不提供整个域的单射与覆盖，也没有从样本最大最小权重构造全域 Δ。主线程 A5 的两个中心×三档噪声充分条件全部未过，噪声半径比约21.77–2344.73；当前全局方向 norm bound 过松，不等于反证局部响应方法。见上级 `audit/A5_INDEPENDENT_REVIEW.json`。

Ordinary Spline、manifold lifting、full66 NLL/SBC/Gamma、独立多链 mixing 均 **NOT_RUN**。RTO 做了实际局部求逆，不声称合格全局密度。没有同可信精度、全覆盖和全 baseline 的成本冠军。

## 成本与停止

共享累计账本（包含 qualification、父线程曲率/接缝 probe、失败和被中断工作）：32705个已保留状态预算，261632 forward RHS，30152 derivative/adjoint RHS，总291784 RHS；CPU315.360秒、受控 wall374.924秒；2次 interruption/failure保留，unfinished actions=0。状态预留有一次在 actual RHS 前中断，因此 states×8 与实际 forward RHS 差8，不回退该预留。上限32768 states、262144 forward RHS、40960 derivative RHS、24小时均未突破，63 states余量不再使用。

`COST_BY_JOB.json` 将 V1/V2/V3、baseline、局部 rank2 与父线程 probe 分账。V1=7685 states，V2=7499，V3=14544；失败费用占比很高，不能只报最后成功节点的成本。真实物理始终一个 CPU worker、BLAS1。没有保存巨型 LU endpoint cache。数组与缓存设计低于1.5GB；OS 峰值 RSS 未单独测量，不写成已测验收。

## 最终门

- 数学归一化、矩与解析分支单元：PASS。
- 所有18份实际数据与原坐标推断尝试：已完成；参考质量逐 case 区分。
- 正常区域 response 积分数值一致性：有通过案例，也有未达到1% SD的案例。
- 混合秩完整质量/尾部、重要 branch certificate、全局 TV与双向 KL：HOLD。
- 同可信精度、完整 setup/失败成本下优于最强 baseline 20%：未建立。
- full66、扩大 campaign：NOT_RUN / HOLD。

值得保留的结果是实际 decoder 的混合秩反例，以及正常 coupling 区中响应积分相对 Gaussian 近似的明确精度差。下一步科学判断由主线程统一审查；本阶段不追加预算、不改 prior、不隐藏旧失败。

Decoder-only geometry、同谱随机几何的独立对照本阶段NOT_RUN，不将三个预声明空间当作这些对照。独立SBC/Gamma与joint proper-score大样本检验NOT_RUN；六个parent仅作机制验证，不报告显著性胜出。
