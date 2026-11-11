# 新增 consensus-reflection toy 独立资格审查

审查对象：`consensus_reflection_toy.py`、`toy_results.json`及四份已保存`.npz`。审查时间：2026-11-08。没有重跑3000步迭代或任何反演优化。只针对“gradient sign / projector / DR sign错误造成假胜出”的具体风险做有限代码资格检查。

**结论：未发现导致当前正例的符号、物理接线或投影错误。** 原物理gradient与独立finite difference一致，两个投影满足对应最优化条件，DR更新是正确的两集合反射。环/开缝壳的空气解是原数据目标的数值支持严格box角点局部极小，direct到空气不由gradient反号造成。审查当时的主脚本仅展示相同初始化下避免direct/AP的不良终点。

**后续补充（同日，主审查完成之后）：** followup_controls.py 已从实际direct空气终点的真实双场重启DR，环/壳相对材料L2分别0.077874和0.172045，相对物理数据残差0.009301和0.009900。该事后机制控制及全部费用保存在toy_controls.json，支持从这两个有限模型坏终点逃逸；不改变原审查的物理资格边界，也不构成等CPU或CIE比较。

本审查新增44次100阶Maxwell LU、336个复RHS；另有固定投影矩阵setup和8次A投影，均为资格检查，不能计入原算法性能。可重现脚本及收据为`toy_qualification_audit.py`、`toy_qualification_audit.json`。

## 1. 物理算子：接线一致，连续精度未认证

采用e^{-iωt}、outgoing H_0^(1)、k=2π，平面波Einc=exp(+ik x·d)。格距h=.1λ，等面积圆盘半径a=h/√π，常数源的圆盘积分给出

G_ij = (iπka/2) J_1(ka) H_0^(1)(kR_ij)，i≠j；

G_ii = (iπka/2) H_1^(1)(ka)−1。

自项中的−1来自H_1原点极限，不应删除。接收矩阵S采用相同体积分因子，故与J=mQ和F=SJ一致，不缺一个k²或cell area。物理求解L=I−diag(m)G，右端diag(m)Einc，矩阵顺序正确。

该模型是等面积disk-pulse/collocation有限VIE。使用同一G/S生成与反演数据属于matched finite model；没有独立网格或连续Helmholtz参考，也不是三维矢量Maxwell资格。面积等价圆盘是网格单元近似，不能把它当精确连续材料几何。现有正例足以用于机制可行性，不能替代论文物理验证。

本次软件检查G^T=G；Hermitian虚部最小特征值约−4.84e−16，处于舍入尺度。这仅排查明显符号/被动性错误，不升级为连续物理证书。

## 2. 实材料物理梯度：符号与共轭正确

令L=I−D_mG，J=L^{-1}D_mEinc，Q=Einc+GJ，r=SJ−Y，f=.5||r||²/||Y||²。则

dJ[h]=L^{-1}D_hQ，λ=L^{-*}S*r，

g_i = Re Σ_s conj(λ_is)Q_is / ||Y||²。

代码`lu_solve(..., trans=2)`求L^{-*}，随后`real(sum(conj(adj)*Q))`完全符合该式，**正号正确**。

四对象各在m=.1及一个不读取GT的随机内点检查一随机实方向，central step 1e−5。八个方向相对误差最大2.41e−9，没有整体sign反转或复共轭错误。

空气角点m=0的原目标f=.5：

|对象|所有材料坐标gradient最小值|最大值|box角点解释|
|---|---:|---:|---|
|环|.0226127|.0390750|全坐标严格正，原目标在合法非负方向一阶上升|
|开缝壳|.00221986|.0131974|全坐标严格正，同上|
|桥|−.0105822|.00560700|存在合法下降方向，空气不是该类KKT极小|
|低对比圆盘|−.0129109|−.00394638|明显能从空气沿非负方向下降|

对可微f，若∇f(0)全部严格正，则由Taylor一阶余项，m=0是box [0,12]^N内严格局部极小；不需Hessian正定。这里数值gradient与单侧finite difference相符，属于有清晰解析判据支持的数值证据，而非向外舍入的严格证书。环/壳direct收敛到空气因此是真实的优化几何现象。

## 3. A集合的graph＋data ball投影正确

代码采用X=J/s、s=6，B=sG、T=sS，等权product norm ||ΔX||²+||ΔQ||²。Aδ={Q=Einc+BX，||TX−Y||≤δ}。

固定graph消去Q，H=I+B*B>0，unconstrained minimizer

X_u=H^{-1}[X_0+B*(Q_0−Einc)]。

若其数据残差超过δ，非负Lagrange multiplier μ给出

X_μ=X_u−μH^{-1}T*(I+μTH^{-1}T*)^{-1}(TX_u−Y)。

代码通过固定port Gram K=TH^{-1}T*的Hermitian eigendecomposition求全source共同μ，正是Frobenius球投影；并未错误地每源分开设δ。球约束是凸的，满足feasibility、μ≥0、complementarity及stationarity即可确定全局最近点。

四对象各检查物理初始点和一个不相容随机product-space点。8次检查中：

- 实际||TX−Y||/δ与1偏差最大约7.1e−14；
- graph等式缺陷约1.2e−15以内；
- 从输出独立恢复μ后的relative KKT stationarity最大4.94e−12。

所以当前正例不来自球投影放松过大、错误白化或错误ball radius。

δ定义也正确：CN噪声实/虚方差各σ²/2，2||ε||²/σ²服从χ²_(2Ndata)，代码δ=σ√(χ²_0.99/2)。σ是模拟中公开已知的噪声尺度；真实实测中必须另由噪声/标定估计。

## 4. C集合的联合材料投影正确

Cbox={X_i=α_iQ_i，α_i=m_i/s∈[0,2]}，每cell单一实α共享6个复source分量。固定α最近Q为

Q(α)=(Q_0+αX_0)/(1+α²)，X(α)=αQ(α)。

条件距离为(aa−2ccα+bbα²)/(1+α²)，aa=||X_0||²、bb=||Q_0||²、cc=Re<X_0,Q_0>。驻点方程

cc α²+(bb−aa)α−cc=0。

代码的两个根、两个box端点及合法性过滤完全一致。cc=0时端点足够，包括常数目标、零场等退化情况；没有漏掉会改变最小值的必要内点。候选通过真实投影距离排名，不读材料真值。

独立用2×2实对称Rayleigh矩阵[[aa,−cc],[−cc,bb]]的最小eigenvector和box端点作reference，对100组多source rows（含零X、零Q、二者全零、纯虚相位、极小cross和反向cross）检查：最大相对cost差4.98e−16，membership residual3.40e−16。

小的工程局限：直接二次公式在|cc|极小时可能发生一根消去误差，建议正式实现用stable quadratic root/eigenvector版本。当前检查没有发现其影响实验结论；不要为了此可逆小风险重跑整套toy。

## 5. DR更新顺序/反射符号正确

实际顺序是C先、A后：

c=P_C z，a=P_A(2c−z)，z^+=z+a−c。

因此T_DR=I+P_A R_C−P_C=.5(I+R_A R_C)，代码`X+=ax−cx; Q+=aq−cq`没有反号。若原报告写A先/C后，应明确只是交换两集合的另一变体，不要将实现顺序抄错。

当前输出材料来自**C shadow** cm，因此每个被检查候选都满足box和共享实材料。A shadow负责graph＋data ball。gap=||a−c||正是完整product-space固定点残差||z^+−z||。准确固定点给c=a∈A∩C。AP则可停在非零gap的最近点对，不能据AP材料已不动断言物理相容。

DR中间的product-space变量不是实际物理解；每50步完整forward审计是必要的。本代码做到了，没有把辅助data fit当物理fit。它没有记录全迭代状态，桥的非收敛可以报告，但不能据60个标量采样点命名为周期/混沌。

## 6. 候选选择没有利用GT，但实验不是盲测

`best`仅按`Physics.forward(cm)`的原数据目标ff更新；合法根由projection距离选择；所有停止上限、scale和box固定。四对象的保存输出material重新做forward，与保存history及初始点中的最小数据目标完全一致。

不过truth确实传进`run_case`，并在每50步history内计算材料误差，用于诊断。当前代码没有任何使用这个值进行更新/排序/停止的路径，所以可以说“算法更新和选图不使用GT”，**不能说在线进程完全未访问GT或预测先封存再评分**。这是四个已知形状的探索性机制试验，不是新盲测结果。单凭当前source和JSON也不能独立认证“在第一次运行前已冻结所有参数”；需保留原执行日志/commit才可作更强协议声明。

## 7. 成本口径必须纠正解释

JSON已正确标注no equal-cost performance claim。AP和DR各同3000次A、3000次C、60次每50步审计，加一次物理初始化，共61 Maxwell LU/366复RHS。seconds计时包括Projectors构造和全部投影，因而两者的实际CPU量级可比较。

但`lu_count/rhs_count`只计原Maxwell系统：每个投影臂还必须单列

- 一次100阶固定H的Cholesky；
- H^{-1}T*的24个setup RHS；
- 3000×6=18000个固定H的迭代RHS；
- 一次24阶Hermitian eigendecomposition；
- 每轮G/G*及端口矩阵作用、求μ标量root等CPU工作。

这些已在seconds内，并不是未付算力，但不能把“366 RHS”解释成全部线性求解费用，也不能据其远少于direct RHS宣称加速。main中的G/S公共setup、truth生成、保存时间不在每方法seconds，属toy公共开销。

direct给600-step上限但环/壳2–3步便到KKT，实际数毫秒；不是与DR等实际CPU或等实际总线性代数预算的全局对照。正确说法是“同一初值下普通局部优化到严格错误角点，而反射迭代在约2秒的额外搜索中找到数据合格候选”。论文必须补成熟CIE及同总预算多起点/direct/CSI/WRI等控制。

200步polish只施于AP/DR各自所选起点，同一函数/上限并单独收费。实际迭代数不同，给定同cap不等实际费用。对于direct自身已终止角点，直接续跑仍会保持角点，这不自动使其成为等预算非局部搜索控制。

## 8. 正例、失败及指标解释

|对象|DR选中材料相对L2|DR原物理相对data residual|CN99相对界|解释|
|---|---:|---:|---:|---|
|环|.0870|.00939|.01098|明显材料改善并进噪声球|
|开缝壳|.17775|.00994|.01097|明显改善且进噪声球，但材料仍有约17.8%相对误差|
|桥|.60491|.18760|.01097|仍失败；最终gap大，不能称已收敛物理解|
|低对比圆盘|.12863|.00945|.01097|与direct .12506材料误差接近；无独立优越性|

此处材料指标是||m−mGT||/||mGT||，不是R11的amplitude-normalized RMSE，不能照搬R11 .15联合门。没有IoU/拓扑资格，因此“开缝壳成功恢复”的强表述应谨慎，用“显著改善且噪声可行”更准确。

桥polish后data residual降至.03473，但材料误差从.6049恶化到.7200，且仍在噪声球外。这保留了用户原研究的核心反例：更低残差不自动给更好材料。环/壳polish改善分别至.05156/.16121，是额外付费局部修正的结果，不能全部归给DR。

所有任务n=100、6平面波、24全孔径接收、同网格、同noise draw/方法、无材料prior；未比较成熟CIE/G19，未验证cross-grid/3D/实测。结果支持一个有明确机制、可进一步研究的候选结构，尚不支持通用算法或论文性能结论。

