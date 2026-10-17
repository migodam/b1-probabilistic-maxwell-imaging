# A30 E0 / E1：材料覆盖、测量选择与合法补偿

2026-10-11。全部证据来自已暴露 development；A29 原冻结 NO-GO 保持不变。

## 结论先行

**E0通过原协议的离线容量门槛。E1未获得进入E2的资格，物理实验停止；E2没有执行。**

完整候选 convex hull 在12个复杂parent上相对N3改善 **13.12%**，探索性95% parent区间 **[11.51%, 14.65%]**。相同bank的真实L-weighted mean却恶化 **3.32%**。这把“候选不存在”与“已有测量权重找不到候选”分开了。

64个实参数的固定支持多尺度稀疏表示达到 **35.69%** 的离线改善，TRAIN-PCA单chart达到 **21.91%**。这些结果使用truth拟合参数，只证明所实施表示的容量，**没有实现相应在线重建器**。

修正后的合法补偿测试获得29/32个候选对。局部补偿确实缩短距离，但最小实际full-F白化距离仍为 **8.84**；复杂条件最小值 **9.18**。没有找到原协议要求的“nominal≥3而compensated≤1”的parent。保守uniform误差界也没有提供有效判别。该结果不排除更大合法集合、其他结构编辑或更好的优化器能产生补偿；它否定的是**用当前已实施局部见证继续推进该路线的资格**。

## 1. 物理问题与历史背景

材料图像与散射数据的关系不是一一对应的线性读出。不同几何和复contrast可能改变内部场，也可能互相补偿。判断一个成像方法失败，需要分别询问：表示里有没有更好的材料；仅依据测量能否选中它；合法结构补偿是否使数据本来就难以区分。

A29原硬候选选择的负结果并不自动否定同一bank的所有连续混合。A30首先零新增全波检查这一点，然后只在容量存在后支付有限的Maxwell补偿验证。研究来源为[A29固定commit](https://github.com/migodam/a29-material-hypothesis-coverage/tree/5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4)、[A28 E1](https://github.com/migodam/a28-maxwell-material-information/blob/c7819fa142bb70f9920835e8e91ba145f15ed8e3/E1_DECISION.md)及随仓库保留的[A30最小协议](theory/package/A30_MINIMAL_PROTOCOL.md)。没有重做理论，也没有把解析恒等式当作成像优越性。

## 2. 数据、材料尺度和信息边界

从[A29 Release](https://github.com/migodam/a29-material-hypothesis-coverage/releases/tag/a29-v0.1.0)取得真实原始数组及复现输入。E0完整包含250065–250080共16个parent，每个weak/ring两条件，合计32；shell、asymmetric_piecewise、mixed_scale共12个复杂parent，4个smooth另列。TRAIN为250001–250064，只用于表示、模板、范围和离散选择规则。

冻结后端是 **24×24二维标量TM**、8源/32接收器、逐源Re/Im packing、complex128/float64、cell mass=1/576。材料约束Reχ∈[−0.5,4]、Imχ∈[0,0.5]。这些是历史二维平台的真实全波结果，**不是3D Vector-Maxwell证据**。标签具有固定复增量相位，不能据此宣称任意独立实虚分布泛化。

在线候选、chart/support及补偿端点不读取development truth。truth仅进入独立offline容量拟合与终点评分。没有读取locked ID/OOD，没有新标签，没有训练网络。测量条件、背景、噪声尺度和已有N3来自同一历史资产；没有用truth选择最有利的洞或补偿方向。

## 3. E0：完整候选凸包与真实测量权重

每个候选池包含58–91个合法材料；小池12个包含N3 no-op。凸组合仍在共同材料box内。full/small均完成32条件优化；独立直接梯度计算的最大FW gap为 **6.33e-09**，满足原1e−8门槛。它是浮点数值收敛证据，不是区间算术证明。

|方法|复杂parent相对N3改善|探索性95%区间|
|---|---:|---:|
|full_hard|2.15%|[1.31%, 3.06%]|
|small_hard|2.09%|[1.22%, 3.02%]|
|full_hull|13.12%|[11.51%, 14.65%]|
|small_hull|7.64%|[6.39%, 8.95%]|
|full_L_mean|-3.32%|[-5.70%, -1.08%]|
|small_L_mean|-3.32%|[-5.70%, -1.08%]|

先平均每parent的weak/ring误差，再计算相对改善；12个复杂parent等权，2000次parent bootstrap，seed291011。不能将32条件或噪声重复当作32独立对象。

full hull的shell/asymmetric/mixed平均改善分别为13.68%、13.47%、12.21%。edge F1平均增加 **0.0607**，support IoU增加 **0.0544**；components/holes计数差的区间跨零，没有据此宣称拓扑恢复已认证。0.15/0.25/0.35阈值的原始数字均保留。

因此此次容量信号没有靠总体牺牲edge/IoU换取。smooth full hull改善很大，但不进入复杂主终点。实际full/small L-weighted mean的ESS中位数均为1，基本退化成hard decision；它并没有有效利用convex组合的容量。真实L权重依赖同一data-conditioned候选bank，称weighted estimator，**不称校准Bayesian posterior**。这不是对所有线性权重估计器的否定。oracle权重不能用于在线算法或作为NN标签偷偷延续实验。

## 4. E0：同16/32/64实参数的连续表示

所有离散起点、支持、TRAIN模板在读取评分truth前确定。每case/表示至多8起点，每起点最多100材料目标迭代；32×4×3=384个请求全部获得精确box可行的已达到拟合。参数预算计入几何以及每区域独立实/虚contrast。

|表示|16实参数|32实参数|64实参数|
|---|---:|---:|---:|
|atlas|-11.82%|-16.98%|-36.40%|
|level_set|-9.04%|-8.00%|-7.69%|
|sparse|-6.98%|11.68%|35.69%|
|single_chart|-10.41%|3.91%|21.91%|

稀疏表示为24→12→6→3正交Haar加3×3粗层DCT；支持来自N3/TRAIN，随后offline拟合系数。另存的truth-selected best-s tail只是无约束上限诊断，**没有替代实际固定支持拟合**。单chart是TRAIN-only非中心化PCA。192个选中的affine固定chart拟合具有数值最优性检查；有限起点之间的truth择优仍是offline操作。

几何atlas使用TRAIN区域连接关系、法向partition-of-unity边界系数和0.5pixel soft分区；没有冒称B-spline。level-set使用DCT signed-distance形变。两者的已达到结果是非凸误差上界，不能写成全局容量失败。尤其atlas在16/32/64预算下只有62/31/5个TRAIN chart满足实现的rank/拓扑条件，atlas32/64的大部分起点达到100迭代上限，保护digital topology还拒绝大量更低误差但改变拓扑的轨迹。它的负结果包含**参数化、栅格化及优化限制**。

稀疏64的95%改善区间为[29.50%, 42.81%]，PCA64为[15.25%, 28.05%]。当前证据支持保留更简单的材料表示；不支持必须采用topology atlas。

![同参数预算的容量](results/a30/figures/fig04_capacity_parameter_budget.png)

![固定四家族材料图](results/a30/figures/fig05_fixed_material_images.png)

图中best-64明确为每个固定case的OFFLINE ORACLE，不是预测图；色标共用0–1.2。原始实部、虚部、参数与轨迹保存。全部16个parent的结果见下图，smooth单列。

![逐parent结果](results/a30/figures/fig06_every_parent_gain.png)

## 5. E1：bounded compensation的实际结果

nominal pair采用measurement-conditioned N3与原TRAIN规则的cavity候选，优先r2/s0，缺失才按冻结r1回退；3条件没有注册候选，保留NOT_APPLICABLE，未按truth补成功样本。原候选还包含固定phase projection，因此这不是单独隔离cavity因果效应的实验。

两端允许相同的4个在线区域的独立复contrast补偿及4个一阶boundary方向，每端12实坐标；contrast半径0.30，每边界系数0.35pixel，联合变化范围与逐系数范围分别记录。实际cavity内格固定，不能通过直接填洞消掉结构区别。F-anchor tangent、background L使用相同材料chart、box与最多100次优化。F实际验证补偿端点；L端点的额外F审核限于预先固定每family首个parent，不根据truth选取。

所有29个合法对都支付了真实F终点验证。下表是有效条件的**描述性中位数**，不是丢掉缺行后的完整parent主统计。

|family|合法条件|nominal full-F距离|补偿后full-F距离|补偿后L模型距离|
|---|---:|---:|---:|---:|
|shell|8/8|70.29|19.18|15.11|
|asymmetric_piecewise|6/8|44.92|20.92|19.90|
|mixed_scale|8/8|49.61|20.21|16.14|
|smooth|7/8|48.80|11.79|11.75|

距离使用复噪声每Re/Im方差σ²/2对应的sqrt(2)/σ白化。独立Gaussian detection draws每对1000个只用于固定对的二元判别；没有将它写成重新生成候选的端到端性能。即使最低距离8.84，等prior固定二元Gaussian错误率也远低于“距离≤1”的模糊范围；这不构成全局材料可识别性证明。

补偿的F优化有6条件、L有6条件达到100迭代上限，联合至少一方未收敛的条件为9个。全部保留可行点与停止原因，**不能用这些局部点证明全局补偿最小距离**。原协议要求至少6/12复杂parent出现清晰→模糊的合法对，本次找到0个；因此当前补偿机制没有取得继续资格。

single-anchor calibrated L真实定义为F(N3)+L(χ)−L(N3)，在同一pair的secant里该偏移恰好抵消，所以本测试中的L/calibrated-L距离相同。这只说明**本次距离任务**不能测偏移校准价值，不表示校准对绝对likelihood ranking无用。已验证的7个L补偿对在真实F下与F-tangent补偿见证相近，样本与优化限制不足以建立full-F的独立材料选择优势。

![实际结构补偿](results/a30/figures/fig08_e1_compensation_outcome.png)

## 6. 必须保留的失败与误差界

首次实现错误把phase-projection改变的全部576格都当成protected cavity，产生零补偿列。15个条件的补偿结果作废；其75次全波调用仍计费，其中15次nominal响应有效，另60次是冗余调用。修正后只保护真实移除的cavity格，并复用首次有效的nominal距离，避免重复支付已算结果。为守住192上限，仅固定四family首parent做额外L-endpoint full-F审核。**总硬上限未超，但60次冗余调用多于原32次失败/复查预留；因此不能声称完全复现原5-state-per-condition加32预留的矩阵。**这是执行偏差，保留计费，不修改科学门槛；实际正常端点计算与额外冗余调用分开公开。这个实现修正与审核范围差异不覆盖A29历史结论。

只检查一个基础uniform Neumann/Taylor误差界，没有额外证书campaign或物理fixture。其q包含整个允许材料参数box；本数据中没有得到可用于3/1判别的有效窄区间。**q≥1表示保守充分条件失效，不表示Maxwell矩阵真的奇异或物理不可区分。**实际全波仍通过原backward residual验收。局部终点余项保存，未把pointwise误差当作uniform强散射证书。

uniform证书gate没有通过；没有错误地用“基础代数检查通过”替代它。unknown背景/仪器nuisance robustness未运行，不归咎于当前已知背景问题。此处停止付款和后续实验，不为获得正结果继续加rank、换大网络或扩表示。

## 7. 成本、资源与复现

累计 **161 /192** 次新增真实full-wave调用，包括首次无效补偿、全部重算和审核；3048个complex LU RHS=1288个full-forward +352个cache-LU审核 +1408个receiver-adjoint。新LU 161次，Green/receiver/full-L residual各1288个源列；聚合计数不重复加总。

E0新增全波为0；没有新truth label，NN更新0，locked访问0。采用CPU缓存计算；GPU占用0，没有为了使用4060额外搬运缓存或建立CUDA设施。到报告生成时已计费CPU **146.97s**；每job按native wall/process CPU较大者，包含下载、失败、汇总和出图。交互式代码编辑/部分读文件未全部计时，数字不冒充所有host CPU消耗的严格上界；最终发布receipt另留账。没有做hash/SHA256校验；用户最新澄清后未追加探针或可追溯备份。

复现输入取A29两份Release资产；容量、补偿端点、所有fit轨迹作为本仓库Release原始证据。运行命令及不可重置同一campaign预算的约束见[README](README.md)。新一次真实物理复现须独立声明预算，不能删除当前账本后假装本campaign尚未付款。

## 8. 科学判断与可证伪预测

**Coverage：**现有bank的凸组合与简单sparse/PCA连续材料表示已具有实质headroom。不能继续将所有复杂材料错误归因为候选数量不足。atlas则尚不能把优化失败与表示局限彻底分开。

**Selection：**已实现的L-weighted estimator没有利用同bank的改善空间；当前值得优先研究的是measurement-conditioned参数/权重与结构先验的共同决策，而非自动换成昂贵full-feedback。没有在线拟合sparse64，因此其35.69%不能被写成算法效果。

**Compensation：**本次局部合法补偿能显著缩短数据距离，却没有找到足以跨过3→1目标的复杂parent。既没有全局求解，也没有有用uniform强散射下界，故不宣称补偿绝不存在。现有证据不足以把它确认为A29主要遗漏机制。

**预测而非已验证结论：**固定支持稀疏/PCA的测量条件化拟合，比继续扩大当前topology atlas更可能以小成本改善材料；但32参数稀疏的11.68%改善是否能由测量找到、64参数是否会吸收噪声，都需新协议验证。另一个明确反证目标是：同表示、同prior下single-anchor calibrated L若已达到full F质量，就放弃feedback-specific主张。

**决定：保留材料容量/选择路线；当前bounded-compensation及certificate-guided F路线NO-GO。E2留待下一次明确授权，不在这里执行。**
