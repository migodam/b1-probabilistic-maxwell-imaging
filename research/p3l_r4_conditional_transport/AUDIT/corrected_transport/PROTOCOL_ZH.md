# 待父调度器释放的显式q完整概率校正

新增src/p3r4/corrected_transport.py及tools/r4_corrected_transport.py。只实现、纯mock核查，无神经采样、付费Full、GT/locked读取或训练启动。后续是否运行、是否进一步冻结production及全部概率资格由父调度器决定，不自动扩大pilot。

freeze spec示例：

```json
{
  "checkpoint":"/absolute/FROZEN_EXPLICIT_Q/BEST.pt",
  "loader":"anchor",
  "parents":["R4DEV_00000","R4DEV_00001","R4DEV_00002","R4DEV_00003","R4DEV_00004","R4DEV_00005","R4DEV_00006","R4DEV_00007"],
  "warmup":0,
  "production":128,
  "seed":20261025141
}
```

`<python> tools/r4_corrected_transport.py freeze --spec SPEC.json --out NEW_PREPARED_DIRECTORY`只接受4或8个明确DEV父样本，拒绝locked；模型只接受raw/anchor加载器的full66 Spline或diagGaussian，拒绝CFM、数值ODE密度、低维补prior和其他未冻结家族。raw模型必须只条件rawy512。保存模型checkpoint和所有在线来源byte/mtime，不做SHA。

每parent只读原nativecomplex观测和canonical StageI truth-free DEV context；验证source sidecar GT_read=false及x0与nativeq0精确probit一致。将条件导出为显式context_x、必要时context_J(已付nativeJ)，不让可能含GT标签的任意packed x冒充online anchor。anchor q-law整个过程中固定条件于同一原Y/x0/J，绝不随当前链state更新条件。没有新J、projection或corrector，原prior/noise不变。

三种固定核：

1. prior-pCN beta=.002，完整x66提案sqrt(1−beta²)x+betaν，原Gaussianized prior可逆，logα=min(0,Φcur−Φnew)。
2. independence q，每次完整66维独立提案q(x'|Y)，显式logα=min[0,Φcur−Φnew+.5(||xcur||²−||xnew||²)+logq(xcur|Y)−logq(xnew|Y)]。
3. 固定.5混合，先独立于state抽selector，按对应component的有效MH比值校正；这是两个各自π不变核的常数混合，不能用错误的“混合proposal density”代替component ratio。

每parent4条链：第一条从online anchor开始，另3条从q独立抽样；同一组4个x初值用于三个核，不挑最低残差或最好GT。**包括anchor在内所有4个初始状态都在每个核中分别实际付费Full**。默认8parents×4chains×128steps×3kernels共12288次候选Full，加96次初始化Full，总计划99072复RHS(每次8RHS/LU)。单parent q初始化采样共同准备成本单独记录在整体model调用费用，不能把既有q训练或online context/J称为冷启动免费。

固定seeds生成独立角色流：每个parent/chain/step的pCN innovation、logU、常数selector、q draw seed相互分开；核间对应component使用共同innovation/U/qseed。没有用state或residual挑随机数。当前logq只在状态改变后更新并缓存；mix的pCN接受也更新该x的logq以备后续q步，拒绝保持原cachedlogq。模型sample/logdensity调用次数及wall另记，物理只由Target的Full8源触发，没有gradient/J/adjoint查询。

后续运行必须有匹配release：

```json
{"released":true,"experiment":"corrected_transport_full66_v1","campaign_id":"P3LR4-20261025"}
```

`<python> tools/r4_corrected_transport.py run --prepared PREPARED_DIRECTORY --release RELEASE_PATH --stage III --job UNIQUE_JOB --max-rhs PARENT_CAP --max-wall PARENT_SECONDS --out NEW_RESULTS`。stage由父调度器选择，III仅示例。max-rhs必须预先覆盖整个冻结名单，且每次初始化/候选前再次preflight8RHS，实际callback也受硬上限；Book/physics.lock、wall和3GiB限制保留。freeze的warmup/production可预先配置，run时不可重选或自动加步数。费用不足、数值异常或model密度失效即保留证据并停止，不换seed或best residual救运行。

Target立即持久化每个实际Full endpoint的q/χ/current/field/prediction，runner只保存独立的x/随机数attempt及endpoint引用，不重复大field数组。每步history和CHAIN保留全部拒绝重复的精确x与同一endpoint引用，包括warmup；production标志单独存。没有删掉“没移动”的样本。每parent每核的4链独立初始化诊断包括材料函数、Φ、prior、accept/component、split-Rhat及每链初始正相关pair ESS/MCSE；不跨parent计算Rhat。常数/停滞函数不伪报ESS=N或MCSE0，多链停在不同常数时Rhat为不合格无穷大。所有指标都FINITE_RUN_UNQUALIFIED，缺共同Full reference，不能据它们自动授予FULL后验资格。

显式q的logdensity包括原normalizer、正交frame和anchor坐标变换Jacobian；神经Spline/diagGaussian变换主要float32，虽有代数密度，仍不是机器级exactMH证书。这里验证的是完整形式正确的比值和有记录的数值实现，不替代父调度器对精度、support、残差/SBC及有限混合的复核。

KERNEL_MOCK_CHECK.json：200对纯mock66维状态，independence detailed-balance误差≤2.84e−14，pCN prior抵消≤6.70e−13，常数.5有效核混合DB误差≤6.83e−13；常数selector、共同角色流、拒绝重复保留与停滞链不伪报ESS/MCSE均PASS。该核查0神经采样、0物理RHS/状态、0训练，不是实际校正实验结果。
