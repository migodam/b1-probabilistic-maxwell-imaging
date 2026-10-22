# Research W P0：可复用资产与后端审计

日期：2026-10-22。范围：只读源码与仓库元数据审计；没有导入后端、运行物理解、计算哈希或改动共享资产。执行协议：`source_package/P0_P1_PROTOCOL_ZH.txt`。以下路径均相对于 Inv_SLAM 根目录。

结论：首选复用 `Gaussian/code/a2/physics.py` 的二维出射 Helmholtz VIE 数学与 FFT/GMRES实现，结合 A24 的直接多 RHS、共轭伴随与费用接口。现有代码不是开箱即用的 W 后端：平方网格、环绕天线、材料上界和单频配置都必须在 W 独立适配。未找到现成可直接接入的二维 FDFD/FEM+PML 独立参考。P0-B仍未运行。

|资产/源码位置|接口与证据|可复用范围及限制|
|---|---|---|
|`Gaussian/code/a2/physics.py:46,68,145`|`Geometry`, `FFTGreen`, `VIE`；`system`, `solve_linear`, `solve`, `forward`, `material_tangent`, `dense_reference`|e^{-iωt}，g=iH0^(1)/4，`j=χ(E+Dj)`，任意有限网格材料，不是 Born。FFT 线性卷积用3n−2补零，避免周期绕回。出射 Green 函数直接表达无限均匀外部域，无需人为吸收边界；这不是执行了PML核验。|
|同文件`:20,25`|`self_cell_green` 等面积圆自项；`square_self_cell_green` 方格自项，`cell_integrated=True` 使用 Gauss–Legendre源单元积分|可复用积分策略。方形单元假定写入构造器；矩形单元必须重新积分。默认中点近似和圆自项不能冒充独立细化参考。|
|`Gaussian/code/a2/matrix_free_tangent.py:20`|`TangentState`, `MatrixFreeTangent.matvec/rmatvec/damped_step`，实参数到复数据的精确离散 Jv/J*|跨频状态列表；切向/伴随 GMRES RHS计数。`partials[p,N]`可接全材料/level-set参数导数。现有`σ`仅标量，不能替代W完整协方差和nuisance合同；须显式约定实复打包。|
|`Gaussian/A24/a24-source-anchored-om-feedback/vendor/a2/physics.py`|A2快照，含同样VIE、切向和dense reference|用于固定来源而非新物理模型。|
|`Gaussian/A24/a24-source-anchored-om-feedback/src/a24/physics.py:70`|`Physics2D.full_state`, `factor_matrix`, `solve_factor(adjoint=True)`, `material_pullback`; `FullState`保存χ、场、数据、L、LU、残差|强直接多 RHS与伴随模板，显式费用统计。硬编码24×24默认、边长1、8发32收、圆周布置、单频k=2π；χ实部[-.5,4]而W要求[0,8]，不能直接沿用。数据是散射场，需要明确直达项/校准合同。|
|`Gaussian/A24/a24-source-anchored-om-feedback/src/a24/relative_feedback.py:344`|`RelativeFeedbackModel.forward/jvp/adjoint`|模型种类必须锁定全反馈/F路径；低秩或局部模型不能作为W最终full-wave reference。|
|`Gaussian/A24/a24-source-anchored-om-feedback/tests/test_physics.py:52,85,120`|原A2 dense比对、G/S/背景共轭伴随、实伴随dot、四步长FD测试代码|是已存在的核验方案。本次未运行；同离散核的不同求解路径不是连续误差证书，也不是独立物理模型。|
|`Gaussian/A25/a25-adaptive-om-feedback/src/a25/torch_physics.py:141,403`|`TorchReference`, `TorchRelativeModel.forward/jvp/adjoint/vjp_complex`，复128，缓存/来源检查|备选导数模板；来自同一Green资产，需重建新几何；CPU优先，无需引入GPU基础设施。|
|`Gaussian/A30/a30-maxwell-material-manifold/src/a30/physics.py:157`|`FrozenMaxwellAdapter` + vendored A29/A24；N3 manifest显式指向此文件|它只接受已付费背景缓存，24×24既有几何；不构建新几何、不生成W观测。不得通过换标签将既有缓存当成W验证。|
|`Gaussian/A17/N3_APPLICATION_R1/src/core.py` 与 `SOURCE_MANIFEST.json`|材料/观测real packing、完整Jacobian、Schur/GN线性代数|可复用代数与验证模式；该core自称FP64代数，既有N3缓存不是W墙体实验。|
|`Gaussian/B1/b1-p1-material-target/src/b1/physics.py:35,64,88`|3D vector multilayer Mie `Forward.run(implementation='matching'/'recursive')`|两套层状球参考，不是任意有限二维墙；可复用独立实现对比的实践，不能给W提供物理参考。|
|`Gaussian/C1/c1-physics-state-generative-imaging/src/c1/problems.py:14`|`DipoleProblem.forward/adjoint_gradient/jvp`，13点向量DDA，4个强度数据|源码明确“not continuum-Maxwell convergence validation”。不是TM墙体或复近场独立参考。|

## W几何适配与分辨率

32×48可以作为1536个逆表示自由度，而不是已合格的前向分辨率。若按(nx,ny)=(32,48)，4λ0×3λ0给出(dx,dy)=(0.125,0.0625)λ0；目前A2统一`n`和`h`的构造器不能表示它。若交换轴，仍是矩形单元。最小材料波长在f/f0=1.2、Reε=9处约0.278λ0；按前一种轴顺序仅约2.22与4.44格/最短波长，不能凭网格名称声称精度。

低改动路线是**逆表示32×48 + 独立更细等方格前向域**：例如长宽比4:3采用64×48方格只是几何可行候选，是否足够须按数值误差要求确定；也可在包围方域零对比padding上复用现有方格FFT，保留真实有限矩形未知支撑。不可把域外zero padding误当成已知域内墙mask。另一条路线是扩展FFT为(nx,ny,hx,hy)，同时重写矩形自项、卷积切片和积分权重；需要新的dot/FD与独立核组装验证。

发射/接收位置必须改为y=-0.3、孔径3λ0、4发12收，并明确端点排列；每频重建k相关E、S、D和直达项；不得使用现有圆周天线。被动号Imε≥0与e^{-iωt}相容。频散、损耗范围、相关长度、最小特征、T、探针集、噪声协方差与跨频共享nuisance仍须另行冻结。

## 独立参考与数值误差限制

A2 `dense_direct`独立于FFT执行路径，用同一Hankel与自项组装；在`cell_integrated=True`时它仍按中点off-diagonal组装，不能当成同一高阶积分kernel的精确dense reference。应将“FFT/直接同核一致性”“积分阶变化”“前向网格细化”“另一离散后端”分开报告。A2 certificate模块只提供有条件的**离散耗散界**；zero contrast需要处理，极小损耗可能无正裕度；不提供连续模型误差证书，且当前核校验函数内部会计算哈希，不应原样调用到本任务。

本次在Gaussian排除vendor的Python文件中搜索PML/FDFD/FEM相关实现，未找到可用于W的独立二维PML后端；不是整个磁盘不存在此类软件的断言。已有VIE可用解析出射条件免除截断反射，但必须披露“analytic outgoing boundary; PML NOT_APPLICABLE/NOT_RUN”，并由主线程明确是否满足协议“PML、网格核验”的操作解释。最终声明按协议限制于最细已核验离散模型；网格差并非严格epsilon_F/epsilon_H上界。

## 费用与预算风险（算术估计，非benchmark）

- n=1536时单个complex128密集N×N数组约37.75 MB；L、G、LU及复制会成倍增加。五频缓存不能称免费。
- 两轴加密2倍使N变4倍：单矩阵约604 MB，密集LU理论算术量约64倍；数组量16倍。FFT/GMRES更省存储，但高对比、小损耗、多次散射可能迭代数大；失败与预条件构造均计费。
- 一材料5频、每频4发的完整端点最少5材料-频率分解+20complex RHS（直接LU路线）；每次完整Jv或伴随遍历另20 RHS。CAL/JOINT和每个任务材料需分别付费，只有**相同材料、频率、网格**可共享LU。
- 2000分解允许的粗略上限仅400个五频完整材料状态，尚未扣除独立细化、任务验证和失败复查。12位置×2对比=24任务若每任务都是新完整材料，一端点仅前向已是120分解/480 RHS；一pair为240/960，独立细化另算。需在运行前明确构建/验证probe与最终任务评估如何分配。
- 2 CPU小时是否可行必须通过协议允许的**一次公开solver benchmark**决定；本审计没有benchmark、没有实测时间，不提出增预算。

## 复用公共仓库

首选既有 `https://github.com/migodam/b1-probabilistic-maxwell-imaging` 下 `research/research_w_p0_20261022/`。本地checkout：`Gaussian/B1/b1-probabilistic-maxwell-imaging`，git origin与此一致；2026-10-22 GitHub API读回public、main、has_pages=true，本次status为空。该仓库已有`research/research_t_20261015`、`research_f3_quotient_exploration`、`p2_material_fiber_20261019`等独立研究子目录，适合W延续。以上W路径只是建议；未创建、未push、未宣称发布。A17仓库更适合计算分支方法证据，B1更适合W墙/目标可辨识性主报告。
