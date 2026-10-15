# 统一材料搜索与archive：冻结设计

物理权威对象是24×24 complex128 χ，Re∈[-.5,4]、Im∈[0,.5]，不是孔数标签。已知background与8源共享材料不变。全局重启、连续位置/尺度的局部加减材料、多尺度分布编辑是共用动作；方法只选择动作/parent，不能更换corrector。

材料距离为 `sqrt(mean((ReΔ/4.5)^2+(ImΔ/.5)^2))`，δM=.03。结构descriptor取`abs(χ−knownbackground)/1.2`，gaussian尺度0、1、2格、边界constant0，occupancy=`clip((smoothed_amp−.20)/.10,0,1)`；dS为三视图RMS差，δS=.08，且原幅度>.25的mask至少改变6格。敏感性dS=.04/.08/.16；6格约束固定。负实/复对比也被计入，不以正Re定义support。posthoc数字拓扑用前景4/背景8互补连通，仅报告≥6格特征，不输入算法。

active材料packing要求每对满足material separation或有效structural separation；结构代表单独pack。合法full端点包含duplicates都进入不可变success ledger。first模式保留最早代表；quality模式只有full residual更低且替换后仍与所有其他代表分离才替换，旧world保留，重插不重复奖励。frontier不合格端点不能进入正式候选。

Top1用最小full residual。TopK保留合法base（若仍在archive）后按固定field/structure覆盖选择，少于K保持短集合，不填invalid。方法不读truth、Best-of-K、孔数、真实位置或未来结局。候选权重none，posterior/mass均null。

四臂：独立base-centered correlated全域multistart、连续尺度/有符号幅度static、固定emitter混合的unstructured QD、轻量Beta bandit。多尺度场重启不是原prior均匀采样；有限Haar/DCT/局部mask库仍有表示与带宽偏好，不冒充全材料空间完备支持。QD/bandit每三步的第一步强制global，比例1/3≥20%，不被筛掉。所有生成提案直接full求值，无proxy筛选，screened-out audit不适用。

原历史32real Haar正确器固定两次有界GN、trust=.35、拒绝后×.3；source函数不改。此有限固定步regression corrector可能过度移动已经相容的提案，这是本比较的限制，不冒称T的proposal-retentive corrector。每轨迹最终最低残差材料再独立freshfull验证，正式安装只看冻结数据/物理规则。

两轮first/quality使用相同动作配置与配对随机种子；archive导致后续parent不同是真实在线差异，不把离线重排当作该反事实。物理states匹配，正常receiver-transfer cache复用保留；实际RHS独立计数，不为凑相同RHS额外求解。
