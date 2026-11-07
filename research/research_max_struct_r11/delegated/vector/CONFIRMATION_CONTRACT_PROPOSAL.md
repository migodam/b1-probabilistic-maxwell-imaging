# 真三维新确认合同提案：待root科学冻结

**NOT_GENERATED / NOT_RUN。** 仅写合同与几何参数，不生成材料GT、Yclean、噪声或新观测；DEV队列照常运行。机器可读合同位于 `configs/vector_confirmation_contract_PROPOSED.json`，生成器专用几何位于 `configs/vector_confirmation_material_design_PROPOSED.json`。两者均明确禁止执行，待root复核。

## 8个不同材料与4主臂

四拓扑分别为非球形长方体包覆、向外开口矩形空腔、相邻分离双体、有窄颈的连桥。每拓扑有high与mixed两个**不同支撑、尺寸和位置**的材料，不是同一几何换contrast。包覆high使用χ=8/4（保证确有包覆界面），mixed使用χ=1/8；空腔high为χ8，mixed为χ1侧壁/χ8底；双体high为χ8/8，mixed为χ8/1；桥high全χ8，mixed为χ8双叶/χ1颈。所有最小几何特征≥.1λ=2个反演格；不是亚格单元逃逸。具体CSG盒项均已写出；新材料没有调用过生成器。

统计单位是8个不同、预先指定的设计，32次反演不是32个独立样本。这里不是从总体iid采样的统计试验。8个高低几何分别不同，但仍只覆盖4个规定家族。

候选共同冻结为sourceCIE、sourceDual、fullGreen、contractedLagged。前者源内PR与β_current=[6,3,1]；sourceDual只用源内PR配β_current=[6,4,3]；后二者均用GlobalDual的global PR和β_current=[6,4,3]，材料β均[6,3,1]。fullGreen为新完整材料依赖metric的8步PG，GN方向曲率初始α，≤8 Armijo试探；contractedLagged为固定T的z坐标凸R0子问题、8步投影方向精确一维线搜索。候选不得等DEV GT后删选。

fresh GlobalDual的额外价值是：在新8对象上把global/sourcewise PR与矩阵分母读出分离。若论文要宣称**新确认对象上的增量唯一由材料分母造成**，应预先给全部8对象加GlobalDual，额外8次、约1.3–2.7小时计划费用；若主张只是四种完整部署的有限确认，DEV已有匹配GlobalDual控制，fresh可保留4主臂，并明确新确认的因果归因范围。此选择须在任何新Y生成前由root决定，不能依据DEV GT挑选。

## 物理、噪声与目标

继承12³、1728实各向同性χ、每源3N=5184交错复current、60真实外部偶极子、60接收位置×3分量的dyadic CM+RR DDA。ROI宽.6λ，k=2π，CM/RR自项不变；m0=.25，box[0,8]，R0，全部cell未知。最大εr9时约6.67反演格/内部波长，只是有限DDA模型，不能冒充连续Maxwell收敛认证或原2019作者复现。

source-driven `orth([Einc,G Einc])`以背景和实际60源独立构造，固定relative thin-SVD阈值1e-10、实际rank120，Galerkin K=PGP。不再扫描rank，不用新Y或GT选basis。该basis只捕获前两次公开背景源传播；已有DEV晚期defect作用遗漏仍大，不把它称全strong或全G近似资格。

每观测使用独立圆复高斯噪声，σ=.03||Yclean||F/√(180·60)，E|n|²=σ²，不重新归一化噪声实现。新噪声根seed提案2026110713，各观测独立SeedSequence子流；seed只用于生成端，inverse只读Y及公开σ，不读clean field、拓扑、regime、支撑、设计参数或seed。

原物理选择Phi=mean_s[||SJs−Ys||²/||Ys||²]，没有1/2，R0；材料metric为0.5ρβ²∑ws||w_s||²，ρ=1/60、ws=1/||Einc_s||²。它们都不是同一标量σ的iid ML白化目标。真实σ白化残差和全数据CN99资格单独报告。DPCS只用公开噪声球，不将跨网格模型误差塞进σ。

## 冻结自然预算与公平费用前缀

每臂MF2/3/4，各max500/min20，原系数RMS阈值.001；共同corrector320 accepted/800 strict requests、maxls12、ftol1e-12、gtol1e-8。保留80/160/240/320 accepted材料及全部已付费original-Phi，不因prefix预算截断优化器或延长自然stopping。

建议每臂CPU与wall的480/960秒前缀，以及physical-LU=80/160/240/320前缀。对每个前缀，仅从**已完成且计费不超过预算**的物理状态中按公开Phi取最小，不插值、不追加forward、不读GT。选中的material要在预测封存包中明确保存。若没有可用状态则null，不伪造零费用结果。原始wall包含暂停；减去已记载协作等待的wall只能另列诊断，不能偷偷替换原wall。

确认runner需要增加一项轻量存档：每个unique paid physical state保存N维m与state_id（320个m仅约4.4MB未压缩），为每个paid事件记completed-arm CPU/wall累计坐标；然后在线按上述公开费用/Φ规则封存prefix材料。这不需要新物理solve。现DEV继承R10保存每个paid Phi与费用事件、80/160/240/320材料；不能把其中未存material的任意480秒点伪称已完成GT材质评分。

共享G/E/S/SVD/basis按unique campaign只收一次，另给单方法冷启动部署视图（共同setup与仅矩阵方法需要的basis费用都明示），重叠视图不能相加。物理5184阶LU、rank120小LU、primal/adjoint60 RHS、辅助FFT、所有line trials、setup、failed starts、helper子进程、传输与评分都保留真实费用。

## 费用预估与可选跨网格组

当前新DEV的v04_sourceCIE已完成320 accepted，实付CPU676.875秒、327个physical LU。以此作量级依据，暂按每inverse600–1200 CPU/wall秒规划：matched32次约5.3–10.7小时；若加2对象×4臂crossgrid，共40次约6.7–13.3小时。此为计划区间，不是wall硬上限；达到800请求或1500个辅助迭代的对象可能更贵，不能将现有一个对象的速度当普适承诺。R10的2个matched12+2个cross18生成合计36.44CPU/37.70wall秒，说明生成相对inverse小，但新几何的GMRES成本和失败仍实记。

可选crossgrid预选c13开腔high、c18连桥mixed，用相同连续CSG几何18³生成、12³反演，另独立噪声并运行同4臂；与matched成功统计分开，新增独立材料数为0。GT评分仍对精确coarse有效cell均值；18→12的模型/几何误差不是测量σ，不是连续Maxwell认证。

全候选、费用prefix和选择规则共同封存后，只有root可读取GT。材料χ相对L2、range归一RMSE、paper Err、双方χ>.5的effective-cell IoU、固定有效contrast分区误差共同报告；不按GT在线挑迭代、rank、算法、家族或根，不设普适胜利宣称。
