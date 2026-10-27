# R6 3D 材料商 prior 独立审计

审计只提取 `stratified.py` 的常量与六个 prior/decoder 函数执行，没有导入 paid physics adapter，没有创建 Book/Target，没有正演、历史 GT 或新训练。初次使用256000个独立 prior 样本；修订复查复用该记录，新增随机样本为0。所有工作CPU单线程，首轮2.33秒、边界修订轮1.66秒，远低于2分钟。

## 代数结果

令 Q=H(1.3z)，D=(Dr,Di)，V=QD。z为[-4,4]截断标准正态；原属性的 probit 变换把两项截断正态精确变为区间[Φ(-4),Φ(4)]上的独立均匀变量。因此 Dr∈[.8+Φ(-4),.8+Φ(4)]，Di∈[.03+.15Φ(-4),.03+.15Φ(4)]，联合常密度的分母为 `.15 Z0² = .14998099785674063`；z自身密度还独立除以Z0。这三个常数没有漏乘。

给定Q>0时，D→V的二维行列式为Q²，所以transition density对z积分包含Q^-2。材料密度还须加上纯相源标签的 `PH*r(V)`；同一个V可能同时收到transition与pure标签的概率，当前代码分别保留两种标签并共享likelihood，处理正确。材料商代表元只是该材料的一种合法原坐标表示，不是纤维条件后验抽样。

Q=0的整个host纤维应被积为一个原子，质量PH=.22086048232212066。它不应作为二维连续密度或一个任意面积单元处理。由于对称截断，纯相质量同为PH，transition质量=.5582790353557587。Q>.5等价z>0，其prior概率恰为.5；host不计该事件，pure全计，transition numerator必须仅积z>0的那一段。当前fibre已正确按z=0拆分。

R=t³、I=θt³的二维行列式为3t^5。它消去host附近连续密度的R^-5/3主奇异阶。改变坐标不删除host原子，也不把3D模型提升为Full66边缘后验。

## 已发现并修复的数值缺口

初版纯相使用径向矩形indicator，导致prior总质量非单调：1.051605、.989680、1.014005、1.007143、1.001093（nr=8,16,28,48,80）。初版source总质量、均值及二阶矩准确，但latent phase indicator没有在z=0拆分：S3 nz10为.507848、nz18为.492416，真实值为.5。`QA_initial_unsplit.json`保留这些未通过的记录，不能覆盖为成功。

主线程在任何A物理查询前修订：纯相独立tensor GL、transition按属性矩形裁切界分径向panel、source增加z=0边界。独立复查保存在 `QA_boundary_patch.json`，被审计源码快照为 `SOURCE_USED_BOUNDARY_PATCH.py`。

|商空间网格nr/nt/nf|prior总质量|transition质量|pure质量|Q>.5 numerator|
|---|---:|---:|---:|---:|
|12/6/48|1.0000164523|.5582954877|.2208604823|.4999919120|
|20/8/64|1.0000043335|.5582833689|.2208604823|.5000372899|
|32/16/96|1.0000012558|.5582802912|.2208604823|.4999919303|
|64/32/128|1.0000002000|.5582792354|.2208604823|.4999986233|

修订source所有计划级别的总质量、三branch质量及Q>.5 numerator均与解析值一致至约1e-15。以上支持继续执行预算内、逐任务检查收敛的实验；**没有建立posterior reference证书**。尤其商空间latent-event numerator误差仍非单调。

初版fibre在5个远离/靠近host的材料点与独立adaptive quadrature一致，相对误差小于4e-15；200个随机pushforward材料的商代表元decoder replay最大绝对误差2.23e-15。256k MC均值与解析(.65,.0525)相容：实部偏差约-1.64个MC标准误，虚部约-.76个标准误，二阶绝对矩偏差约-1.05个标准误。MC只核查prior，未产生posterior样本。

## 仍需明确的界与任务测试

`prior_mass_quadrature`是连续prior积分的数值诊断，**不是被遗漏posterior质量界**；`1-prior_mass_quadrature`甚至可为负，不能取绝对值后充作证书。当前A没有提供严格的未覆盖likelihood或normalizer区间。有限正节点也不表示节点之外的连续空间被证明没有后验质量。

若将来真实地界定某遗漏域prior质量ε，并证明该域likelihood上界U以及已保留证据下界Zlo>0，才可得遗漏后验质量≤εU/(Zlo+εU)。无这些条件，小prior质量可能被likelihood放大。对未移位Gaussian likelihood可用U≤1，但Zlo仍不可省略；shifted evidence需要一致移动U。任务均值还须界定遗漏域任务幅度。

后续必要检查为：每个observation/noise分别比较evidence、各branch posterior质量、材料均值/方差、latent Q>.5与material R>.02事件；source与quotient互相核对；避免把prior normalization通过当成尖锐likelihood下收敛。material R>.02事件还有径向indicator，需要在R=.02分panel，或明确保留该事件的任务级收敛不确定性；latent numerator可在.5倍属性边界处补径向panel以改善收敛。严格误差或概率证书仍需可信区间/余项界，普通GL级别接近不足以自动提供。

文件：`QA_initial_unsplit.json`、`QA_boundary_patch.json`、`MATERIAL_EVENT_PRIOR_QA.json`、三个可重复执行脚本，以及修订源码快照。没有编辑main A代码。

material事件的独立prior数值证据：解析一维积分 P(R>.02)=.7181396482583506；修订商空间12/6误差+.00242868、20/8误差−.00207727、64/32误差+.00014984。修订source S3 10/8误差+.00114240、16/9误差−.00010862。这是一个明确反例：prior总质量接近1仍不足以认证某个indicator任务的积分。此项只用确定性prior积分，新增MC样本仍为0。
