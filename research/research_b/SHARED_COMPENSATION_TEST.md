# 共享补偿：物理见证与边界

每个世界可以有自己的ε和gain；该世界必须对整段history沿用同一个有界gain。首先堆叠旧新Jacobian再消去nuisance，不能分别消元后把信息当作joint信息。非线性远分支则优化合法补偿并验证full端点。

## 已运行的三类机制

**A：新增频率破坏部分合法补偿。** 固定central r=.25R，旧频优化proxy距离.89991；加入2.8后的joint full距离3.22530，3.2为3.17618。大孔3.2在同一ε端点的shared-vs-independent gain差的平方约.58907。但ε频谱本身也有贡献，不能将整个距离提升全部归给gain。

**B：散射变化而材料仍相近。** central r=.08R，分别新增1.2/2.0/2.8/3.2的joint full距离.03910/.06859/.14027/.14227。分频得到的优化材料不同，因此这些数值本身不证明三步也无效。另行同世界审计如下。

**C：间接校准价值。** known-reference校准不直接依赖未知材料，raw material差为0，但可通过限定shared gain改变joint Fisher与材料风险。例如E0 solid校准预测风险下降.019930、repeat .016637；2.4仍更好.032561。并非校准总优于频率，raw conflict单指标全选校准反而表现差。

共12个joint端点完整full-Mie验证。bounded优化是best-found，不是全局最小；独立gain/只新频等额外优化对照为proxy，未宣称全部fullvalidated。详见results/compensation/COMPENSATION_MECHANISMS.json/.csv及原优化trial。29次callback请求、缓存复用后22次真实新attempt，失败0。

## 固定同一小孔世界的三步上限见证

独立6次full审计拟合一个共同ε/g的r=.08R世界，对全六频的白化距离约[.02590,.05705,.07193,.09519,.13138,.13521]。calibration与独立geometry动作也列入同一有限菜单。

对这两个固定世界、已声明的相同Gaussian协方差以及只依赖历史的策略，KL链式法则给：

`KL(P_history || Q_history) ≤ (d_old² + T max_a d_a²)/2`。

取T≤3，Pinsker给等两世界**设计先验**最小分类错误率下界`.5*(1-min(1,sqrt(D²)/2)) = .4397423865`。不把design prior叫posterior。该见证是固定两世界与菜单的数值结论，有限Mie误差和所声明噪声/constant ε/gain模型均是条件；不是一般材料空间、现实仪器或所有可能测量的不可能定理。

完整原始材料、距离、full字段与预算：results/SMALL_CAVITY_HORIZON_WITNESS.json/.npz。体积MAE对微孔很不敏感，结构判别错误率应单列。
