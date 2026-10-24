# 失败、修复与未运行模块

本页保留科学失败与工程失败的区别。原CLOSED/ERROR、source snapshots、partial数组及失败账单均不改。

|对象|实际失败/局限|采取的不同机制|最终状态|
|---|---|---|---|
|点阈值材料矩|求积偏差、有限head风险差|精确h质量bin，再anchor mass初始化|三类均不胜强CNN，停止扩大|
|完整条件生成|相对NLL改善但样本Full残差很大|完整物理力+成对监督对照、三角和几何CFM|覆盖/response冲突未解决，FULL未合格|
|CFM显式密度|ODE容差和离散map不一致|双精度全trace adaptiveDOP853|895秒后partial，未完成项NOT_RUN；不用于exact-MH|
|active-row精确反馈|变化集大，setup不划算|与Full及普通goal/Krylov公平对照|加速NO-GO|
|首轮区间链|RHS少但wall更慢|收紧逆范数包络、同核重放|新locked8数值加速通过；原失败不改写|
|Euler commutator|二阶自身场加速度污染|Heun及直接bracket，距离门复核|Heun仍不能保持任务距离/成本，停止扩大|
|localDINO-like/J学习|误差门全部失败|value/JVP、动作ranker、普通fieldupdate/Broyden|不宣称文献复现或净加速|
|普通SLSQP代表|callback只记录major trial，未覆盖最佳有效query|零新物理best-query；普通Full-J TR双目标|原预算停事实保留，data-only更拟合却risk变差|
|fusion checkpoint元数据|闭账6份moment_K被参数数目覆盖|只读compat从封存TRAIN/VAL门恢复K32/0|权重/原元数据不改，canonical推断完整|
|fusion批次复播|2e−6技术容差原FAIL|另列1e−5只读技术核验；fresh统一batch64|原FAIL保留，不靠容差改写主评分|
|首次代表离线评分|Python list矩阵二维索引TypeError|附加wrapper仅转等价ndarray、新输出目录|新RHS0，原失败收据及partial保留|
|PPC/后拟合release接口|producer状态名/绑定键错误，物理前停止|独立source-only绑定修复|无伪造物理结果，真实失败wall仍计费|

未执行：连续PDE严格interval certificate、可信FULL posterior teacher/credible intervals、整个LazyDINO文献复现、完整动态全秩状态相关RL策略、48/64二维扩展、分布外noise/采集变化、3D和新Mie。任何尚未运行的原R3建议不标PASS。

没有因概率reference不合格停止监督或几何路线，也没有在测试GT公开后继续调参。现实限制来自高信噪比下完整law的response calibration、过渡相的非线性平均、局部谱的不稳定性、参数/field事件耦合和当前网络/优化机制；不把这些概括成“Maxwell无多解”或“posterior mean不如MAP”。
