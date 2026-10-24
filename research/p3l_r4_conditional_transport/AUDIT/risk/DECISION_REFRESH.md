# StageII固定星形决策刷新原型：实施与数学边界

本次只新增src/p3r4/decision_refresh.py及tools/run_r4_decision_refresh.py，完成小矩阵数学检查；没有运行物理pilot。工具在读取父调度器的明确release前，不能构造Book/Target或调用付费物理。StageII目前的物理几何作业不受修改。

令Lχ=I−diag(χ)G，L0为原已付精确基点的矩阵。R0由backend.solve_anchor(I576,anchor.state)一次取得，回调要求实际账本恰增加576个RHS，禁止未计费np.linalg.inv。数值R0无需被假装为精确逆：定义

ε0=||I−R0L0||F，η=ε0+||R0||F·||diag(χ−χ0)G||F。

因为||I−R0Lχ||2≤η，当η<1时，R0Lχ可逆，因而Lχ也可逆，并且

||Lχ⁻¹||2 = ||(R0Lχ)⁻¹R0||2 ≤ ||R0||F/(1−η)。

这是近似左逆的Neumann级数界。ε0=0时退化为原理想精确R0公式。它不是把一个数值谱范数估计当作已认证上界；这里明确使用Frobenius范数及乘积上界。浮点范数、矩阵乘和ε0本身没有定向舍入，因此实际实现仅是**精确算术条件下的上界原型**，不是机器级区间证书。

给定近似电流jhat，Ehat=e_inc+Gjhat，ρ=diag(χ)e_inc−Lχjhat=diag(χ)Ehat−jhat。真电流误差满足||j−jhat||F≤||Lχ⁻¹||2||ρ||F；白化观测误差上界为

e=(√2/σ)||S||F·||R0||F·||ρ||F/(1−η)。

设rhat白化范数为a=(√2/σ)||Sjhat−y||F，反三角和三角不等式给出

Φχ∈[0.5·max(0,a−e)², 0.5·(a+e)²]。

η≥1时没有这个界，记录UNKNOWN并需Full fallback。Richardson从j0+R0diag(Δχ)E0开始，使用jhat←jhat+R0ρ，最多4次修正；每次都沿用同一个logU，不因区间不确定重新抽U。

提案是完整66维Gaussianized pCN：x'=sqrt(1−β²)x0+βν，β=.002，u'=NormalCDF(η')。完整Gaussian prior可逆，所以固定基点一步MH的log ratio为Φ0−Φχ。若logU<min(0,Φ0−Φupper)则接受；若logU≥min(0,Φ0−Φlower)则拒绝；否则继续精化或付费Full。基点Φ0来自同一已付真实基点。该推导不把prior项错误丢弃于任意非pCN提案。

计划固定首8个新DEV(R4DEV_00000至00007)，每个128个星形提案，seed20261025089，beta.002，基点来自新96+64在线context及其完整已付cache。不读取任何GT。所有1024个候选均另付8-source Full作为研究审计；每个Full点由Target立即保存，Richardson每次近似迭代和GMRES最终current/field/prediction/residual也保存。基点从不沿接受事件更新，因此这不是pCN链，不评价posterior MCSE。

普通Krylov强对照是8源GMRES，固定restart8、maxiter1、rtol1e−10，使用同一个基点LU作预条件。每次实际预条件solve_anchor(vector)必须由回调核对账本增加1RHS；不能把这些LU-RHS或密集R0 setup摊成免费。每个源的迭代数、matvec列数、preconditioner RHS、计算时间及未知fallback均保存。GMRES区间也依赖已付R0提供的逆界，其standalone区间成本须计576 setup；只作普通残差精度对照时不能偷称同等认证。GMRES达到不了容差的info保持原样，不追加迭代救结果。

工具只接受匹配的release JSON：

```json
{"released":true,"experiment":"decision_refresh_star_v1","campaign_id":"P3LR4-20261025"}
```

父调度器释放后才可运行`<python> tools/run_r4_decision_refresh.py --release ABSOLUTE_RELEASE_PATH`。物理Book固定StageII、250000RHS、1200s及3GiB上限；使用现有physics.lock，忙时拒绝执行。不得由worker自己创建release。所有失败、未知、fallback、Full审计与实际callback费用进入结果；不自动重复付费作业。

成本界限：R0 setup共4608RHS，1024个Full研究审计至少8192RHS，GMRES预条件实际RHS额外全计。Richardson/G/S密集matvec和范数虽然不计为求解器RHS，仍有实际时间、内存和矩阵动作计数。每一独立部署arm的fallback费用单列为8RHS/未知候选；研究审计的Full只实际支付一次，不能再把counterfactual费用加进实际账本，也不能从研发成本里删掉这批Full。

Full的target.evaluate时间包括压缩端点持久化，而近似helper计算时间不包括其数组写盘。结果分别保存Full物理时间、含持久化Full时间、近似计算和近似持久化时间，整作业wall才是实际端到端费用；不能直接用不匹配的小计宣称加速。

任何区间接受都不意味着获得新精确cache、场或逆作用。原型若通过全部Full包含检查且实际成本有希望，仍须由父调度器另行判断是否授权真正链实验；本次不实现状态推进、双端势区间链或路径等价采样。经典Neumann/残差界也不是本项目原创定理。

DECISION_REFRESH_MATH_CHECK.json：100组任意18×18复线性系统候选，刻意扰动R0使左逆缺陷2.64e−4；所有200个Richardson/GMRES已判断事件与直接线性代数核对一致，所有区间包含真Φ，η≥1正确fallback，少计setup/preconditioner回调被拒绝。这仅验证数学/接口，不是Maxwell pilot或成本成功。新物理RHS/状态均为0。
