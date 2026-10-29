# 理论到实验：哪些命题真正被检验

输入报告全文、源代码、协议与历史审计见inputs及audit/input_evidence、history_contract。原两cell对称、branchprior及monodromy仅回归，不算新发现。新实验编号按用户任务：E1非对称连续图册，E2独立离散分支，E3跨表示，E4概率；与附件原建议编号差异明确保留。

## Maxwell与材料约束

采用真实disk Green：自项 iπka H1(ka)/2−1；互项 iπka J1(ka)H0(kd)/2。完整状态L=I−diag(m)G，j=L⁻¹diag(m)Einc，y=Sj。不是Born。Full-J为(SL⁻¹)_i E_is；transpose解支付一个receiver RHS，source-major实化与sigma一致。

在接收方程Sj=u下，j=j0+Nc；shared-material minors为0只是非节点条件的一部分。实际算法保留实m并解j−mE=0，所以同时处理所有源、实值支持和节点。当前构造没有节点，节点一般仍须保留m。单独dark-current已经得到不可实现反例。

## 维数、局部证书与全域覆盖

n=6,k=4且合法内点、邻域满行秩时，标准IFT给二维局部纤维；这不是原创维数定理。interval_chart使用储存系数为精确常数的复矩形区间算术，R近似逆的Neumann余项包围L⁻¹，再包围J。T(d)=d−R_J(F(d,a)−u)将dependent box映到自身且带权导数范数<1时，每个自由a有唯一dependent根。三参考半宽为[0.00625, 0.00625, 0.0015625]。不包含真实Hankel系数的严格区间误差、不认证整个材料域或continuum。

所有738端点的后付伴随J均满行秩，最大dependent条件数7.403；这只是节点诊断。原在线corrector未执行逐节点rank-stop，事后无越过阈值不追认为原算法在线证书。两个方法都没有认证跨过奇异chart、整个连通分量或全部原prior质量。

## 测度

材料图册m(a)=[d(a),a]的coarea密度p(m)/J_F * sqrt(det(DmᵀDm))，与p(m)/|det J_dep|相等；current恢复需要的是材料切向metric，不是高维c的平坦密度。E4只检验同一原prior的局部核。固定精确(p,q)时的0.75，与有限噪声的0.6751644709不同。均匀branch不能由群的存在推出。

## 原任务差异与未实现

原默认32/64双层E4的解析响应预算不适合保守30k合并账；本轮只新32并明确缓存64。最强普通基线具有最小接收伴随RHS；成本后补控制挑战而非挽救current结果。复杂monodromy因有限覆盖/判别式前提缺失NOT_RUN。完整Fourier/Wavelet成功只记换基。新模型不是P1，原Full66后验仍UNQUALIFIED。
