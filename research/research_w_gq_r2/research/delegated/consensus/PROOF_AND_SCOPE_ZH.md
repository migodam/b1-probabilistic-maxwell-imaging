# E3 条件共识的范围与证明

此单元只做旧 118 墙 × 25 位置 counting prior 的受控逻辑回放，标签为 `REPLAY_NOT_NEW_PDE_SPEED_EVIDENCE`。目标固有参数已知；不验证未知材料、孔洞、连续目标或冷 PDE 加速。

设条件目标分布为 K_w，背景后验权重为 beta。精确目标分布为 beta K。若 delta=max TV(K_w,K_v)，则 TV(beta K,beta' K) <= delta TV(beta,beta')。证明：令 t=TV(beta,beta')。t=0 时直接成立；否则把 beta-beta' 的正负部分分别除以 t，得到两个概率权重 u,v。其条件混合之差的 TV 不超过 sum_wv u_w v_v TV(K_w,K_v) <= delta，因此原差不超过 t delta。此为一般 Dobrushin 概率工具，不主张原创。所有条件距共同 q 的 TV <= e 时，任意背景混合距 q <= e，无须求 beta。

平方位置损失下，共同动作 a 的条件 regret 为 ||a-E_w X||²。若对所有 w 都 <= e，则任何混合下 regret <= e；只检查采样墙的条件均值，不能把未见墙纳入。后验共识与共同动作均不保证真实目标误差小，必须另报 GT RMSE。

Gaussian 白化模型的精确式为 D_w(x)=<r_w0,Delta_w(x)>-||Delta_w(x)||²/2，其中 r_w0=y-F(w,x0)，Delta_w=F(w,x)-F(w,x0)。若 U 包含全族增量，则只需 P_U r_w0 与完整增量函数；相同增量本身不够。两墙差满足 D_1-D_2=<rbar_x,deltaDelta>-<delta f0,barDelta>。本实现保留全部残差投影项；离线全表增量空间与距离不用于选组或免费界。

背景逐行尺度反例：等 prior，L1=(1,eps)，L2=(eps,1)。将 L2 乘 eta，不改变该行内部 contrast，却改变目标混合。`ALGEBRA_CHECKS.json` 保存数值实例。缺失原子反例固定所有墙 anchor 为 t，并完全查询第一行 (t,t²)；未查询行可为 (t,t²)，也可为 (t,1)。两张表的查询记录完全一致，而目标后验可相差接近 1。因此采样条件一致不能覆盖未查询条件。此是抽象概率表，不冒充 Maxwell 物理反例。

本回放没有可独立获得的未查询物理界；未查询 likelihood 默认 [0,1]。在线的全体条件 TV 半径上界为 1，停止判定始终 NOT_CERTIFIED；即使离线发现接近，也不能追溯升级为免费证书。仅当独立物理界覆盖所有正 prior 背景时，才可重审冷 PDE 入口。
