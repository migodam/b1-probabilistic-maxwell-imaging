# 只读跨会话公式复用：2026-10-29 17:54 UTC 快照

本轮按用户要求读取了同项目活跃会话，没有向其发送新任务、修改其源码、启动其物理作业，亦不把正在变化的结果当作R3已运行实验。

- `推进 A18 Route B 科研` 当前实际工作目录为 `Gaussian/B1/b1-p3l-r10-bilateral-structure/research/p3l_r10_bilateral_structure/`。读取了 `BILATERAL_OPERATOR_THEORY.md` 和 `SHAPE_AMPLITUDE_TRANSFER.md`。
- `Explore algebraic Maxwell solver` 当前在 `research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r3/`。只读到其当前运行摘要及代码片段；不把尚未定位的完整理论称为已验证。
- `科研复盘与R7提示詞` 最新assistant内容仅返回content-reference，无法从这次工具输出读取具体正文；只使用可见用户问题及已定位的本地R10推导，不假装读完不可访问回复。

## 可复用的确切结构

R10在固定宿主与形状下写出 $D=\delta A$、$P^2=A$，从完整反馈得到

$$F(\delta)-F_h=\delta L(I-\delta B)^{-1}R,\qquad
\partial_\delta F=L(I-\delta B)^{-2}R.$$

这给出真正非线性的形状—材料幅度分离；不是Born线性幅度拟合。原文已明确普通Schur/resolvent不是新定理，profile也不等于marginalization。其现有256阶缓存条件实验报告，形状固定后需约21–58个幅度节点摊回Schur预付；这是邻会话的DEV范围，**不是本轮FEM或未知目标已获得的速度**。

向当前FEM迁移须保留参数依赖右端。固定墙W与目标形状S，令

$$L(\delta)=L_W-\delta M_S,\quad b(\delta)=b_W+\delta b_S,\quad u_W=L_W^{-1}b_W.$$

则直接相减并移项得到

$$u(\delta)-u_W=\delta(I-\delta G_WM_S)^{-1}G_W(b_S+M_Su_W).$$

它可以支撑固定几何内的低维材料profile；每次墙或形状改变，相关算子/载荷也改变。不能把其他会话的256维VIE缓存当作当前FEM的免费目标缓存，不能删去$b_S$，也不能把复$delta$按频率独立拟合。跨频仍只有共同实介电率与参考损耗两个实材料变量。

R10的双边导数 $dF=O_D\operatorname{diag}(dd)E_D$ 及Khatri–Rao Gram提醒：收端谱需要与激发共同看。这与本轮公式集和A17X一致，不作为重新发现。R10噪声定义与本轮不同，不复制其$\sqrt2/\sigma$白化因子。

代数会话当前强调“小块子观测的合法根未必可达，子观测换支必须由完整数据接受”。该观点可用于将来离散/连续形状分支搜索；本轮不把它写成已经改善当前目标成像的实证。

## 当前决定

优先把上述结构交给Astra的FFT/k-space及幅度profile审阅。保持R3已冻结主比较不变；若提出新优化器，应另列开发方法和费用，不能替换原失败行。没有必要再运行同一个resolvent恒等式机器精度测试。
