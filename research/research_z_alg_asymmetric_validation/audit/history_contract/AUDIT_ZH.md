# Z-ALG-E1 历史合同与结论防火墙审计

审计日期：2026-10-29。范围为继承的 R6/R7 冻结文档与 Z-ALG 输入包。此审计只核对历史对象、声明、费用和外推界限，不作本轮最终科学裁决。没有运行物理、训练、历史重放、网络检索或 truth 读取；没有修改旧来源。机器表见 `PROVENANCE.json`，每项给出仓库相对路径及一基行号。这里的“已完成”均指来源所报告的历史完成，不代表本次重新验证数值账本。

## 1. 两份合同必须分别命名

|项目|原 P1 / Full66（R6、R7）|新 Z-ALG-E1 小型 Green VIE|
|---|---|---|
|材料来源|Gaussian64×Uniform2，经精确 probit 的 Gaussian66；quintic 平台、host、被动复材料|本轮另行冻结的 6–8 个实材料自由度及明确新 prior/support|
|材料维数|latent 66；decoder 像处于 576 复网格材料空间，不能以环境空间维数替代 support 维数|n 为独立实材料坐标数；不是 Full66 的边缘、重参数化或自动等价限制|
|离散物理|24×24 scalar-TM/VIE，576 全网格，256 ROI|小型真实 Green 系数 VIE；几何、pulse/self 项、频率、标签、材料域须自行冻结|
|采集|8 源×32 复接收，256 复坐标=512 实坐标；source-major|须在运行前冻结全部源/接收标签；k 必须是独立实观测坐标数，若有 Ns×Nr 复通道，原始上限为 2NsNr|
|噪声|sigma=.00119211798514519；R7 source-major sqrt2/sigma 白化|新 covariance/实化约定；不可沿用 sigma 名称却更换 complex variance 含义|
|目标|Full66 材料弱响应、有限噪声竞争性及实际方法比较|合法内点 n>k 的正则局部纤维及 current/material chart 可计算性|
|可转移证据|历史边界和普通强基线的必要性|新模型上的结果只属于新合同；不能回填旧 GO/NO-GO|

输入报告建议的是非对称 **4–9 cell** 新系统；本轮“6–8 实材料自由度”是另行缩定，不能称输入包已经执行。cell 数与自由材料维数必须逐项写清：若每 cell 有复自由属性，实维数通常翻倍；若固定相位仅有实幅度，则每 cell 一个实自由度。冻结 n>k 时不能只数复通道而把 k 少算一半。

原 Full66 与输入报告另述的 Full576 同样不同：576 固定复相位实幅度、512 实观测得到的条件维数至少 64，属于**新的 support/prior**，不是原 P1 的 66 维 decoder 像结论。CLAIM_LEDGER C7 明确是条件维数定理，未执行新的 576 实验。

## 2. 四层物理/几何对象不能互证

1. **同材料 decoder gauge**：不同 latent 映到同一材料。R6 的 V=QD 双属性补偿达到材料 RMS 约 1.48e−16；R7 特定无 pure、transition 行可调整 strata 的同材料路径最大约 2.42e−16。这些重复表示不会生成新材料世界。R7 同时存在材料 Jacobian 满 66 的 regular 区域；没有通用一维尺度商。数值 SVD 不能无条件替代 exact rank。
2. **真正不同材料**：必须材料距离非零、满足同一支持、被动性、共享材料构成关系和稳定 resolvent。receiver-dark 电流只是接收方程的解；多照明须有共同标量材料，field node 须单独处理。CLAIM_LEDGER C3 的反例只证明暗电流可不合法，不证明任意暗方向可生成材料纤维。
3. **有限噪声歧义**：R7 Full-J 的 0.02 材料 RMS 端点 8/8 likelihood 竞争，其中仅 7/8 另满足 pair≤1；两个计数不合并。五节点弦/路径通过只是 sampled-path 资格，既非节点间全区间保证，也非精确 F(m)=F(m0)。近噪声邻接不传递，不能称群轨道或等价类。
4. **精确连续不同材料纤维**：新 n>k 小系统若在合法内点、无 pole 的邻域具有常秩 r，则常秩定理给局部维数 n−r≥n−k>0；若证实满行秩 r=k，则一个可逆 k×k minor 在邻域保持非零可给 IFT chart。单点小奇异值、浮点残差阈值、有限候选或单点 rank 计数都不是该邻域条件的严格证书。奇异点或边界可能截断纤维，正维不推出显著幅度、所有分量覆盖或 posterior mass。

R7 的八中心材料有效响应数值秩 66、最小奇异值约 .0308–.0956，所选单位材料速度一阶响应 norm 约 8.33–16.53，并非 exact null。条件数学满列秩给局部单射，不能由此排除远处分离离散多解；也不能把新低维欠定采集得到的连续纤维当成推翻历史满秩诊断。

## 3. 概率、防重与 coarea 资格

R6 的材料商密度包含 Q^−2；host 原子、pure 与 transition 的重叠像应加权相加。局部 GH15→21 收敛不证明全域 tails/外部分支质量；二维/三维活动 Gaussian 截断于 [−4,4]，不是 Full66 marginal。低维 host 约 22%不能替换 Full66 prior host 约 .0688%，Full66 posterior 原子质量仍未知。

R6 固定材料条件联合的独立边缘乘积给 v E(Q)E(1/Q)，会制造错误材料。R7 四个非均匀 Full66 受控等权离散纤维也验证同类 mapping 效应，但等权四点不是原 prior 条件核；不能称已恢复 Full66 posterior/copula 质量。R7 两点 KL/TV/风险界仅为未来固定候选对的条件两点子实验，不是同数据搜索所得的 posterior credible interval 或原 prior Bayes risk。

新 chart 上须从**新声明的原材料 prior** 推导条件测度。若实秩 r 的 quotient π 在光滑材料层为 submersion，则条件密度相对于 H^(n−r) 为 p_M/Jπ，再按整条目标纤维归一化；Jπ=√det(DπDπᵀ)。若原始观测维数大于秩，不能把奇异 det(DF DFᵀ) 当分母；须用真实像坐标或适当的非零奇异值乘积及对应像测度。以图册参数化时还须面积 Jacobian、重叠去重/partition 与全部分支质量。仅在一个 chart 上归一化给局部条件分布。

有限 regular 根的权重为 p(m_j)/|det Dπ(m_j)| 再归一化，既非一根一票，也非按 solver 命中频率。跨不同维 strata、原子和奇异 prior 应使用原测度 disintegration，不能硬塞同一 Lebesgue density。输入 C9 是标准理论应用，不是新数学发现。

只有 likelihood 沿 π 纤维常值，才能复用原 prior 条件核 K；否则须条件 posterior K_y。有限噪声的完整输出还要对 quotient posterior ν^y 积分；固定一个 noiseless F(m0) 的精确纤维积分不是整个 noisy posterior。输入 C4/C5 的两 cell E4 是有限噪声、明确新 prior 的数值积分阳性/错误权重对照，非严格区间积分证书，不能替代本轮非对称 E4。

## 4. 历史 gate 保留

|历史声明|冻结状态与准确边界|
|---|---|
|R6 概率 GO-A|REDIRECT；Q1/Q2 局部证据，无 Full66 posterior、连续 TV 或任务认证成本分离；4D 未扩大|
|R6 成本 GO-C|原 global8 两噪声通过 1% NI/20% wall；强化 Full-q 低噪 NI 未确立，同窗口低噪 20% 节时下界仅 16.82%；不改历史通过或升级强化结果|
|R6 GO-TAP|REDIRECT；普通 direct Probit 有竞争力|
|R7 GO_HIGHER_ORDER|NO_GO；C Jets 6/8 覆盖而 B Full-J 8/8，C另有 setup；只针对本实现，不否定所有高阶理论|
|R7 GO_PHYSICAL_FIBER|NOT_ESTABLISHED；非零 c1 与有限节点不满足 exact fiber 条件|
|R7 GO_FINITE_NOISE_AMBIGUITY|SCOPED_GO；0.02 RMS、B8/8竞争且7/8 pair≤1|
|R7 FULL66_POSTERIOR|UNQUALIFIED；并未以 12 个 DEV/新上下文认证全后验|
|R7 TASK_IDENTIFIABILITY / GO_TAP|LOCAL_OR_EMPIRICAL_ONLY / REDIRECT|
|输入 C8|非对称一般分支覆盖及加速 NOT ESTABLISHED，下一步实验而已|

## 5. 实际冻结费用与建议预算分开

|来源|已报告实际费用|不可混写为|
|---|---|---|
|R6|6,860,984 complex RHS / 8,000,000 cap；210,355 state/LU；33 新训练阶段；A158,144，B0，C2,864,200，D3,838,640 RHS|不能把历史 R5-B 1,999,504 再计入，不能作为本轮免费数据生成费用|
|R6 时间|物理 job wall 5147.36s，solver kernel 3310.51s；保守 MPS 训练 .2414 device-hour|不是全部编辑、文献和 I/O 的完整 profiling；不同优化配置的 1.73s 不能与 R7直接换加速比|
|R7|37,456 complex RHS / 50,000 cap；10,600 forward + 26,856 derivative；1,325 state/LU；19 物理 job；物理动作失败0；无网络训练|不是 50,000 已耗用，也不是科学负结果等于运行失败|
|R7 时间|物理 job wall72.21s、动作wall14.77s、CPU67.47s，峰值约190MB；零物理材料秩 E1另7.18s、490MB|不是所有人工/检索交付时间；两类计时不可冒充同一求解器 wall|
|R7 单位计费|新 state8 RHS，首次 receiver-J32 RHS，每阶 jet8 RHS；C setup24 RHS，四档摊销每档6|拒绝、恢复与失败路径必须计费；缺端点而费用低不等于加速|
|输入两 cell witness|8 次各2复RHS前向，1次2RHS灵敏度，另1次证书矩阵逆、1次实2×2商求解、2次独立Green求积；优化/训练0|18个显式 solve RHS并未包含证书逆矩阵的隐含列求解，勿虚构合计统一账本|
|输入两 cell E4|15,360 likelihood；direct解析响应10,245，quotient有理响应5,125；另5次2×2 factorization、10复RHS；CPU .01107099s、wall .01107962s|解析 likelihood 不是免费；其毫秒时间不能宣称一般 Maxwell 加速|
|输入后续建议|每run≤10,000复RHS、2CPU小时、8GB；整轮≤30,000复RHS；建议20%效率门|这是建议预算/门槛，不是本轮已付费、自动继承的最终冻结门或成功收据|

本审计自身新增物理 RHS=0，训练=0。以上费用为历史文档/JSON抽取，没有重跑中央账本逐项对账。输入包已有摘要字段中的 SHA 字样仅属来源元数据，本审计没有执行摘要校验。

## 6. 对新结果的合规表述边界

若新系统找到合法不同材料和数值 chart，可写“新合同的浮点局部纤维候选/图册”；若另有常秩/IFT 条件证明或区间证书，才能按证明范围写精确局部纤维。若只完成几何，概率与净费用优势仍为 NOT_RUN/NOT_ESTABLISHED。新系统成功不能改写 R7 的 exact fiber 未建立、Full66 posterior 未合格和高阶 NO_GO。

新 current chart 与直接材料 Full-J/IFT chart 应在同一中心、prior、材料域、观测、覆盖范围和任务误差下比较；结构构建、校正、失败、重复覆盖、认证与积分均计入。n>k 保证的局部维数本身不是新算法贡献；最终科学新颖性、证据是否支持主张与 GO/NO-GO 由主 Codex 判定。
