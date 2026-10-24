# 冻结几何的概率价值审计：只准备，尚未运行

新增src/p3r4/frozen_geometry_pcn.py及tools/r4_geometry_probability.py，保留全部封账StageIII结果。源资产为dynamic_geometry_first/R4DEV_00000…00007/d66_r4的BASE_FRAME.npz和BRACKET_h0.001.npz，另用其原已付context cache及原在线观测。源账本须已CLOSED COMPLETE，q与cache相同，Dchi须与原P1 decoder在Gaussianized66坐标的纯导数相同；核对decoder不产生新Maxwell/J。缺失资产明确MISSING并停止，绝不重算教师。

当前canonical准备目录为AUDIT/geometry_probability/PREPARED_V1_VALIDATED；先前PREPARED_V1是同矩阵的开发准备副本，正式调度使用VALIDATED入口。该目录已保存8×4个固定Q/rho/sqrtQ、来源byte/mtime身份、source derivative核对和完整FROZEN_PROTOCOL.json。新增物理RHS、J及训练均0。

四个形状为Ciso=I、Cweak=.05I+VwVwᵀ、Cbracket=Cweak+bbᵀ、Crandom=Cweak+rrᵀ；b是已付full66 bracket的Euclidean单位向量，r是固定seed单位Gaussian向量，二者rank1相对权重均1，不把bracket强行投到weak正交补，也不让随机对照使用结果依赖方向。Vw为已付rank4正交弱frame。每个C都满秩，不把其余62维补成零。

令Q=αC，rho=sqrt(I−Q)，x'=rho x+sqrt(Q)ν，ν~N(0,I66)。Q、rho和sqrtQ都在每parent开始前冻结且不随当前state改变。全部66个Q特征值严格为正且小于1；约束最大值≤.05。因为rho对称，rho²+Q=I，联合Gaussian(X,X')的协方差为[[I,rho],[rho,I]]，交换两端不变，故核对原标准Gaussian66 prior可逆。显式full-rank Gaussian提案密度满足log q(x|x')−log q(x'|x)=0.5(||x'||²−||x||²)，与target先验项严格抵消，MH只需原Full Φ差。该恒等式在运行时每一步也核对；没有伪行列式、退化Euler提案或动态frame适应。

归一化不能只配噪声trace。以源D=concat(Re Dchi,Im Dchi)/24，在原在线anchor x0处配

E||D(x'−x0)||² = ||D(rho−I)x0||² + tr(DQDᵀ) = .002²。

标量α由纯代数求根；若λmaxQ=.05仍达不到目标，固定cap、明确cap_binds和未匹配，不声称同RMS。当前32个实际准备核均匹配.002且未触cap。这只是anchor的一阶期望尺度，不能说沿整条链每一步finite材料RMS都被匹配；实际提案和接受移动RMS另存。

父调度器后续才可明确release：

```json
{"released":true,"experiment":"frozen_geometry_probability_v1","campaign_id":"P3LR4-20261025"}
```

运行入口：`<python> tools/r4_geometry_probability.py run --prepared CANONICAL_PREPARED_DIRECTORY --release RELEASE_PATH --max-wall PARENT_BOUND_SECONDS --out NEW_RESULT_DIRECTORY`。max-wall必须显式给出，BookIII控制；physics.lock及3GiB峰值上限保留。固定8parents×4kernels×1024steps，每候选付8复源Full RHS，共硬上限262144RHS及32768次Full LU，不扩次数救结果。源frame/bracket/J不重新查询。已付教师459712RHS及其StageIII所有权单列，不把它们称为免费、不再重复收费。

四核从同一已付在线anchor开始，各自重置SeedSequence([20261025111,parent_index])；每步normal66后uniform1，故同parent的创新与U逐项相同。因Q不同而路径不同是预期，不作路径相等断言。不读GT，不以residual挑候选，不改变原prior/noise/decoder。每个Full候选立即由Target及本实验保存q/χ/current/field/prediction/Φ/prior/accept/U，接受后或拒绝重复状态立即保存相同全字段和材料函数。CHAIN含初值x路径、全部reject重复、Φ、prior能量、提案及实际接受材料移动RMS。

材料函数固定为ROI实/虚均值、fraction≥.5体积比、离初始anchor的全网格RMS及4个固定材料投影。各parent各核仅一条1024步轨迹；保存初始正相关pair ESS、MCSE及单链前后两半split-Rhat，但一律FINITE_RUN_UNQUALIFIED。常数函数ESS/MCSE不评估而非报零成功；不同parent是不同target，绝不跨parent造Rhat；不同核用共同创新也不是独立链。没有共同posterior reference，不能由这些短链指标认证FULL质量、MCSE或采样混合。

KERNEL_MATH_CHECK.json验证32个准备核及320对显式Gaussian提案比：最大prior抵消误差约5.39e−11；所有Q满秩、rho²+Q=I、anchor材料期望RMS匹配，人工cap例子正确披露未匹配，退化Q被拒绝。该数学/接口通过不代表几何核有概率或费用收益。动态bracket非零只解释候选几何，不能证明任何核会更好；实际结论由父调度器在全费用和负结果下判断。
