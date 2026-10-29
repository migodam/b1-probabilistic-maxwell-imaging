# R3-DIGITS：十个有用公式的来源与实施边界

状态：2026-10-29 19:38 UTC 静态独立审查。不超过十项；历史文献沿用相邻 R3 已核验材料，不把本表称为新一轮穷尽检索。当前实现与最终成绩分开：代码存在不等于盲成像目标已完成。

|公式及作用|来源与必要条件|当前实际路径与成本/限制|
|---|---|---|
|1. `J=(I-DG)^-1 D E0`、`Q=E0+GJ`、`F=SJ`|`../research_z_alg_r2/FORMULA_COMPENDIUM.md` F06；`../research_z_alg_r3/MATHEMATICAL_THEORY.md` §3–4；固定几何、非极点。|`src/physics.py::Model.forward`；一次576阶 LU＋12源 RHS；完整多次散射，无 Born/对角交换。零参考复用不构成物理降阶。|
|2. `W E0=Y`、`W=(S+WG)D`；转置为 CSI/SOM 状态/数据方程|`../research_z_alg_r3/theory/ALTERNATIVE_STRUCTURE.md` A1–A2；`../research_z_alg_r3/audit/PAPER_NOVELTY.md` §2；E0满列秩、物理可逆。历史来源：[CSI 1997](https://doi.org/10.1088/0266-5611/13/6/013)、[SOM 2010](https://doi.org/10.4310/MAA.2010.v17.n4.a6)。|`MATHEMATICAL_SCALING.md` D1的推导来源；本轮不显式存 W、C、Z，没有每步576×564 QR。接收对偶不需要 Green 对称或物理互易性，也不是独立首创。|
|3. `psi=tr[R (Q*Q)^-1 R*]`|相邻 R3 A8–A9；本轮 `MATHEMATICAL_SCALING.md` D1；源满列秩和非极点。标准 variable projection：[Golub/Pereyra 1973](https://doi.org/10.1137/0710036)；参数 Gram 度量直接先例：[van Leeuwen/Yang 2025](https://doi.org/10.1088/1361-6420/adab89)。|`src/physics.py::Model.smallgram`；12×12 Gram 及 Cholesky。与硬数据消元目标精确同值，但不是原 likelihood；近极点可虚降。|
|4. `Q=UT`、`psi=||R T^-1||²`|`MATHEMATICAL_SCALING.md` D2；Q满列秩。标准薄 QR/三角求解。|`Model.smallgram` 的条件/残差触发薄 QR回退，尺寸576×12；主路径是小 Gram，不是每次薄 QR。无 jitter；物理 LU 仍付费。|
|5. `M*Λ=S*A-G*(Q A*A)`、`g_i=2ReΣ conj(Λ_ia)Q_ia`，其中 `A=R H^-1`|本轮 D3–D5，从上述精确消元和标准隐式/伴随微分得到；实材料、固定参考归一化。|`Model.smallgram`；复用物理 LU，12伴随 RHS；保留 Gram 导数项。省去 Full-J存储，不保证比普通数据伴随更便宜。|
|6. `T=SM^-1`、`∂F_ra/∂m_i=T_ri Q_ia`|本轮 D7；相邻 R3公式复用表 R9 F-OP-010/F-OPT-002；相同保留 LU。|`Model.jacobian`；48普通转置 RHS，1152×576实 Jacobian。`Run.gn` 另支付阶段法方程与 LM 解；不是逐材料有限差分。|
|7. `m=P_h c`、`g_c=P_h^T g_m`；倒数系数 `c=1/x` 给 `J_x=(J_m P_h)diag(-c²)`|本轮 §5和D8；P正、行和1，最终P=I；box正。多尺度与坐标链式法则是标准构件。|`src/multiscale.py::material_basis`、`resample_material`；`src/solver.py::Run.gn`。两普通臂使用同一个物理粗空间；不是 `m=1/(P_h x)`。576物理系统从未缩小。|
|8. `R(m)=λ_s mean(d²)+λ_TV mean(sqrt(d²+ε²)-ε)`|`src/multiscale.py::Regularizer.evaluate`；固定物理差分 `d=Dm`、固定权重和epsilon。标准平滑TV/二次正则。|每次在24×24物理图像计算；三臂共同排名和GN使用。GlobalDual提案乘固定参考场能量a/(2||Y||²)，a来自初始均匀材料1的已付费总场，不随阶段重置；源方向Gram度量仍区别于原数据目标。|
|9. `argmin_c [.5||F(m_c)-Y||²/||Y||²+R(m_c)]`；独立验算 `B E=E0, F=SDE`|相邻 R3关于原目标验收及近极点反例；本轮 `MATHEMATICAL_SCALING.md` §6、`BLIND_PROTOCOL.md`；仅付费、合法候选，GT禁止参与在线选择。|`Run.score/finish`、`experiments/score_predictions.py::independent`。top1与GT oracle best@K分列；验算每端点额外一次576阶LU＋12源RHS。是同模型独立方程求解，不是连续/跨网格认证。|
|10. 普通归因对照 `M*Λ=S*R`, `g_i=ReΣconj(Λ_is)Q_is/||Y||²+reg_gradient`|`MATHEMATICAL_SCALING.md` D6与半平方归一化；标准数据伴随，实材料、相同物理前向。|`data_adjoint_control.py::DataAdjointRun.dual`；12源+12伴随RHS、零Gram；同L-BFGS/Full-J家族，后生成前GT补充。与budget-filled GlobalDual比较仍有重启/计时差别，非严格单因素逐位对照。|

没有在576材料上复用 Q6/Q20 商根枚举，也没有获得全材料根覆盖。小 Gram 是已有硬数据 CSI/SOM式变量投影的计算重排；本轮可检验的增量是规模、数值实现和冻结图像任务表现。最终是否有材料恢复收益须读全九数字 top1、失败、预算和 oracle差距，不能由本表代替。

19:48 UTC更正：旧表关于阶段初始psi归一化的文字错误，已依据封存V1和当前V2源码修订。V2共同native粗路径、fine初值一致性和主实验前精确计算重排的静态审阅见 `audit/SCALING_REVIEW.md` 第7节；不增加本表公式项数。

主实验源/配置冻结见 `audit/HOLDOUT_FREEZE.json`（19:47:38 UTC，16文件）；独立静态审阅确认快照与当前文件16/16字节相同、三臂各90秒。性能证据须等待27份冻结预测及九例评分。

后生成、前GT的budget-fill与ordinary data-adjoint补充源码静态审阅见 `audit/SCALING_REVIEW.md` 第8节；不改变原27份主证据身份。补充累计root成本与新增incremental成本必须分开。

父线程补充一次576状态旧宽QR与小Gram表示内核计时，见 `audit/SCALING_KERNEL_BENCHMARK.json` 和审查第9节：单次wall比2.7005，固定/逐目标setup另列，不升级为成像或普通伴随普遍加速。

核心63份已全部封存后评分：公式1–9在原主/补齐臂实际运行，公式10在ordinary九例运行且Gram=0。实际576细层导数/公共fine初值一致性通过；数字9有额外恢复，但其余六例图形失败，32个存档均只有1个RMS≥1贪心代表。因此公式复用不支持全域分支/九数字全恢复/新代数定理。读 `audit/FINAL_SCIENTIFIC_REVIEW.md` 第9节；源方向白化解释D9及固定材料噪声放大D10不新增本表项数。

最终科学冻结（20:41:43 UTC）后的解释边界：三TRF诊断及坏端点Jv/独立总场验证未揭示所检查点的简单物理/Jacobian不一致，但没有恢复原高对比失败；同公开fine初值TRF残差0.39096仍远高于原Gram9形0.0002283。两低幅度控制可拟合，fine各1J/0接受步，属于改变Y的事后敏感性实验，不是原任务补救。噪声14完成/13NOT_RUN，无部分均值。公式实现PASS不升级为成熟普通基线支配或新定理；最终裁决见 `audit/FINAL_SCIENTIFIC_REVIEW.md` §12及 `audit/FULLJ_CROSS_SESSION_REVIEW.md` §6。
