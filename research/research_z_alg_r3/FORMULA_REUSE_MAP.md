# R3 公式复用与拒绝清单

本表是设计审计，不是运行收据。状态截至本理论工作交付：R2与R9全部只读；本工作没有编辑它们、没有新增物理RHS、没有运行数值实验、没有读取私有GT。主线程已授权一次全部接收global-dual/variable-projection DEV pivot，其明确矩阵/梯度接口已交付R3 solver；实现/运行结果尚未由本工作核验。R3实际采用状态须根据源代码与成本收据确认，不能由“写入本表”推断已实现。

## 来源与记号

- R2：同级 `../research_z_alg_r2/FORMULA_COMPENDIUM.md`，辅以 `src/general_quotient_solver.py` 与公开 `experiments/GENERAL6_CONTRACT.json`。`R2Fxx`对应原文`Fxx`，不重新赋予原创性。
- R9公式库（只读）：`Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/formula_collection/` 的01、02、03、05、06、08章与`SOURCE_FORMULA_INDEX.csv`。
- R9复核（只读）：上述研究目录的`audit/theory/FORMULA_COLLECTION_FINAL_REVIEW.md`。
- 汇总atlas（只读）：`Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/maxwell_formula_atlas/FORMULA_ATLAS.md`。已核对其同名公式与接收投影/材料合法性段；不是另一个实验来源。
- A17X源码（只读，由索引定位）：`Gaussian/B1/b1-p3l-r4-conditional-transport/src/a17x/finite.py`，`q_operator`、`exchange`、`source_chart`与`row_fit`。`row_fit`明确返回zero-field的undetermined标记。

**索引限制：** 实际读取的`SOURCE_FORMULA_INDEX.csv`为F-MAT系列提供历史源码路径/函数/身份；本表所列F-OP、F-OPT、F-ERR与F-POST不能伪装成由该CSV逐项列出的源码实现。这些项直接引用对应公式章节，推导项不冒充历史已实现算法。

R3统一记号：N=全部独立实材料数，s=源数，r=复接收数；D=Δ=diag(δ)，固定核为(K,q,L,Y0)，响应为Y0+L(I−ΔK)^−1Δq。R9的全网格N=576、ROI q=256与其latent66是不同符号/支持；对应到R3时，R9的C→K、E_R→q、O→L、F_h→Y0、d→δ。当前初始DEV是N12/s3/r3、18实观测、box[.5,30]，不得沿用R9的σ、576归一化或Gaussian66先验。R3材料RMS应使用自己的N分母。

## 必选公式逐项审计

|公式|成立条件与R3映射|真正节省或改变的计算|采用/拒绝与实现资格|
|---|---|---|---|
|R2F06 固定活跃材料有理闭包|参考M0可逆；U是独立格点选择；这里U覆盖全部N未知格点。K、q、L、Y0固定，Δ变化；无需Δ^-1。|复用几何与参考因子；零参考时M0=I可直接用原G/S/Einc。所有N均活跃时没有N→k降维。|`ADOPTED_IN_THEORY_COMMON_BACKEND`；R3代码实际采用待主线程核对。拒绝把低维参数流形误当低秩对角更新。|
|R2F09 测量先消元|q满列秩；原文单接收，R3扩为Wq=T、W=(L+WK)Δ，全部r行共享Δ；W=W0+CZ，Zq=0。|缓存Z、ZK；新观测只更新W0/V0；替代每源完整电流未知量。C维数r(N−s)复数，不自动低于N。|`GENERALIZED / SELECTED_FOR_ONE_DEV_PIVOT`；实现结果尚未核验。拒绝单接收部分纤维当全部观测解。|
|R2F11 切片Schur复用|补集δB在一次块提案中固定，TB=I−ΔB KBB可逆；Keff/qeff/Leff/Ybase按原公式更新。|避免每个小块重新分解原大参考系统；TB及变化补集仍需实际计价。|`VALID_CONDITIONAL_BLOCK_TOOL`；是否用于ABAR由主线程核对。拒绝把条件块全根当N维全材料域遍历。|
|R2F12 单坐标有理最小化|固定其它材料、一个实增量、原纯二次数据误差/固定白化；Y(t)=Yc+a t/(1−κt)。|将多次forward线搜索缩为二次驻点求根与端点比较；逆更新仍有O(N²)及重分解成本。|`RECOMMENDED_STRONG_CONTROL_OR_POLISH`；本工作未实现。若额外σmin硬门产生可行边界，单纯过滤根不再保证受限最优。循环不保证全局。|
|R2F14 6维商代数|4未知、2独立源、1接收；精确共轭配对；leading T可逆；根/分离/回代与极点资格另查。|4×4归约+6×6特征解代替6条数值同伦路径；仍支付系数、候选过滤和完整验证。|`RETAINED_BLOCK_MECHANISM`；不能直接套N12/s3/r3。未将数值小交换子当完整根证书。|
|R2F15 一般平衡商族|2s未知、s独立源、1接收；leading block可逆及同类共轭/根资格。原R2实际s3/N6，s4仅理论。|小规模代数候选生成；商维数binom(2s,s)指数增长。|`REJECTED_AS_CURRENT_GLOBAL_N12_FORMULA`；当前3源3接收不满足结构。924/12870仅分别属于假设N12/s6/r1与N16/s8/r1，不是当前根界。|
|R9F-OP-008 精确宿主/ROI消元|宿主固定、完整ROI选择、背景可逆；C→K、E_R→q、O→L、F_h→Y0。必须含宿主再散射观测项。|固定背景一次消元，后续只解完整未知区域；与F06同一固定核结构。|`ADOPTED_AS_EQUIVALENCE_AND_FAIR_BASELINE_RULE`。本R3若无已知外区，不虚构576→256收益；精确核必须各臂共享。|
|R9F-OP-010 固定端口导数|L_D=I−ΔK；E=q+Kz；Tport=L L_D^-1；dY=Tport diag(dδ) E。转置/共轭求解约定一致。|完整N材料Jacobian由s个前向RHS+r个端口RHS构造，不用N×s条tangent。|`ADOPTED_AS_STRONG_FULLJ_REQUIREMENT`；当前s3/r3为每状态6个reduced RHS，共享一因子，而非36。不表示本工作已运行。|
|R9F-OPT-002 梯度/Full-J/矩阵自由费用|与上一行同一材料、同一factorization；scalar梯度用s adjoint，full-J用r端口；目标与实化固定。|当前全J传播成本s+r=6；矩阵自由J*J作用约2s=6，在已有场后按每次迭代支付。|`ADOPTED_ACCOUNTING_RULE`；拒绝把“只用少方向”自动算作省物理费。LU、收缩、重启与失败另计。|
|R9F-ERR-001 固定对角化障碍|全部独立对角材料族、跨格Green耦合；[Δ,K]ij=(δi−δj)Kij。|没有可合法记账的通用一次对角化收益；只允许复用固定K。|`ADOPTED_NEGATIVE_RESULT`；拒绝把均匀幅度Schur公式推广到N12独立未知材料。|
|R9F-ERR-003 原方程缺陷|近似状态z~，rphys=Δq−(I−ΔK)z~；真状态误差=(I−ΔK)^−1rphys，数据误差=L乘该式。|可先检查实际缺陷，避免盲目精求；完整Kz~与误差传播仍需成本。|`ADOPTED_VALIDATION_PRINCIPLE`。本备选的对偶对应为E Qδ=Yobs−Y(δ)，不是把E本身当数据误差。没有低秩核实施。|
|R9F-ERR-004 缺陷校正|声明近似预条件器；若有‖Tδ‖≤γ<1才有对应几何误差界；否则Krylov真实残差或完整回退。|可能减少重复精确分解，必须包含迭代与fallback；没有必胜定理。|`NOT_ADOPTED_AS_NEW_R3_ARM`，仅保留回退原则。N12直接求解已很小，不为本工作新增rank sweep。|
|R9F-ERR-005 双边/dual误差|线性任务、正确adjoint；余项界另需‖L_D^-1‖≤β；近似dual随材料变化。|少量任务可用dual针对性评估误差；不能免费为全部格点给证书。|`NOT_IMPLEMENTED / REJECTED_AS_NEEDED_SPRINT_ADDON`；原完整数据验证足够，未经β不声明严格误差区间。|
|R9F-ERR-006 likelihood包络|ε必须是原噪声单位下的白化预测误差上界；包含原packing/σ，不能借R9数值。|可在门槛已经可裁决时停止精求；仅有数值残差不能免费成为上界。|`CONDITIONAL_FUTURE_RULE`；R3实际使用其冻结残差门，未新增概率证书。ψ小不能据此获放行。|
|R9F-POST-001 径向角向覆盖|同一材料欧氏/RMS嵌入；误差²=R²+ρ²−2Rρcosθ。R3实材料应以√N归一化。|不节省forward；防止把较大移动/候选数量当覆盖。GT方向仅封存后解释。|`ADOPTED_INTERPRETATION_BOUNDARY`；本工作不读GT、不作在线oracle选根，也不沿用R9经验数字。|
|R9F-POST-002 随机目标|原公式需要线性Gaussian或局部Gaussianized prior；完整方向随机化data/prior中心并正确保留梯度项。R3产品box先验与非线性散射不符合精确Gaussian posterior资格。|固定局部H可共享分解；继续非线性优化需要真实每状态成本。|`REJECTED_AS_DIRECT_POSTERIOR_TRANSFER`；可另作成熟随机优化基线，但本备选不新增该臂、不把Sobol起点说成RTO/rMAP。|

## A17X/P3L接缝与本轮新增组合

|来源/式|本轮用途|限制|
|---|---|---|
|F-MAT-029与A17X `row_fit`（源码58–72行）|将多源共享材料的逐行LS思想对偶为多接收共享实δ的逐列box拟合；显式保留v=0节点。|本R3为独立实box，不采用P1耦合decoder。若未来换回Full66，逐格clip不再保证decoder像合法。|
|F-MAT-031与A17X `q_operator`（9–11行）、`exchange`（18–25行）|确认(I−ΔK)^−1Δ的端口关系与材料/电流空间不可混同。|没有采用Galerkin exchange或声称它提供材料全域坐标；core/Schur可逆条件不省略。|
|R9第03章F-MAT-015/017/020–026|提醒材料归一化基、Gaussian prior度量、anchor平移与受限校正是不同对象。|R3不借用weak4/weak8限制；所有N材料允许变化。没有把换基当条件数改善或全域支持覆盖。|
|A1–A4（见theory/ALTERNATIVE_STRUCTURE.md）|F09全接收观测仿射空间与共同实材料的精确可行性。|非节点、可逆性与节点单独处理；所有接收必须同一δ。|
|A5–A7|全材料兼容性ALS/variable projection及精确包络梯度。|`SELECTED_FOR_ONE_DEV_PIVOT / IMPLEMENTATION_IN_PROGRESS_AT_HANDOFF`；标准CSI式/可分离最小二乘构件，没有已证速度或全局收敛。|
|A8–A9|精确揭示消元目标为材料相关的加权观测误差，给近极点解析反例。|这是拒绝仅凭兼容性残差接受候选的依据；不把辅助目标称为原likelihood。|

## 主线程可据此作出的有限决定

1. 共同精确后端、正确Full-J成本和原观测验收应作为所有臂共同要求；这是公平实现条件，不是额外实验臂。
2. 如果预算已用于ABAR，保留本备选为`NOT_RUN`是完整且诚实的理论结论，无需为了本表增加实验。
3. 若选择执行备选，使用全部3接收的共同材料varpro/ALS，公开同一起点，原目标校正，完整失败与费用记录；不把初始单接收平衡讨论误当当前合同。
4. 小块代数、全材料辅助目标、F12精确坐标步与全材料全根覆盖是四个不同声明。前三者的可实施性不证明第四个。

历史工具先例：本构残差路线已有[CSI（van den Berg/Kleinman, 1997）](https://doi.org/10.1088/0266-5611/13/6/013)，线性内层消元已有[variable projection（Golub/Pereyra, 1973）](https://doi.org/10.1137/0710036)。本次仅核验出版社元数据/摘要；本表不作新颖性裁决。

## 主研究者最终实现核对（科学封账后）

本节更新上文交付时的状态，不改写当时尚未得到的实验结果。

|公式/接口|最终实际采用|证据与边界|
|---|---|---|
|R2 F06/F11|IMPLEMENTED|src/physics.py、src/solver.py：固定核、精确Schur、完整观测约化评分；全部费用保留。|
|R2 F12|IMPLEMENTED|coordinate及ABAR预扫描；二次驻点和端点；有限坐标对照已运行。|
|R2 F14/F15|IMPLEMENTED_SMALL_BLOCKS|DEV运行4/6活动材料；V6只运行4材料Q6。没有全12材料商环。|
|R2 F09全接收扩展、A5–A9|IMPLEMENTED_AND_BLIND_EVALUATED|globaldual及V7保护组合；全12材料、3源3接收、QR/梯度、原物理评分与修正。属于CSI/SOM与变量投影的专门化。|
|R9 F-OP-010/F-OPT-002|IMPLEMENTED_COMMON_BASELINE|完整N12系统每状态3源RHS加3接收端口RHS，共用LU；Jacobian并非免费。|
|R9 F-ERR-003原则|IMPLEMENTED_VALIDATION|独立完整矩阵验收；没有低秩传播近似或严格误差界。|
|R9 F-POST边界|ADOPTED_IN_SCORING|在线top1只看原数据；GT仅在冻结后用于RMS/best@K。无posterior声称。|
|latent/characteristic连接L1–L4|THEORY_ONLY_NOT_IMPLEMENTED|见theory/LATENT_AND_CHARACTERISTIC_SUBSPACES.md：decoder链式梯度、更新秩反例及嵌套子空间目标不等式；无Full66实测。|

阶段结果分别见results/HOLDOUT_V6_SUMMARY.md与results/CONFIRM_V7_SUMMARY.md。实现完成不把HOLD的新颖性或覆盖判断改成PASS。
