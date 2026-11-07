# R10 实际三维合同：18臂，GT仍封存

当前实际是四个材料的 matched-grid 三臂，加两个材料的18³→12³ cross-grid三臂，共18。12个设计、36条主反演只是未执行的初始草案，不能算作已完成样本。可执行配置为 `configs/vector_actual_contract.json`、任务表 `configs/vector_all18.json`；初始提案另存 `HISTORICAL_PROPOSED_CONTRACT.json`。旧R1–R9不修改。

|材料|几何与介电常数|matched平台|cross18→12|
|---|---|---|---|
|v02|平移实心块，边.3λ，ε9|Mac|Mac|
|v09|外边.4λ ε2，芯边.2λ ε9|Mac|Mac|
|v04|双块各边.15λ，ε9/ε2，gap.1λ|Mac CPU-dense（旧Windows仅engineering）|无|
|v08|立方壳外边.4λ、腔边.2λ，ε5|Windows CUDA retained-dense（旧FFT仅engineering）|无|

v04低介电体在任何v04生成/评分前由ε1.5改为ε2，理由是明确high/low任务和固定support χ>.5，不是按GT救结果。已有v02/v09完全不改。未来v11/v12含χ=.5的设计仍是未执行草案。

所有反演箱宽.6λ、12³=1728个独立实χ未知量、3N=5184个每照明复矢量场分量，60源完整J/Q数组各含311040复分量；χ∈[0,8]，m0=.25，R=0。60个外部电偶极，三圈各20个，法向/极化分别z/x/y；固定radius=.9λ，中心(.04,0,-.02)、(.02,0,-.03)、(-.01,.02,0)。同60位置接收全部三个散射电场分量，Y为180×60复矩阵。细网格18³使用完全相同的绝对源/接收几何，不随网格缩放。

主DPCS是实际S的economy SVD，公开iid CN3%标量σ，E|n|²=σ²；完整噪声球radius²=σ²χ²(.99,21600)/2，取最小数值支持的连续左奇异前缀使discarded Y能量不超过该球。一个full-noise事件同时适用所有投影，不做逐rank显著性/GT调rank。历史.1阈值只留rank/能量诊断。噪声球不界定cross-grid模型误差。噪声实现不按sample norm重标定。

CIEsource、GlobalDual、A+global三臂同MF2/3/4、各max500，至少20步后才允许sqrt(mean_source ||Δγ||²/||γold||²)<1e-3停止；min20为R10公开保护条件，不冒称论文要求。原生projected-IDFT、频率嵌套transport，无γ白化。material β6/3/1；sourceCIE current同6/3/1、sourcewise PR；GlobalDual/A+ current均6/4/3、global PR。**CIEsource与A+的直接差异同时含current β与聚合；GlobalDual才是匹配current规则的材料A+消融。** 不加NIE/Twofold/新3D算法。

A+为实各向同性、同一标量self的条件Mobius弧，两轮8色块坐标，并按完整固定current的A块loss接受CIE seed与候选；seed/root/candidate/bound全部一致χ上限8。没有复材料、张量材料或全Maxwell秩一结论。每阶段保存真实联合m/J/Q/γ，不能用physical J替换辅助J。

共同原Phi物理L-BFGS-B最多80 accepted、strict200 objective requests、maxls12，ftol1e-12、gtol1e-8，存40/80接受步检查点。选择最低已付公开corrector Phi，GT不选终点。maxiter是预算停止，不是驻点证书。三臂只共享immutable G/FFT/E/S/SVD和同publicY的DPCS，材料相关LU/J/gradient缓存和优化历史每臂重置；臂内同材料的value/gradient复用LU，全部60源与60伴随为同LU的多RHS。

## 实际物理约定

时间因子exp(-iωt)，λ0=1、k=2π、μr=1、Z0=1。令n=r/|r|：

`g_k(r)=exp(ikr)/(4πr) [k²(I−nnᵀ)+(ik/r−1/r²)(I−3nnᵀ)]`。

非self `G_ij=V g_k(r_i−r_j)`，`G_ii=s I3`，`s=−1/3+i k³V/(6π)`，V=h³。`J=χ Q`、`Q=Einc+GJ`，J为归一化对比极化密度/算法contrast-current坐标，实际极化电流为`j_phys=−iω ε_b J`，偶极矩`p=ε_b VJ`；本归一化ε_b=1才简写p=VJ；exciting field排除self，`Eexc=Q−sJ`。对应CM+RR极化率 `α=Vχ/(1−sχ)=3Vχ/(χ+3−i k³Vχ/(2π))`。源Einc为中心点电偶极场；receiver块为`V g_k(receiver−cell)`，只接收散射场。没有Born/标量TM/单分量替代。

实际CPU物理来自R9 `paid_density_state` retained SciPy LU，K=I−DχG，所有源共享一次分解；adjoint为K†，不是另一个LU。辅助传播沿用R9 padded dyadic FFT。CPU的dense primal适配明确只支持χEinc，不伪装任意rhs求解器。Windows只占自有P12/14、BelowNormal、BLAS1，connected foreground JobObject；GPU在用户释放并检查空闲后使用。Mac保持单numeric worker/BLAS1。最终v08 matched整组三臂为Windows CUDA retained-dense，v04整组三臂为Mac CPU-dense，v02/v09 matched及cross仍Mac CPU-dense。旧Windows CPU与CUDA FFT/GMRES仅为独立engineering attempts，不进入primary材料配对表。

## 电尺寸与物理限制

inverse h=.05λ，ε9内部λ=1/3λ，只有6.67格/内部λ；fine18为10格。这对高介电仍偏粗，不称充分连续分辨。v02边.3λ=.9内部λ；v09高芯.2λ=.6内部λ、低壳厚.1λ=2inverse cells；v04高体边.15λ=.45内部λ、gap.1λ=2cells；v08 ε5外边.4λ=.894内部λ、壳厚.1λ=2cells、腔边.2λ=4cells。它们不是原2019的3λ箱/.8λ块/30³→60³物理复现。

两例cross-grid独立组装18³并重新求解，不是同LU/同预测加噪声。体素本构为精确轴对齐box体积分数的线性χ混合；不混用旧Mie球的MG。保留R9 ε9球CM+RR误差大、IGT未修复的阴性证据，不重复失败Mie系列。新cross-grid只能量化本组模型/离散差异，不是连续Maxwell误差界。

## 范数能量与评分隔离

S右模式与G模式作用于每照明complex3N电流列，完整60源矩阵有311040复元素（这不是独立优化DOF），材料是realN，禁止直接把标量χ投到dyadic3N模式。均匀体积的电流范数能量为h³||J||F²；它不是完整电磁储能/Poynting功率。记录receiver retained与正交补空间、实际projected-IDFT APCS范数、未投影γ范数及实际affine重建误差，不能把γ范数当projected current能量。

G背景-only部分谱计算请求top32，CPU2.126454s，无物理forward/GT。G含上述CM/RR self，不是Goff。返回σ约1.05696–1.23476，最大Gram相对残差5.15e-8，正交误差2.55e-14。只标top32 captured与5152维orthogonal remainder/uncaptured；没有full5184谱、.1全strong或严格kernel结论，截断简并空间也不保证完整。该数值子空间不作为新反演算法。

评分合同预先固定 `configs/vector_scoring_contract.json`：χ>.5同时用于effective-cell GT和预测，Otsu仅诊断；连续χ相对L2、amplitude-normalizedRMSE、paperErr及低/中/高有效χ-bin误差均须报告。严格阈值不是连续几何support，部分填充体素可被排除。所有18个终态封存后仍须root明确授权才能读新GT评分；当前没有材料成功/通用胜出结论。

## 资源与账目

GPU方案原1.5GB是条件GPU工作集政策，**不是CPU内存限制或实测显存**。自己的GPU探针在CUDA初始化前等待外部视频；778.84秒等待后转CPU，0CUDA context/0科学RHS。保留4.03125 CPU秒主进程、2.0 CPU秒77次query子进程及传输/停止费用，非算法失败。Mac公开m0=.25完整CPU后端探针为wait4 4.416365CPU/4.591230wall、峰RSS1.903GB；其中sharedLU+60forward约2.533CPU、sameLU60adj升级约.184CPU。实际case峰值/全成本以各process回执为准，不能把benchmark峰值称所有任务上限。

首次CPU探针在forward完成后的回执组装处失败，保留1LU/60RHS已付、CPU未知，不能伪记0。wholeprocess wait4/GetProcessTimes优先于nested setup/perarm组件；跨host不合成单平台公平速度，墙时间不相加当makespan。暂时外置盘I/O停顿保留在wall，不为此重启/重投任务。所有GT仍封存。

当前与历史资源分开：778.84秒/0CUDA context是已终止的历史等待，不代表现在没有GPU任务。用户释放GPU后，完整公开uniform VG探针4.953125主CPU秒/8.375137wall，warm forward+adj4.54050wall，原残差1.43e-8，tensor峰reserved128MB；m0弱点不外推高χ全反演速度。root据此授权v08完整三臂CUDA，以首完整CIE实际任务判断部署成本。CUDA-event记录默认stream时间戳区间（含host/IO/pause间隙），不是active kernel或SM time；主CPU、query子CPU、GPU-event、峰内存和实际算量分账。

DOF补注：共享材料只有1728个real未知量，不随60源重复。MF2/3/4的γ数组名义列数为192/648/1536每源，数组大小分别11520/38880/92160复元素；只有actual projected-basis数值rank×60才是自由APCS current复DOF。该rank在RESULTS/vector/representation_rank用公开S及实际native Fourier投影单独核算，不把存储元素数当rank，也不增加任何forward。该background算子核算不读取目标GT。

已核算的实际representation rank：对公开S完整180维receiver rowspace，MF2/3/4的(I−UFUF*)V最小奇异值分别.053731/.023985/.010418，远高数值cut1.15e−12，故其所有DPCS前缀与UF均无数值intersection。实际B rank=192/648/1536，×60后为11520/38880/92160 complex自由APCS DOF（real为两倍）；四个已在Mac存在publicY的cohort逐一同结果。计算1.121519 CPU秒，0physical RHS/GT。这是数值representation计数，不是场信息/材料可辨识性定理。

最终工程迁移修订：在未读GT下，root按公开m0与完整CIE成本将v04整组三臂重跑Mac CPU-dense，v08整组三臂重跑CUDA retained-dense，原publicY和σ字节完全复用，绝不重新造噪声。旧Windows CPU v04、旧CUDA FFT/GMRES v08的完整/partial轨迹、已保存stage/40/80、全部CPU与中止尾部费用保留在RESULTS/vector_engineering_attempts。低于40且未序列化的paid尾部只保留TRACE/算量下界，不能冒称已保存最新材料。新primary从同m0=.25独立重跑完整500×3/min20＋80，不把旧partial当额外对象或新的primary。

API与counts：CPU物理为R9 paid_density_state的SciPy lu_factor/lu_solve(trans=2)；CUDA物理为相同G的torch.linalg.lu_factor_ex(check_errors=False)、lu_solve、lu_solve(adjoint=True)，同材料正向/伴随共享LU。辅助current均仍为原packed/vector FFT。physical_LU是一材料系统分解一次，primal_RHS/physical_source_RHS是同一计数的别名，adjoint_RHS/physical_adjoint_RHS也是别名，不能相加。Green/FFT/coefficient计数是传播动作，不是额外独立LU或physical RHS；每RHS维度5184，不把all60数组311040当LU维度。fine18生成的每RHS维度17496另列。

CUDA dense公开m0单点实测：LU2.008776wall、共享LU全60正向.113644、全60伴随.078230，warm完整value/gradient+原残差2.727353wall，原primal/adj残差9.61e-16/1.00e-15，tensor峰reserved1.319GB。当前dense生产允许2.0GB tensor-reserved（包含活跃Fourier/current与allocator余量）；历史FFT1.5GB政策不是CPU内存限值。完整新v08 CIE实际194.09375CPU/200.067416wall，physical166.323302wall，80accepted/81requests，峰reserved1.797GB；同Y/预算旧FFT CIE882.213975wall，因此实际完整任务4.41×工程改善。此为backend成本对照，不等于A+成像收益或材料成功，GT仍未读。

最终执行状态：18/18 primary已封存，所有数值worker退出；source seal及local byte验证位于PRIMARY_REMOTE_SEAL_VERIFICATION.json和PRIMARY_LOCAL_BYTE_VERIFICATION.json。新v04三臂整进程1292.304098 CPU/1326.884030 wall；新v08三臂604.031250 CPU/626.103118 wall。最终费用只按ACCOUNTING_FINAL.json的additive rows累计，旧慢队列、被中断尾部和必要探针单独保留。

公开backend端点核对（0新forward/未读GT）：旧FFT v08 CIE Phi=.000904347649、全局Frobenius数据残差.0300510562、PG∞8.416e-7；新dense Phi=.000903497170、残差.0300371488、PG∞4.949e-7。两者native69/171/63、物理80accepted/81requests，整任务882.213975→200.067416 wall，4.40958×。有限容差导致优化路径略异，不能声称相同精确端点或该速度等于材料收益。完整公开证据见V08_BACKEND_ENDPOINT_COMPARISON.json。

原物理目标精确归一化：`Phi3=(1/60) sum_s ||S J_s−Y_s||_2² / ||Y_s||_2²`，无1/2、R0。它继承3D逐源归一化，不是`||SJ−Y||_F²/||Y||_F²`的全Y归一化，也不是iid同σ的白化高斯ML；各源信号能量不同会形成不同权重。三臂及v04 posthoc continuation完全保持同Phi3，不能将其raw值跨2D/3D合并。公开CN99 DPCS另用真实σ的全噪声球事件，不能把Phi3直接称该球的规范化统计。V08_BACKEND_ENDPOINT_COMPARISON另报两个不同残差：sqrt(Phi3)是逐源相对残差RMS，global_Frobenius_relative_data_residual是整体数据残差，二者明确区分。
