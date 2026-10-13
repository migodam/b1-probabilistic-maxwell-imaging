# B1 研究日志

## 2026-10-13 — 输入审计与首轮假说

用户指定方向在本地统一称 B1；A32 仅作为来源标签。ZIP完整性检查通过，10文件安全提取至 inputs/source_package/，原件保持不变。完整阅读 README_ZH、NEXT_DECISION、THEORY_REPORT、检查脚本及检查/资源JSON，并核对原A31脚本。历史74次二维解析调用、A31的2个slab及1个微型dipole状态分别保留，不把它们当新增三维成像。

已有 A32 缓存连续MAP导致复杂材料MSE恶化39.46%，不能直接沿用为概率成像成功。A30 Sparse/PCA容量是offline truth-fit，不是posterior支持保证。A17冻结论文不修改。

假说：低维、明确的三分支proper prior可用可信积分输出多张材料代表及分支概率。先核验3D向量Mie的电/磁多极、切向E/H连续、角向/能量/极限，再通过3/5/7阶积分比较概率及log Bayes factor。不同模板不是连续材料单点的非零概率。

初始计划：物理审计约20次，三分支连续积分最多249次（有效缓存减少重复），数据/控制约12次，48次density-corrected importance；预留未分配预算至384。所有新2D/3D实际调用共同计入384，不使用附录80上限重置历史成本。只读已公开代码，不打开任何旧locked/truth数组。默认CPU。

概率可靠性优先：未通过reference门槛的case报告REFERENCE_UNRESOLVED，不进行算法优越性/连续校准宣称。Top2展示保留第三分支mass；模板之外质量为未建模未知，绝不伪造百分比。

文献：ScholarQA三条query的6个operation均HTTP429，空候选不是文献不存在。保留原始错误；原包已识别论文的官方全文/元数据可直接核对，检索完整性与全局原创性保持未认证。

## 第一轮实际结果与归因

19项3D forward审计、16项代数/概率测试通过。原生3/5/7阶连续积分共243次新状态（中心重复有效复用）。总356/384，零GPU/NN/旧数据读取。竞争案例P=(0.08932,0.50642,0.40427)，5/7阶ΔP=2.72e-6、ΔlogBF=1.10e-5；信息不足案例近先验。24个新连续材料parent×30噪声，全部5/7参考门槛通过；720次噪声不是720个独立物体。

失败：强信号ΔP极小但ΔlogBF=.366；gain/表示错配ΔP=.031/.046且ΔlogBF=1.71/3.10，参考不足。不能放宽原门槛。48次密度修正importance最大概率误差=.02565，未达到.02；高ESS不证明概率准确。暂停扩张随机采样。

归因：本阶段主要的reference失败是likelihood积分分辨率，forward本身通过有限审计；绝不是证实Bayesian方法或Maxwell不可行。gain错配与表示错配的PPC约.00156/.00176，提示模型不合适，即使概率集中也不能当正确成像。full-versus-calibrated-L在竞争例的ΔP=.01232；A17签名在2D参数上rank2，没有压缩收益，不能声称native feedback加速。

## 第二轮设计 — 冻结剩余28次预算

先用已支付5/7阶网格构造tensor-Legendre forward插值；再在同一个continuous prior上以17/25阶廉价积分分辨sharp likelihood。此时target是明确的近似forward，不能偷称全精度认证。新付费验证：每拓扑4个box corner+4个独立均匀点，共24；另4次使用已存在A15/treams多层矢量Mie参考（独立radial kernel，14多极），避免重复开发已有物理工具。所有28次先reservation，失败同样计费，不允许越过384。

比较：5/7阶forward插值之间、17/25阶积分之间、24个新full-forward点、独立treams。原生full-grid reference失败仍保留；若替代诊断通过，只标EMPIRICALLY_REFINED_REFERENCE，不标严格uniform remainder/概率认证。新点是自适应development核验，不冒称locked/blind validation。


## 第二轮结果与最终交付

总计384/384次新材料–频率物理评估，预算已封闭。24个新插值验证点的最大相对场误差为5.9995e-7；4次独立treams矢量Mie核对的最大差异为1.8521e-15。17/25阶积分与5/7阶forward插值对照支持EMPIRICALLY_REFINED_REFERENCE，原始四项REFERENCE_UNRESOLVED保留，未升级为统一误差或严格概率认证。

竞争案例的三事件质量为：中央空腔50.6418%，中央与环形空腔40.4265%，实心8.9317%。Top2保留第三事件质量。24个新连续材料parent各30次噪声共720个development观察，不能视为720个独立物体或真实仪器概率校准。材料区间仍有积分/CDF分辨率与漏模态限制。高概率但低predictive consistency的gain/representation负控制说明模型错配不能靠概率集中掩盖。

同一暴露panel的冻结importance控制在等75次名义cold forward预算下明显逊于5阶积分，停止当前48样本方案的精度宣称。A17/A30/A32的小参数签名提供局部解释与反例核查，未建立native A17压缩或加速价值。结论为GO_LIMITED_TOPK_PROTOTYPE；一般几何、严格概率认证与强TAP原创性尚未建立，不启动大NN或扩张campaign。

用户要求停用DeepSeek后，已停止相关进程，后续概率复核、8张图与文献ledger均由原生Codex subagents完成。复核修正了H边缘材料coverage标签以及gain nuisance元数据丢失；历史错误文件和成本保留。最终24项测试通过，19项真实3D forward审计通过。局部HTML浏览器实测因工具拒绝file URL而NOT_RUN，未绕过限制；静态图已视觉检查。4次独立oracle原始数组缺失如实记载，不越预算重跑。

中文报告、理论、结果、裁决、可执行CLI、配置、概率数组、Top-K图、交互HTML及完整资源账本已交付。结束文献检索与实验；未访问旧locked数据，未自动公开发布，未修改A17冻结稿。最终清单只记录文件路径与大小，不执行SHA核验。

## 2026-10-13 用户授权公开发布

用户要求直接publish，并明确不新增gate、不核验SHA-256、不增加smoke或零碎代码测试。按此要求发布既有代码、科学报告、材料图、概率结果与成本/失败记录；不重跑实验。公开状态与已有科学限制分别记录。

## B1 E0 — 连续posterior摘要修复（独立本地分支）

用户授权新增物理0。主线在 research/b1-e0-posterior-audit 上执行；未发现额外variable-radius阶段，未借用C1预算。对六旧development案例、三条件事件与材料边际实施25/41/65/97同一density积分及独立二维分区自适应参考。全部指定连续摘要经验条件通过；四项历史full-grid失败仍保留。

关键发现：原生node moments/bin-CDF组合违反必要方差下界，72历史条件坐标中5项不相容。原生与细化q05在修复density中的实际CDF约0.006096/0.030012，不能视为5%。修复后strong H0 Re SD=.027064934288、q05=2.467483321080。只修复summary，不把该问题称evidence归一化失败。

新CLI显式携带gain_sd，实际gain观测重放概率与E0一致，旧缺失metadata不能静默改变likelihood。提供216个直接CDF probes、连续数组、六TopK图和两诊断图。保留slogdet警告、独立Cholesky检查、旧treams数组缺口及三次原生报告worker容量失败。初版图文字重叠已修，旧图/绘图成本保留。

新的同一polynomial target是经验参考，无uniform/full-Mie/现实校准认证；原先720数据材料coverage不自动沿用。E1 full-Mie secant/epsilon+gain补偿、多频、统一implicit field、active sensing、NN与新发布均未执行。科研裁决：继续修复概率benchmark；机制解释力未测。不放宽门槛、不追加物理、不重置历史账本。完整证据见 research/b1_next/e0/。

## 2026-10-13 — 用户要求默认公开、复用仓库并整理主页

E0发布并入当前B1 main，最新入口改为同一项目网页。新增研究总览、发布维护规则，历史A系列19仓库归档但不改名/删除，A17论文保持独立。主页只突出B1/C1/A17，其他项目与历史折叠。此工作没有新物理、训练、SHA或额外测试。
