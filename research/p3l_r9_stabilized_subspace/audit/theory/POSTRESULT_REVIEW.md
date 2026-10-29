# R9-S 已运行结果后的独立审查

2026-10-29；读取当前 SUMMARY.json、全部B/C在.02的CORRECTION.json、probability_integrate.py及已存概率结果。新增Maxwell RHS=0。这里的E3修订是开发后提出的新变体，不能回写E2或算作新独立确认。

## 1. E2现在实际说明什么

.02材料RMS下：A普通Full-J 8/8竞争，B固定weak 7/8，C动态weak 7/8，Frozen-J 8/8；pair≤1分别7/6/6/6。B/C各7个对象在第一trial之后就因 NONPOSITIVE_PREDICTED_REDUCTION停止。B/C的首步全forward并没有验证新的corrected材料，因为预测门已拒绝。两者多数退回同一个初始材料；因此当前结果没有识别“冻结子空间转动失效”机制。

Frozen-J是必要强普通对照，其.02 operational RHS192、wall约.6182秒；B为360、1.0965秒，C为368、.9940秒。此为已暴露8对象与本窗口描述性成本；完整research账和诊断RHS继续另列，不从成功子集挑速度。未校正U在.01已8/8，进一步削弱该半径上的复杂方法独立主张。

源码对B/C叠加了 Qwᵀδ=0、材料signed-progress一阶约束及材料半径切向约束，然后用仅在Qs中的decoder Newton再强制精确半径。这与A/F的全位移radial retraction不同，A-v-B不是单纯换basis的因果消融。

## 2. 不需要转动也会出现的限制机制

最简单的线性例子已经足够。令decoder材料映射M=I，初始位移d=αe₁，weak方向e₁固定，要求weak幅度α不变；强变量s⊥e₁。若目标材料半径ρ=|α|，则

‖αe₁+s‖²=ρ² ⇒ ‖s‖²=0。

任何非零强修正都不可行。可是半径的一阶切向条件 dᵀδ=0 对所有强δ自动成立，故切向GN可提出非零步，再被有限半径重投影取消或扭曲。这个失效与weak子空间转角完全无关。非线性、非等距材料metric中，Qs上的半径梯度也可能很小，Newton分母变小，投影修正被J放大后使预测下降转负。

该例证明一种合理机制，不单凭它断言全部R9实例已被因果解释。现有日志只保存投影后pred、投影改变量及最终stop；建议从已付J/材料导数重算投影前pred作零物理核查，若正→负则直接定位retraction交互。没有必要为此重跑Maxwell。

## 3. 固定weak方向、自由幅度的合法修订

以固定Q=[Qs,Qw]，选v₀∈range(Qw)，考虑线性子空间

x−x₀=αv₀+Qs s，α可变。

它固定weak方向而不固定weak幅度。等价地，只将weak子空间内与Qwᵀv₀正交的另外三维置零；总允许子空间维数63，再加材料半径条件。此时线性例子可取 α(s)=sqrt(ρ²−‖s‖²)，局部有非平凡强修正。

对非线性decoder，以 h(α,s)=‖G(x₀+αv₀+Qs s)−G(x₀)‖−ρ。若在当前合法点 ∂αh=nᵀDG v₀≠0，则IFT给局部唯一α(s)。这是局部可行充分条件，不能保证整个半径、跨seam或大s都有解。原quintic C²足以做此一阶IFT；不需要新增高阶Maxwell结构。

**普通radial retraction也是合法实现，但不是α-only solve的同一个映射。** 对可行候选d=αv₀+Qs s，用x′=x₀+t d匹配真实材料半径，则

Qwᵀ(x′−x₀)=tα Qwᵀv₀。

因此弱方向及其余三weak零坐标精确保留，幅度和强坐标共同乘t。全forward/合法性、signed progress、latent半径和信任步仍须重新验收。不能把它报告为“weak坐标值全部固定”。沿ray的decoder距离不保证全局单调/唯一；只有原实际bracket/距离门成功才返回。

为最小修复，可以保持原GN tangent约束，但允许这一radial retraction改变weak幅度；也可直接在63维允许子空间内算相同GN步。二者不是同一算法，须预先明确是哪一种。后者较符合“只约束另外三weak坐标”的字面定义；前者只是最小retraction消融。二者都不是新posterior prior，prior目标仍在原x上以‖x‖²/2评价。

## 4. 动态C不能省掉的定义

更新Qw,t后，当前d_t=x_t−x₀在该weak frame中的向量 a_t=Qw,tᵀd_t 一般四维都非零，即使原frame只有一个非零坐标。Procrustes对齐只选basis表示，不把a_t自动变成(α,0,0,0)。

一个明确、便宜、公平的动态定义是：在本步固定当前weak投影方向a_t/‖a_t‖，允许其幅度变；强修正δ∈range(Qs,t)，然后全位移radial retraction。于是a_{t,new}=t a_t。下步再刷新frame，重新以真实当前投影定义方向。若将三维约束写成B_tᵀQw,tᵀδ=0，其中B_t列张成a_t的正交补，当前点就在允许子空间中；不需额外投影。

若方法一定要让weak投影沿重新选择的首向量而非a_t，则必须显式transport/projection当前总位移到该方向，支付产生的真实新trial及其验收。不能只旋转basis就声称另外三坐标仍为零。a_t≈0时方向不适定，应保留退化状态/事前规则，不能任意选有利方向。

以上仅建立定义和局部可行性，完全没有证明修订会赢；不能事后把E2的失败改称成功。四DEV先行再扩到八的门、费用和动作应由主线程冻结。

## 5. probability_integrate.py逐项审查

**正确项。** 非零host transfer中的O保留；完整支持在[0,RHI]×[0,IHI]被矩形分区包住；各未保留cell只要整个cell满足Φ−Φ_mode≥50，就用整个排除union prior≤1得到scaled evidence≤exp(−min gap)，没有对每个cell重复使用单位prior质量。triangle majorant连续两次逆作用给输出导数的代数正确；host原子明确留在source积分中，尾界重复包住它只会保守。当前两个retained矩形严格处于合法bounding rectangle内部，覆盖几何没有倒序。

**已发现且worker确认修正。** 原reference qualifier的eventerr只取finest source-v-GH，漏掉source和GH各自的末两阶事件/层质量变化。PURE_OVERLAP真实最大source层变化约2.551132942e−5，原字段仅9.02119e−8；仍小于reference事件门.0004，不改变pass，但必须用v2/最终ledger记录正确数值。源代码与原输出留存，不回写掩盖。

**证据强度。** 原`.999*min(finestZ)`不是可证明的下界，只是数值denominator proxy。source粗→细Z的相对变化在REGULAR约.58%、overlap约1.3%，因此不能将.1%系数描述为经过该序列认证。最终source-v-GH已经很近，但需同时报告证据与任务的细化差异。worker计划v3采用包含全部规则的更保守proxy，并明确不是保证下界；适合数值敏感性审查，仍不产生严格upper mass。

**无需新物理的改进。** 保留整个cover的cell/endpoints/majorant诊断并记录重放算术成本；加初始cut排序/包含性断言，避免未来mode窗口越界产生负面积cell；各事件和源层分别检查两序列refinement；把strict bounds继续设null。数值遗漏质量字段应叫proxy/estimate而不是严格证书。对现有硬事件，retained矩形远离阈值；不能据此宣称一般硬事件积分成本已解决。

当前完成cover节点数REGULAR748、overlap1582；minimum excluded Φ gap约50.0484/50.0227。记录的遗漏posterior proxy约2.18e−18/8.22e−18，只说明在该浮点模型及数值分母下的外域贡献非常小；Schur/系数/舍入无外包，strict certificate始终未建立。

## 6. 本次审查的决策建议

允许把weak方向自由幅度修订作为新的限定开发实验；先保留旧变体负结果与明确几何机制，再与Frozen-J、原Full-J和uncorrected强对照比较。若修订仍无收益，停止把子空间冻结当作独立算法贡献。

概率结果可有NUMERICAL_REFERENCE_QUALIFIED与物理复用PASS，同时保留CERTIFIED_CONDITIONAL未达成。数值参考修正后即使PASS，也不推出R9探索候选具有posterior质量。Full66、连续PDE及全域exact材料fiber继续未认证。
