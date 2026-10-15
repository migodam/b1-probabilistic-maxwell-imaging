# 原材料先验约束下的 Maxwell 竞争解释：本轮负结果

本轮已完成 E0 和受限 E1；**未建立原 A25/A30 先验内的竞争材料见证，停止 E2，预算封账。** 这不证明唯一性，也不证明更充分的原先验搜索不可能找到竞争解释。

1. **是否找到了竞争材料解释？** 没有。缓存中 8,064 个同 acquisition、同 amplitude 的 TRAIN parent 配对，1,632 对通过事先冻结的结构条件；其最小预测距离为 44.826717731720585。新物理搜索中最小为 49.223622812870374，远大于本轮优先寻找的 2。对四组各自最接近的结构对，两个端点轮流作为 truth、各 32 次固定独立噪声，256 条重放中同时满足 noise-ball 和主似然条件的数量为 0。
2. **旧 T/G 候选败在哪里？** 首先是原材料先验不一致。绝对空气核不是原背景上的非负固定相位增量；历史 corrector 还允许全域实虚部分别调整。735 条有真实 full prediction 的历史记录，没有一个同时通过原先验的必要 PHASE 条件。G 的 18 条 noise-ball 通过记录相对其历史 base 的最小 Δell 仍为 13.370668，高于 log(10)，甚至高于 log(100)。这项比较只是历史诊断，因为该 base 本身也没有原先验资格。
3. **N3 是否仍值得保留？** 是，保留强 Top-1 基线。既有固定 8 条件的材料 complex MAE 为 0.0417555943；历史最小改动的相容修正版为 0.0507565790，绝对增加约 0.009001。可是后者没有 P0_EXACT 资格。本轮没有把“较小 MAE 损失地恢复 noise-ball”升级成“在原先验中恢复相容”，也没有运行新的 N3 先验修正实验。
4. **新方法收益属于哪里？** 本轮的有效收益是识别并纠正先验与似然解释错误；尚无生成策略收益。共同 C_P、四臂比较、Best-of-4 改善均 NOT_RUN。只实现并运行了原生成器内的有限参数 secant / paid exact-J 排序见证搜索。
5. **下一步选择？** 不推进本 target 上的多解生成器竞赛或更大网络。保留 N3，下一阶段若获授权，优先研究先验保持的最小必要修正及其模型失配；将本轮 target 作为辨识性对照，而不是预设它必须多解的 benchmark。没有理由为得到正结果放宽噪声、材料相位或结构条件。

## 问题为什么改变

生成不同图像很容易；原问题需要的是原材料分布允许、且对同一测量有接近数据支持的不同解释。历史 T 已证明某些 box 合法端点能进入固定噪声球，但这个集合远宽于 A25/A30 的原生成分布。两个端点各自通过 95% noise-ball，也不意味着二者点似然接近，更不意味着它们各占显著 posterior mass。

本轮分别检查 BOX、PHASE、P0_EXACT、noise-ball、点似然与 posterior 区域质量。最后一项没有积分或采样证据，所有 posterior_mass 均为 null。PHASE 通过只是必要条件；生成器回放与参数范围同时成立，才标 P0_EXACT。没有拟合回生成器的端点不会被武断标为不属于 P0；但明确违反必要相位、非负或峰值规则的端点确实与原 P0 冲突。

## E0：真实 prior 和零物理审计

从实际生成器与配置恢复：

\[
\chi=\chi_b+a\,p/\max(p)\,v,\quad
v=(1+0.2i)/\sqrt{1.04},\quad a_1=1.2.
\]

原 p 来自 shell、asymmetric_piecewise、smooth、mixed_scale 四类 TRAIN 混合；没有新增 hole-count 桶。weak 背景为 0.1+0.04i；ring 背景在原固定环区为 2+0.1i、其他位置仍为 weak。rho=0 表示背景，不是空气。原配置的 a1 是索引，不是幅度 1。

64 个独立 TRAIN parents 的全部 512 条条件材料由原生成器逐元重放，最大误差 0。新增 256 个新 seed 的 material-only parents，未新增 paired Maxwell 训练数据。另有 256 个 TRAIN 图像的半格插值诊断变换，只校准 descriptor 稳健尺度，不当作 P0_EXACT 候选。

审计共 2,835 条历史材料记录；735 条有精确 full prediction，2,100 条 pool 记录只有背景线性 L，因此后者的 full ell 保持空值。记录含 trial / final / best_valid 别名，不是 2,835 个独立材料世界。同观测 30,068 对中，230 对满足 |Δell|≤log(10)，但没有一对的两端同时满足原先验的必要 PHASE 条件。所有 phase 合格 pool 的线性预测不能被偷换成真实 Maxwell 似然。

旧 32 条件属于 exposed development。历史 truth 的 4/32 个噪声实现落在 95% noise-ball 外，保留为预期可能发生的随机事件，没有改变 sigma 或 tau，也没有据此宣布材料物理非法。未打开旧 locked truth / geometry / clean response。

## 同一个物理合同和数据支持

原 24×24 scalar TM/VIE，8 源、32 接收，complex128，所有源共享同一材料，gain=1；复 RMS 噪声为 sigma，实/虚通道各以 sqrt(2)/sigma 白化，共 512 实通道：

\[
\ell(x)=\tfrac12\|W(F(x)-y)\|^2,\quad
\tau=\sqrt{\chi^2_{512}(.95)}=23.785449004057032.
\]

每个 pair 固定同一个 sigma。缓存主筛查用较小 TRAIN parent ID 端的原 sigma，并用于双方；另存 median-development-sigma 的 4,032 对替代诊断，绝不混成主表。新 E1 用原公式从事先固定的 anchor differential response RMS 生成 sigma，随后不再随候选改变。

若 δ=W(F_b−F_a)，d=||δ||，且观测来自 a，在本模型下：

\[
\ell_b-\ell_a=\tfrac12 d^2-\xi^T\delta,
\qquad \xi\sim N(0,I),
\]

因此其均值为 d²/2、标准差为 d。噪声球相容与似然竞争是不同条件。这里的 2 是本轮预声明的见证搜索尺度，不是 Maxwell 物理常数；完整距离和 log(100) 敏感性均保存。没有将最小 residual 或 profile minimum 当 evidence。

## E1：实际执行了什么

物理资产可复原，因而在 E0 完成后进入 E1。四个新独立 anchor parents 覆盖四个原 family，背景组合依次为 shell/ring、asymmetric/weak、smooth/weak、mixed/ring；amplitude 固定 1.2。它们是离线 benchmark construction，不是免费在线发现。E1 的材料 reservoir 和所有见证参数从未来 E2 online view 隔离。

每个 anchor 运行两轮原生成器有限参数搜索，每轮生成 128 个合法 parameter secants。候选保持原离散 rasterizer、归一化、参数范围及增量相位。第一轮付费 exact material Jacobian 对 local 与独立 global reservoir 提案排序；实际 full 求值后再刷新 Jacobian 搜第二轮。全部端点由原 A30 full solver 验证，未用任意全材料补偿、可变 gain、truth 排序或新神经网络。

结构条件来自 TRAIN nearest-neighbor 尺度与半格 descriptor 扰动尺度的较大者，排除 DC 的 centered 多尺度 RMS≥0.08796801446、occupancy 平均差≥0.02549189815、boundary 平均差≥0.03399884259，并要求至少 3/5 excursion thresholds 的跨尺度 occupancy 差≥0.05，以及至少 6 格的归一化 rho 差≥0.1。它是本轮较保守的结构尺度，不能覆盖一切微小孔洞或细节；未按结果降低。

|固定 anchor family / background|新计算中最近结构合格 P0_EXACT pair 的 d_pred|A 作 truth 的竞争频率|B 作 truth 的竞争频率|
|---|---:|---:|---:|
|shell / ring|49.223622812870374|0/32|0/32|
|asymmetric_piecewise / weak|310.3676883669929|0/32|0/32|
|smooth / weak|257.2295870642063|0/32|0/32|
|mixed_scale / ring|173.2215987936159|0/32|0/32|

表中 family 指 anchor；候选可来自原四类混合中的另一类。32 个预冻结噪声 seed 对每个固定 pair 的两个端点分别轮流重放，y 始终来自真实端点加正常噪声，没有使用两响应中点。log(100) 敏感性也未得到竞争见证。噪声 rows 不能当作新的几何 parents。

### 执行偏差与失败

两处实现故障分别发生在标量参数转换和 metadata 重复字段写入，已修复。失败消耗的 3 states、88 RHS 与数组、原日志全部保留；继续运行使用同一账本，没有重置额度。

预注册文字要求第二轮中心来自首轮 local；实际实现从首轮 local 和 global 的已验证结构合格端点中取最近者。三个 anchor 首轮没有结构合格 local，因此切换到 global。这一偏差在 EXECUTION_DEVIATIONS.json 单列，未改写冻结协议，未补跑。实际结论只适用于这一允许跨 family 重启的有限搜索，不宣称严格预注册的 PC-SP 优越性。完整 PC-SP 广义特征 / C_P 四臂算法未激活。

## 成本、数值范围与停止

E0 新 full state / RHS 均为 0。E1 累计 **45/64 full-state attempts、680/2048 complex full-solver RHS**，包括失败、重复 anchor、10 次 receiver transfer 构造和所有端点。完成的主执行消耗 42 states；39 个结构 pair 比较中有一次 paid-cache 复用，缓存没有被重新计成新增物理。未比较未经缓存披露的 standalone 算法成本。

成功主执行 wall=2.792059416 秒；两次失败缺单独进程 timer，保留动作级计时，并以包含 Git 开销的 command elapsed 给出 1.621224458 秒保守上界。physics-controlled 三次执行 wall 上界约 4.413284 秒；这不是整个科研会话耗时。E0 审计、材料生成、画图、来源记录及发布后处理另列，不能把它们统称免费。材料生成保守合计 1,536/2,048。GPU 训练 0，旧 T/G/A30/P2 额度使用 0。

FULL_VALIDATED_NUMERICAL 仅指原有限网格 full solver 与数值残差核验；没有 uniform forward-output error bound，没有连续 Maxwell 证书。本轮距离远离 2 的边界，无 BORDERLINE_NUMERICAL 正结果；这不让数值校验变成全域物理误差认证。

没有两组不共享 parent 的合格竞争 pair，E2 四臂全部 NOT_RUN，共同 C_P 仅为冻结规格。原 N3 保留不变。未新增训练、未开启 3D、未访问旧 locked 数据。本轮已封账，不继续消耗剩余额度寻找成功。

## 可复核入口与科学边界

- [完整图库](gallery.html)：同观测 N3 / 最小改动 base / 历史 T 对照，以及四个新 pair 的实部、损耗、rho 与差分。
- SOURCE_MANIFEST.json：实际 commit、Git blob、数组路径与字段；没有新增 SHA-256 检查。
- PROTOCOL_FROZEN.json、E1_SEARCH.csv、PAIR_NOISE_REPLAY.csv、WITNESS_PAIRS.json、METRICS.json、COST_LEDGER.jsonl：精确数值与实际执行。
- E0_AUDIT.csv、PRIOR_SUPPORT_AUDIT.md、LIKELIHOOD_AUDIT.md、CORRECTION_ATTRITION.md：旧端点为什么不能直接作为本轮竞争解释。
- FOUR_ARM_RESULTS.csv：四臂明确 NOT_RUN，而不是编造空的比较成绩。

本轮建立的是先验审计和一个可复现的有限负结果。P0 的低维生成约束可能使当前 8 源全孔径观测比自由 1,152 维材料场更可辨识，但当前搜索不足以证明这一点。没有 posterior mass、全局唯一性、生成器增益或强 IEEE TAP 原创算法证据。应先保护合理材料分布与 N3 的成像质量，再决定是否需要多解生成。
