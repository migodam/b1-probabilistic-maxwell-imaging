# E4 概率与材料任务：结果和后续交接

本页只整理已经完成的E4固定模型。新增物理调用、likelihood计算、benchmark、fresh32统计执行及GT读取全部为0。原始negative、初次报告、修复重放、append-only费用账均保留；不修改正在运行的fresh32冻结源码或config。

## 可保留的结论

原R6 REGULAR与PURE_OVERLAP在相同prior、原复噪声、原观测、完整TM24下，精确幅度Schur传递通过数值实现门；强普通GN适应求积在direct与共享物理两侧都通过旧数值任务门。物理复用带来了实测组件节省。其依据是成熟的消元/Schur恒等式，不是新的概率原理。

全支持任务中心化包络没有显示总体概率算法优势。R6均值门在冻结32²局部cells后仍未闭合；一维原P1条件模型虽然得到数值任务通过，强普通adaptive quad仍更便宜。所有数值区间与reference不确定性代理都没有完整物理系数、算术或prior任务矩外包，故strict certificate仍未建立，Full66仍未合格。

材料pushforward/coarea、条件期望与平方风险分解、Gaussian probit矩、Laplace/GN、RRR、Schur/Woodbury、任务区间的归一化记账均有成熟数学基础。公式集的作用是把假设、原prior、噪声口径与实现位置写清楚，不把标准恒等式改名为原创。任何最终novelty判断由root独立完成。

## 结果分层

|模型/任务|已有结果|资格与边界|
|---|---|---|
|R6物理实现|16个固定幅度核对，白化最大差4.90349e-12|NUMERICAL_PHYSICS；没有外向系数/舍入资格|
|R6正确prior普通GH9/11|direct与Schur均过均值绝对1e-5且.01SD、概率/源层.002旧门|NUMERICAL_TASK_PASS，依赖显式数值reference资格|
|R6reference|GH21/31、原source坐标两级细化、完整支持数值审查|数值交叉与遗漏代理；不是严格Z下界或证明半宽|
|R6 task envelope均值|REGULAR半宽.00808，PURE半宽.01149|旧1e-5门未闭合；不宣布达同误差更快|
|R6 256-ROI硬事件|32²cells下数值半宽通过.002|两例峰远离uniform阈值，不能代表一般hard event|
|原P1一维offset条件模型|两坐标adaptive mean一致；512cells均值数值半宽2.6922e-4|原未截断Gaussian尾部按平台原子保留；不是非均匀P1或Full66|
|完整posterior/连续Maxwell/一般图册|没有建立|不由条件模型成功或候选集合数量解除|

## 数值reference与未闭合的证明项

R6 S3是活动[0,64,65]三维Gaussianized坐标独立截断[-4,4]、其余固定0的材料子模型。V=QD；host=.1+.04i、ROI256。transition材料density保留Q^-2，pure/transition在同材料像上的原prior贡献相加，host单独成原子。REGULAR移除fibre prior会产生.111977SD位移；PURE_OVERLAP源贡献保持.568896345/.431103655。

|案例|均值数值不确定性代理|代理/材料SD|事件/源层数值代理|support审查cells|
|---|---:|---:|---:|---:|
|REGULAR|5.875664e-07|0.0004585259|1.445466e-15|748|
|PURE_OVERLAP|2.116191e-07|0.0001446292|2.551133e-05|1582|

完整V支持外矩形分为保留的mode±12 GN边缘SD矩形和其补域。补域递归cell用整个cell的响应三角majorant推出残差下界，全部cell闭合后，排除并集的总prior质量≤1给缩放evidence上界。host与两连续源层都在该支持账中，没有用有限优化重复收敛充当global审查。

但Schur系数/基列/Q非酉/rounding没有严格外向包围，分母采用有限细化的正数代理而非证明下界；prior task moments仍是数值积分。以上值只能叫NUMERICAL_REFERENCE及numerical proxies。source粗→细evidence约变化.58%/1.31%；最终fine-source对GH31约差1.53e-6/3.78e-6；这些normalizer差异保留，没有被mean通过掩盖。

canonical v3修正了早期event qualification漏计source粗细差的问题；PURE最终proxy为2.55113e-5，仍小于reference要求.0004。早期v1/v2文件没有删除。

## 普通求积与任务包络：公平的负结果

普通direct与Schur的概率规则、GH节点、prior与任务完全相同，故该对照只支持物理复用收益。表中setup+integral只是实测组件和，共同fresh profile、验证、reference、读文件和imports另计；不能把它称完整cold start。单次计时不代表统计速度保证。

|案例/同一GH规则|direct积分wall(s)|Schur积分wall(s)|Schur一次setup+积分(s)|
|---|---:|---:|---:|
|REGULAR/GH9|1.051957|0.107782|0.186043|
|REGULAR/GH11|1.666307|0.245037|0.323299|
|PURE_OVERLAP/GH9|0.489621|0.078790|0.157052|
|PURE_OVERLAP/GH11|0.708039|0.124214|0.202476|

任务包络使用相同Schur evaluator，另花cell majorant、prior矩与根求解费用。R6均值门还未闭合，所以没有为其补跑昂贵direct全概率规则；二因子矩阵该格明确NOT_RUN，不能把普通积分的物理收益据为新包络规则收益。R6 event数值通过与mean未通过分别保存。

在一维条件模型中，原z0与prior CDF坐标的完整普通积分路线各约.4s，数值mean=.2189064273064、real ROI mean variance=8.05967e-7。任务包络首次mean通过在512cells，累计约1.57s；2048cells约6.98s。这里两方同属数值资格，普通求积的error estimate也不是严格证书。观察结果不能支持相对强普通积分的总费用优势，更不能把多花算术换来的数值interval称为strict认证优势。

一维event在普通浮点中下溢为0。首轮报告保留，完整重放以positive upper likelihood floor1e-300与atom上界保留未知量；512cells数值event区间[0,6.45969e-298]。零point estimate不是实际零概率，极小上界也不是已认证的有效位数。该event很稀有，因此其早停不能推广为一般不连续任务总比平滑任务更容易。

## 费用总账与机器可读交付

本分支总计3624新增完整complex RHS、421个576阶LU、1个256阶complex Schur；181,551个八源reduced节点即1,452,408个256阶complex RHS。CPU累计70.337s，作业wall累计72.847s，单进程RSS峰值150.75MiB。部分reduced-only作业短暂重叠，因此wall累计不是端到端elapsed；imports未计。剩余376完整RHS未消费。

- `PROBABILITY_AND_TASK_RESULTS.csv`：43行，含数值reference、匹配direct/shared普通积分、错误prior控制、source细化、任务包络、旧失败复核与总账。reference/component/trajectory/total行相互重叠，不能把所有费用行相加。
- `formula_collection/04_PROBABILITY_COST_TWO_FACTOR.csv`：56行组件费用，区分setup、online、reference/qualification、重放、未运行格及非加和subtotal。
- `results/probability/final_numerical_comparison_v3.json`：R6 canonical；`original1d_results_v2.json`：一维修复后canonical；`task_envelope_results.json`：R6停止轨迹。
- `results/probability/aggregate_costs.json`与`audit/probability/ledger.jsonl`：权威费用；`failed_attempts.jsonl`：报告修复、缓存限制与方法未达门。
- `formula_collection/04_PROBABILITY_AND_MATERIAL_TASKS.md`：35个稳定F-PROB条目，逐项给条件/推导/源码位置。

## 与当前fresh32可能材料集合目标的关系

当前新实验允许approximate ensemble weights，目标是有用的可能材料集合：在有限集合及epsilon-tube的宽度/多样性/完整费用下，评价fresh-parent经验覆盖与可部署动作风险。真实posterior准确率不是这个新目标的硬门。E4严格证书未建立并不禁止这一新目标；反过来，新目标的经验覆盖成功也不能追溯升级E4为posterior质量、Full66全域质量或严格概率证书。两套冻结合同与旧负结果都保留。

主集合只由已封存、原绝对95%noise-ball与decoder/radius/progress合格的实际材料组成；shared anchor也须满足自己的噪声/合法门才进入主集合。best-of-K与nearest GT只用于offline capacity；不能让GT反向选择online集合、权重或模型。finite cloud、tube与逐像素范围的区别，在统计草案中明确。

## 统计CLI：只列命令，尚未执行

必须先得到全部32 parent×3arm×16 direction slot的`PREDICTIONS_SEALED.json`，然后由root检查并明确授权GT评分。以下变量只说明路径，不读取任何私有GT：

```sh
P="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace"
PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
```

全32 seal之后的prediction-only验证：

```sh
"$PY" "$P/src/possibility_statistics.py" validate_predictions \
  --output "$P/results/possibility_volume" \
  --statistics "$P/config/possibility_statistics_frozen.json" \
  --protocol "$P/config/possibility_volume_protocol.json"
```

验证stdout预期：status=`SEALED_PREDICTIONS_VALIDATED_NO_TRUTH_READ`、parents=32、groups=96。它不创建评分输出、不读GT。

仅在root给出seal之后的评分批准文件后：

```sh
"$PY" "$P/src/possibility_statistics.py" score \
  --output "$P/results/possibility_volume" \
  --statistics "$P/config/possibility_statistics_frozen.json" \
  --protocol "$P/config/possibility_volume_protocol.json" \
  --root-approval "$P/config/possibility_volume_offline_score_approval.json"
```

评分批准文件由root创建，不由本文代建。必须包含approved_by_root=true、action=`OFFLINE_TRUTH_SCORING_AFTER_ALL32_SEALED`、statistics_id=`R9_POSSIBILITY_VOLUME_STATISTICS_V1`、protocol_id=`R9_POSSIBILITY_VOLUME_V1`、utc不早于seal.utc。code还核验冻结统计config/方法/parent/顺序/实际资格/文件seal与timing rotation。

私有truth磁盘读取只在`score_sealed_predictions`内部。若`offline_statistics/`已存在，代码拒绝静默覆盖；若评分失败，保留TRUTH_ACCESS_BEGIN和现场，不得以修复名义更改冻结在线模型或事后重选endpoint。

## 评分输出预期结构

目录为`results/possibility_volume/offline_statistics/`：

|文件|报告所需keys/含义|
|---|---|
|TRUTH_ACCESS_BEGIN.json|utc、prediction_seal、statistics_id、protocol_id、root_approval、role；只在实际允许读GT时生成|
|SUMMARY.json|status、parent_count、primary_four_contrasts、secondary_yield_noninferiority、coverage、deployable_approximate_actions、secondary_paired_action_differences、material_task_coverage、cloud_width_summary、descriptive_width_coverage_bins、probable_GT_empirical_claim、anchor_unqualified_parent_count、primary_empty_cloud_counts、claim_boundaries|
|PARENT_PREFIX_SCORES.json|每parent/arm/K0/1/4/8/16的qualified unique/yield、best RMS/MSE、epsilon hits、tube geometry、ROI/threshold/soft tasks、uniform/likelihood action损失与fallback标记|
|PARENT_RADIUS_SCORES.json|四个预声明radius的secondary材料/任务结果；不把radius当独立replicate|
|RAW_NOISE_QUALIFIED_DIAGNOSTIC.json|全部noise-qualified visited states的variableK诊断；不替代固定16方向主集合|
|FULL_SYSTEM_COSTS.json|每parent/arm的全部尝试、anchor/frame/runtime/overhead分配、historical setup、失败数、first-returned-qualified/right-censor、timing窗口|
|WIDTH_COVERAGE_COST_POINTS.json|parent×arm的pair diameter、qualified_unique、epsilon.075 hit、nearest MSE、standalone wall、动作损失，供非空洞tradeoff图|
|COST_ACCOUNTING.json|authoritative_campaign_resource_snapshot、generation_receipt、parent_end_to_end_wall_seconds、cost_scope；区分生成费用与部署视图|
|SCORING_COMPLETE.json|status=COMPLETE、statistics_id、new_Maxwell_RHS=0、new_forward_calls=0、no_post_truth_model_changes=true|

四主contrast的精确key：`coverage_fixed_minus_dynamic`、`coverage_fixed_minus_random`、`offline_nearest_MSE_fixed_minus_dynamic`、`all_attempt_standalone_wall_fixed_minus_dynamic`。每项有status、parents、difference、confidence=.9875、percentile_interval、bootstrap_resamples=20000；若空cloud导致nearest error未定义或费用缺失，则该项为NOT_ESTIMABLE，不能删parent重新计算。

coverage行提供successes/parents/rate、Wilson95、Clopper_Pearson95、exact_one_sided95_lower。bootstrap单位始终是parent；directions/seeds/sources不增加n。uniform与likelihood-weighted动作分别保存全32 parent_losses与bootstrap95，不使用GT选择哪种权重。`probable_GT_empirical_claim`只对应预声明epsilon/delta下fresh-parent总体命中频率，不是每个测量的posterior置信度。

本文完成时没有调用validate或score，也没有查看当前fresh32 GT或新评分结果。
