# R5 measurement 复核报告（deepseek-pro reviewer，2026-11-18）

复核对象：`research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r5_measurement`
复核方式：只读。未运行、未改动任何实验或历史文件；未触碰 `private/`；未改任何线上/冻结决策。
因推理已冻结，允许读 GT（`private/truth.npz`、`private/response_bank.npz`）用于事后核对；GT 只用于验证，
不用于任何重选、重调或重跑。

## 总判定

1. **证据链闭合，无结果级错误，无伪造，无预算/盲测破坏，无科学过度声明。** 本包自身用词极其克制，
   所有未证明项均显式标 UNKNOWN，与 `EXPERIMENT_PROTOCOL.md` 的 status vocabulary 一致。
2. **区间证书的“序列化盒”处理正确**：qualification 用“向内一 ULP 的安全框 + 向外扩宽的 K + 严格包含”
   重建了一条可复核的证明链（详见 C 节），我独立复算 16 个 CERTIFIED 根全部通过，菜单响应包围 26/26。
3. **局部 vs 全局覆盖区分正确**：12/12 case-arm 的全局覆盖均为 UNKNOWN，没有任何一处把局部证书写成全局结论。
4. **被拒匹配方向已显式分析**：菜单对超指向（feed norm>4）的拒绝有落盘审计；但 T0 的拒绝量级偏弱（见 E 节）。
5. **Gate2 不支持** → 按协议应为 NOT_RUN、NO-GO for expansion。理由见 G 节。

## A. 正确性与实现核对（通过）

- 公式/实现自查全部通过（只读复看既有收据）：
  - `results/debug.json`：非零宿主消元相对误差 3.70e-16；解析 Jacobian 对中心差分 5.43e-10；
    dual 行 vs 完整正演最大误差 1.29e-15；直接补基恢复材料误差 8.9e-14；噪声移动 w0 检查 4.26e-17。
  - `results/formula_tests.json`：多项式残差 <1.6e-16、poly-Jacobian 差分误差 <1.1e-11；
    power Gram Q 与 2048 角积分相差 2.44e-15；node（G=0,S 分量=0）与 pole 反例处理正确。
- `validate()` 用独立 full numpy solve（响应相对 ≤1e-8、状态残差 ≤1e-10、σmin ≥1e-6、box 边界）过滤根；
  Q/TRF64/reciprocal64 三方法的 accepted roots 与 `candidate_registry.csv` 的 GT_RMS 与 `private/truth.npz` 一致。
- direct 补基公式 `m=w/(R+w^H H)` 推导正确（`THEORY_DERIVATIONS.md` §5；debug 直接误差 8.9e-14）。
- TRF 的 reciprocal 链式法则 `-1/x²`、64 Sobol+中心起点、max_nfev120、tol 1e-12 与协议一致
  （`src/experiment.py` solve_trf）。

## B. GT 顺序、盲测与公平性（通过，含 2 个 LOW 事项）

- GT 与靶定义一致：T0 距历史锚 RMS 0.0743（max 0.112，符合 ±0.15）；T1∈[0.6,1.5]；
  T2/T3/T4∈[0.5,30]；T5∈[0.5,0.7]。truth.npz 键序与 scorer 的 case 循环一致，无错位。
- clean 臂的点预测过滤在全部 6 个目标命中真根（T0→root0、T2→root1、T3→root1 等，
  `results/source_design_ledger.csv` hard_retained），证明预测模型正确。
- `score.py`→`score_v2.py` 仅一行差异（qualification 的 interval 字典不再被当作复数解码），
  `inputs/scoring_repair_receipt.json` 完整记录：GT 未在修复前暴露、原文件保留、sha 已存。
  修复发生在评分码冻结之后，但属于协议允许并记录的“engineering repair”，不影响结果。LOW。
- 公平性：同一源向量在不同策略得到相同噪声（broker 用 case+arm+a 的哈希），跨策略可比；
  direct completion 用 n-s=3 次单位功率测量并公开列出（`score_table.csv` direct_new_measurements=3）；
  TRF64/reciprocal64 用各自发现集与各自设计做归属对照；same_active_TRF64 用同源 64 起点做 multistart 对照。
  无事后挑源或挑靶。
- T0 是锚定靶（历史 GENERAL6_00 根 + iid ±0.15），协议已声明非完全盲（`EXPERIMENT_PROTOCOL.md` 末段）。INFO。

## C. 区间证书与“序列化盒”处理（重点复核：结论为正确）

证书本身（`src/certification.py`）：
- 输入 G/E/S/y 作为精确 binary64 常数；mpmath.iv 定向区间 64→128→256 bit；
  Krawczyk `K=x-CF(x)+(I-CJ(X))(X-x)`；要求 K ⊂ int X 严格（`_strictly_inside`，certification.py:406-424）
  且行和收缩 ||I-CJ(X)||∞<1（certification.py:599-604）；区间 Gauss–Jordan 每个 pivot 严格排除 0。
  该组合是 Krawczyk 存在唯一性定理的保守充分条件，实现正确。

序列化审计（`src/qualify_predictions.py:19-20`）：
- 证书里的 X 端点用 `float()` 最近舍入落盘（不是向外扩宽）；K 经 `_report_interval` 向外各扩一 ULP。
- qualification 对每个序列化 X 端点再 `nextafter` 向内各收一 ULP 得 `safe` 盒，然后检查
  `safe.lo < K.lo 且 K.hi < safe.hi` 全坐标严格成立。
- 方向性正确：K_true ⊂ K_serialized（向外）⊂ safe（向内）⊂ X_true，故原收缩界仍适用，
  且唯一根 ∈ K_true ⊂ safe。这使 `certified_safe_binary64_box` 是可由落盘数据精确重建的合法域。
- 我用同一 venv 只读复算：16 个 CERTIFIED 根（含 clean/noisy）serialized strict inclusion 全部 True；
  T0 两臂两根基的 26/26 菜单响应区间包围实际点响应；safe 盒包含根点。
  safe 盒不交 ⇒ 两根不同 ⇒ CERTIFIED_AMBIGUOUS 成立（根唯一性保证根 ∈ safe 盒）。

结论：**证书有效**，其有效域严格限于“离散代数映射的局部实根小盒”，排除项
（无全局覆盖、无复数根、无 Hankel/几何/PDE/离散化误差证书、无扰动输入证书）
在 `certification.py` SCOPE/EXCLUSIONS 与 qualification 的 `interval_prediction_scope` 中如实写明。

LOW（仅序列化类，可修）：
- 原始 `Q.json` 中 `certificates[].certificate.intervals.X` 的 float 端点不是保证包围，
  只能以 `qualification.json` 的 `certified_safe_binary64_box` 作为正式域。
  建议（纯序列化修复，不重跑）：certify() 落盘 X 时也走 `_report_interval` 式向外加宽，
  或直接持久化向内 safe box，避免他人误用原始 X。此修复不影响任何已得结论。

## D. 局部 vs 全局覆盖（区分正确，一处 LOW 措辞）

- 每个证书只证明“小盒内唯一实根”；qualification 与 summary 全量返回
  `global_material_box_coverage=UNKNOWN`、`noise_feasible_set_coverage=UNKNOWN`。
  CERTIFIED_AMBIGUOUS 仅表示“≥2 个不交（safe 盒层面）认证局部根”，status_scope 明写 no exhaustiveness。
  没有任何 residual-based completeness 或噪声可行集完备性声明。与协议 §噪声/§状态词汇一致。
- LOW 措辞：`material_boxes_disjoint`（qualify_predictions.py:40）是在**向内 safe 盒**上判的不交。
  这足以证明“两根不同”（根 ∈ safe 盒），但**不等价于**完整认证盒 X_i、X_j 不相交
  （子集不交推不出超集不交）。若将来需要“盒级分离”，须改用向外盒的交点测试。对当前
  CERTIFIED_AMBIGUOUS 声明无影响，仅防误读。

## E. 被拒的匹配方向与源菜单（显式分析，一处量级注记）

- 菜单规则冻结：2 old + 12 primitives + 12 random + ≤18 pair matched + 1 fisher，≤48；
  pair 方向 `a=Q†d*`（Cauchy–Schwarz/广义 Rayleigh 最优），归一化功率=1，硬限制 ||a||₂≤4，
  超限即拒绝并记录省略（`src/experiment.py` Physics.menu；`EXPERIMENT_PROTOCOL.md` 源与功率段）。
- 审计重建全部被拒 pair 方向并落盘 `omitted_or_attempted_pair_directions`：
  T0 clean/noisy feed=6.4（仅 1.6× 超限，**拒绝量级偏弱**）；T3 clean feed=109.5；
  T4 noisy 三对 feed=452–721（明确超指向）。即“不可行的超指向匹配方向”这一叙事在唯一完全认证的
  二根情形（T0）只是边界性成立，强烈超指向只在 T3/T4 出现。这是科学量级注记，不是代码缺陷；
  阈值 4 是冻结设计规则，本轮拒绝合规，且协议已声明 max-min 只是有限菜单枚举、非约束连续全局最优。
- T0 选定源 random_06：点间隙 0.004772、区间分离下界 0.0040682（clean）/0.0056466（noisy-central）、
  chosen_noise_model_margin 0.00383/0.00541。我复算点间隙与下界方向一致（下界 ≤ 点间隙，正确保守）。

## F. 噪声比较与模型失配（公平，两处限制）

- 噪声：固定 |noise|=1e-4/观测；旧数据相位来自私有种子，新测量相位来自 sha256(case+arm+a) 的确定性哈希。
  后者与协议“deterministic random phase from private seed”措辞有轻微出入，但确定性、与 GT 无关、
  同源跨策略相同（公平），不构成盲测风险。LOW。
- 模型失配（噪声臂半径 .04004 vs 重建 .04）实际超过声明 allowance 2e-5 于 **4/6** 目标
  （T0 1.09e-4、T2 3.9e-5、T3 3.6e-5、T4 5.8e-5），`score_table.csv` 的 old_model_bound_ok=False 如实呈现。
  硬过滤与 tube/margin 用 eps+allowance=1.2e-4：该值对“噪声+失配”只在部分源方向勉强覆盖，
  因此过滤仅保证保守拒绝（可假阴性），`chosen_noise_model_margin` 不是已证噪声分离；
  qualification 的 scope 已声明“无不确定数据管、无物理核误差证书”。这是科学限制，非缺陷。
- 结果层面：noisy 臂 maxmin 与 fixed 几乎重合（T1/T2/T3/T4/T5 noisy 均退化到 primitive_00，
  仅 T0 选 random_06），设计优势在噪声臂基本消失；fisher 仅在 T1(0.0228)/T2(1.29)/T5(0.1729) 局部更好，
  且后两者仍 >0.1 失败。阈值 0.1 下 noisy 成功仅 T0/T1/T4 三例。

## G. Gate2 判定：不支持（NOT_RUN / NO-GO for expansion）

协议 Gate2 条件（`EXPERIMENT_PROTOCOL.md` 预算与停止段）：“strict coverage 且
physically useful noise-margin 测量优势在 direct completion 与普通 multistart 归属下仍存活”。

- 条件一（strict coverage）：12/12 case-arm 全局根覆盖与噪声可行集覆盖均 UNKNOWN，不满足。
- 条件二（测量优势存活）：
  - maxmin vs fixed：noisy 臂无优势，T0 反而略差（0.0724 vs 0.0680），其余同源相同；
  - fisher 的不一致局部优势不足以成为“physically useful”且 T2/T5 仍失败；
  - direct completion（3 次测量）系统性压过单次测量设计：T2 noisy 0.288（仍 >0.1）、
    T3 noisy 0.0869（唯一成功）、T4 noisy 0.0579；
  - multistart 归属对照（TRF64_active/reciprocal64_active、same_active_TRF64）无优势残留。
- 结论：**Gate2 NOT supported**。建议状态 NOT_RUN、NO-GO for expansion。该结论与本包
  `THEORY_DERIVATIONS.md` §7/§8（噪声可行集覆盖未证、失配超 δ 于 4/6、测量优势未成立）一致，
  也与“未成立命题不得改写为成功”的自述一致。

## H. 严重度清单

无 HIGH（无结果级错误、无过度声明、无预算/盲测破坏、无证书失效）。

- MEDIUM（科学限制，非缺陷，不建议改动靶/调参）：全局与噪声可行集覆盖未证；
  模型失配 allowance 在 4/6 目标被真实失配超过；噪声臂设计测量优势消失（Gate2 的直接原因）。
- LOW：score_v2 冻结后一行修复（已记录、GT 未暴露）；certify 原始 X 序列化非向外
  （有 safe box 兜底，建议向外加宽落盘）；material_boxes_disjoint 为 safe 盒级（当前声明仍成立）；
  新测量噪声相位用哈希而非私有种子（公平、无 GT 依赖）；T0 匹配方向拒绝偏边界；
  Q 方法无显式 30s 墙钟（实测 ≤0.035s，受 RLIMIT 1100s 兜底）。
- 仅记账/序列化类的具体修复（可选，不重跑、不换靶）：
  1) certification.py 落盘 `intervals.X` 时改为向外加宽（或直接持久化 safe box）；
  2) 若未来需要盒级不相交声明，改用向外盒交点测试；
  3) 将 score_v2 的 sha 并入 freeze 清单（sha 已存于 scoring_repair_receipt）。

## I. 主要证据索引

- `EXPERIMENT_PROTOCOL.md`：Gate2 条件与 NOT_RUN/NO-GO（:27）、噪声与 allowance（:11）、
  源菜单与 feed 限制（:15）、状态词汇（:23）、预算（:26）。
- `THEORY_DERIVATIONS.md`：序列化审计原理（:72）、根管与噪声可行集（:74、:82）、未成立命题（§8）。
- `src/certification.py`：_report_interval 向外扩宽（:380-386）、_strictly_inside（:406-424）、
  Krawczyk 判定（:599-604）、intervals 落盘（:670-678）。
- `src/qualify_predictions.py`：safe box 与严格包含（:19-20）、不交/分离下界/margin（:33-40）、
  状态与 UNKNOWN scope（:46-47）。
- `src/score_v2.py`：GT 只读评分与 model bound 检查（:31-40）。
- `results/qualification.json`、`results/score_table.csv`、`results/source_design_ledger.csv`、
  `results/candidate_registry.csv`、`results/executor_receipt.json`、`inputs/scoring_repair_receipt.json`、
  `results/qualification_failure_01.json`。
- 复核实测：16/16 证书序列化严格包含；T0 两臂 26/26 响应包围；预算 CPU 43.76s/7200s、
  RSS ≤178MB/8GiB；证书墙钟 ≤0.43s/target（上限 600s）；discovery 各方法墙钟 ≤1.3s（上限 30s/臂）。
