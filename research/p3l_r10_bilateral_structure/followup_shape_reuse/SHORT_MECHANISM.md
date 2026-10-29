# Support activation and failure of the frozen weighted-state basis

**中文结论：** 24个路径形状中20个（7/8个parent）激活旧占比严格为零的格点，最多27格；新支持上的输入能量占比最高3.079%。旧基的源投影相对缺陷为0.652%—24.288%。但另外4个没有新支持的形状也有32/32节点失败，所以“支持激活”不是全部失败的唯一原因。投影缺陷与每形状平均实际白化响应误差的描述性相关为0.506，新支持输入能量与误差相关仅0.147。

该条件障碍针对加权状态u的anchor基。总场形式的注入E_R固定，因此不能扩大为“物理要求每步更新材料弱谱”，更不能否定所有TSOM/固定状态表示。原共享66维材料约束下的固定总场/电流TSOM对照仍为 `NOT_RUN`。

Zero-Maxwell-RHS analysis of saved arrays only; no truth or forward solver read/import. Exact-zero occupancy is defined by a0 == 0; activation by anew > 0. ROI indexing is read from the historical setup.

## Observed mechanism
- 20/24 path shapes (7/8 parents) activate previously zero occupancy; newly active cell counts span [0, 27].
- New-support source-energy fractions span [0.0, 0.030793448762061474]; relative source-projection defects span [0.006519866278933605, 0.24287980538404727].
- The largest saved anchor-basis row norm on exactly zero anchor occupancy is 8.233e-13. The source-energy lower-bound check holds for all 24 shapes: True.
- 4/24 shapes have NO new support, yet all 32 corresponding reuse nodes fail the gate. Their source defects span [0.006519866278933605, 0.04430860085321456]; shape-mean actual white errors span [0.061931530203421854, 1.7989323941491169].
- Descriptive Pearson correlation across 24 shape means: source defect versus actual response error = 0.505926; new-support source-energy fraction versus error = 0.147124. The 24 shapes are grouped within 8 parents; these are not independent samples or causal tests.

## Conditional algebraic obstruction
With B0=P0 C P0 and R0=P0 E_R, an anchor cell with a0_i=0 has row i identically zero in R0 and every B0^k R0. Hence every exact vector in the anchor Krylov space is zero on those cells. If newly active cells carry nonzero Rnew, that injection cannot belong to the frozen anchor space. For S={i:a0_i=0,anew_i>0}, the exact-space lower bound is ||(I-V0 V0*)Rnew||F/||Rnew||F >= ||Rnew[S,:]||F/||Rnew||F. Saved finite-precision rows are separately reported in the CSV. This proves an injection-containment obstruction, not a lower bound on port error or proof that support activation is the unique failure cause. No-new-support failures show it cannot explain all failures.

## Coordinate and scientific boundary
The failure is specifically for frozen V0 in the weighted state u, where h=P u and u=delta P e. It does not refute TSOM or every fixed state representation. The total internal field obeys e=(I-delta C A)^(-1)E_R and F-Fh=delta O A e, whose injection E_R is fixed. Thus the new-injection obstruction is partly representation/preconditioning dependent; it is not proof that physics requires a fresh material weak eigenspace at every step. Operator/port changes and nonnormal amplification can matter even without support activation. A fixed total-field or current TSOM construction retaining the original shared 66-dimensional material constraint is a minimal future control: NOT_RUN.

Reuse failed 192/192 nodes in the saved experiment. Rebuild passed 184/192; its remaining eight nodes had an unavailable rank-64 model because the saved basis rank was 61. Their placeholder approximate-output discrepancy is not the error of an executed rank-64 reduced solve. No material-coverage improvement, general impossibility, or unique causal explanation follows.

Files: SUPPORT_ACTIVATION.csv (24 rows), SUPPORT_ACTIVATION_SUMMARY.json, SUPPORT_ANALYSIS_RECEIPT.json.
