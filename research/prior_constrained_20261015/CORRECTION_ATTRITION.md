# Cached raw → corrected prior attrition

T has 40 frozen raw edits. All 40 fail the original increment phase/range/peak necessary conditions before correction: strict absolute-air cores set chi=0 even though the original background is nonzero, so they are not P0-preserving edits. All saved T corrected endpoint rows remain off-phase. Corrector success in the historical noise ball does not restore prior membership. `CORRECTION_ATTRITION.csv` retains core preservation and outside-core compensation per endpoint. Final/best_valid aliases are explicit, not independent discoveries.

G has 72 raw attempts and two saved correction steps per attempt. All 72 raw and all correction/endpoints fail the same necessary original-prior check. 18/72 final validated endpoints meet the data noise ball; none is P0_EXACT or PHASE-qualified. G edits were global in many cases; outside-raw-support compensation is then zero by definition and must not be interpreted as local correction. DC and low-frequency changes remain separately reported.

No prior-accepted raw trajectory with an actual full endpoint was found in these T/G histories, hence prior-retention success is NOT_ESTABLISHED. Historical T/G experiments tested box-and-noise-ball feasibility, not prior-preserving competitive explanations. PHASE failures here are algebraic exclusions from the declared P0 form; passing PHASE would still not prove original-generator membership.

Structure uses continuous increment rho and a shared multithreshold/multiscale descriptor, with 4-connected foreground, 8-connected complement and >=4-cell component/hole persistence. Absolute near-air is retained only as a historical diagnostic, never used to label original increment topology.

## 分母与归因边界

T corrected104条记录含best_valid/new_valid/final等别名，不是104个独立raw；G的72 fresh validation只有36个精确材料、其中9个noise-valid。原corrector为32实参数Haar GN/SLSQP，最多两步，允许全域DC/低频补偿，不是1152自由参数full-space GN。普通exact Jacobian包含resolvent feedback。新共同C_P尚未激活，不能据此宣称其改善，也不能将历史corrector的收益归于PC-SP。

本轮需报告的先验合格raw到corrected保留率没有可用轨迹，记NOT_ESTABLISHED而非0%或100%。详情见delegated/cached_audit/CORRECTION_ATTRITION.csv与delegated/g_cache/G_CACHE_AUDIT.md。
