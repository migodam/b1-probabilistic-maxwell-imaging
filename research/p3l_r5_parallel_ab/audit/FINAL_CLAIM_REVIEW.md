# Final independent claim review

Status: SCIENTIFIC_SCOPE_CONSISTENT; TWO_PRESENTATION_CORRECTIONS_REQUIRED. No new physics, training or numerical experiment was run. Reviewed root RESEARCH_REPORT_ZH.md, GO_NO_GO.md, START_HERE.md, A final research report, B REPORT_ZH.md/B_SUMMARY_ZH.md, and retained compact ledgers/replay/bootstrap receipts. Raw packaging was not reviewed.

## Required corrections

1. Add the timing regime to overview cost claims: “CPU单线程、batch64、每parent”. Root START_HERE.md's2.54%/54.57% line, GO_NO_GO.md's2.54% table cell, and root report's opening point3 currently omit the batch denominator. These are batch64 ratios; single-parent provider/full ratios are9.82%/47.65%, respectively. Both provider regimes pass1/3, but their numbers cannot be interchanged. Root§4.3 already gives the correct batch64 context.
2. Replace “24个线性RRR拟合” in root START_HERE.md and root report§6 with “24个最终选定的RRR模型（12个原控制、12个嵌套控制；另含rank/ridge网格和OOF重复拟合）”. Twenty-four counts final arm×seed models, not every covariance/SVD regression fitting call inside ridge/fold search. Recorded elapsed costs include the search; this is a count-label correction, not a missing cost or a new arm. The42 neural fits are correctly6 matched heads+9 compact providers+9 compact heads+18 nested baselines.

Suggested concise headline clarification: START_HERE.md “统计压缩通过开发门” may say “开发集点估计通过压缩门”. The full report already states that the98.333% retention interval has lower endpoint below80%, so no new scientific decision is requested.

## Verified numbers and boundaries

- A:32,705 state reservations/attempts;261,632 actual forward RHS =32,704×8. These RHS counts are reservations before attempted solver actions, not a successful-endpoint count. Exactly one reserved state attempt was interrupted before its forward RHS; the8RHS difference is explicitly disclosed and not refunded. Derivative/adjoint30,152; total291,784; two failures; no unfinished actions. Values agree across reports and ledger. No unused-budget transfer or new budget is claimed.
-18 observations are6 independent parents×3 paired noise levels across3 declared2D submodels, with other64 coordinates fixed. Reports do not call them18 independent subjects or full66 posterior marginals.
- A charts are numerical moment/normalization comparisons, not certified global TV. Gaussian TV/KL are explicitly against qualified numerical references, and invalid old-reference metrics are excluded. RTO is local with unknown global success/branch normalization. Spline/manifold lifting/Full66/SBC remainNOT_RUN.
- S1 lower and upper prior strata each have mass0.2208604823. Final lower-host posterior estimates0.926520/0.941277/0.953144 come from direct retained-node readback; stage A explicitly discloses remaining convergence failure and pure-stratum Float64 underflow. The old0.70–0.78 values are not promoted. Root report correctly labels these estimates and withholds a branch certificate.
- A5 is fixed-unit-W, fixed-weak-center line-segment evidence. Root§3.4 explicitly says this, and does not turn legacy L_on_rho_ball into an ambient2D/66-ball claim. All six R3 sufficient conditions fail; Float64 norms are not directed-rounding enclosures.
- B raw retention92.9591% is mean parent loss over3seeds, not ensemble-image MSE. Its98.333% exploratory parent-bootstrap interval[71.60%,135.63%] is conditional on reused VAL selection and has lower endpoint below80%. Reports preserve this and do not call it fresh confirmation. Clipped91.72% has separate numerator/denominator.
- Canonical batch64 additional-provider ratio2.5446% and full provider/features/head ratio54.5672% are separate. The1/3 primary gate applies to additional provider. Old receipts remain preserved; main timing table uses canonical raw30-repeat values.
- CPU/MPS original3e-6 replay isFAIL:3 violations, maximum3.945e-6, mean-risk difference≤9.31e-10. Root and stage reports retainFAIL without widening tolerance. Separate CPU-exact double-load replayPASS is described as a different test.
- B's final outcome is conditional development progress with no new256+256 same-anchor confirmation. Budget-basedRESOURCE_LIMITED/NOT_RUN and original96+64 anchor semantics remain explicit. Maxwell-specific linear-feature increment is not established; the negative result is not generalized to all nonlinear contexts.
- Historical budget/negative-result edits are not claimed. Sources and assets are explicitly reused; no history rewrite or SHA256 verification is asserted. No contradiction was found between these textual claims and the audit/provenance records reviewed. File-stat/Git provenance does not independently certify every historical payload byte, and the reports do not claim that it does.

No scientific blocker was found beyond the already declaredHOLD/NOT_RUN/FAIL gates. Apply the two count/denominator wording corrections before publishing the overview, then retain the current conditional conclusions.
