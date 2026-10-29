# Statistical analysis: canonical frozen Fresh10 family

This document uses the root statistical plan frozen before fresh generation and the independent QA bootstrap. No bootstrap seed, confidence interval, or comparison family was reselected. Canonical plan: `config/candidate_fresh10_statistics.json`. Canonical bootstrap: `audit/independent_qa/INDEPENDENT_PAIRED_BOOTSTRAP.json`; its byte-identical convenience copy is `results/candidates_fresh10/CANONICAL_PAIRED_BOOTSTRAP.json`. The plan copy and copy receipts are in that directory.

The independent unit is the parent, N=10, not the480 proposals. Each arm uses K8. Differences are candidate minus reference: positive coverage and negative MSE favor the candidate. The20,000 paired parent resamples use SeedSequence([20261030,101010,999]); one shared resampling index matrix serves all four primary comparisons. Percentiles0.00625 and0.99375 give98.75% two-sided intervals, the frozen Bonferroni-four descriptive family. N10 remains an exploratory mechanism check.

Coverage means a legal, deduplicated, truly noise-qualified material cloud member lies within full576 material RMS0.075 of truth. The shared anchor is included when qualified. Empty qualified clouds count as coverage=false; actions explicitly fall back to the anchor. Old arms additionally preserve their radius/progress gate. All10 parents participate. In this run all10 anchors qualified, with no action fallback or right censoring.

| Candidate minus reference | Metric | Mean parent difference |98.75% paired interval |Wins/losses/ties |
|---|---|---:|---|---|
|BILATERAL_TARGET_RML_BLOCK - FULL66_RML|coverage_eps075|-0.1|[-0.4, 0]|0/1/9|
|BILATERAL_TARGET_RML_BLOCK - FULL66_RML|uniform_action_MSE|3.682555405e-05|[-0.0002343646566, 0.0003616069675]|5/5/0|
|FULL66_RML - FIXED_Q8_RAY|coverage_eps075|0|[0, 0]|0/0/10|
|FULL66_RML - FIXED_Q8_RAY|uniform_action_MSE|-0.0003854600132|[-0.001482911377, 0.0007900552987]|6/4/0|

Bilateral has no observed primary coverage advantage over Full66:2/10 versus3/10, one parent loss and mean difference-0.1. Its uniform-action MSE interval crosses zero. Full66 and FixedQ8 have identical coverage on these10 sampled parents; the[0,0] bootstrap interval describes these all-zero paired outcomes and does not prove population equivalence. The Full66 uniform MSE mean is lower than FixedQ8, but its interval also crosses zero.

| Arm, K8 |Coverage.075|Qualified fraction|Offline oracle MSE|Uniform action MSE|exp(-Phi) heuristic action MSE|
|---|---:|---:|---:|---:|---:|
|FULL66_RML|3/10|1.0000|0.008144385509|0.008352861149|0.008401172875|
|SHAPE_TARGET_RML_BLOCK|3/10|1.0000|0.008221914432|0.008426370497|0.008400523964|
|BILATERAL_TARGET_RML_BLOCK|2/10|1.0000|0.008323248787|0.008389686703|0.008400329071|
|FIXED_Q8_RAY|3/10|0.6125|0.007245361091|0.008738321162|0.008496688068|
|DYNAMIC_Q8_RAY|3/10|0.6125|0.007262005488|0.008591982134|0.008477812226|
|RANDOM66_FULL|3/10|0.2500|0.008130255492|0.008388748303|0.00840674572|

Oracle uses truth offline and is not an online selector. Uniform and exp(-Phi) means are actual finite-cloud actions, without prior/proposal/Jacobian corrections; neither is a posterior mean. Remaining metrics and K1/2/4/8 cost prefixes are descriptive, with no added post-hoc significance tests.

EqualK across six arms is not equal cost, constraints, or correction freedom. Only the three new arms share full66 freedom,12F8J, initialization and common eta/xi draws. Old Fixed/Dynamic/Random retain6F4J, radius/progress constraints, latentcap4 and their own correction geometry. Cross-family differences cannot be attributed solely to target distribution or called matched-cost causal advantages.

Independent QA recomputation: `audit/independent_qa/FRESH10_SCORE_RECOMPUTATION.json`, PASS,0newRHS. QA guarded GT access on the complete480-row prediction seal and GT_SCORING_OPENED. Persisted producer `audit/independent_qa/recompute_bootstrap.py` uses only the independent parent-score artifact and frozen plan, with no physics or GT IO. Its default replay file preserves canonical JSON; QA reports replay equality apart from UTC. This documentation task copied canonical JSON without recomputing it.

Statistics-only replay command (no physics):

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/p3l_r10_bilateral_structure/audit/independent_qa/recompute_bootstrap.py --output /tmp/r10-bootstrap-replay.json
```

DEV8, heldout-shape response kernels and the DEV836-node conditional replay are separate scopes, not additional independent material parents. See REPRODUCE.md for staged replay and provenance.
