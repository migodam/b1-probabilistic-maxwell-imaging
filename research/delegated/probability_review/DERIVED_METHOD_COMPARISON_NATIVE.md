# B1 derived frozen-method comparison

This is saved-array reanalysis of the exact existing calibration panel, with zero new physics. There are 24 continuous-material parents, stratified 8 per topology, and 30 previously recorded noise realizations each. Replayed parent parameters agree exactly with the manifest; all Gauss7 probability records agree within 1.45e-15.

The 48 importance points and densities remain frozen from one competitive y0 proposal. Its cold setup is 27 order3 states plus 48 samples=75; Gauss5 is 3*25=75. Gauss7 is 3*49=147. These are nominal standalone setup counts, not new budget charges. Actual original Gauss5 pipeline new calls were 72 because 3 centers were already cached. Reference and solver-audit costs are additional shared study costs.

| Method | Nominal cold states | Mean / maximum event-mass error vs Gauss7 | Both P/logBF diagnostics passed | Brier / log score |
|---|---:|---|---:|---|
| gauss3 | 27 | 0.00314929 / 0.0211448 | 712/720 | 0.518486 / 0.852173 |
| gauss5 | 75 | 3.43196e-05 / 0.000290358 | 720/720 | 0.518211 / 0.851616 |
| gauss7 | 147 | 0 / 0 | 720/720 | 0.518204 / 0.851607 |
| frozen_importance48 | 75 | 0.0387663 / 0.232232 | 124/720 | 0.530183 / 0.871099 |
| grid7_min_residual_softmax_HEURISTIC | 147 | 0.0733014 / 0.200956 | not a posterior estimator | 0.530219 / 0.878695 |
| grid7_hard_MAP_ONE_HOT_DECISION | 147 | 0.460127 / 0.816755 | not a posterior estimator | 0.883333 / positive infinity |

Both-diagnostic success requires max event-probability difference <=.02 and max log Bayes-factor difference <=.1 against the exposed Gauss7 quadrature reference. This is numerical agreement, not a uniform certificate. Gauss3 fails 8 of 720 under the combined diagnostic despite failing only 1 under the probability-only diagnostic. Frozen IS passes only 124 of 720 under the combined diagnostic.

The grid7 min-residual softmax is a normalized peak-fit heuristic. Hard MAP is encoded as an explicit one-hot decision only for diagnostic scoring: 318 wrong one-hot decisions make its exact mean log score positive infinity. JSON uses null plus infinity flags for these entries and does not floor probabilities.

Frozen importance48 minus Gauss5 paired changes at equal nominal cold cost are Brier .0119718 [95% .0033808,.0208755] and log score .0194827 [.0039134,.0349092]. These intervals resample 24 material clusters within true-topology strata (2,000 resamples, 8 parents per stratum). They apply only to this frozen proposal, prior/noise design, and exposed panel. No general superiority claim or newly adapted per-datum importance experiment is implied.

Two instrumented comparison passes used 1.162910s wall and 0.889905s process CPU; both kept the physics ledger 384/384 unchanged.

Data and receipts: [research/delegated/probability_review/DERIVED_METHOD_COMPARISON_NATIVE.json](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging/research/delegated/probability_review/DERIVED_METHOD_COMPARISON_NATIVE.json>), [research/delegated/probability_review/DERIVED_METHOD_COMPARISON_RECEIPT_NATIVE.json](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging/research/delegated/probability_review/DERIVED_METHOD_COMPARISON_RECEIPT_NATIVE.json>), [research/delegated/probability_review/DERIVED_METHOD_COMPARISON_RUNS_NATIVE.json](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging/research/delegated/probability_review/DERIVED_METHOD_COMPARISON_RUNS_NATIVE.json>), [research/delegated/probability_review/DERIVED_CALIBRATION_OBSERVATIONS_NATIVE.npz](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging/research/delegated/probability_review/DERIVED_CALIBRATION_OBSERVATIONS_NATIVE.npz>), [research/delegated/probability_review/derived_method_comparison_native.py](</Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging/research/delegated/probability_review/derived_method_comparison_native.py>).
