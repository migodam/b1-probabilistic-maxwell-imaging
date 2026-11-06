# Post-result direct-physics descriptive cost contrasts

These ten contrasts were added after seeing the primary results. They are descriptive post-result analyses, not new algorithms, confirmation tests, or preregistered significance tests. No training, solver, physics, interpolation, GT selection, or raw truth access was performed. Existing scored records alone were read.

All 42 primary objects are retained: ID 24, held-out family 12, cross-grid 6. Cohorts and actual hardware remain separate. Positive differences mean physical L-BFGS / NN is larger; lower error, Phi and cost are favorable, while higher IoU is favorable.

CPU is measured same-hardware deployment CPU, including frozen preprocessing and inference; solver-only CPU is additionally reported. Training and cold loading are excluded from this deployment comparison. LU and complex RHS are separate axes. Final CPU is full paid run cost, including failed or ineffective work.

For each object and axis the common allowance is min(final actual cost A, final actual cost B). Each endpoint is the minimum original-Phi paid checkpoint at or below that allowance, with latest-index tie break. No GT chooses a checkpoint. Actual checkpoint costs can differ, so unused allowance, next checkpoint costs, bracketing gaps, and incumbent ages are reported. A stopped branch receives no invented work.

Intervals are 95% percentile intervals from 10,000 independent-object paired bootstrap resamples within cohort/hardware; they are descriptive, unadjusted for multiple comparisons, and may be unstable with small hardware strata. Missing pairs remain explicit.

A quality gain against a uniform direct method that failed or stopped quickly cannot establish acceleration: the NN physical solver may spend more CPU to attain its better final result. Use both final costs and common-allowance endpoints. No cost-to-GT-target claim is inferred.

Endpoint status counts: `{"PAIRED_PAID_INCUMBENTS": 1680}`.

| Cohort / hardware | Reference | n | final CPU A / B (s) | final NRMSE delta [95%] | common CPU NRMSE delta [95%] |
|---|---|---:|---:|---:|---:|
| CROSS_GRID / Mac | a_plus / nn | 6 | 30.823 / 66.837 | -0.24501 [-0.53389, 6.0546e-05] | -0.24197 [-0.53559, 0.0036466] |
| CROSS_GRID / Mac | a_plus / uniform | 6 | 30.823 / 66.516 | -0.24401 [-0.532, 0.00043314] | -0.24102 [-0.53207, 0.0044916] |
| CROSS_GRID / Mac | a_plus_sourcewise / nn | 6 | 30.823 / 67.007 | -0.24379 [-0.53083, 0.00010542] | -0.24182 [-0.53377, 0.0037228] |
| CROSS_GRID / Mac | a_plus_sourcewise / uniform | 6 | 30.823 / 67.424 | -0.24399 [-0.54398, 0.00056206] | -0.24144 [-0.53492, 0.0040682] |
| CROSS_GRID / Mac | cie_dual / nn | 6 | 30.823 / 61.57 | -0.24404 [-0.53146, 0.0001744] | -0.24302 [-0.53574, 0.0042682] |
| CROSS_GRID / Mac | cie_dual / uniform | 6 | 30.823 / 62.647 | -0.24389 [-0.53174, -1.2654e-05] | -0.24312 [-0.53487, 0.0042525] |
| CROSS_GRID / Mac | cie_single / nn | 6 | 30.823 / 62.04 | -0.244 [-0.53169, -5.7667e-05] | -0.24333 [-0.53413, 0.002877] |
| CROSS_GRID / Mac | cie_single / uniform | 6 | 30.823 / 60.038 | -0.24314 [-0.531, 0.0013393] | -0.24319 [-0.53443, 0.0036888] |
| CROSS_GRID / Mac | global_dual / nn | 6 | 30.823 / 62.783 | -0.24402 [-0.53154, 0.00047437] | -0.24345 [-0.53806, 0.0048913] |
| CROSS_GRID / Mac | global_dual / uniform | 6 | 30.823 / 61.404 | -0.24442 [-0.53189, -0.00043921] | -0.2436 [-0.53673, 0.0043507] |
| HELDOUT_FAMILY / Mac | a_plus / nn | 6 | 39.123 / 74.633 | -0.66341 [-1.1629, -0.2734] | -0.65055 [-1.1459, -0.27247] |
| HELDOUT_FAMILY / Mac | a_plus / uniform | 6 | 39.123 / 80.772 | -0.74939 [-1.194, -0.43385] | -0.71943 [-1.1402, -0.43075] |
| HELDOUT_FAMILY / Mac | a_plus_sourcewise / nn | 6 | 39.123 / 69.807 | -0.65865 [-1.148, -0.25693] | -0.64797 [-1.1276, -0.26075] |
| HELDOUT_FAMILY / Mac | a_plus_sourcewise / uniform | 6 | 39.123 / 82.574 | -0.73702 [-1.183, -0.411] | -0.7119 [-1.1417, -0.41707] |
| HELDOUT_FAMILY / Mac | cie_dual / nn | 6 | 39.123 / 63.743 | -0.61105 [-1.1044, -0.1562] | -0.60508 [-1.1246, -0.15723] |
| HELDOUT_FAMILY / Mac | cie_dual / uniform | 6 | 39.123 / 69.893 | -0.75332 [-1.1894, -0.43589] | -0.74929 [-1.1823, -0.43475] |
| HELDOUT_FAMILY / Mac | cie_single / nn | 6 | 39.123 / 72.574 | -0.68989 [-1.1715, -0.29096] | -0.69331 [-1.1852, -0.30278] |
| HELDOUT_FAMILY / Mac | cie_single / uniform | 6 | 39.123 / 74.696 | -0.69761 [-1.1864, -0.30835] | -0.69362 [-1.1623, -0.29858] |
| HELDOUT_FAMILY / Mac | global_dual / nn | 6 | 39.123 / 70.065 | -0.685 [-1.1667, -0.28702] | -0.67056 [-1.1383, -0.29146] |
| HELDOUT_FAMILY / Mac | global_dual / uniform | 6 | 39.123 / 74.015 | -0.74726 [-1.1828, -0.44624] | -0.74791 [-1.1656, -0.44929] |
| HELDOUT_FAMILY / XINAN | a_plus / nn | 6 | 128.52 / 177.93 | -0.19956 [-0.49735, 0.00058554] | -0.19888 [-0.49662, 0.0029072] |
| HELDOUT_FAMILY / XINAN | a_plus / uniform | 6 | 128.52 / 183.67 | -0.20068 [-0.50219, 0.0017106] | -0.19854 [-0.49347, 0.0012674] |
| HELDOUT_FAMILY / XINAN | a_plus_sourcewise / nn | 6 | 128.52 / 182.25 | -0.19973 [-0.49799, 0.00082647] | -0.19829 [-0.4944, 0.00351] |
| HELDOUT_FAMILY / XINAN | a_plus_sourcewise / uniform | 6 | 128.52 / 181.46 | -0.2 [-0.49943, 0.0011247] | -0.19552 [-0.48987, 0.0047431] |
| HELDOUT_FAMILY / XINAN | cie_dual / nn | 6 | 128.52 / 157.71 | -0.19964 [-0.49743, 0.00045387] | -0.19908 [-0.49318, -0.00023885] |
| HELDOUT_FAMILY / XINAN | cie_dual / uniform | 6 | 128.52 / 159.74 | -0.1983 [-0.49641, 0.0027549] | -0.19741 [-0.49503, 0.0041335] |
| HELDOUT_FAMILY / XINAN | cie_single / nn | 6 | 128.52 / 157.75 | -0.19906 [-0.49743, 0.0018872] | -0.20164 [-0.49459, -0.0018847] |
| HELDOUT_FAMILY / XINAN | cie_single / uniform | 6 | 128.52 / 158.09 | -0.19776 [-0.4962, 0.0039178] | -0.19903 [-0.4938, 0.0023211] |
| HELDOUT_FAMILY / XINAN | global_dual / nn | 6 | 128.52 / 161.53 | -0.19973 [-0.49736, 0.00033983] | -0.19925 [-0.49426, -0.00022678] |
| HELDOUT_FAMILY / XINAN | global_dual / uniform | 6 | 128.52 / 138.28 | -0.1786 [-0.47657, 0.040775] | -0.18163 [-0.49814, 0.036547] |
| IN_DISTRIBUTION / Mac | a_plus / nn | 10 | 32.034 / 63.901 | -0.4269 [-0.79167, -0.12548] | -0.35719 [-0.72602, -0.072953] |
| IN_DISTRIBUTION / Mac | a_plus / uniform | 10 | 32.034 / 67.341 | -0.53278 [-0.90226, -0.21851] | -0.52195 [-0.88142, -0.20716] |
| IN_DISTRIBUTION / Mac | a_plus_sourcewise / nn | 10 | 32.034 / 67.396 | -0.53552 [-0.89524, -0.22303] | -0.5229 [-0.88977, -0.2102] |
| IN_DISTRIBUTION / Mac | a_plus_sourcewise / uniform | 10 | 32.034 / 64.958 | -0.53156 [-0.90132, -0.22516] | -0.51921 [-0.88129, -0.20351] |
| IN_DISTRIBUTION / Mac | cie_dual / nn | 10 | 32.034 / 61.853 | -0.45861 [-0.82025, -0.14072] | -0.45341 [-0.82202, -0.13642] |
| IN_DISTRIBUTION / Mac | cie_dual / uniform | 10 | 32.034 / 62.119 | -0.55313 [-0.92115, -0.23767] | -0.55194 [-0.92446, -0.23608] |
| IN_DISTRIBUTION / Mac | cie_single / nn | 10 | 32.034 / 62.628 | -0.40062 [-0.78419, -0.080519] | -0.39433 [-0.78122, -0.074334] |
| IN_DISTRIBUTION / Mac | cie_single / uniform | 10 | 32.034 / 62.114 | -0.49234 [-0.87539, -0.16037] | -0.48991 [-0.87344, -0.1648] |
| IN_DISTRIBUTION / Mac | global_dual / nn | 10 | 32.034 / 60.908 | -0.45867 [-0.83276, -0.14103] | -0.45397 [-0.83131, -0.1363] |
| IN_DISTRIBUTION / Mac | global_dual / uniform | 10 | 32.034 / 60.296 | -0.57102 [-0.93909, -0.24516] | -0.56852 [-0.9422, -0.24618] |
| IN_DISTRIBUTION / XINAN | a_plus / nn | 14 | 69.765 / 145.23 | -0.037269 [-0.13114, 0.049197] | -0.015435 [-0.07802, 0.044326] |
| IN_DISTRIBUTION / XINAN | a_plus / uniform | 14 | 69.765 / 168.81 | -0.037979 [-0.13311, 0.048812] | -0.044978 [-0.13244, 0.031992] |
| IN_DISTRIBUTION / XINAN | a_plus_sourcewise / nn | 14 | 69.765 / 156.44 | -0.0037994 [-0.073391, 0.065061] | -0.016864 [-0.079001, 0.041837] |
| IN_DISTRIBUTION / XINAN | a_plus_sourcewise / uniform | 14 | 69.765 / 159.38 | -0.044093 [-0.15388, 0.04859] | -0.049445 [-0.14728, 0.031721] |
| IN_DISTRIBUTION / XINAN | cie_dual / nn | 14 | 69.765 / 153.52 | -0.069709 [-0.17255, 0.016545] | -0.065851 [-0.16702, 0.016752] |
| IN_DISTRIBUTION / XINAN | cie_dual / uniform | 14 | 69.765 / 144.84 | -0.071121 [-0.17601, 0.01435] | -0.065702 [-0.16656, 0.016536] |
| IN_DISTRIBUTION / XINAN | cie_single / nn | 14 | 69.765 / 155.37 | -0.073387 [-0.17585, 0.015969] | -0.068467 [-0.17694, 0.015696] |
| IN_DISTRIBUTION / XINAN | cie_single / uniform | 14 | 69.765 / 151.92 | -0.073129 [-0.17974, 0.015498] | -0.067145 [-0.17107, 0.016615] |
| IN_DISTRIBUTION / XINAN | global_dual / nn | 14 | 69.765 / 143.85 | -0.0027453 [-0.071076, 0.065695] | 0.0028113 [-0.060692, 0.064009] |
| IN_DISTRIBUTION / XINAN | global_dual / uniform | 14 | 69.765 / 149.52 | -0.038558 [-0.13915, 0.049293] | -0.033618 [-0.12937, 0.048214] |

Machine-readable tables retain original Phi, NRMSE, IoU, final residual and final costs where available; the copied endpoint logic reports Phi, NRMSE and IoU. `final_object_pairs.csv` includes individual cost ratios. No cross-hardware absolute CPU pooling is used.

Artifacts: `final_object_pairs.csv`, `final_paired_effects.csv`, `object_matched_at_most_budget_endpoints.csv`, `object_matched_budget_paired_bootstrap.csv`, `ANALYSIS.json`, `COST.json`.
