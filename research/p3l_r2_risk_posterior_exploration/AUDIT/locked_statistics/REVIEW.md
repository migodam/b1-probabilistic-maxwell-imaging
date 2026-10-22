# Independent locked statistical audit

REVIEW_REQUIRED

Both fresh scoring COMPLETE markers were required before reading scored risks. All256parents and77registered ablations are retained.

## Two fixed primary comparisons

| Comparison | Relative risk reduction | Two-sided97.5% CI |
|---|---:|---|
| fused_vs_probit | 0.2526728 | [0.22336827707953513, 0.2798811522568588] |
| curvature_arc_vs_probit | -0.01956624 | [-0.04325680786138226, 0.002138474682012736] |

## Secondary fixed-three-seed mean risks

Intervals bootstrap256parents, conditional on3fixed seeds. These are mean risks, not prediction ensembles or n768. Secondary multiplicity is unadjusted.

| Comparison | Relative reduction | 95% CI | Improved seeds |
|---|---:|---|---:|
| local_geometry_cnn_vs_local_anchor_cnn | 0.03453091 | [0.026156880107943814, 0.042169821303324306] | 3/3 |
| local_geometry_cnn_vs_local_capacity_cnn | 0.03487741 | [0.026104426950610382, 0.043316602840164545] | 3/3 |
| local_geometry_cnn_vs_local_permuted_geometry_cnn | 0.03469538 | [0.02583243360148816, 0.043207557126405484] | 3/3 |
| weak_strong_decoder_vs_random_decoder | 0.008597694 | [0.0008363551006305709, 0.01684412283483116] | 3/3 |
| weak_strong_decoder_vs_all_decoder_matched | 0.03220728 | [0.018506843480982763, 0.046125684133373236] | 3/3 |
| weak_strong_decoder_vs_u2_decoder | 0.04646716 | [0.032250991070789964, 0.06036963880978014] | 3/3 |

## Metric coherence

- fullgrid_ROI_MSE: True
- fullgrid_ROI_MAE: False
- point_primary_reproduced: True
- same_primary_baseline: True
- trajectory_primary_reproduced: True
- frozen_bootstrap_protocol: True
- two_frozen_primary_comparisons: True

## Density coverage

| Frame | Macro90% coverage | Parent-bootstrap95% CI | Fallbacks |
|---|---:|---|---:|
| density | 0.8902344 | [0.876953125, 0.9031250000000002] | 0 |
| weak_density | 0.8847656 | [0.872265625, 0.8968750000000001] | 0 |

Per-component Wilson95% intervals and individual training-seed ranges/counts are in RESULTS.json.

- All-constant MCMC chains can have MCSE0 and remain wrong/unmixed. MCSE0 is not accuracy, mass or posterior-fidelity evidence; no chain data are read by this audit.
- Informative8 and risk-weak8 densities target different fixed-coordinate marginals. Their raw NLL values are not directly rankable across frames; compare each with its own unit-prior reference.
- These are supervised material-risk and fixed-coordinate density scores. No Bayes-mean identity, full posterior calibration, novelty or GO decision follows.
- Two co-primary comparisons use two-sided97.5% bootstrap intervals (Bonferroni for planned family size2). Secondary intervals do not provide simultaneous coverage.
- The method family, three training seeds, n512 and all pairings are fixed in this audit source; no held-out model/seed/size reselection is performed.
- Uncertainty is conditional on fixed trained checkpoints and1024interval draws; training-seed and finite density-sampling variability are not resampled.
