# Frozen V6 holdout result

Four fresh targets; five frozen methods; maximum 3 s numerical wall per method/case. Top1 is selected by all-data residual before truth scoring. Oracle coverage is a separate offline diagnostic.

| Method | Top1 data-fit cases (residual ≤ 1e-6) | Mean top1 RMS | Median top1 RMS | Mean oracle best RMS | Mean numerical s |
|---|---:|---:|---:|---:|---:|
| ABAR4_coordpre | 1/4 | 8.68048 | 10.5601 | 4.57169 | 1.0643 |
| Coordinate | 1/4 | 8.18233 | 10.2709 | 5.44329 | 0.192883 |
| FullJ1024 | 1/4 | 7.91041 | 10.4076 | 6.30546 | 2.25324 |
| FullJReciprocal1024 | 2/4 | 4.8076 | 4.61388 | 3.76623 | 1.58415 |
| GlobalDual | 2/4 | 4.24797 | 4.15804 | 3.05578 | 1.80982 |

GlobalDual succeeds on cases 00 and 03, fails on 01 and 02; reciprocal FullJ succeeds on 00 and 01. Native FullJ, coordinate and ABAR4 succeed only on 01. GlobalDual therefore does not dominate the reciprocal control. On case 02 its lower residual and better oracle coverage do not produce a recovered image. These four finite targets support complementarity observations, not global enumeration, an accepted general advantage, or a novelty claim. The V6 scientific disposition remains HOLD.

All 3,500 archived candidates passed independent fresh NumPy full-VIE validation, at a cost of 10,500 additional full RHS. This validates physical evaluation and stored residuals; it does not upgrade unsuccessful reconstructions. Forty-nine raw dual endpoints that had merged with nearby archive candidates were additionally validated at 147 full RHS for exact phase analysis.

## Pure dual endpoints and full-J corrections

| Case | Pure endpoint residual top1 | Pure endpoint material RMS | Final residual top1 | Final top1 RMS | Endpoint full RHS | Correction full RHS | Inner QR / LS RHS |
|---|---:|---:|---:|---:|---:|---:|---:|
| HOLDOUT12_00 | 1.1319e-06 | 9.45654e-05 | 5.21095e-11 | 6.33262e-09 | 36 | 2562 | 1106 / 3318 |
| HOLDOUT12_01 | 0.259419 | 8.26475 | 0.250984 | 8.31609 | 426 | 50001 | 13738 / 41214 |
| HOLDOUT12_02 | 0.0801225 | 8.65379 | 0.0588328 | 8.67579 | 468 | 42789 | 15987 / 47961 |
| HOLDOUT12_03 | 3.60095e-06 | 0.000249047 | 1.10886e-15 | 1.51708e-13 | 156 | 16392 | 5211 / 15633 |

Phase counts are exact for these frozen receipts: completed correction fit counters plus raw endpoint fees reproduce the total full cost. QR/LS work is charged separately. The dual objective is a proposal metric; all rankings in this table use the independently evaluated original-data residual.

## Public-receipt complementarity audit

These rows inspect existing public prediction receipts only. They do not run a guarded pipeline or reuse V6 targets as new blind evidence.

| Case | Native FullJ first-start residual | First fit s | First fit available s | GlobalDual final top1 discovered s | GlobalDual goal reached |
|---|---:|---:|---:|---:|---|
| HOLDOUT12_00 | 0.400719 | 0.0101887 | 0.0226375 | 0.182953 | yes |
| HOLDOUT12_01 | 2.05891e-13 | 0.00425275 | 0.0116406 | 0.0394568 | no |
| HOLDOUT12_02 | 0.451407 | 0.00641021 | 0.0130833 | 0.622541 | no |
| HOLDOUT12_03 | 0.676565 | 0.0108665 | 0.0202905 | 1.05485 | yes |

No method was retuned after V6 truth scoring. `HOLDOUT_V6_METRICS.csv` preserves every case/method metric; `HOLDOUT_V6_SUMMARY.json` retains this summary plus first-fit operation counters. Complete paid exploration, setup, full/reduced/inner-LS work, independent validation and end-to-end subprocess time remain in `COST_LEDGER.json`.
