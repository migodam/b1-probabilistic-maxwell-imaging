# Final fresh V7 confirm result

The separately frozen guarded comparison fits one additional fresh target at the declared all-data residual threshold: GuardedGlobalDual 3/4 versus GuardedReciprocal 2/4. The sample is four targets. The independent-paper scientific disposition remains HOLD.

The material RMS values below are descriptive. No material-RMS success gate was predeclared; none is introduced retrospectively. A data fit means relative residual ≤1e-6. The stricter optimizer stopping goal is 1e-10; both methods meet that goal on 2/4 cases.

| Case | GuardedGlobalDual top1 RMS | GuardedGlobalDual residual | GuardedReciprocal top1 RMS | GuardedReciprocal residual |
|---|---:|---:|---:|---:|
| NEWCONFIRM12_00 | 2.11855e-10 | 5.9145e-12 | 10.9402 | 0.0970395 |
| NEWCONFIRM12_01 | 11.1686 | 0.107983 | 16.0621 | 0.132418 |
| NEWCONFIRM12_02 | 1.7067e-07 | 3.75382e-10 | 4.91879e-10 | 1.32611e-12 |
| NEWCONFIRM12_03 | 9.57434e-10 | 4.00144e-11 | 6.28527e-15 | 7.21951e-16 |

Case 00 supplies the additional data fit. Both methods fail on 01. Both fit 02 and 03; GuardedGlobalDual case 02 reaches residual 3.75e-10, which passes the declared data-fit criterion but misses its stricter search stop goal. None of these four fresh cases stops at the common native guard; all invoke a continuation. The paired guard materials, residuals and operation counts are identical, and every root cost ledger equals the exact sum of its phases.

| Method | Data-fit cases | Strict search-goal cases | Mean top1 RMS | Median top1 RMS | Mean outer s | Full RHS | Inner QR | Inner LS RHS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| guarded_globaldual | 3/4 | 2/4 | 2.79214 | 8.58136e-08 | 1.7187 | 116994 | 45389 | 136167 |
| guarded_reciprocal | 2/4 | 2/4 | 6.75056 | 5.47008 | 1.69244 | 207213 | 0 | 0 |

The total outer budget is 3 s including guard/setup/continuation/wrapper accounting. Actual terminal/bookkeeping overrun is reported: maximum 0.0128075 s. All 968 returned candidates pass independent full-VIE validation, charged 2,904 full RHS. Online top1 remains residual-selected, while best@K is labeled as offline oracle coverage.

V6 predictions/scores remain unchanged and serve as development evidence for this V7 revision. The four NEWCONFIRM12 targets were generated after the V7 source/config freeze. Their final eight predictions froze before validation or GT scoring. Target material arrays were disclosed only after scoring in `truth_after_freeze/`; random target-generation seeds stay local. Public disclosed-truth scoring reproduces every row exactly and is explicitly not a new blind trial.

See `CONFIRM_V7_METRICS.csv`, `CONFIRM_V7_SUMMARY.json`, `../COST_LEDGER.json`, and `../EXPERIMENT_TECHNICAL_NOTES.md` for the exact phase costs, original receipts and replay boundary. No physical experiments continue after `audit/EXPERIMENT_SCIENCE_FREEZE.json`.
