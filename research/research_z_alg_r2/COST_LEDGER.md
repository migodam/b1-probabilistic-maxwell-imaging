# Resource and cost audit

**Distinct actual full complex RHS: 53,951.** This preserves the verified prior R2 subtotal **53,903** and adds exactly **48** for the fresh six-material extension. Discrepancies: 0. No physics was run for this accounting update.

Each canonical event book is charged once. Model counters, target logs, pool summaries, and reused anchors are not additional charges. Full-source and material-J RHS are separate from small algebraic RHS.

| Phase | Full complex RHS | Known receipt wall (s) | Known receipt CPU (s) |
|---|---:|---:|---:|
| blind_engineer | 9,974 | 0.591672 | 0.566886 |
| main_branch | 9,710 | 30.546853 | 30.130062 |
| adaptive_NEW | 23,591 | 7.397774 | 6.977569 |
| quotient_prototype | 10 | 0.011729 | 0.005668 |
| quotient_benchmark_NEW | 10,514 | 3.501823 | 3.438385 |
| control64_NEW | 22 | 0.348282 | 0.319985 |
| holdout_solver_NEW | 64 | 2.149258 | 2.079704 |
| holdout_generation_NEW | 8 | 0.012935 | 0.002806 |
| CLI_reproduction | 10 | 0.005434 | 0.004435 |
| general6_extension_NEW | 48 | 1.612868 | 1.588242 |

## Full RHS breakdown

| Category | RHS |
|---|---:|
| forward_source_RHS | 32,656 |
| independent_validation_source_RHS | 748 |
| J_adjoint_RHS | 14,621 |
| J_tangent_RHS | 5,088 |
| reference_setup_RHS | 768 |
| generation_source_RHS | 70 |
| other_full_RHS | 0 |

Setup source and active-column subcounts are included in setup RHS, not added again. Historical adaptive anchor 551 RHS is charged only in the original blind phase. The unsaved 28-RHS cost remains reported with limited saved evidence; its timings/algebraic details are not reconstructed.

## Fresh six-material extension

- Total 48 new RHS =18 shared reference setup +21 fresh candidate validation +9 isolated generation. A separate pre-generation geometry audit used 0 full RHS.
- Three sources, one receiver, and 6 independently unknown real material cells;15 quartic minors,15×15 leading block, and 20dimensional quotient.
- The quotient found 2/1/1 qualified legal roots; the frozen64-start reduced Full-J found 1/1/1. All 60 quotient endpoints were numerically distinct,with 0 Newton polish steps.
- The extra same-observation material pair has RMS separation 1.1587903; fresh nonlinear response residuals are 1.94e-14 and 1.71e-14.
- Leading-block condition maximum 106.64; pairwise relative commutator maximum 3.38e-16; eigenvector condition maximum 17840; backsubstitution condition maximum 8.72. These are numerical diagnostics, not exact certification.
- Single-run total wall including setup, validation, and target logging: quotient 0.10442s versus reduced FullJ64 1.51803s. No timing repeats or post-score tuning were performed.
- The64-start baseline has 24 paid optimizer nonconvergences among 192 starts. Its missing root is a failure of this frozen budget, not a proof that stronger local searches cannot recover the root.

Evidence: [solver freeze](results/general6_SOLVER_FREEZE_RECEIPT.json),[prediction freeze](results/general6_PREDICTION_FREEZE_RECEIPT.json),[numeric audit](results/general6_numeric_audit.json),[extension resources](results/general6_RESOURCE_SUMMARY.json).

## Small algebraic work (not full RHS)

| Saved counter | Count |
|---|---:|
| reduced_complex_RHS | 143,150 |
| reduced_LU | 36,461 |
| algebraic_linear_solves | 291,792 |
| quotient_linear_solves | 140 |
| quotient_eigenproblems | 140 |
| quotient_recovery_least_squares | 876 |
| quotient_polish_linear_solves | 0 |
| algebraic_determinants | 3,477 |
| polynomial_evaluations | 366,957 |
| polynomial_J_evaluations | 365,475 |
| tracked_paths | 1,482 |
| determinant_permutation_terms | 1,080 |
| polynomial_linear_multiplications | 4,320 |

These selected counters describe saved reduced systems and algebraic operations. They exclude incomplete unsaved-development details. Small RHS are never added to the full-RHS total. General6 baseline counts actual cached small solves: 3 source RHS per new state and 1 receiver-transpose RHS per Jacobian.

## Negative evidence and limits

- Initial six-path studies did not achieve a total-wall GO against FullJ8 despite their full-RHS reduction.
- In the exposed four-material control,both eight-start misses were recovered by 64-start reduced Full-J.
- Fresh four-material holdout coverage was identical across all four methods: root counts 1/1/1/2.
- The old 24-path batch has 8 slices with endpoint collisions: 672 terminal paths represent 662 distinct endpoints.
- Every paid failed/nonconverged attempt is retained. The full inventory now contains 48 optimizer nonconvergences among 1920 inspected attempts;no saved paid physics-operation error was found.
- Numerical finite-dimensional root evidence does not establish exact completeness,continuous-PDE claims,connected components,global material coverage, or posterior inference.

## Hardware,timing gaps, and publication scope

Apple M4, arm64, 10 logical CPU cores, 24 GiB RAM; local CPU with one declared BLAS thread. Python 3.13.13, NumPy 2.5.3, SciPy 1.18.1. GPU jobs 0, cloud jobs 0, neural training runs 0. External compute spend reported 0; local hardware/energy cost was not metered.

Known wall/CPU sums are partial additive receipt measurements, not campaign elapsed time. The 28 unsaved development RHS and 45 initial aggregate generation RHS have no audited timings here. The separate general6 geometry audit has known wall but no CPU timing.

Canonical raw receipts and non-overlapping rows are linked in [COST_LEDGER.csv](COST_LEDGER.csv),[COST_LEDGER.json](COST_LEDGER.json), and [RESOURCE_SUMMARY.json](results/RESOURCE_SUMMARY.json). These retain the original phase totals and all known algebraic counters.

Publication/privacy review remains with the top-level reviewer. This ledger copies no private truth arrays,seeds, or score values and does not authorize publication.
