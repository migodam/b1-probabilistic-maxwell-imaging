# Fresh10 candidate mechanism check

Fresh10 bounded mechanism check, small N10; not paper qualification; all6arms fixed before any freshdata or DEVGT review by candidate worker

| Arm | coverage .075 | qualification | oracle MSE | uniform action MSE | likelihood heuristic MSE | diameter | rank | RHS256 | LU256 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| FULL66_RML | 3/10 | 1.000 | 0.00814439 | 0.00835286 | 0.00840117 | 0.1462 | 8.00 | 21504 | 948 |
| SHAPE_TARGET_RML_BLOCK | 3/10 | 1.000 | 0.00822191 | 0.00842637 | 0.00840052 | 0.1468 | 8.00 | 21672 | 941 |
| BILATERAL_TARGET_RML_BLOCK | 2/10 | 1.000 | 0.00832325 | 0.00838969 | 0.00840033 | 0.1526 | 8.00 | 22152 | 941 |
| FIXED_Q8_RAY | 3/10 | 0.613 | 0.00724536 | 0.00873832 | 0.00849669 | 0.1379 | 4.90 | 11392 | 332 |
| DYNAMIC_Q8_RAY | 3/10 | 0.613 | 0.00726201 | 0.00859198 | 0.00847781 | 0.1454 | 4.90 | 11384 | 331 |
| RANDOM66_FULL | 3/10 | 0.250 | 0.00813026 | 0.00838875 | 0.00840675 | 0.0394 | 2.00 | 12520 | 313 |

sameK8; optimization caps, prior-target distributions, radiusconstraints and correctorfreedoms differ across old/new arms and are explicitly structural; current3newarms matched; costprefix diagnostics, no costmatched whole-six claim

The shared anchor is a legal true-noise-qualified baseline whenever available. Oracle errors use private truth offline; uniform and exp(-Phi) actions use finite qualified clouds without posterior weighting claims. K1/2/4/8 cost-prefix rows are in METHOD_COMPARISON.csv. Fresh parents are independent, but N10 supports only a bounded mechanism check. All six methods and parameters were frozen before fresh generation; DEVGT did not select arms or retune them. All480onlineprediction records preceded fresh truth scoring.

Amplitude block Schur solves are ordinary joint GN, not nonlinear profiling. Coupled nonlinear profile remains NOT_RUN. New/old latent caps and radius/progress qualification differ explicitly; equalK does not imply equal physics budget, geometry constraints, or posterior reference.

Actual inclusive fresh cost=118760 RHS256/4313 LU256;0full576. See CANDIDATE_COST_AND_PROVENANCE.md for exact proposal/anchor/generation/wasted-adapter reconciliation and corrected fresh-J provenance.
