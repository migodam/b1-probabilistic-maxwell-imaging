# R3 implementation and mechanism audit

This audit reads exactly the 20 public holdout solver receipts listed in `MECHANISM_STATS.json` and the frozen V6 implementation. It reads no private material targets, scoring files, or independent validation receipt, and runs no new physics. All values below concern the online unweighted complex all-observation residual. They do not establish material reconstruction accuracy.

## Result and causal interpretation

The Q6 branch generator produces legal and often distant materials, but its direct proposals never improve the global best data residual in any holdout. ABAR improves the fixed eight-sweep coordinate arm on the three nontrivial cases only through subsequent ordinary FullJ correction. The all-receiver variable-projection arm changes the proposal metric and finds different basins: it succeeds on two cases, improves a third, and fails on the case that ordinary initialization solves immediately. Neither mechanism supports a uniform advantage.

| Case | Native FullJ | Reciprocal FullJ | Coordinate | ABAR+coordinate pre | GlobalDual+FullJ |
|---|---:|---:|---:|---:|---:|
| HOLDOUT12_00 | 0.0680012232 | 5.42032163e-13 | 0.159169949 | 0.093355582 | 5.21094867e-11 |
| HOLDOUT12_01 | 2.05890969e-13 | 3.69445431e-16 | 2.05890969e-13 | 2.05890969e-13 | 0.250984297 |
| HOLDOUT12_02 | 0.126533227 | 0.0916969995 | 0.257922942 | 0.169747046 | 0.0588327909 |
| HOLDOUT12_03 | 0.0457520823 | 0.0457520823 | 0.306735126 | 0.0579594687 | 1.1088615e-15 |

Native and reciprocal FullJ use the same physical start prefix and exact common forward/Jacobian implementation. All methods share cached public `G,E0,S`. A 3-second cap is a maximum search budget, not a requirement to spend the entire budget: coordinate and ABAR can finish their frozen finite schedules early.

## Exact coordinate and Q6 branch generation

The one-cell restriction is an exact fixed-complement rational VIE response. Its all-data squared loss is a real quadratic divided by a real quadratic; stationary candidates come from a degree-at-most-two derivative numerator. Endpoints and incumbent are included, poles rejected, and the all-data incumbent is retained. Eight sweeps are a finite coordinate procedure; they are not convergence certification.

| Case | Coordinate steps | Unique coordinate outputs | Coordinate improving steps | Q calls | New legal Q roots | Q roots RMS ≥1 | Direct Q global-best updates | Q GN corrections |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HOLDOUT12_00 | 96 | 287 | 68 | 288 | 118 | 103 | 0 | 32 |
| HOLDOUT12_01 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| HOLDOUT12_02 | 96 | 263 | 54 | 288 | 183 | 172 | 0 | 52 |
| HOLDOUT12_03 | 96 | 286 | 65 | 288 | 237 | 187 | 0 | 60 |

“New legal” means a physically filtered root that created a new archive entry; “RMS ≥1” is absolute material RMS distance from that root’s own incumbent. It is distinct from the beam separation rule and from the box-normalized archive merge tolerance. Repeated/unchanged roots are excluded from correction slots. Worse initial full-data proposals remain eligible for charged FullJ correction, while the incumbent stays archived.

| Case | Coordinate-pre best | Best new raw Q root from coordinate/initial incumbent | Best raw Q root anywhere later | ABAR final |
|---|---:|---:|---:|---:|
| HOLDOUT12_00 | 0.159169949 | 0.305063309 | 0.0945335908 | 0.093355582 |
| HOLDOUT12_01 | 2.05890969e-13 | — | — | 2.05890969e-13 |
| HOLDOUT12_02 | 0.257922942 | 0.262284975 | 0.178221983 | 0.169747046 |
| HOLDOUT12_03 | 0.306735126 | 0.321279088 | 0.0729200033 | 0.0579594687 |

No new root produced directly from a coordinate-pre/initial incumbent beats the coordinate-pre best. Some later raw Q archives have lower residual than the original coordinate result, but they occur after GN has already found better basins. The exact best-so-far trace has zero updates with a `q6_` origin across all four cases. Therefore those later raw archives cannot be used as evidence that direct Q solving caused the improvement. The observed causal route is legal branch proposal → FullJ basin correction. A comparison against other equally charged proposal generators would be needed to isolate a Q-specific benefit.

The Q6 quotient uses one receiver and two source channels at a time, current-complement exact Schur restriction, target τ∈{1,.5,.25,0}, physical conjugacy/material/box/pole filters, and all-observation scoring via a four-dimensional reduced solve. Numerical resolution of six endpoints is explicitly not a rigorous root-completeness certificate. Source rank, quotient conditioning and masked material-J regularity are different requirements.

## All-receiver variable projection

At the zero reference, q=E0, K=G and L=S. The observed-data constraint Wq=Y is represented as W=W0+CZ with W0=Y pinv(q), Zq=0. For material m, the shared receiver solve minimizes ||W−(L+WK)diag(m)||² over C. Each objective/gradient call shares one rectangular QR across all receiver right-hand sides. The envelope gradient is −2 Re Σ_receiver conjugate(error)·(L+WK). The box L-BFGS-B variables are all material coordinates.

This dual norm is a material-dependent metric, not the original data likelihood. It can distort the landscape, including near poles. Each completed endpoint is scored with a fresh full VIE forward, archive ranking uses the original data residual, and a new endpoint can receive the frozen FullJ polish budget. A low dual norm alone never qualifies an output.

| Case | Dual starts attempted | Scored endpoints | Best exact unmerged pure endpoint | Final after charged search/polish | QR | Rank failures |
|---|---:|---:|---:|---:|---:|---:|
| HOLDOUT12_00 | 12 | 12 | 1.13190109e-06 | 5.21094867e-11 | 1106 | 0 |
| HOLDOUT12_01 | 143 | 142 | 0.259419151 | 0.250984297 | 13738 | 0 |
| HOLDOUT12_02 | 157 | 156 | 0.0801226391 | 0.0588327909 | 15987 | 0 |
| HOLDOUT12_03 | 52 | 52 | 3.60094697e-06 | 1.1088615e-15 | 5211 | 0 |

The per-start records in `MECHANISM_STATS.json` give dual objective/termination, exact pure-proposal residual when retained unmerged, the archival representative otherwise, correction accounting, and exact global-improvement traces. Deduplication can merge a pure endpoint or corrected endpoint into a nearby existing material; the receipts then do not separately retain that endpoint’s true full residual. Such scores are marked unavailable rather than reconstructed by an uncharged forward. “Best exact unmerged pure endpoint” is therefore an observed retained minimum, not a guaranteed minimum over every paid endpoint.

The case01 failure is decisive negative evidence against universal improvement: GlobalDual uses its full 3-second budget and returns residual ≈0.251, while both ordinary FullJ parameterizations solve from their first shared low-material start. Case00 is also not a unique GlobalDual success: reciprocal FullJ reaches the data goal. The native/reciprocal difference demonstrates that an ordinary coordinate change can substantially alter the landscape.

## Generation, scoring and correction costs

| Case | GlobalDual QR / algebraic RHS | Endpoint-score full LU / RHS | All GN full LU / RHS | GN wall (s) | Total wall (s) |
|---|---:|---:|---:|---:|---:|
| HOLDOUT12_00 | 1106 / 3318 | 12 / 36 | 451 / 2562 | 0.098580 | 0.183145 |
| HOLDOUT12_01 | 13738 / 41214 | 142 / 426 | 8734 / 50001 | 1.889675 | 3.000532 |
| HOLDOUT12_02 | 15987 / 47961 | 156 / 468 | 7488 / 42789 | 1.698732 | 3.000360 |
| HOLDOUT12_03 | 5211 / 15633 | 52 / 156 | 2869 / 16392 | 0.652042 | 1.055258 |

These four GlobalDual receipts have no rank failures, no failed starts and no budget interruption inside a FullJ correction. Cases01/02 stop during the dual objective. Hence the sum of completed fit counters is the exact total GN cost for these receipts, and the remainder is exactly one 3-RHS endpoint forward per recorded scored endpoint. There is no unallocated full LU/RHS in this table.

Frozen V6 has a real accounting limitation: if budget exhaustion occurs inside a fit, that terminal fit row omits its counters. Total solver operation counters still remain exact. The audit preserves this limitation rather than claiming completed-fit sums are always sufficient. With a terminal GN interruption, endpoint-event evidence and failed-forward branches must be audited before assigning the residual work. No such ambiguity occurs in these four GlobalDual holdouts.

The dual-event wall includes proposal optimization, endpoint scoring and its subsequent correction. Subtracting recorded GN fit wall gives proposal+endpoint+bookkeeping time, not isolated pure-dual time. Exact QR work and separate operation timing are retained in the JSON. Endpoint scoring inside the solver and later independent frozen-candidate verification are separate costs; this audit does not read or count the latter receipt.

| Case | ABAR full LU / RHS | ABAR reduced LU / RHS | Q-correction full LU / RHS | Coordinate full RHS / reduced RHS |
|---|---:|---:|---:|---:|
| HOLDOUT12_00 | 1173 / 6930 | 885 / 3039 | 1142 / 6735 | 195 / 1698 |
| HOLDOUT12_01 | 13 / 78 | 0 / 0 | 0 / 0 | 78 / 0 |
| HOLDOUT12_02 | 2030 / 11721 | 876 / 3012 | 2002 / 11544 | 177 / 1656 |
| HOLDOUT12_03 | 2528 / 14463 | 913 / 3123 | 2492 / 14238 | 225 / 1692 |

ABAR full cost includes initial GN and the one all-N kernel setup in addition to the listed Q corrections. Reduced cost includes eight coordinate-pre sweeps, complementary Schur restrictions, and all-data rational scoring of Q roots. Frozen receipts have aggregate reduced counters without stage snapshots, so a precise Q-only versus coordinate-pre reduced split is not asserted. The all-N setup is charged as a full factorization with ns+N complex RHS; the method does not hide that setup behind a “small solve” label.

## Scope and implementation judgment

The implementation is consistent with the intended finite discrete nonlinear scalar-TM VIE model and online observation-only selection. It does not establish continuous-PDE claims, global inverse-image enumeration, uniqueness, posterior coverage, general imaging performance, or mathematical novelty. Q is an exact conditional slice representation with numerical enumeration; variable projection is an all-receiver alternative search metric with mandatory original-data scoring.

The safe claim is that these public finite holdouts show distinct optimization basins and an auditable proposal/correction mechanism. They also contain strong failures and ordinary-control successes. Any claim of improved material recovery must use the separate blinded scorer and independent verification evidence, neither of which this audit accessed.
