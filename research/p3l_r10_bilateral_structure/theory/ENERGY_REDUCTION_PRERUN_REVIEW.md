# Independent pre-run review: energy and residual-qualified reduction

Reviewed 2026-10-29 18:02 UTC. Read-only review by the R10 theory worker. Reviewed `PHYSICS_ENERGY_AND_RECIPROCITY.md`, `config/reduced_response_frozen.json`, `src/reduced_response.py`, and the shape construction in `src/transfer_campaign.py`. No execution, no Maxwell RHS, no numerical identity test, no edits to parent code.

**Decision: mathematics PASS; execution can proceed under the frozen budget. Comparative setup accounting requires the clarification below before any speed claim.** This review does not certify unseen shapes, numerical rounding, continuous Maxwell, or strict minimal realization.

## Mathematical checks

| Item | Finding | Status |
|---|---|---|
| TM energy convention | For \(e^{-i\omega t}\), the H expression, Poynting flux, Helmholtz source sign, absorption sign and source-work sign agree. | PASS |
| Discrete current identity | From \(e=E-Gj\), \(j=D_mE\), the sign of \(j^*\operatorname{Im}_H(G)j\) is positive on the right. Interpretation as nonnegative radiation still needs the declared PSD/normalization assumptions. | PASS |
| T dissipative identity | \(\operatorname{Im}_H A^{-1}=-A^{-*}\operatorname{Im}_H(A)A^{-1}\) yields the stated T identity with the correct order. Zero contrast is appropriately excluded from inverse-D reasoning. | PASS |
| Coercivity | \(-\operatorname{Im}(v^*K_\delta v)\ge b_\delta\|v\|^2\) implies \(\|K_\delta v\|\ge b_\delta\|v\|\), regardless of normality. | PASS |
| Rectangle bound | On Re delta in [.8,1.8], Im delta in [.03,.18], the minimum of Im delta / abs(delta)^2 occurs at (1.8,.03). Negative beta is retained, not silently clipped. | PASS |
| Reduced solve | Schur of \(V^*BV\), solving \((\delta^{-1}I-T_r)u=Q_r^*V^*R\), gives physical approximate state \(VQ_ru\). | PASS |
| Full residual | `R-(VrQ@u)/delta+BVQ@u` is exactly \(R-K_\delta VQ_ru\) in exact arithmetic. Cached BVQ does not change the residual's mathematical meaning. | PASS |
| Port bound | Code uses Frobenius norm of L, while E05 states spectral norm. Since \(\|L\|_2\le\|L\|_F\), the implementation is a conservative valid bound. It must be labeled as Frobenius-based in reporting. | PASS with labeling |
| Whitening | \(\sqrt2\|\Delta F\|_F/\sigma\) matches the original realified complex residual convention. | PASS |
| Likelihood potential | \(\|\widetilde r_y\|e+e^2/2\) is the valid absolute Gaussian-potential error bound. It uses unperturbed observation y. | PASS |
| Acceptance logic | Both declared bounds and positive numerical coercivity are required. Actual full-reference error is logged, not used as the acceptance gate. | PASS |
| Source/support | Transfer construction uses original nonnegative quintic occupancy, \(P=\sqrt a\), and original complex amplitude rectangle. This conditional family is legal for current P1. | PASS |

Eigenvalue evaluation, residual arithmetic and norms are floating-point quantities without interval enclosure. The record correctly calls the result a numerical discrete residual bound. The analytic estimate is uniform in the rectangle given an exact certified beta, but current sampled acceptance outcomes do not establish a certified all-amplitude or unseen-shape result.

## Findings for execution and reporting

### R-COST-01 — Include exact-fallback setup in each method total

`setup_conservative_seconds` adds energy setup, the full basis-to-128 setup and the selected reduced model setup. `online_with_fallback_seconds` adds full-Schur node evaluation when a reduced node fails. However, the first quantity does not add the full-Schur **setup** for a parent that needs fallback. The full setup is available separately as `full_setup_seconds`/`SETUP_COSTS.csv`, and the physical experiment ledger has already paid it; the remaining risk is an understated *standalone method comparison*.

For fixed rank r, the standalone total must be
\[
C_r=\sum_p(C_{energy,p}+C_{basis128,p}+C_{model,r,p})
+\sum_{p:\exists\,failed\ node}C_{fullsetup,p}
+\sum_{p,j}(C_{solve,r,pj}+C_{guard,r,pj}+1_{fail}C_{fullnode,pj}).
\]
For the current adaptive path, all precomputed rank-model setups should be charged unless a distinct on-demand construction is actually implemented. Add the exact full setup for every parent whose adaptive path falls back. The exact reference comparator includes its full setup plus all exact node evaluations. If full Schur is declared a pre-existing cache instead, show both cache-reuse and standalone views consistently for every arm.

**Status: REPORTING_HOLD until these totals are included.** This does not block a physical run because all necessary raw costs are already captured; it blocks an unqualified speed claim from `setup_conservative_seconds + online_with_fallback_seconds` alone.

### R-ROBUST-02 — Early Krylov breakdown below the first listed order

`basis` may legally stop before rank 32. The current code then leaves `models` empty, followed by `max` on empty lists and `rr[0]` access. A low-dimensional invariant reachable space would be a meaningful result rather than evidence of failure. Before execution, prefer handling this case by an explicitly recorded available-rank model or exact fallback, without silently inventing a new accepted rank sweep. At minimum, preserve a clear FAILED receipt with all spent costs and do not interpret incomplete outputs.

**Status: ROBUSTNESS_GAP, not a demonstrated mathematical flaw in these eight inputs.** This review does not assume or measure their achieved Krylov rank.

## Scope of a possible result

The first sampled rank satisfying the two numerical bounds is an **empirically useful order for the stated nodes and basis family**. It is not the minimal realization order, nor even the smallest possible accurate rank over all reduced bases. Failure to pass at rank 128 would limit this Krylov/coercivity combination, not prove that all compression is impossible. A successful guard with no end-to-end savings would establish an accurate approximation without computational benefit. Preserve all three distinctions in the paper.

The conservative RHS cap is arithmetically consistent: 8 parents × 30 nodes × 8 sources × 5 orders = 9,600 triangular RHS, comprising 1,920 at order 256 and 7,680 at lower orders if every model is present. Orders must continue to be listed separately; they are not cost-equivalent merely because the conservative budget counts them together.

## Completion addendum, 2026-10-29 UTC

The preceding findings are preserved as the pre-run record. Current code and result readback closes R-COST-01 for the reported **fixed-rank cold totals**: `fallback_Schur_setup_seconds` is separately recorded and added to the cold cost. The order-32 result charges the full Schur setup for all eight fallback parents. The independent compression audit also reports agreement. Adaptive online rows must still not be presented alone as setup-inclusive cost. No such adaptive speed claim is used in the finalized narratives.

R-ROBUST-02 was addressed in the DEV probe by an explicit `KRYLOV_RANK_BELOW32` failure with preserved cost, rather than an empty-model crash. The separately frozen heldout-shape probe handles all-host occupancy analytically and unavailable rank 64 through exact fallback. Read-only code review confirms these branches; they were not triggered in the reported data, and no extra physics test was run for them.

Completed outcomes: order 64 passes 240/240 DEV nodes and 320/320 heldout-shape nodes, with zero fallback in the latter. The heldout M32 workload pays its own formation, coercivity, basis-to-64, projected Schur, transforms, solves, and full-residual guards. It costs 0.115014 seconds versus 0.317346 seconds for full Schur and 0.238701 seconds for direct LU. These are cold **kernel** costs; complete job time, shared host loading, and artifact writes are reported separately. Other M values are extrapolated, not additional executed workloads. This closes the fixed-rank cost qualification while retaining the numerical-bound, empirical-order, and non-blinded-material limitations above.
