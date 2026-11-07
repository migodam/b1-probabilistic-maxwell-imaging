# R4 scientific freeze — 2026-11-07 09:40 UTC

Numerical science is complete; publication and artifact transfer are in progress. Six-hour deadline remains 10:44:27 UTC. Earlier checkpoints are retained as historical records.

- 504 TRAIN / 78 VAL / 108 independent TEST mothers, paired clean/noisy: all 1,380 strict GN anchors complete. Noisy strict endpoints are mostly budget stops, not claimed converged.
- All 1,116 VAL/TEST traditional control endpoints complete: weak regularization, original final-only discrepancy, and corrected all-scales discrepancy stopping.
- E0 historical disclosed-data audit: 3,264 method/observation cases, separately labelled and excluded from fresh TEST statistics.
- 18 from-scratch models plus two matched fine-tunes, 20 weights. Soft-alpha VAL choice was 1; this aliases J-full rather than inventing a separately trained sixth model. G-weak was the selected constrained family.
- Methods/checkpoints frozen at 08:03:14 UTC. All 26 runs x 216 TEST observations = 5,616 canonical physical candidate sets sealed before first GT scoring at approximately 09:02 UTC. No TEST-based selection followed.
- 5,510 COMPLETED and 106 BUDGET_STOP; the latter are noisy continue-TRF. Zero NOT_RUN canonical TEST pairs. All have data-qualified top1; this is distinct from joint material success.
- Local M4 handled 36 TEST mothers; i9 handled 72, using four own P-core workers, single-thread BLAS. The finite remote command returned and released its resources. RTX4060 was not taken from the other session.
- TEST true-current energy audit complete: 9,921 LU / 119,052 complex RHS, separate offline cost. Final energy, two-fold, loss, generalization, seed and cost figures complete.
- Main four image panels and complete 1,296-page lossless WebP atlas complete, with image/mother/class/style exposure labels. All planned failures retained.
- Three recurrent families trained at two passes. Separate VAL depth-16 diagnostic worsened true physical residual. No whole-network contraction or fixed-point success claim.
- Optional paired full-VIE fine-tune ran after the declared VAL trigger. Mixed TEST outcome retained; no universal loss/feedback advantage claimed.

Final scientific judgment: clean material-learning benefit supported, current-space superiority and fixed-point algorithm claims HOLD. No high-contrast generalization claim, no probability posterior, no full-root enumeration, no whole-system speedup. Original frozen physics unchanged.
