# MAX-STRUCT-R1 baseline and dynamics audit

This audit reads the frozen FBG-R3 reports, mathematical audit, actual CPU/CUDA CIE source, E3 scorer and selected raw receipts. It does not rerun a 576-cell inverse task. D1 and D3 are new **OFFLINE_DIAGNOSTIC** calculations on the existing genuine 4×4 Hankel–Green operator/data geometry. They establish no new imaging or blind-test result.

## Historical result that must remain intact

- R2's **8/12** is a historical result, not a fresh blind result of this study.
- The disclosed R3-E mature dual-beta CIE plus common original-data corrector has **10/12** joint passes: all nine particular digit rasters 1–9 and `two_islands`; **`ring` and `curved_bar` fail**. “All nine” describes this finite disclosed set, not arbitrary digits, fonts, noise or amplitudes.
- This E3 scorer explicitly says `POST_CLOSURE_SCORED` and `DISCLOSED_HISTORICAL_MATURE_BASELINE_DIAGNOSTIC`. Historical dataset disclosure and prior development are not reversed by delayed final scoring. The readout identity and R3 derivative mathematics remain attributed to R3, not to MAX-STRUCT-R1.
- E3's material denominator is amplitude **12**, verified by row metadata and scorer implementation. A stale generic definition mentioning “8 or18” cannot replace that executed contract.

On these 12 objects the gate is amplitude-normalized NRMSE≤.15, reconstruction-Otsu IoU≥.7, original physical relative residual≤.001. The selected estimate is the online original-Maxwell objective incumbent, not GT-ranked `best@K`. There are 576 independent real material cells at the final stage; Fourier restriction acts on auxiliary current, not material labels or a digit decoder.

## Faithfulness and adaptations of the mature baseline

The actual code uses receiver SVD DPCS, projected Fourier APCS in native coefficient coordinates, fixed positive source-normalized residual weights, one sourcewise PR+ current step with exact quadratic line minimization and non-descent reset, shared-real material alternation, zero-padding transport at mode expansion, and stagewise beta schedules. PR memory resets at stage boundaries. Dual beta restores one physical `m` but uses different current/material metrics; it is not a two-block descent algorithm for one common scalar objective.

The material readout clips the modified contrast in the legal monotone image of the material box before inverse transport. When the common smooth/TV prior is active, an actual bounded material inner solve replaces the uncoupled closed-form readout; the analytic derivative includes its off-manifold term. This is a disclosed real-material adaptation, not literal complex-material/3D reproduction. Dense explicit Fourier columns do not establish FFT complexity or speed. General correlated joint whitening was absent from the historical implementation.

The original paper coefficient stop and the separately frozen adequate-duration stop must stay separate. At 4×4, `mf=1` and `mf=2` both span the same rank-four APCS after receiver rank12 extraction; the schedule changes beta/coordinates but not independent capacity. Calibration alone was therefore not a mask-growth imaging result. E3 later supplies a separate 576-cell budget-point experiment, still wall-limited.

## D1: local exact-current material dynamics

Source: `src/dynamics_audit.py`; evidence: `research/delegated/audit/D1.json`. Fix the full current space, six source columns, beta_current=1 and one disclosed nonconstant real material root on the loaded genuine geometry. The current block is solved exactly by small Hermitian positive-definite Cholesky systems. With no prior and an interior material root, the composed material map has analytic derivative `M`; the derivative of the real shared-material readout at a physical root is independent of beta_material. This is a consequence/application of the existing R3 structure, not a new readout identity.

For beta_material .2/1/4, root defects are about 1e-15. Central directional differences agree with the same analytic `M` to maximum errors 5.15e-8, 5.16e-10 and 1.09e-10 at perturbations 1e-3/1e-4/1e-5. The local 16×16 material derivative has spectral radius **.8366359** and 2-norm **.8374044**. This is a local statement on one small exact-current root; it does not locate a basin or imply global recovery.

Away from the root, beta choices separate. For perturbations .01/.05/.2, the maximum material-map difference against beta_material=1 is 2.21e-6/5.45e-5/8.30e-4; division by perturbation squared gives .02208/.02179/.02076. Thus equal root first derivatives coexist with finite-step differences.

Three single-assumption controls locate the boundary:

1. The same common material prior makes the unregularized physical root cease to be a fixed point. Root defects for beta_material .2/1/4 are 8.97e-5/6.87e-6/1.89e-6. There is then no justification for applying the prior-free root contraction directly.
2. A material root on the box's lower face has a nonsmooth clipped extension. Opposite directional responses have sum norm .72253. Their first-order beta coincidence survives in this control; the box does **not** automatically refute beta invariance. The negative test direction leaves the feasible input box and is explicitly an external nonsmoothness diagnostic, not a permitted solver step. Within the box one must use a tangent-cone analysis.
3. Replacing exact current minimization with one freshly reset PR/steepest-descent step changes the dynamics. The directional additivity defect tends to **.05358**, and its tested directional derivative differs from `M h` by **.42141**. Exact line search has a direction-dependent Rayleigh quotient at zero gradient: this control cannot be assigned the exact-current map's Fréchet Jacobian or contraction spectrum. Carried historical PR memory would add further state; it is not reconstructed here.

## D2: actual historical successes and failures

Evidence: `research/delegated/audit/D2_historical_evidence.json` extracts actual receipts for dual-CIE `digit3`, `ring`, and `curved_bar` without opening raw truth. Their scalar iteration histories contain 1191/973/1246 entries including stage ends. Each current trajectory terminates at **WALL_BUDGET**. `digit3` reaches mf10 and later passes (residual .0004178, NRMSE .11764, IoU .76238); `ring` reaches only mf6 and fails (residual .16740, NRMSE .42020, IoU .36975); `curved_bar` reaches mf10 and fails (residual .11188, NRMSE .39312, IoU .40777). The final physical corrector is separate from the current trajectory.

The scalar logs record residual terms, relative coefficient changes and inner statuses, but omit per-iterate coefficients, current gradients, search directions, PR coefficients, line coefficients and full material states. The full dynamical state and PR recurrence are therefore unavailable. We cannot reconstruct periods, estimate trajectory Lyapunov exponents, classify fixed points, or claim chaos. A wall cutoff is not a basin fixed point; a repeated scalar loss or small coefficient change would not establish recurrence of the full state either. No new expensive replay was justified by these logs.

## D3: source-basis covariance of complete iterations

Evidence: `research/delegated/audit/D3.json` and `D3_trajectories.npz`. Mix six incident columns by a fixed complex unitary `U`: `E'=E U`, `Y'=Y U`, `J'=J U`. For the quadratic metric `Re tr(R* R Omega)`, transport **Omega'=U* Omega U**. Both data and state metrics are transported. Every arm includes the actual coupled residual, analytic current gradient, legal shared-real material readout, 16 fixed iterations and the same initial/stage memory reset. beta_current=1 and beta_material=.2 are fixed, so the recorded current metric objective is not asserted monotone through the material step.

| Arm | Initial objective original / mixed | Final material difference after undoing U |
|---|---:|---:|
| Native per-column normalization recomputed |17.999975 / 18.153001|.00296891|
| Full transported metrics + sourcewise search/PR |17.999975 / 17.999975|.00324540|
| Full transported metrics + global PR and global exact line search |17.999975 / 17.999975|8.47e-16|

The first arm changes the objective itself. The second preserves that objective but retains a source-coordinate-dependent iterative rule. Its line coefficients are solved as a **coupled real 6×6 subspace minimization** under offdiagonal source metrics; this is a diagnostic extension of the independent-column historical line step. It already differs at iteration1, so the residual basis dependence belongs to sourcewise search/line minimization as well as subsequent PR memory. It cannot be blamed on PR memory alone. The third transports the complete global iteration and agrees to roundoff. This is a dynamics/covariance result, not an imaging advantage or a proposed replacement baseline.

## Decision

**PASS within stated diagnostic scope:** exact-current interior root derivative/spectrum and complete source-covariance controls are implemented on genuine small Maxwell operators. **HOLD for historical dynamical classification:** absent PR-state history and wall truncation block fixed-point/cycle/chaos attribution. **NOT_RUN:** new blind imaging, 576-cell replay, trajectory stability study, global-basin claim or a new method comparison. Costs and provenance limits are in `audit/data_cost_review.md`.
