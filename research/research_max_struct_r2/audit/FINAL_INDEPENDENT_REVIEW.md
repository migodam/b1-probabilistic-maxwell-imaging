# R2 independent scientific/data/cost review

Final scope: all12 original blind predictions, post-closure material scoring, two no-GT diagnostics, one explicitly post-hoc capacity control, and the separate restricted-family3D acquisition demonstration. No inverse or forward solver was launched by this final review/scorer. Generation was performed once by the separate authorized generator stage and charged separately. R1 source/evidence is unchanged.

## Verdict

**The final numerical claims are supported within their stated finite scope. No new imaging superiority over mature CIE is established. PAPER_HOLD is appropriate.** Native/global both recover the same two of three fresh objects; global provides no extra joint coverage. Known full-current profile recovers only one, and direct physical control none. The failed curved bar remains failed after the one enlarged-capacity development run. These results do not establish universal failure, global convergence, basin boundaries, continuous uniqueness or a new theorem.

The final `RESEARCH_REPORT_ZH.md` correctly preserves finite work, prior/source-metric differences, all failures, M/T/ROM qualification stops, historical attribution and active-family scope. Its rendering-inclusive measured time must use the final cost ledger value **103.807448s**, including .025044s of diagnostics-summary assembly. The physical/current operation totals are unaffected by that metadata work.

## Blindness and score independence

The generator wrote `configs/blind_cases_frozen.json` before running once. It draws fresh private128-bit entropy, realizes three new geometries, uses independent opaque IDs and emits only public operators/acquisition, noisyY, calibrated sigma/radius and known common priors. No alpha, mask, material, cleanY, family label or seed enters public NPZ or the inverse runner. The independent R1 read audit opens only the R1 tree.

`experiments/run_blind.py` checks an explicit public-key allowlist. Its original-physical objective selection includes the common initial material, paid auxiliary candidates and common corrector; it does not import the generator/scorer or open private truth. Source inspection and recorded `truth_access=NONE` support this dataflow statement; it is not a claim of operating-system isolation against a malicious worker.

The scorer independently checked all12 completed receipts and finite576-cell predictions and wrote `audit/PREDICTIONS_SEALED.json` at **06:09:07.838542UTC**. First private2D truth access occurred at **06:09:12.846938UTC**, after seal/path/size/mtime verification. No hashes were computed. Scores are for the sealed online top1; no GT-ranked replacement or oracle candidate is substituted. The12 original predictions and their denominator remain unchanged after scoring.

Physical residual/objective were independently recomputed from the already-paid stored `Y_prediction`, the public observedY and unchanged prior. Differences from raw receipts are at floating roundoff (physical objective at most about5.6e-17). This adds no forward solve and does not independently qualify the continuum model. All predicted materials satisfy the common box; no candidate was rescued by clipping during scoring.

## Actual gates, not inferred success

The noise is iid circular complex Gaussian, real/imag variance sigma²/2; it differs from historical fixed-norm noise. The marginal99% per-object whitened radius is25.165422. Primary joint success requires this data gate, amplitude12-normalized NRMSE≤.15 and reconstruction-only256-bin Otsu IoU≥.7. Historical noisy relative residual≤.003 is separately reported, not silently substituted. The prior is exactly historical normalized-data smooth/TV: in a unit-weight whitened objective its multiplier would be ||Y||²/sigma², not the reciprocal.

|Original fresh case|Native joint|Global joint|Profile joint|Direct joint|
|---|---|---|---|---|
|Annulus|PASS|PASS|PASS|FAIL|
|Curved bar|FAIL|FAIL|FAIL|FAIL|
|Two-island control|PASS|PASS|FAIL|FAIL|

For native/global, annulus NRMSE is .056567/.056045 with IoU .908163 for both; two-island NRMSE is .077698/.077222 with IoU .925373/.878788. Thus a tiny median NRMSE change does not establish dominance: the control's IoU worsens and joint coverage is unchanged. On the failed curved bar global improves NRMSE .510695→.388216 and IoU .245902→.418367, while original residual worsens .100735→.107831. No noise/material/shape gate passes there. Profile's two-island NRMSE1.206119 is retained in the table and fixed-scale figure.

The actual GT noise realizations all lie inside the frozen radius. Post-closure GT self-Otsu IoUs are .98913/.97561/1.0, so the image gates are not impossible due to the earlier coarse-grid self-IoU problem. These are diagnostics only; the generator did not reject or regenerate a realization using them. Three particular objects support no population success-rate, significance, basin-volume or universal morphology claim.

The3×5 `figures/BLIND_RECONSTRUCTIONS.png/.pdf` uses shared [.02,20] contrast limits, identical raw array orientation, all four methods and all failures. Visual inspection confirms readable labels and no per-panel contrast normalization. The numerical labels are taken directly from independent scored rows.

## Strong baseline and implementation audit

Native/global each execute all1500 scheduled current steps (mf3→6→10,500 each, coefficient_tol0) without a short wall cap, with the same box, initialization, historical prior, DPCS cutoff, material inner solver and1000-iteration physical corrector. The native solver is the frozen mature CPU implementation with a forwarding-only snapshot wrapper. Its mask basis and all material variables remain genuine; no shape decoder or known support is introduced.

The global extension holds full positive source metrics fixed, transports them under a unitary change, uses the coupled gradient, Cholesky-whitened shared readout and a global PR/line-search/stop. The new16-step wiring check agrees to1.43e-14 in material trajectories. This checks implementation covariance, not image superiority. Main native/global source normalizations agree; changing only the old scalar CG flag is not asserted to support arbitrary full correlated whitening.

The known profile arm uses a full576-current Cholesky and an envelope derivative retaining `w' e = -w²e`; its one disclosed-root gradient check has absolute directional error6.58e-13. The code's total-energy scalar source normalization differs from native per-source normalization. Consequently profile-versus-native differences cannot be attributed solely to exact current elimination. Profile has250 outer material iterations, a different inner work unit and common physical correction. Its negative result is not a theorem against all variable projection/WRI/CSI implementations.

The direct runner performs **two physical L-BFGS passes each capped at1000 iterations**, since the direct auxiliary pass is followed by the shared corrector. It is not1000 total steps. Actual early stops are retained; this stronger allowance is disclosed rather than retrospectively rerun. Equal scheduled current steps, different profile factorization work and early direct stopping are not equal CPU/FLOPs/search opportunity.

## Failure explanation remains conditional

The no-GT endpoint diagnostic finds all six native/global corrected endpoints above the optimizer's1e-8 projected-gradient threshold: approximately1.4e-8–2.7e-6. Five correctors hit1000 iterations; global curved bar stopped by relative objective change with projected gradient2.73e-6. It is therefore unjustified to call these certified wrong local minima or true basins, despite active box faces and high residual.

At the failed native/global endpoints, ordinary current-L2 projection error is25.76%/30.20% in the mf10 affine space. Freeing all DPCS coefficients plus mf12 spans the complete current space; retaining the actual noisy-data-fixed DPCS leaves4.59%/2.64% mismatch. This is a representation diagnostic of final predicted physical currents, not a projection of GT, a feasible PDE/material correction, or an H(curl) statement. Observed one-residual amplification is not an operator-norm resolvent certificate. M's compatible H(curl) contract and T's usable uniform exception/error bounds remain unfilled; ROM remains NOT_RUN.

The parent's one further mf12 native experiment was selected after GT disclosure. Its receipt correctly says **POST_HOC_DEVELOPMENT_AFTER_GT_RELEASE**. The runner did not read GT, but that alone does not make the selected experiment fresh blind. Separate scoring in `results/POSTHOC_CAPACITY_SCORE.json` gives NRMSE .558744, IoU .257813 and residual .097693: lower data residual than the original native .100735, worse NRMSE than .510695, all gates still FAIL. Whitened residual2344.64 is far above25.165; projected gradient3.13e-7 is not converged to the requested tolerance. More complete APCS capacity did not rescue this particular finite-budget run, without proving unlimited work cannot help. No blind table or denominator changes.

## Active3D: extra information with a strong family prior

The active example reuses the qualified R1 finite vector DDA and a public one-dimensional family `epsilon=I+t*b*zz^T`, t∈[0,3], with known smooth bump/geometry. Planning uses only public endpoint transfer differences. Acquisition, fitting and post-fit score are separate invocations; fitting reads public observations and not `active_truth`. The acquisition seed is fixed in shared source, so this should be called a scripted restricted-family demonstration, **not a strict fresh-blind unknown-material benchmark**. The separate primary2D secret-entropy blind test is unaffected.

All old-family members have identical zero old scattered field, and the old noisy observation fits the whole interval. New observation yields t_hat1.762889 versus t1.798734, absolute error .035844, with whitened residual2.960436<6. The fine-grid same-parameter prediction gap .561382 is a finite-model discretization control, not an independent solver, continuum certificate, global-minimizer or injectivity proof. Source coefficient power is idealized, not an antenna calibration. One coherent new observation is acquired synthetically, while all four library plus old source RHS are computationally charged.

The transverse-component coverage target is a separate disclosed off-family control. Its entire candidate-family best possible old-data residual14.89030 exceeds radius8.48528, correctly rejecting the whole family. It is not a generic inversion failure and does not trigger retrospective prior expansion. Neither the one-parameter fit nor this control establishes superiority over576-variable CIE or arbitrary tensor recovery.

## Final costs and gaps

`experiments/consolidate_costs.py` reads every outer task/process receipt once. `COST_LEDGER.md`, `audit/COST_LEDGER.csv` and `audit/COST_TOTALS.json` retain source references and nested counter details. Original inverse totals are9279 physical LU/222048 complex RHS; generation adds3/36, profile wiring1/6, active16/92, endpoint/source wiring11/174 and post-hoc capacity951/22740. **Final physical total:10261 LU,245096 complex RHS =122996 forward +122100 adjoint.** There are no uncharged new forward validations in this review/scorer.

Profile adds769 dense576×576 current Cholesky/9228 current RHS. Gradient wiring adds3 small16×16 current Cholesky/18 RHS. Hence current totals are **772 Cholesky/9246 current RHS**, separate from physical LU/RHS. Source-metric Cholesky, basis SVD/QR and Green/Gram products have raw counters but their complete internals are not exhaustively aggregated; unknown totals are null with reasons rather than zero. LU dimensions differ across the scalar576, diagnostic16, vector192 and vector1536 problems; counts are not equal FLOPs.

Outer numerical process sums are **102.699278s wall /103.537034s CPU**. Including measured seal, material scoring, plots, post-hoc scoring and diagnostics-summary metadata gives **103.807448s measured wall**. The known CPU component sum is104.344143s; full CPU is null because ACTIVE_PLOT CPU and other overhead were not measured. Nested solver/current/corrector/setup/state timers are not added again. Process-time sums are neither research elapsed time nor a speed comparison.

No original inverse disposition is a runtime failure; MAXITER/relative-function stops and scientific gate failures remain exposed. Failed or repeated tooling/retrieval/startup attempts without a numerical receipt, interpreter imports, source preparation, some JSON/CSV serialization, independent read-only audit, packaging/publication and historical operator assembly have incomplete or absent timers and remain unmetered. Monetary fees are unknown. No new SHA check or hash calculation is introduced. These accounting limits must survive publication.
