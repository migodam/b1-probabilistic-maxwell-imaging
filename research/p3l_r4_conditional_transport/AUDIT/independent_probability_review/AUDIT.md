# Independent probability-contract review

Status: COMPLETE for the source snapshot and CLOSED receipts recorded in `SNAPSHOT_SUMMARY.json`. This is a contract/algebra review, not a posterior-qualification or novelty verdict. No shared source was automatically fixed; PC-001 was reported to the parent, whose explicit exporter was then independently checked.

## Frozen inventory

`CLOSED_FITS.csv` and its JSON equivalent contain every discovered CLOSED density, anchor-density, field, and risk training attempt at snapshot time 1792873309854500000: **72 attempts, 71 completed eligible fits, 1 failed attempt**. Two live jobs are listed separately in `OPEN_JOBS_EXCLUDED.json`; no live metrics were read. `CLOSED_FAMILY_SUMMARY.csv` groups only identical module/kind/data panel/parameter count/frame/criterion/LR configurations and preserves individual seed/update lists. Later completions are outside this snapshot.

| Completed-fit panel | Count | Interpretation |
|---|---:|---|
| Raw joint TRAIN32768 / VAL2048 | 34 | Raw-y density/field fits; a global TRAIN-only frame may be an inductive bias |
| Raw joint TRAIN16384 / VAL2048 | 2 | Density scaling |
| Raw joint TRAIN4096 / VAL2048 | 5 | Scaling and LR controls, including update0 selections |
| Historical anchor TRAIN512 / VAL128 | 28 | Historical weaker anchor pipeline; distinct labels/context budget |
| Strong anchor TRAIN1024 / VAL256 | 2 | Closed strong-context risk fits only at this snapshot |

There are **no CLOSED strong-anchor density fits** in this snapshot. Do not describe the raw 32k density results as a comparison with an equally trained 1024-anchor density family. Strong1024/256 uses the expanded 96+64 context construction; historical512/128 uses the older 48+48 construction. Seed, parameter count, selected update, NLL/energy/MSE criterion, command, output, and CLOSED receipt are retained per row. Historical ROI/component initial losses must not be subtracted from final full576 complex MSE. Frozen R2 CNN has historical update400 but zero new optimizer updates. Classical controls' zero neural parameter count does not imply zero fitted coefficients. The original failed MPS watermark startup remains in the inventory beside its completed runtime-fix retry.

## Verified probability conventions

* `core.py:51,237` packs each source's 32 real then32 imaginary measurements. Its whitening is sqrt(2)/sigma, so `Phi=||complex residual||²/sigma²=.5||white residual||²`; independent relative error is 1.62e-16.
* Explicit Gaussian/spline/triangular families cover all66 Gaussianized coordinates. The reported density is with respect to raw Gaussianized `dx66`: TRAIN affine Jacobian, actual fixed-frame determinant, and conditional-coordinate determinant are included. Native uniform properties are `u=Phi(eta)`; a density with respect to native `dq` would additionally subtract the two normal log-pdfs. Current IS and pCN use Gaussianized x throughout, so that conversion is not missing from their target.
* Triangular factorization is `q(a|y)q(b|a,y)` with full innovation in all b coordinates. **Default first16 is a fixed coordinate split, not a Maxwell weak/strong frame.** The separately supplied TRAIN Fisher frame is the ascending eigensystem of the mean prior-relative likelihood Gram from256 TRAIN anchors, with no VAL/DEV/GT direction selection. Its16 weakest directions are global TRAIN directions, not object-specific material-risk directions.
* GN coordinates use positive-diagonal QR of `[Jx;I]`, `U.T U=I+Jx.T Jx`, and `r=U(x-x0)`. Raw log density adds `log|det U|`; inverse sampling uses a triangular solve. Random alignment is `U R.T`, with fixed TRAIN-independent orthogonal R. It rotates the raw induced covariance while preserving its spectrum. It is a different control from a global random fixed frame. Coordinate-isotropic base and explicitly transformed raw base have different intended covariance meanings. None is itself a GN posterior claim. Double-precision algebra is sound; the retained extreme float32 stress limitation must not be concealed by a numerical H-floor.
* Explicit importance weights are `exp(-Phi-.5||x||²-logq_x)`; normalized weights and ESS are correct. The independent flat-likelihood fixture gives equal weights and ESS9. **`log_mean_unnormalized_weight` is not normalized Bayesian evidence**: the Gaussian prior constant is omitted (33log(2pi)=60.649943), as is the normalized likelihood constant. Preserve the existing unnormalized label.
* Frozen full66 pCN uses SPD Q and the symmetric matrix `sqrt(I-Q)` from the same eigensystem. Its standard-Gaussian prior reversibility is consistent with the implemented full proposal ratio; cancellation error is 1.56e-13. Thus likelihood-only acceptance is valid for the fixed per-y kernel. This says nothing about finite-chain mixing. One chain per parent, finite ESS/MCSE, caps, and dependence remain descriptive; no cross-parent R-hat or common-reference qualification follows.

## Actionable findings and resolved blocker

Machine findings with status, evidence, and recommended action are in `FINDINGS.json`.

**PC-001, fixed for DEV128:** original canonical online anchor caches had `x`, not explicit `context_x`; the strict reader correctly rejected that ambiguous schema because supervised joint NPZ uses x for GT. The parent's exporter now verifies q/x, pure decoder chi, prediction/residual, canonical role, and source identity, then writes explicit context_x/context_J. All128 source/export byte+mtime identities and completion receipt were independently checked; first/last strict-reader shapes passed. Use this exported panel for online anchor inference. The original packed online_train companion still has ambiguous x; density_train separately has explicit context_x. Future panels require the same explicit export gate.

**PC-002, CFM numerical-density gate failed:** first4 already exposed DEV parents, complete66-coordinate divergence, NFEs16/32/64, frozen .05nat agreement threshold. Absolute32/64 NLL deltas are .114944458, .037582397, .058952332, .119895935: **3/4 fail**. All12 evaluations being finite does not repair this. CFM density remains a numerical continuous-ODE approximation, not the exact finite-RK sampling-map density, not an exact-MH proposal, and not a qualified posterior. Evidence: `../cfm_density_trace_audit_1792871479806924000/RESULT.json`. No tolerance relaxation or extra NFE rescue was performed.

**PC-004, optional hardening:** the paid audit checks finite sealed log_q_x for a declared explicit group but does not independently whitelist the producer kind. Current evaluator never exports CFM log_q_x, so current misuse fails closed. A future external producer should additionally be checked against algebraic Gaussian/spline/triangular kinds.

## Calibration and field-moment labels

Material energy score is the proper distinct-pair sample estimator under full-grid complex RMS (realification divided by24); projection CRPS uses the same convention. Randomized inclusive ranks handle plateau atoms using `(less+U*(ties+1))/(K+1)`. These marginal coordinate/projection rank diagnostics are distinct from ordinary empirical interval coverage and are **not simultaneous material credible bands**. Host-background coverage can be trivial; ROI diagnostics are separate. Paired bootstrap seed is20261025099 with4000 resamples; seeds are averaged within parent before parent resampling.

Gamma256(shape256,scale1) is a **prior-predictive marginal** benchmark across independent parents for the fixed sample index0, not a claim that Phi has that law conditional on one y. Paid likelihood-rank SBC is a data-dependent necessary diagnostic. Repeated/rejected chain states stay in the roster; dependence limits IID rank and importance-ESS interpretation. Sampling/NLL/energy alone cannot mark Full Maxwell likelihood, Gamma or data-dependent SBC evaluated. Those require the separate paid predictions and closed online gate. The paid driver preserves fixed indices/failures and separates online physical work from different-PID truth scoring; no residual-based candidate selection is allowed.

Field heads regress conditional per-pixel raw moments E[Re chi], E[Im chi], E[Re²], E[Im²], E[Re Im]. Their2x2 covariance is a plug-in pointwise action and can be indefinite; PSD projection changes that estimate and is reported. There is no cross-pixel covariance, joint spatial sampler, latent law, or ambient Gaussian posterior density. The exact original-P1 prior artifact is unconditional only: prior independence permits E(delta)E(h), not conditional factorization after observing y. Its source analytic moments and no-training audit are separate from the fit CSV.

The5842-parameter local CNN is the reoptimized R2 function class, as verified by `../risk/CNN_CLASS_EQUIVALENCE.json` (max replay error1.20e-7); do not label it a new architecture or established Maxwell mechanism. No Full posterior qualification is established by this review, improved NLL, energy, or risk fit.

## Cost and provenance

Final collector receipt: `../../I/jobs/independent_probability_review_data_1792873309751008000/CLOSED.json`: COMPLETE, wall0.255344s, CPU0.102572s, peakRSS238370816B, no GPU, no new RHS/states/LU/optimizer updates. `SNAPSHOT_SUMMARY.json` records source byte/mtime identities and earlier scaffold receipt transparently. These are collector measurements; static human/agent review time is not a physical-solver measurement. The separate CFM audit used wall5.068494s, CPU4.852681s, RSS236896256B, zero physics and updates. No training, new physical query, model modification, or shared-source autofix was performed in this review.
