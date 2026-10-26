# Confirmation figure captions

**FIXED_PARENT_MATERIAL_A / B.** Fixed within-batch parent indices 0 and 1, chosen before looking at losses. Truth and raw predictions at frozen training seed 2026102601; real and imaginary parts have separate scales shared across both batches, all four selected parents and all methods. These images do not display the seed-mean loss estimator or an ensemble.

**RISK_AND_PAIRED_BATCH_EFFECTS.** Left: pooled raw full-grid complex MSE, averaging each parent's three frozen seed losses first; method intervals are descriptive. Right: paired improvements, positive favoring the named candidate, with Bonferroni four-comparison 98.75% bootstrap intervals. Pooled raw inference is primary; A/B intervals describe replication and cannot be selected as alternative primary tests.

**ANCHOR_PHASE_RISK_AND_GAIN.** Phase is fixed from the online anchor A q0 / halfwidth, using the original prior's ROI. Conditional phase-cell risk divides summed squared error by summed phase-cell count. Phase contributions divide each parent's phase squared-error sum by576 and average; contributions add to total full-grid gain versus StrongCNN. Transition is an anchor-defined prior saturation diagnostic, not a Maxwell-specific mechanism claim. Empty phase groups have no estimated conditional risk.

**RETENTION_AND_PAIRED_CONTRASTS.** The retention contrast is (Strong−Distilled)−0.8(Strong−K32). Intervals are two-sided98.75%; the triangle separately marks the predeclared one-sided98.75% lower contrast bound. The ratio is secondary and is displayed only with the frozen stable-positive-denominator qualification. Raw is primary, clipped secondary. This figure does not change frozen statistical decisions.

**RAW_CLIPPED_SEEDMEAN_AND_ENSEMBLE.** Mean of seed losses and loss of mean prediction are separate estimands. Raw and clipped losses occupy separate panels. Ensemble values are descriptive and never replace primary seed-mean losses.

**OFFLINE_ERROR_SELECTED_DIAGNOSTICS.** The two largest parent raw seed-mean MSE excesses of distilled versus StrongCNN, selected after offline scoring, with stable parent-order tie breaking. Images show mean predictions across seeds for explanation; the selection statistic is mean seed loss. This is explicitly OFFLINE_ERROR_SELECTED, never an online choice rule or a representative sample; excess may be nonpositive if no parent is worse.

**CPU_PROVIDER_AND_PREDICTION_COST, if present.** Measured CPU1 median and p10–p90 repeated-execution timings per whole batch, with batch1 and64 separated. Provider, full prediction and fresh-action load plus first prediction remain distinct. These timings exclude physical acquisition and anchor construction; filesystem/import caches may be warm. They do not establish anchor-inclusive end-to-end runtime or independent training-seed timing variability. Missing timing data produces no substitute cost chart.
