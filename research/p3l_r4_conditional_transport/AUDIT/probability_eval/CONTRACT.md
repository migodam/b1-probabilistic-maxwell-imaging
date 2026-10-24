# Independent full66 generator evaluation

Implemented in `src/p3r4/probability_eval.py` and `tools/evaluate_r4_density.py`. No physics backend, fitting, importance correction, filtering or candidate selection is imported/launched. Bootstrap uses the protocol seed **20261025099**,4000 resamples of independent parents, and98.333333% pointwise intervals. Training-seed losses are averaged within each parent before family inference; samples/images are not ensemble-averaged across training seeds.

The interface is three separate processes:

```text
python tools/evaluate_r4_density.py freeze --methods SPEC.json --out FROZEN.json
python tools/evaluate_r4_density.py infer --freeze FROZEN.json --contexts research/p3l_r4_conditional_transport/inputs/online --out NEW_ONLINE --device cpu --batch-size 16
python tools/evaluate_r4_density.py score --freeze FROZEN.json --online NEW_ONLINE --truth research/p3l_r4_conditional_transport/inputs/offline_only --out NEW_SCORES --device cpu
```

SPEC has `parents`, and `methods` entries `{id,family,seed,checkpoint}`. Checkpoint format must be R4 full66; seed must match the checkpoint. Optional comparisons predeclare `candidate_family`, `reference_family`, `metric` (`raw_sample_mean_full_mse`, `material_energy_score`, `full66_nll`). CFM comparison families carry `__nfe16`, `__nfe32`, `__nfe64` suffixes. Checkpoint path/size/mtime, decoder, projections, parent IDs, sample seed, K128 and index0 are frozen before inference. A separate PID and all sealed artifact identities are checked before any truth access. No SHA256 is used.

ONLINE reads only explicit `observation`, identity `parent`, optional `target_id` and, for an anchor-base model, y-only `base_mean/base_chol`. It never materializes q/x/chi/truth/noise/prediction label fields. Source-major packing preserves each source's Re32 followed by Im32. The synthetic audit NPZ deliberately contains forbidden object label arrays that would fail with allow_pickle=False if accessed; inference succeeds without reading them.

Inference streams one parent at a time, limiting effective network batches to<=16. Every parent gets128 draws from each declared model; CFM gets16/32/64 NFE variants using identical base seeds. Sample0 is frozen before truth. Each sealed per-parent NPZ contains `x_samples[128,66]`, original-coordinate `q_samples[128,66]`, `chi_samples[128,576]`, full complex empirical mean and per-pixel variance plus real/imag variances, phase-event probabilities, sample indices and observation-only features. Explicit Gaussian/spline/triangular models additionally save `log_q_x[128]`. CFM sample log densities are not saved as exact quantities. Failed parents remain in the gate; there is no fallback probability law or silent removal. Aggregate/paired uncertainty is not evaluated if required parent values are incomplete.

SCORE computes the raw128-sample material mean MSE/MAE, full-grid RMS material energy score using distinct sample pairs, full66 proper raw dx66 NLL for algebraic models,20 fixed projection CRPS values (16 real random +4 spatial), phase threshold.25/.5/.75 event Brier and reliability arrays, and50/80/90/95 marginal coverage. Coverage includes both randomized inclusive-rank handling of exact atoms/ties and ordinary empirical central intervals, separately for latent66, full real material1152, ROI real material512, and20 projections. It is never whole-image simultaneous coverage. Material realification uses concat(Re,Im)/24; interval widths in those realified coordinates carry that scaling.

Raw event probabilities can give infinite Bernoulli logscore. Raw zero-support counts and raw logscore are saved without clipping. A separately labelled finite Jeffreys predictive logscore uses(count+.5)/(K+1); it is not the exact sample-law probability. Variance is empirical ddof0, not a qualified posterior covariance. Pure original C2 decoder replay must agree with online decoded worlds and offline supplied truth material.

CFM NLL remains NOT_RUN unless each method's SPEC includes `density_audit: {parents:[...], nfes:[16,32,64], agreement_tolerance:0.01}` before freeze, with at most16 declared audit parents. Then the complete trace gives labelled **numeric continuous-flow NLL** at those NFEs, with32/64 convergence difference and the frozen tolerance. It is not the finite RK4 map determinant, exact MH density, full-test proper-NLL table entry, or qualification certificate. Pure inference never performs this trace audit. Shared CFM audits are cached across sample-NFE variants.

FULL Maxwell likelihood, Gamma256 residual law and data-dependent SBC are always `NOT_RUN_AWAITING_PARENT_PAID_PREDICTIONS` in this tool. The parent-owned physical checker can use sealed `q_samples[0]`/sample IDs and observation features to compute paid diagnostics separately. Cheap latent/material marginal ranks are clearly named as diagnostics, and do not activate those physical flags or establish posterior qualification.

`PIPELINE_AUDIT.json` records three synthetic parents, one update0 Gaussian and one update0 CFM checkpoint, four generated variants, no failures, no forbidden-label reads, separate online/scoring PIDs, same-PID rejection, and CFM/physical NOT_RUN flags. No campaign truth, optimizer update, or Maxwell call was used. These synthetic fixture scores make no scientific calibration claim.
