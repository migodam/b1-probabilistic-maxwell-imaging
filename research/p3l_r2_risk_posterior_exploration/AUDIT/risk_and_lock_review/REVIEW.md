# Supervised risk and lock design review

No physics, training, locked labels, locked estimates, or module edits. Review scope is specified source and completed train/validation summaries.

Primary actionable lock issue: matched176 is currently discovered from live file presence and reread after truth access. Freeze inclusion and score saved inference estimates.

## L01 — P1_BEFORE_TRUTH: Freeze matched176 inclusion and score the saved truth-free estimate

evaluate_p3r2_locked.py101-102 determines method inclusion by any live MATCHED file;123-125 reads raw MATCHED after truth opened at111. infer_p3r2_locked.py234-237 already saves matched_point176 when present.

Require an explicit manifest include_matched flag and frozen all-parent completion receipts before label access. Score matched_point176 from each saved estimate archive. If absent, apply frozen no-exclusion fallback with counts; do not discover methods from runtime file presence.

Does not directly alter fused-vs-probit primary statistic, but current optional matched comparison has weaker isolation and unplanned-subset behavior.

## L02 — P2: Label and count baseline fallbacks separately

infer_p3r2_locked.py178-186 substitutes cache/native/prior if context parsing fails and initializes anchor/probit to that fallback; evaluator174 always calls baseline strong native96+probit64. Per-parent failures are retained but no summary counts/qualification for point baselines.

Before truth access validate the strong baseline estimate contract or explicitly freeze fallbacks. In results summarize fallback counts for anchor/probit/native and each point method, qualify the strong-baseline label if substitutions occurred. Keep every parent.

Potential labeling/fairness failure; no observed locked failure asserted, no locked outputs read.

## S01 — CLAIM_SCOPE: Geometry CNN increment changes inputs and parameter count

residual_risk.py150-158/190-196: local_anchor has4 channels and5266 parameters; local_geometry has8 channels and5842 (+576,+10.94%). Same16-channel/3block convolution family, different contextual inputs.

Describe as geometry-feature-augmented local CNN, not a parameter-matched geometry-only ablation. State added weak/strong sensitivity and residual-backprojection channels. Existing results support combined feature/capacity comparison only; no required retraining to report this accurately.

Mean validation MSE0.0177651 vs0.0171738 (~3.33% extra reduction); selected min/max seeds are not independent test evidence.

## S02 — CLAIM_SCOPE: Most validation gain already appears with small anchor-only local filters

Anchor MSE0.02297525; local_anchor mean0.01776506; local_geometry mean0.01717383. physics_cnn322992 parameters mean0.02297066. Local models5266/5842 parameters; both selected400 updates and hit100epoch limit; large physicsCNN selects8-12updates then validation patience.

Attribute comparison to architecture/context/optimization package. Report anchor-only local model alongside geometry model and training stopping differences. Do not call the larger physicsCNN a capacity-matched control or infer geometry is the sole cause of gain.

All are supervised validation-selected point estimators, not posterior mean/calibration demonstrations.

## S03 — CLAIM_SCOPE: All-coordinate decoder control has different correction radius

residual_risk.py204-208: all_decoder aggregate norm<3; weak_strong204ff uses16 coefficients each bounded3 plus12 bounded0.75, max combined norm sqrt(150.75)=12.278. Head dimensions66 vs30; parameters185154 vs180510.

Call all_decoder a bounded all-coordinate inductive-bias control, not a matched full-span replacement proving low-rank geometry superiority. Record both radius and parameter differences. random_decoder is cleaner: same context, parameter count180510 and coefficient radii, fixed random frame.

Existing random control only tests output correction frame; it retains geometry inputs, so it is not no-geometry control. One frozen random frame supports a conditional comparison, not all random subspaces.

## S04 — REPORTING: Validation reuse and context-budget shift constrain claims

select_p3r2_fusion.py loads validation labels and selects family/projection/alpha; residual selection uses the same128 validation parents over models/seeds/epochs. Fusion metadata explicitly states trained48+48 contexts and future96+64 shift.

Keep validation risk descriptive and report all searched alternatives. Freeze primary fusion before new locked observations; no locked reselection. Label the anchor-optimization covariate shift and use prespecified locked baseline96+64. Density scores remain fixed10 marginal diagnostics, separate from material-risk estimators.

No independent validation significance/GO claim from these summaries; fresh locked evaluation remains necessary.

## Existing safeguards

- Inference opens a strict allowlist of online fields and never opens offline_only labels; material methods consume observations or observation/anchor/J only.
- Freeze verifier checks source revision, tracked differences, untracked Python files, artifact path/size/mtime registrations, artifact time before freeze, and256 unique fresh parent inputs. This is a concrete ordinary mutation guard; no hash checks requested.
- Evaluator requires inference COMPLETE with exact parents/revision/time, completed generator exact locked cohort, context closure jobs, all256 files, and96/64 truth-free sidecars before its first label archive read.
- Material/density exceptions use explicit anchor/unit-prior fallbacks and preserve failure strings; all256 parents remain in final metrics. Source failure prevents projected/fused branch from silently using a changed learned method.
- Bootstrap resamples paired independent parents with shared indices;2planned co-primary comparisons use97.5% two-sided intervals. Other method comparisons remain labeled secondary.
- Point-estimator and density claims are separated: residual checkpoints say calibrated_density NOT_CLAIMED; fusion says no density/posterior claim and no legal-single-decoder-material claim; evaluation labels fixed10 NLL/coverage without full posterior calibration.
- Random decoder preserves context, width, coefficient radii and total parameters relative to weak_strong; u2-only fixes z exactly and preserves property-odds map. Controls are now available for all3 seeds.

## Completed summaries

| Model | Parameters | Mean validation ROI MSE | Selected updates |
|---|---:|---:|---|
| anchor | 0 | 0.022975250 | [0] |
| physics_cnn | 322992 | 0.022970656 | [12, 12, 8] |
| weak_strong_decoder | 180510 | 0.021965917 | [28, 40, 24] |
| local_anchor_cnn | 5266 | 0.017765063 | [400, 400, 400] |
| local_geometry_cnn | 5842 | 0.017173825 | [400, 400, 400] |
| all_decoder | 185154 | 0.022689854 | [8, 8, 12] |
| random_decoder | 180510 | 0.022232314 | [32, 32, 28] |
| u2_decoder | 176898 | 0.022939599 | [20, 20, 28] |

Validation comparisons are descriptive after model/seed/epoch/fusion selection. No fresh locked result, posterior claim or novelty decision is asserted. Source/artifact changes made by parent after this snapshot require their normal final freeze check.
