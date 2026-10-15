# Research G: executable next-round task

Date: 2026-10-15. This is a prospective protocol, not a report of executed experiments.

## Objective and authorization

Build a topology-class-free, physics-corrected search for a small set of distinct material fields explaining the same electromagnetic observation. Measure discovery AFTER correction and independent full-forward validation. Do not train PPO, GRPO, diffusion, or a larger MDN in this round.

When the user explicitly starts this prompt, execute E0 and the bounded E1 development pilot below. E2-E4 are protocol-only and remain NOT_RUN. Do not reuse historical remaining budgets. Do not modify Research T, Research B, historical results, their ledgers, or existing candidate schemas. Work in a new `research/research_g_20261015/` directory or an isolated new repository.

First read the sources, write a short implementation plan, and execute within the specified gates. Do not replace scientific work with repeated approval requests, excessive hashing, large generic audits, or publication housekeeping. No GPU training. Missing indispensable inputs must be reported as BLOCKED, not silently synthesized.

## 1. Read and freeze provenance

Research A repository: `migodam/b1-probabilistic-maxwell-imaging`
Commit: `f1b674eb4893e39126701ed48f3d3dcef3098c3c`
Read:
- `research/research_a_20261015/RESEARCH_REPORT_ZH.md`
- `research/research_a_20261015/results/SUMMARY_E1.json`
- `research/research_a_20261015/results/STRATIFICATION_ATTRIBUTION.json`
- `research/research_a_20261015/results/PAIRWISE_AMBIGUITY.json`
- `research/research_a_20261015/results/E2/E2_FAILURE_ANALYSIS.md`
- `research/research_a_20261015/results/E2/RESULT_SUMMARY.json`
- actual solver adapter, corrector, raw proposal histories and append-only ledgers referenced by these files
- actual original and execution-overlay MaterialHypothesisPackage schemas; use a separate G sidecar, not incompatible additional fields.

A30 repository: `migodam/a30-maxwell-observable-material-manifold`
Commit: `717d3d26a290248bd2b9aeb55f3d03694b407b84`
Read `A30_E0_E1_DECISION.md`, actual configuration and capacity results. The 35.69% Sparse-64 improvement is an OFFLINE truth-fit capacity result, not online reconstruction or full-data feasibility.

C1 repository: `migodam/c1-physics-state-generative-imaging`
Read `START_HERE.md`, `p1/C1_P1_RESEARCH_REPORT_ZH.md`, `p2/C1_P2_RESEARCH_REPORT_ZH.md`, accepted continuous comparison, and actual frozen proposal implementation. Record the revision actually read. P2 models are radial four-dimensional proposals and must not be represented as trained A30 image generators.

Read Research A's upstream-history summary for A17 provenance. Do not invent an A17 Git revision. Its local frozen-GN results do not certify nonlinear material endpoints.

## 2. Immutable scientific boundaries

1. Preserve the A30 24x24 scalar TM/VIE backend, eight sources and 32 receivers unless a separate experiment is explicitly registered. This is NOT general three-dimensional vector Maxwell.
2. Preserve complex128/float64, source ordering, real/imaginary packing, quadrature weights, material bounds, background and noise convention. For the inherited A30 setting, Re/Im noise variances are sigma^2/2 and whitening uses sqrt(2)/sigma. Gain remains 1. Do not import the radial Mie noise/gain convention.
3. Preserve the inherited E2 acceptance threshold sqrt(chi-square(512,0.95)) = 23.785449004 where that exact target identity applies. Do not silently substitute the radial threshold or retune noise. Check the actual target identity before reuse.
4. Every published candidate must pass the frozen full-discrete forward rule, material legality, shared material across sources, and available solver/constitutive diagnostics. Store unavailable diagnostics as null with a reason, never as zero.
5. Formal candidate counts exclude invalid and unverified states. Keep invalid stepping stones only in a separate frontier/failure ledger.
6. No test truth, hole count, true boundary, true material error, oracle representation coefficients, or future outcomes enter proposals, correction, archive selection, rewards, feature normalization or model selection.
7. Offline scoring runs in a separate process which cannot be imported by the online search package. Freeze and serialize online outputs before starting it. Never rewrite outputs after scoring.
8. Proposal density and search frequency are NOT posterior probability. Store posterior_weight, represented_mass and unmodeled_mass as null/unquantified.
9. A true material difference with nearly identical receiver response is NOT a duplicate. Deduplicate decoded physical fields, not latent coordinates or observation distance.
10. Do not name material-field clusters disconnected solution branches or posterior modes without additional evidence.

## 3. E0: zero-new-physics audit

Do not instantiate a solver to recompute quantities already stored. Use cached predictions, canonical fields, trajectories, timestamps, actions and recorded costs.

### 3.1 Required outputs

Produce `E0_DIAGNOSTICS.json`, `E0_REPORT_ZH.md`, and a compact figure directory. Report missing arrays as missing, not reconstructed guesses.

- Separate raw proposal diversity, corrected endpoint diversity, full-valid diversity and archive-retained diversity.
- Separate material-distinct counts from geometry-distinct counts.
- Compute duplicates under the original metric AND a predeclared multiscale structural metric. Do not overwrite historical counts.
- Repeat threshold diagnostics at 0.5x, 1x and 2x the new primary threshold. These are sensitivity analyses, not a route to select the most favorable threshold.
- Inspect first-member deduplication: how many within-neighborhood quality improvements were discarded? Use only truth-free quality for any proposed future replacement rule; GT-based reranking is diagnostic only.
- Build a parent/proposal-to-endpoint transition table. Quantify correction collapse into the same material neighborhood.
- Report invalid->valid, valid->invalid, structural retention, bound saturation, correction magnitude, trust-region exits and failures by operator and spatial scale.
- For A30, distinguish invalid base, missing raw hole proposal, erased hole, retained-but-invalid hole, and valid-but-not-archived hole.
- Preserve the actual finding that protected edit cells in Research A E2 were unchanged; do not assert that this corrector filled the holes.
- Preserve NA for MAE when no valid candidate exists. Report failure rate and all-case Coverage@K separately; an empty valid set has coverage zero.

### 3.2 Does history contain actionable information?

Predict future corrected validity and corrected structural novelty using grouped cross-validation by material parent. Compare:
A. observation summary + current material + action + remaining budget;
B. A plus prefix-only success/failure summaries;
C. B with history shuffled between matched observations;
D. history-only ablation.

Prefer logistic/ridge models or a small tree ensemble; use no new physics and no large neural model. Report calibration/log loss as well as discrimination. Any normalization is fitted on training folds only. Avoid treating noise repeats or search steps as independent objects.

Off-policy performance claims require logged propensities and overlap. Without them, report predictive association, not an identified advantage of a new adaptive policy. Reordering a cached candidate pool diagnoses archive selection only; it cannot recreate the counterfactual proposals an adaptive policy would have made.

### 3.3 Representation and topology audit

Use all legal decoded material fields and the inherited TRAIN-only Sparse/PCA representations. Truth-fit reconstruction is permitted ONLY inside the isolated diagnostic scorer, never as online labels.

Check whether a hole can appear in a continuous coefficient sweep without fixing its count. Show that fixed support, smoothness and decoder restrictions can still exclude certain holes even without an explicit topology class.

Freeze spatial scales in physical/grid units, foreground convention, boundary handling and complementary digital connectivities. For complex/negative contrast, do not define support solely by positive real contrast. Persistent homology is optional: pair amplitude persistence with spatial smoothing/minimum physical feature size; a high-contrast one-cell artifact can be persistent in amplitude.

## 4. Common material domain and proposal interface

Use a unified bounded real/imaginary material grid as the authoritative physical field. Proposals may be parameterized by local smooth masks, multiscale Haar/DCT coefficients or a continuous implicit field. Do not assign one neural expert per hole count or one archive cell per named topology.

Expose a common action library:
- global correlated material-field restart;
- local material removal OR addition with continuous center, scale and amplitude;
- multiscale distributed material change;
- compensated weak-direction proposal, only with valid tangent/box handling;
- topological-sensitivity-guided local proposal where a correct implementation is available.

The operator names are not object classes. The same action library and legal material support must be available to all selection policies. At least 20% of proposals in adaptive/archive methods remain global restarts, chosen before surrogate screening. This is a robustness choice, not a proven optimal mixture.

A local action should initially use approximately 8-32 continuous parameters, rather than an RL action per material pixel. Correction can use a richer shared representation. A fixed low-dimensional chart is not allowed to masquerade as the full material domain.

## 5. Corrector and validator

Before comparing policies, freeze one common corrector. Prefer a feasibility-restoration objective anchored to the RAW proposal, not always to the base reconstruction. Use the same implementation, stopping rule, tolerances, representation expansion and per-trajectory caps for every policy.

Maintain two modes during development only:
- inherited corrector, for regression;
- proposal-retentive corrector, for a separate small factorial diagnostic.

Do not credit a search policy for an improvement caused solely by changing the corrector. The primary E1 policy comparison uses only one frozen corrector.

The proposal-retentive formulation minimizes change from the proposal subject to the discrepancy threshold, material bounds and a local trust region. A reversible edit-retention constraint or penalty may be relaxed along a predeclared schedule. It must not impose a target hole count. Stop reducing residual once feasibility and solver accuracy are adequate; further movement must have an explicit purpose.

Record whether a valid endpoint retained the intended edit. If correction erases an edit, retain the outcome in the ledger but do not count the raw edit as a discovered structure.

If a proxy is used, record its target and error estimates. Independent full validation is mandatory for every public candidate. Reserve a predeclared 10% of validation opportunities for randomized audits of screened-out or uncertain proposals, within the same total budget. Log selection probabilities. These audits estimate missed opportunities; they do not prove global absence of false negatives.

## 6. Unstructured archive

Maintain:
- `success_ledger`: immutable full-validated endpoints, including duplicates;
- `active_archive`: a compact epsilon-separated set of materially/structurally meaningful representatives;
- `frontier`: invalid or uncertain promising stepping stones;
- `failure_ledger`: raw action, endpoint, correction trace, reason and cost.

Use a volume-weighted complex material distance, multiscale material/support/boundary differences and, optionally, scale-robust persistent topology. Freeze descriptor versions before confirmation. Learned embeddings may be auxiliary only in the first round.

Use post-correction novelty. Do not reward image noise, decoder reparameterization, infinitesimal cavities or unrestricted contrast saturation as new structure. Distinct complex contrast can still count as material diversity even when it does not count as geometry diversity.

Do not use single-link clustering that chains distant endpoints into one family. Use explicit pairwise separation/packing checks. If a local representative is replaced, retain the old point in the immutable ledger. Re-inserting an old world must not earn another discovery reward.

Select Top-K without truth: preserve a legal base where required, then use a frozen material/geometry coverage rule. Report this separately from the minimum-full-residual Top-1 choice. A collection with fewer than K valid distinct entries remains short; never pad it with invalid outputs.

## 7. E1 development pilot: six primary arms

Use eight distinct exposed development parents, with one acquisition condition per parent assigned before reading truth. Balance inherited conditions by ID/order, not by favorable outcomes. These are development diagnostics, not fresh confirmation.

Primary arms:
1. IID prior/correlated-field multistart + common corrector.
2. Static field stratification across continuous proposal scale/volume coordinates + common corrector. Do not claim that the old radial eight-interval experiment is directly portable to arbitrary 2D fields. Any explicitly topology-aware control must be labeled biased.
3. Unstructured QD with fixed emitter mixture and post-correction novelty.
4. Lightweight adaptive CEM or a contextual bandit plus the SAME archive/action library. Start diagonal covariance; no full-field covariance matrix.
5. Direct conditional field proposal + common corrector. Start with TRAIN-only conditional retrieval/empirical mixture plus multiscale residual noise; no new large NN. Label this a nonparametric generative baseline, not an MDN transferred from Mie.
6. A legitimate topological-sensitivity + GN proposal baseline, with the SAME common correction/validation budget. If the correct derivative for this backend cannot be implemented, mark NOT_IMPLEMENTED and do not rename a heuristic as topological sensitivity. Then the six-arm comparison is incomplete.

Numerical deflation, differential evolution, population repulsion and beam search are protocol/cheap-development comparators, not additional unbounded full-physics campaigns. Deflating stationary points is not equivalent to enumerating a noisy feasible set.

### Hard E1 caps

Per parent and policy: at most 32 new material states/factorizations, 1,024 complex LU RHS, and 12 full-validated endpoint attempts. Endpoint solves are INCLUDED in the state and RHS counts, not added on top. Stop when ANY cap binds.

Common base per parent: at most 16 new states and 512 RHS, shared across all policies and recorded once. The base can fail; preserve that result. The total E1 cap for eight parents/six policies is 1,664 states and 53,248 RHS, excluding only explicitly available historical predictions. Existing observed measurements require no new truth solve. If new measurement generation is necessary, its calls consume the global reserve; do not hide them.

The RHS/state relation is NOT assumed to be fixed. These are independent caps. Record wall time, LU factorizations, primal/adjoint/JVP/VJP solves, failed attempts, cache hits and proxy calls. Shared setup costs appear once globally and in an amortization table, not once per candidate package.

Randomize method order for wall-time comparisons. Reuse only exactly matching immutable cached states. Cache replay is not a new physical speed result.

## 8. Frozen endpoints and gates

For every arm and budget prefix report:
- legal Top-1 and its selection rule;
- Best-of-2/4/8 real, imaginary and complex material MAE, offline only;
- all-case truth coverage at fixed material AND geometry thresholds;
- material-distinct and structure-distinct full-valid counts;
- robust hole/support/boundary recovery with spatial-scale sensitivity;
- discovery per state, per RHS and per wall time;
- proxy false positive rate and audited false negative rate;
- raw proposal success, correction success, correction collapse, retention;
- support/bound saturation, out-of-training-morphology status;
- complete failures, empty archives, solver-cap failures and archive-shortfall rates.

Do not combine NA MAEs into an apparently successful all-parent mean. Report conditional-on-valid MAE together with all-case coverage and failure rate. Do not replace missing candidates by the base unless the base itself passes validation and is explicitly included by the selection rule.

Use paired material-parent resampling for uncertainty; seeds/noise repeats are not independent parents. Eight exposed parents do not establish publication-grade superiority.

E1 is a mechanism GO only if a hierarchy-free method repeatedly discovers a full-valid, spatially meaningful new structure that is not merely a pre-existing uncorrected proposal, parameter duplicate or archive-threshold artifact. Lower loss alone is NO_GO.

For a later fresh confirmation, pre-register a joint superiority gate: at least 20% more valid structure discoveries at matched physics cost, at least 10% lower Best-of-4 complex material MAE, no meaningful Top-1/legality deterioration, and positive paired uncertainty evidence against the strongest completed baseline. These percentages are engineering decision margins, not theoretical constants. Freeze them before confirmation; do not tune after seeing results.

If the full-vector metric shows only contrast diversification without geometric changes, report that narrower outcome honestly.

## 9. E2 protocol only: hole birth/death and correction attribution

Write a plan for eight fixed witness cases x four modes x 16-state/512-RHS cap = 512 states/16,384 RHS. Do not execute in this authorization.

Modes: raw/no correction; inherited correction; proposal-retentive correction; richer/full-field reference restoration. Count any costs for obtaining witness seeds inside the cap or explicitly outside in a separately authorized setup ledger.

Trace from both solid-like and hole-bearing starting fields when such legal seeds exist. Otherwise state NO_VALID_SEED. A truth-seeded restoration is allowed only as a clearly isolated oracle diagnostic and must never receive discovery credit.

Use continuous local additions/removals and distributed field changes, without target hole counts. Track:
raw topological birth -> retained birth -> full-valid birth -> archived birth.

Vary only one mechanism at a time. If a hole survives but residual stays above threshold, the failure is not mode coverage of the proposal. If full-field restoration succeeds but low-rank correction fails, prioritize correction representation. If valid holes exist but are discarded, fix archive selection. If all tested finite paths fail, report a best-found failure, not uniqueness or impossibility.

## 10. Later stages: do not run

Fresh confirmation suggestion: 24 new parents, two finalists, same 32-state/1024-RHS policy caps, common bases 16-state/512-RHS, plus 24 eight-source measurement states. Maximum 1,944 states and 61,632 RHS. This is a suggested separate authorization, not an entitlement to execute now.

E3 RL: only after E1/E2 success, predictive history signal and an economically plausible deployment volume. Compare identical low-dimensional actions under PPO, original-GRPO-style individual utility, group-aware credit, and no-RL CEM. A possible initial training-physics screen is 4,096 total states, not an unlimited training budget. Include proxy-bank creation, all-fail groups and confirmation calls. Retain NO_GO if this cannot plausibly amortize.

E4: general, nonconcentric three-dimensional vector Maxwell and active sensing only after nontrivial 2D success. Radial Mie is not a substitute. No information gain claim from uncalibrated candidate weights.

## 11. Deliverables

Write:
- `RESEARCH_G_REPORT_ZH.md`
- `E0_DIAGNOSTICS.json`, `E0_REPORT_ZH.md`
- `E1_PROTOCOL_FROZEN.json`, `E1_RESULTS.json`, `E1_DECISION.md`
- `CORRECTOR_ATTRIBUTION.md`
- `SEARCH_HISTORY_SIGNAL.md`
- `GSearchEpisode.schema.json` as a versioned sidecar, with actual schema tests
- `COST_LEDGER.jsonl`, `FAILURE_LEDGER.jsonl`
- canonical before/after complex fields and saved full predictions
- `NEXT_STAGE_PROTOCOL.md` with E2/E3/E4 explicitly NOT_RUN
- figures showing every pilot parent, including failures, at shared color scales

End with one of: GO_FEASIBLE_DISCOVERY, GO_BASE_RESTORATION_ONLY, NO_GO_SEARCH, or BLOCKED_INPUTS. State exactly what is supported, what failed, and which resource cap stopped computation. Do not claim a strong IEEE TAP algorithm from interface tests, a decreasing residual, a learned reward curve, a few radial images, or a larger raw candidate count alone.
