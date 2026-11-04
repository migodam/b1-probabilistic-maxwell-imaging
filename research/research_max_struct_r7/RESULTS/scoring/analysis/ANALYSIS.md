# Independent R7 primary analysis

30 independent object realizations, six paired arms; acquisition strata are object attributes. Final images use all576 material variables. No forward/adjoint calls are made by scoring or analysis.

|Arm|Material successes /30|Same-grid material+noise /24|Median NRMSE|Median IoU|CPU total s|LU total|Complex RHS total|
|---|---:|---:|---:|---:|---:|---:|---:|
|cie_single|21/30|15/24|0.08722|0.8103|1543.60|54873|1154196|
|cie_dual|20/30|14/24|0.08584|0.8051|1539.15|53408|1137060|
|global_dual|20/30|14/24|0.08653|0.8141|1512.97|53923|1137684|
|a_material|21/30|14/24|0.08177|0.8266|2495.14|54384|1148820|
|a_plus|21/30|14/24|0.08233|0.8297|1711.34|52156|1103040|
|physical_lbfgs|10/30|9/24|0.3149|0.2186|838.18|47754|1009764|

## Paired object bootstrap95% confidence intervals

Percentile bootstrap10000 paired object resamples, fixed seed2026110421. Mean deltas are left minus right. Two prespecified contrasts and one separately labeled secondary contrast; no multiplicity-adjusted inference. Failed tasks remain failures in success denominators; continuous metrics state their finite-pair count.

|Comparison|Metric|Pairs|Mean delta|95% CI|
|---|---|---:|---:|---|
|a_material vs global_dual|material_nrmse|30|-0.043816|[-0.10088, -0.0056782]|
|a_material vs global_dual|material_image_success|30|0.033333|[-0.066667, 0.13333]|
|a_plus vs a_material|material_nrmse|30|0.00060392|[-0.00092329, 0.0028723]|
|a_plus vs a_material|material_image_success|30|0|[0, 0]|
|a_material vs cie_dual|material_nrmse|30|-0.038383|[-0.093229, -0.0033848]|
|a_material vs cie_dual|material_image_success|30|0.033333|[-0.066667, 0.13333]|

## Scope and cost matching

The24 same-grid objects permit the frozen independent CN99 noise discrepancy interpretation. For six36-to24 cross-grid objects, raw residual and noise-only pass are shown separately: physical model mismatch is unbudgeted, and sigma is not a complete data-error budget. Material recovery is always separate from data fit.

Equal-cost grading uses each object’s common budget=min(maximum saved cost across six arms), separately for process-prefix CPU, physical LU, and total complex RHS. It takes the last saved original-Phi incumbent at or below budget; no interpolation or GT-based iterate selection. Snapshots every ten corrector observation events plus sparse auxiliary stage endpoints, and float32 material storage, make this a checkpoint-resolution approximation. Auxiliary CPU gaps can be large; unused-budget distributions are reported. CPU includes auxiliary work; equal LU/RHS do not imply equal CPU. Final serialization after last budget snapshot is outside that last checkpoint; BUDGET_AUDIT.json quantifies the tail.

Pairwise A-versus-globalDual and A+-versus-A frontiers are primary cost comparisons; A-versus-sourcewiseDual is secondary. For each object and pair, the observed common window is [max(first saved cost),min(last saved cost)] on CPU/LU/RHS separately, sampled at11 fixed fractions. No stopped method is extrapolated. The all-six common-budget table is retained but can collapse when a method stops early at a poor local point, so it is not the sole cost comparison. Mean paired NRMSE differences and95% object-bootstrap intervals accompany exact per-object budgets and unused gaps. All30×6 images are shown in six fixed groups of five, without score-based selection. Continuous distributions, condition/family summaries and matched-budget material scores are retained. Acquisition conditions jointly vary contrast, noise, sources, aperture and generation grid; their descriptive differences do not isolate one factor. Bootstrap intervals across this fixed heterogeneous acquisition mixture do not establish population-wide universality.

Historical R6 results may be quoted separately as already-public development/old-confirmation evidence; they are not added to this30-object denominator and are not rerun.
