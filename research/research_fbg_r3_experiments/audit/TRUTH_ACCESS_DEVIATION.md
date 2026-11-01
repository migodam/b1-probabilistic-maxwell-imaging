# Truth access deviation: eager post-E2 loading

Status: **REAL PROTOCOL DEVIATION; E1 IS NOT A FRESH BLIND CONFIRMATION SET**.

On 2026-11-01, the scoring worker received explicit parent authorization to score the already sealed E2 mechanism batch (`e2_mechanism_results.json`, 240 runs). Before opening any private file, it communicated its label-free scientific inference: the 33 legal conditional selected-channel roots did not fit all physical data; algebraic data passes came from empty-root anchor fallback, and the matched seed-count control did not support a basin advantage. That inference was not changed after GT scoring.

The initial scorer implementation read `mapping.json` and eagerly loaded the material/mask NPZ arrays for **all 12 objects** into process memory, rather than limiting truth loading to the **3 objects actually represented by the 4 E2 observations**. E1 predictions were not yet closed. The worker discovered this immediately after the first E2 scoring pass and disclosed it to the parent. No E1 metrics, truth values, support arrays, or reconstruction figures were computed or printed by that pass. Parent, N1, and N2 did not receive E1 truth values, supports, or scores. Nonetheless, opening the arrays is actual access and violates the intended batch-specific holdout boundary; absence of printed values does not reverse the access.

The E1 algorithm/configuration and data-free design decisions had been frozen before E2 scoring. The existing E1 batch is therefore classified as a **controlled mechanism set with disclosed scorer exposure**, not a fresh blind confirmation set. No E1 algorithm, hyperparameter, support selector, proposal rule, or acceptance threshold may be tuned from this exposure. This record does not reclassify the E1 object family as previously unseen or repair the holdout retrospectively.

The scorer was corrected to build its truth-load list from only the case IDs in the closed result batch. The explicit parent authorization and closed-artifact status guards remain before all private access. It performs no new Maxwell solve. No further E1 truth access or scoring is authorized until the parent closes E1 predictions. E2 will not be repeated merely to obtain a cleaner access story; its first scoring outputs remain part of the evidence.

A future fresh confirmation would require newly generated objects after algorithm freeze, with every truth array deferred until the complete prediction batch closes. Such a study is **NOT RUN** here and is not implied by this correction.

Affected first-pass outputs: `results/e2_scores.json`, `results/e2_scores.csv`, and E2 rows in `results/COST_LEDGER.csv`. Material metrics were computed only for the 3 E2 objects. No additional scoring physics, error-floor checks, or new optimization was run.
