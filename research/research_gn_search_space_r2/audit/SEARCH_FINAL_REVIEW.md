# Search final technical review

Scope: independent source/receipt review of candidate-set semantics and transverse reference points. No experiment or scorer was rerun; only `src/score_search.py` truth-directory literals were updated and this review was added.

## Released truth path

The scorer now reads `research_z_alg_r4_diversity/truth_after_freeze/`, including its retained `mapping_PRIVATE.json` filename (`score_search.py:37,45`). All arrays/keys in the 24 observation truth files used here are exactly equal to the former archived files. Every mapping record's scoring fields (observation ID, shape, amplitude, background and noise) also agrees; only the mapping's top-level provenance/entropy-policy metadata differs. Existing scores therefore remain applicable and were not regenerated.

## Frozen execution record

- Closed/scored predictions: 36 frozen E2, 96 frozen E3 and 12 separately labeled tighter-stopping diagnostics; all 144 finite.
- Twelve clean E3 continuation predictions reuse the exact E2 full-regularization run: 132 unique successful physical runs.
- No successful task reached its 45 s / 500 LU / 30,000 complex-RHS cap. Eight corrector failures remain in receipts and measured costs. One native worker-death partial attempt has explicitly UNKNOWN costs in `results/e3/PROCESS_FAILURES.json`.
- Main predictions, tighter-stopping diagnostics and the separately charged stationarity check remain distinct. The review does not reinterpret early stopping as exact stationarity or a global optimum.

## Oracle candidate set and online ranking

**Every before-seed included by the formal oracle was physically evaluated and online-ranked.** `search.py:131` calls `evaluate(seed, ...)`; `evaluate` performs a full nonlinear forward solve, appends the scalar trace and updates the common-objective global incumbent (`search.py:60–67`). Corrector evaluations use the same path (`search.py:71–76`). This was verified for all 576 before-seed rows and all 576 saved after-material rows: their labels/objectives match evaluated trace entries. All 144 selected top1 objectives equal their trace minimum exactly. Eight of the 72 directional-search task top1 labels are before-seeds, which is permitted by the implemented common ranking.

**The formal best@K set is not identical to the complete evaluated trajectory archive.** For each directional-search task, `score_search.py:48–52` uses anchor + selected top1 + eight before-seeds + eight per-seed correction-best materials: 18 slots, including duplicates. The complete scalar trace has 22–74 evaluations per task. Other intermediate TRF materials are not retained/scored unless they become a per-seed best or global top1. Thus the oracle is a disclosed finite subset of physically evaluated/online-ranked points; it adds no unranked before-seeds, but it is not an oracle over every evaluated iterate and its slot count is not a distinct-root count. A per-seed “after” is the minimum-common-objective point observed during that correction, not necessarily its last iterate (`search.py:76,96,133`).

## Transverse reference points

The stored formulas are internally correct and use two different reference points for different quantities:

- Proposal displacement is `h = seed − anchor` (`search.py:128–130`). Across all saved rows the largest component discrepancy between `before − anchor` and `proposal_filtered` is 5.56e-17.
- Saved `correction_weak_norm` / `correction_strong_norm` project **after − seed**, using the fixed anchor-J right basis (`search.py:134–136`). Recalculation from arrays matches all 1,152 saved projection norms exactly.
- Scored `selected_change_*_squared` projects **top1 − anchor** (`score_search.py:57–61`), which is the correct selected displacement from the anchor.

The strong cutoff is the declared fixed-anchor, normalized-J rule `sigma(K) * 0.03 * sqrt(576) > epsilon`; epsilon is 1e-4 clean or 1e-3 noisy. The same basis is used for seed and correction decomposition. It is not recomputed along the nonlinear correction trajectory.

**Correction weak norm is not weak retention.** Retention must compare `||Pweak(after − anchor)|| / ||Pweak(seed − anchor)||`, because `after − anchor = (seed − anchor) + (after − seed)` and correction can cancel the seed. Neither the search receipt nor scorer currently stores this ratio; the existing correction plot labels its axes as correction components, accurately. The saved arrays permit the correct posthoc calculation without a new physical solve:

| Arm | Noncollapsed rows | Median weak retention |
|---|---:|---:|
| Haar unfiltered | 180 | 0.936671 |
| Green pullback | 180 | 0.324600 |
| Soft Jacobian | 180 | 0.106820 |

These review-only ratios exclude 12 collapsed rows per arm using weak seed norm <=1e-12; they are norm ratios rather than signed alignment, error recovery or basin certificates. For example, noisy digit4 soft task `e3_obs_e9c0c99d4be34e95_soft_jacobian`, proposal 3, has weak seed norm 1.6845264, weak correction norm 1.6742773, but weak displacement remaining from the anchor only 0.0277578 (retention 0.0164781). Calling the correction norm “retained weak displacement” would reverse the interpretation.
