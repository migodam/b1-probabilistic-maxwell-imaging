# G archive software contract frozen before tests

2026-10-15. Zero Maxwell/optimization/training. This module accepts external full-validation receipts; it cannot manufacture physical validation or certify a backend. T/root owns the physical validator and accounting.

Canonical material is the complex 24×24 physical χ field, with fixed known background. Bounds are Reχ∈[-.5,4], Imχ∈[0,.5]. No latent, response distance, hole count, truth field, material MAE or future outcome is used to decide insertion or Top-K.

Primary fixed metric `G-multiscale-support-v1`:

- `dM = sqrt(mean((Re Δχ/4.5)^2 + (Im Δχ/.5)^2))`, δM=.03.
- For each candidate, amplitude is `abs(χ-known_background)/1.2`. Gaussian smoothing is σ=0,1,2 cells, constant-zero boundary and scipy truncate=4. Occupancy in each view is `clip((smoothed_amplitude-.20)/.10,0,1)`.
- `dS = sqrt(mean((occupancy_views_A-occupancy_views_B)^2))`. Structural distinction requires BOTH dS≥.08 and XOR≥6 original-grid cells under `abs(χ-background)/1.2>.25`.
- Pairwise separation is `dM≥.03 OR structural_distinct`. Active archive must satisfy this for every pair. No single-link component merging.

Every externally full-validated, solver-validated, legal, finite endpoint within the declared data threshold and with a nonempty receipt is appended to `success_ledger`, including duplicates. Endpoint fields/predictions are copied into immutable backing bytes. Invalid/uncertain endpoints are held separately in `frontier` and `failure_ledger`. `posterior_weight`, `represented_mass`, `unmodeled_mass` remain null.

An endpoint separated from every active member is inserted. Under default `replacement_mode='quality'`, a single nearby representative may be replaced only for lower full residual AND if the new representative remains separated from every other active point. Under `replacement_mode='first'`, the original representative is retained. These are the two declared archive ablation modes. If it overlaps several active representatives, none is deleted. The old success record remains. Discovery reward requires separation from **every prior success record**, so pruning/replacement/reinsertion cannot mint repeat reward; this conservative historical rule may reduce reward relative to a packing-only count.

`material_representatives` and `structural_representatives` are separate greedy subpackings of the current active set (stable candidate-ID order). Their cardinalities are order-dependent packings, not disconnected branches or posterior modes. Main active set uses the combined rule.

Top-K starts with an externally validated designated base if present, then fixed farthest coverage using `sqrt(dM²+dS²)`, enforcing combined separation from already-selected entries. Tie breaks by candidate ID. It returns fewer than K when necessary. Top-1 separately minimizes full residual over the complete success ledger. No truth-dependent reranking.

Optional `max_offers` caps the number of submitted software endpoints, including failures; it is not a physical-state/RHS budget and never resets or refunds an external ledger. Exhaustion retains all previously validated results. No pruning at exhaustion.

Posthoc digital topology uses 4-connected foreground/8-connected background. Component/closed-background sizes below six cells are reported separately; minimum size is an analysis scale, never a requested number of holes or archive cell label. The archive's XOR≥6 requirement does not imply six-cell holes or continuum topology certification.

Necessary tests are limited to ≤10 controlled wall seconds, CPU single thread, no physics calls. Receipt will contain actual process CPU/wall and be merged into G's independent cost ledger by root.
