# Frozen checkpoint reproduction

Status: **PASS**. Existing VAL256 only; CPU1 canonical. No new truth, model training, or Maxwell RHS.

| Action | CPU seed-mean MSE | Sealed MSE | Absolute difference |
|---|---:|---:|---:|
| physics | 0.0108894868145 | N/A | N/A |
| strongCNN | 0.00746135367081 | 0.00746135413646698 | 4.656612873077393e-10 |
| original_K32_frozen | 0.00714549794793 | 0.007145498879253864 | 9.313225746154785e-10 |
| deterministic_matched | 0.00703415414318 | 0.0070341541431844234 | 0.0 |
| distilled_w128_n32768 | 0.00716773793101 | 0.007167737931013107 | 0.0 |
| det_frozen | 0.00704105570912 | 0.0070410557091236115 | 0.0 |

All action/seed pairs use two fresh loads across all 256 VAL parents; exact equality and global torch RNG preservation are individually recorded. MSE reproduction gate is 1e-8; optimized K32 stored-moment equivalence gate remains 2e-6.

The original CPU-versus-MPS replay remains **FAIL**, at the original 3e-6 tolerance. Exact CPU repeatability is a separate claim and does not repair that failure.

Normalization, versions, allowed observation keys, K32 sampling contract, and source dependency paths are in CHECKPOINT_REPRODUCTION.json. BLAKE2b-512 checkpoint identities are in audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json. Compact predictions and per-parent losses are in audit/old_val/.
