# Strong TRAIN/validation context packing: zero physics

New standalone tool: tools/pack_r4_context_training.py. It does not generate contexts, train models, edit existing model code, modify original joint data, or alter the shared ledger.

Availability check (no output written):

`<python> tools/pack_r4_context_training.py --train-count 1024 --validation-count 256 --out research/p3l_r4_conditional_transport/inputs/strong_context_1024_256 --check-only`

After all physical jobs complete, parent can schedule the same command without --check-only under its analysis ledger. Optional `--base-chol-scale 0.25` exports a fixed fullrank0.25*I66 base factor, not a posterior covariance. Omit it to let CFM use an explicitly declared positive --anchor-scale. Defaults read stage/results/contexts/{train,validation}/R2_group_NNNNN.npz and stage/inputs/joint/{train32768,validation}.npz. Counts512/128,1024/256 etc select deterministic firstN identities; maxima remain2048/512. Counts do not trigger new generation.

Each requested context requires its finished JSON sidecar: correct parent/target, GT_read=false and native_cap96/probit_cap64. Missing/incomplete contexts report WAITING with no output; --check-only exits0, actual packing exits3 until ready. No polling/background job is started. Joint parent uniqueness, correct fold prefix, requested roster, TRAIN/validation disjointness and context-parent identity are checked. Context source byte/mtime snapshots are retained without hashes.

Packing reads only original online context fields observation,q0,chi0,J,prediction,residual,phi,x,H plus parent metadata; no GT file/key is consulted there. Original joint label x/material is joined by exact matching parent in the same fold. Context observation must match the original joint source-major y512 exactly after float32 packing; residual and context Gaussianization must agree. This prevents quiet fold swaps or pairing another observation with a label.

Each output directory contains:

- anchor_train.npz / anchor_validation.npz: compact R2-style arrays, chi_GT complex64, q_GT reconstructed from original float32 Gaussianized label, x_GT original float32 label. The x key remains online Gaussianized x0. Existing read_anchor(supervised=True/False) accepts this schema; online readers still materialize only their allowlist.
- density_train.npz / density_validation.npz: x66 truth label, original raw y512, parent, base_mean=contextx0, optional positive fullrank base_chol, context_x/context_q0/context_chi0/context_residual/context_J and original material/s labels. Truth and online coordinates have explicit names. No y augmentation is silently applied.
- online_train.npz / online_validation.npz: truth-free companion containing original online arrays only.
- PACK_COMPLETE.json: exact schemas, dtypes, counts, identities, joint observation checks, resource profile, no-training/no-physics receipt and budget-shift explanation.

Real source observation/material context arrays preserve complex128; J,H,q0,x and other real context arrays usefloat32, with label materialcomplex64 and original float32x. Atomic compressed NPZ writes and final directory rename prevent consumers from seeing a partial successful pack. Existing output directories are never overwritten. Consume only directories with PACK_COMPLETE.json.

Historical512/128 contexts were native48+probit48. These new contexts are native96+probit64 on the same underlying joint parents; they are not new independent label parents and must not be pooled as one unchanged context estimator.

Current density trainer compatibility limit: CFM --base anchor consumes base_mean/base_chol, but spline/Gaussian conditioner still uses raw y512 and ignores context_* fields. Exporting these fields does not by itself create an anchor-conditioned spline. Conditioner augmentation remains a separate parent/density-worker implementation decision.

Synthetic schema test CONTEXT_PACK_SCHEMA_ONLY/SCHEMA_ONLY_RESULT.json passed supervised anchor loading, full66 density loading, fullrank base diagonal and absence of GT keys from online companion. This is not a generated strong-context experiment. At the availability check in this worker turn, new requested1024/256 contexts were0/0 ready; no real pack or training was launched.

Static ledger review: physics/training each have their own exclusive campaign lock, and budget mutations share budget.lock. Thus one physical and one training job may coexist, consistent with current protocol. The wrapper's error text saying one campaign worker is broader than its actual separate-lock behavior; no demonstrated counter/reservation race was found. Analysis jobs do not acquire either worker lock, so parent remains responsible for scheduling large packing alongside memory-heavy work.
