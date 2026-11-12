# Per-run metadata contract wrappers

New standalone script: `src/build_run_contracts.py`. It creates new `RUN_CONTRACTS.json` / `RUN_CONTRACTS.csv` wrappers around frozen metadata; it never edits the original receipt, manifest, seal or result arrays. It performs no physics, factorization, candidate generation or optimization and never opens GT/labels/analysis quality tables.

```text
python -m research.research_max_struct_r13.src.build_run_contracts \
  --root research/research_max_struct_r13 \
  --manifest research/research_max_struct_r13/configs/CAMPAIGN_1.json \
  --manifest research/research_max_struct_r13/configs/CAMPAIGN_2.json \
  --output research/research_max_struct_r13/audit/NEW_RUN_CONTRACTS_DIRECTORY
```

Output must be a new directory outside RESULTS, queues, private sealed truth and seal directories. Existing files are never overwritten. The final study wrapper should be rerun into a new directory after remaining submissions/queues have completed.

## Availability and null policy

Every expected task is retained from the actual manifest task list. Only a locally available explicit terminal queue (`COMPLETED`, `COMPLETED_WITH_ERRORS`, or explicit completed recovery, with no running/pending work) permits reading execution receipts/result metadata. Missing/nonterminal campaigns retain expected rows with null execution fields and a reason; they are not treated as completed zero-cost runs. Even within a terminal queue, missing/unreadable receipts remain null rows.

JSON uses null and CSV writes the literal `null`. Unknown frequency units, missing executed-source evidence, unretained stage settings, optimizer states, physical defects and costs are not manufactured from current code or replaced by zero. Project-declared conventions and manifest-declared contracts are labeled separately from actual frozen execution metadata.

## Source provenance

| Campaign | Shipped/executed source evidence | Algorithm/config reference |
|---|---|---|
| C2 core48, mechanism16, physics_v1 | `audit/CONFIRMATION_FINAL_BUNDLE.json`: shipping `689850d725ddcbd6b65152595ac0ddf2b6eff1bd` | algorithm freeze `3a91b35605634a18a731441693a46152eeb3363e` |
| C3 online controls16 | parent session explicitly declared shipping prefix `4bd`; full commit is not invented | inherited manifest algorithm/config reference `3a91…` remains separately recorded |
| local reference16 | terminal `queues/.../SOURCE_SNAPSHOT.json`: `a1317c738` | local manifest reference, not automatically an i9 shipping bundle |
| local TRF16 | terminal source snapshot `2cd5516cb` | the ready/container commit is not substituted for optimizer source |
| Older/other campaigns | source snapshot when available; otherwise actual execution source null | a manifest source declaration alone is retained as a declaration, not promoted to proved shipped code |

`source_commit`, `shipping_commit`, `algorithm_freeze_commit` and detailed provenance are separate. Current HEAD is never substituted. The C3 prefix is explicitly session-declared, not a cryptographic content attestation. Every per-run and top-level `code_hash` / `config_hash` is null with exact reason `userforbidsSHAchecks`; no SHA checks are performed.

### Physics dependency-recovery composition amendment

`physics_recovery_v1` and composed `physics_v1` additionally reference `PHYSICS_DEPENDENCY_RECOVERY_PLAN.json`, `PHYSICS_C4_PREPARATION_RECEIPT.json` and `PHYSICS_RECOVERY_COMPOSITION.json` when those files exist. The immutable numerical base remains C2 shipping689850… / algorithm3a91…; `runtime_physics_dependency_C4.zip` supplies the unchanged R1 vector_backend/a9_green dependencies. It is not a new numerical algorithm or prior. Documentation/recovery/portability commit `a611` is explicitly excluded as a numerical-source claim.

For a terminal composed queue, `job.attempt_source_campaign` labels each row's actual C2 `physics_v1` or C4 `physics_recovery_v1` attempt. Missing labels remain unknown rather than guessing from method name; missing composition evidence does not assert a completed future composition. The wrapper preserves mixed-attempt provenance instead of calling every composed row an original C2 execution. No pending row gains executed-source/attempt fields.

C3 keeps shipping prefix4bd unless an exact compatible longer value is present in registered manifest/source-snapshot metadata. Such a value can resolve shipping provenance without substituting it for inherited algorithm3a91 or using current HEAD. No Git/content hash check is introduced.

## Public metadata access

G/E/S/m_init/coordinate arrays are read **only as NPZ/NPY headers** for shape/dtype, never as numeric payloads. Y is not loaded. Small whitelisted public metadata covers mesh shape, components, explicit frequency/unit fields, k, bounds, width/spacing, material/time convention and noise qualification. Hz conversion occurs only for explicitly Hz/GHz-labeled metadata; otherwise raw frequency/k may remain present while frequency_Hz is null with a reason. The wrapper does not infer physical Hz from normalized k.

The declared `exp(-i omega t)` convention and real susceptibility `chi=epsilon_r-1` fall back only to the project's explicit THEORY_DERIVATIONS opening declaration, with that source labeled. Numeric initializer summaries remain null: the wrapper records task declaration and m_init header, not material-array values or NN quality.

## Wrapped fields and limits

Each run contains identity/data path, source provenance, material convention/class, mesh/source/receiver/frequency metadata, box/prior/init, actual current-space rank/DOF versus stored ambient scalars, baseline family/qualification, retained stage/stop metadata, selection/candidate IDs, original saved residual/constraint metadata, cost coverage, online truth declaration, available public acceptance checks and artifact coverage.

Existing `audit/{campaign}_PUBLIC_ACCEPTANCE/AUDIT.json` defects/index checks are reused if present; the wrapper does not run a new matrix-product validation. Missing validation is null. Saved physical residual metadata is copied with its original meaning; no new Y-fit or material quality is computed. Auxiliary gap/event metadata may be retained, but full CIE/G19 PR trajectories and exact optimizer continuation remain null with the known coverage gap. Unretained per-stage beta/iteration settings are not silently recreated from current source; source-qualified settings are documented separately in `BASELINE_SOURCE_AUDIT.md`.

Physical aggregate RHS takes precedence over disjoint role counters, and vector source/adjoint roles are preserved where appropriate. Internal and external CPU remain separate alternatives; nested auxiliary physical deltas are not added again. The wrapper does not claim native schedule completion merely because a process wrote `COMPLETE` or a current-space stage was set up.

## Narrow metadata validation

The script was exercised on terminal exposed `development_v3_local` plus the locally unavailable C3 controls manifest. It emitted all54 expected rows: six terminal execution wrappers and48 pending/null execution rows. No C3 result/public payload or fresh quality was read. Two successful implementation checks are retained in `audit/RUN_CONTRACTS_METADATA_VALIDATION/` and final `audit/RUN_CONTRACTS_METADATA_VALIDATION_FINAL/`; no failed metadata build was hidden.

This is a metadata-completeness/provenance wrapper, not a new scientific validation, seal authorization, recovery certificate or proof that all original requested fields were historically saved.
