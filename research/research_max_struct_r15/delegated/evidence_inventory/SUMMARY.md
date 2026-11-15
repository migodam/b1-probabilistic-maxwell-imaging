# R15 mechanical evidence inventory

No scientific validity/novelty judgment. Read-only tables, source definitions and NPZ/NPY headers; no array payloads, truth payloads, physics, scoring, SVD, remote requests or hashes. Old files unchanged.

## Finite exposed cohort and actual final tables

The observable19 follow-up is **24 already exposed r14_new_0000–0023 cases ×6 fixed arms**. ALL_VIEWS has2304 rows over8 view labels, including144 each primary/raw/corrected. The extracted432 rows and18 summary rows retain original final values; no score recomputation.

| Arm | raw recovery /24 | corrected recovery /24 | primary recovery /24 | raw mean NRMSE | corrected mean NRMSE | primary mean NRMSE |
|---|---:|---:|---:|---:|---:|---:|
| g19_anchor_observation_galerkin | 10 | 12 | 12 | 0.172510847 | 0.159953636 | 0.159953636 |
| g19_anchor_observation_right | 6 | 12 | 12 | 0.178751657 | 0.162208756 | 0.162208756 |
| g19_internal_galerkin | 5 | 10 | 10 | 0.207666547 | 0.198432839 | 0.198432839 |
| g19_internal_right | 4 | 10 | 10 | 0.206821156 | 0.19766735 | 0.19766735 |
| g19_observation_galerkin | 9 | 12 | 12 | 0.173226793 | 0.158919149 | 0.158919149 |
| g19_observation_right | 7 | 12 | 12 | 0.176242559 | 0.161073889 | 0.161073889 |

Recovery columns are the inherited **posthoc descriptive** NRMSE≤.15 and IoU≥.70 rule, not GT selection or a new inference criterion. Physical Phi/noise status and raw/corrected separation remain in the compact table. Historical CIE/Dual/A+/11/G19/DR rows are existing R14 references, not new equal-cost six-arm runs.

## Endpoint reuse and separation

- Public inputs: all24 new-case files exist at parent R14 data_public/independent24; NPZ headers include G(576×576), S(48×576), E0/Y, public initialization/noise/bounds. Truth sources remain in sealed_truth/independent24 and were not opened. Their original scoring receipt uses Windows paths; the local original-generation counterpart is separately named.
- Supplement full primary/raw/corrected endpoint entries:432 in existing archive PACKAGE_PLAN, **0 local payloads**. No remote assets fetched. Local final tables alone do not supply current/field arrays.
- Existing R14 historical new24 primary exports:96 files over DR/global_dual/G19/direct; **all are m-only**. Old raw/corrected and DR independent iterator are absent locally in this scope.
- Old R13 confirmation_core48: three relevant methods each have48 primary/raw/corrected endpoints, each with m,j,Q,Y_prediction (432 full-field endpoints total). DR also has48 FINAL_AUXILIARY_STATE files with X,Q. Full-current/full-APCS lack that separate saved iterator file. Headers include physical j/Q shapes, not evaluated values.
- R13 exposed48 public inputs/truth source paths are independently indexed. This is another cohort/acquisition:12 sources,24 half-circle receivers, matched24-grid; new24 is48 full-circle receivers, with matched/crossgrid distinctions retained in COHORT_PUBLIC_TRUTH_SEPARATION.csv. Do not pair r13_fresh_0021 with r14_new_0021.
- Twelve previously curated pair01/11 complete checkpoints for new0012/0021 at128/512/1500 exist and expose J,Q,A/old_grad/old_direction and representation caches. They are not the six-arm observable19 endpoints.

## Existing definitions, without new derivation

- GN-Green / GN-NN energy: fixed complex right-domain orthogonal projections; squared Euclidean/Frobenius norms, not electromagnetic stored energy. Main sigma/sigmamax cutoff .1, numerical complement1e-12, intersections1e-9. Old geometry G19/557/0 and S9/39/528; GN training visible S includes **all48**, giving bright union67 and dark509 with Gstrong19.
- Sequential Pin→Pob and reverse-order energy tables are separate; neither is strict intersection. Numerical intersections retain47 unclassified dimensions in the old geometry. Projected/orthogonalized spans can overlap and exceed100% when separately summed; do not renormalize away duplicate coverage.
- R10 implemented chart is [Vo,strong,orth((I−Po,strong)Vg,strong)]: it spans the **sum/union** of the strong spaces, not their intersection. The second fold projects into joint receiver weak+null while keeping shared internal coefficients; independently parameterizing both pieces changes the span. Frozen .1 / residual rank cutoff1e-10. R10 old ranks19/9/28 and its fresh24 ranks19/11/30 are geometry-specific.
- Supplement basis contract is fixed19 for G, raw S, and O0=S(I−Dm0G)^−1 at public m0=.1; PGP/GP are distinct projections. It is not the old Sstrong9 or all48 visible projection, and O0 is not a material Jacobian.

## Files and unresolved reuse limits

SOURCE_EVIDENCE.csv gives exact source paths/line anchors and scope cautions. FILES_AND_NPZ_HEADERS.csv indexes1047 explicitly selected paths,663 present, with key names and selected NPY shape/dtype headers. AVAILABILITY_COUNTS.csv is the compact file-presence summary. No NPZ header errors were recorded. Saved GN-Green fixed_spectra/error_arrays/cumulative_mode_energies exist; the latter errors are truth-derived diagnostics, not public inverse inputs.

UNRESOLVED_ISSUES.csv records7 mechanical reuse limits. Principal blocker for the **new24 six-arm** offline error decomposition is local absence of original full-field endpoints; old R13 endpoints are available for a separate declared cohort, and selected-m exports cannot substitute for independent current/iterator states. Parent chooses scientific scope and any later truth access.

Inventory pass cost (including24+48 public header scans): CPU 0.194050s; wall 0.253768s. Earlier endpoint-only pass CPU.198044s/wall.450234s; both were read-only. No publication performed.
