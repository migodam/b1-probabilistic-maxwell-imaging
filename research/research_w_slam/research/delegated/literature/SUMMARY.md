# Focused literature audit — live verification BLOCKED

2026-10-25 UTC. Question: what closest primary work distinguishes ordinary incremental smoothing/selective relinearization from unknown-wall–target joint dynamic compensation, and does a directly adjacent method combine unknown wall, target pose/material, and task-risk reopening?

**Outcome:** live evidence retrieval and final citation verification are BLOCKED by Semantic Scholar DNS/network failure. No live primary-paper claim is verified in this run. No absence, originality, or global novelty finding follows. The five papers below are **local-package leads awaiting fresh verification**, not newly verified citations.

## Executed scope and receipts

Applied `/Users/migodam/.codex/skills/scholarqa-research/SKILL.md`, `references/retrieval.md`, and `references/reporting.md` in focused Evidence-QA mode. Read `source_package/Literature_Audit.md` and `source_package/References_Main.md` as seed records.

Subquestions were (1) what update/relinearization capability is already standard; (2) which wall–target compensation papers are physically closest; (3) whether the retrieved primary record supports a task-risk reopening bridge. Two formulations only:

1. `incremental smoothing selective relinearization task risk`
2. `unknown wall target full wave joint inversion`

`scholarqa-cli collect` requested one paper per formulation, no snippets: maximum two discovery records. A single `scholarqa-cli verify` batch requested five seed DOIs, yielding a planned maximum seven unique papers. No expansion, retry loops, generic web fallback, alternate providers, PDE work, or credential output occurred.

Both CLI commands exited 1. Both discovery operations failed with `Semantic Scholar API error: [Errno 8] nodename nor servname provided, or not known`. Batch verification returned the same error. Candidate IDs are empty. This is transport failure, not a successful empty search. Receipts:

- `rawbundle/collect.json`: question, date, two queries, empty candidate list, two operation errors.
- `rawbundle/verify.json`: top-level verify error.
- `rawbundle/verification_request.txt`: intended seed IDs and limits.

## Local-package comparison, explicitly not freshly verified

| Lead | Local audit says | Boundary that must remain | Local evidence pointer |
|---|---|---|---|
| Kaess et al., iSAM2 (2012), DOI 10.1177/0278364911430419 | Fluid/selective relinearization, partial state update, historical correction; author reports inspecting §4.1–4.3. | Ordinary revisit/cache/selective update is a strong inherited baseline. Its conditional-independence/sparsity complexity cannot be assumed for dense material Schur information after eliminating fields. | `source_package/Literature_Audit.md:103-109`; `source_package/References_Main.md`, S1 |
| Arablouei, IGG/SPO (2026), DOI 10.3390/robotics15050087, arXiv:2601.08110 | Information gate chooses local/global GN and retains observations. | Package reports full-text access to preprint but only metadata/abstract for journal version. Its information-only counterexample is not a fresh full-proof critique or a refutation of journal convergence results. | `source_package/Literature_Audit.md:75-99`, `111-119`; `source_package/References_Main.md`, S2 |
| Watson–Andre–Lionheart (2024), DOI 10.1088/1361-6420/ad5b83, arXiv:2403.10354 | Full-wave wall-background Green/BEM, reduced POD model and variable projection jointly handle wall nuisance/target reflectivity. | Package says target model uses Born; full-wave wall propagation is not by itself arbitrary all-order target–wall feedback. Do not describe this as already establishing wide unknown-support material posterior or task-risk reopening. | `source_package/Literature_Audit.md:211-219`; `source_package/References_Main.md`, S12 |
| Wang–Amin (2006), DOI 10.1109/TSP.2006.879325 | Multiple stand-off observations and common target consistency break modeled wall ambiguity. | Multi-observation wall–target consistency is already a direct physical precedent. This alone does not settle finite full-wave nonlinear compensation, pose/material-wide inference, or task-risk reopening. | `source_package/Literature_Audit.md:177-183`; `source_package/References_Main.md`, S21 |
| Yao et al. (journal 2025/DOI online year 2024), DOI 10.1109/TIM.2024.3502875 | Distributed TWI jointly estimates target location and wall parameters using target consistency. | Local author had abstract-level access only; probability, iteration, and revisit details remain unknown. Missing abstract statements are not evidence of absence. | `source_package/Literature_Audit.md:195-201`; `source_package/References_Main.md`, S11 |

The stable identifiers above are supplied seed strings, not successful DOI-resolution receipts from this run. Five-ID verification failed as a batch. The ledger records all five as `UNVERIFIED_LOCAL_PACKAGE_RECORD`; it does **not** relabel them fresh Tier A/B evidence.

## Evidence qualification and remaining gap

Tier A requires an identifiable retrieved full-text passage; Tier B requires retrieved abstract support; Tier C requires retrieved bibliographic metadata. **None of these tiers was newly obtained.** The local author’s F/A/M labels are recorded provenance, not the worker’s fresh A/B/C tier. No cited paper is promoted to verified merely because its DOI appears in a local bibliography.

The package itself frames the possible distinction as innovation-aware task residuals plus Schur coupling plus finite Maxwell feedback/error control at lower total cost, rather than ordinary smoothing or information gating (`source_package/Literature_Audit.md:63-73`, `355-365`). That is a **candidate claim requiring primary-paper comparison and physical/cost evidence**, not a novelty verdict.

The requested direct closest-prior question — whether an actual paper simultaneously imposes unknown wall, target pose/material, and task-risk reopening — remains **UNRESOLVED_DUE_TO_RETRIEVAL_FAILURE**. This run provides no justified yes/no. Parent retains all mathematical and scientific judgment. Further work requires restoring scholarly transport or explicit permission to use an alternate retrieval path.
