# R5 focused prior-art evidence audit

2026-11-08 (Asia/Singapore). Evidence QA only; no global novelty certification, new-theorem review, experiment or training. The parent owns those judgments.

## Retrieval status

Read `scholarqa-research/SKILL.md`, `references/reporting.md` and `references/retrieval.md`. One CLI collection used exactly two targeted formulations; query text and all four `operation_errors` are preserved in `collect.json`. Both paper and snippet endpoints returned HTTP 429. No discovery candidates were returned. One subsequent batch verification of the two supplied arXiv IDs also returned 429, preserved in `verified.json`. No repeated rate-limit retries or generic web search occurred.

The two papers were already explicitly identified in the task. Their official arXiv abstract and full-text HTML pages were inspected directly to fill that documented gap. This is primary-source inspection of identified papers, not successful ScholarQA/Semantic Scholar retrieval. Canonical Semantic Scholar IDs remain **UNAVAILABLE / UNVERIFIED** and are not cited. ArXiv identity, authors and year were checked against the official records below. No third paper was selected: the two-query discovery failed and no directly matched learned-SPD-metric paper was verified within this bounded task.

## Verified primary-source references

1. Pengju Liu, Hongzhi Zhang, Kai Zhang, Liang Lin, Wangmeng Zuo (2018), **Multi-level Wavelet-CNN for Image Restoration**, [arXiv:1805.07071](https://arxiv.org/abs/1805.07071), [full text v2](https://arxiv.org/html/1805.07071v2). Official arXiv comments specify **CVPR NTIRE Workshop 2018**; avoid presenting this record as a main-conference CVPR paper without separate verification.
2. Yu Liu, Hao Zhao, Rencheng Song, Xudong Chen, Chang Li, Xun Chen (2022), **SOM-Net: Unrolling the Subspace-based Optimization for Solving Full-wave Inverse Scattering Problems**, [arXiv:2209.03567](https://arxiv.org/abs/2209.03567), [full text v2](https://arxiv.org/html/2209.03567v2), [related journal DOI](https://doi.org/10.1109/TGRS.2022.3217162). Official record says accepted by IEEE TGRS.

## Claim ledger

Tier A = identified primary full text; Tier C = bibliographic metadata only. Supporting passages below are paraphrases and exact equation identifiers, not verbatim quotations. Source sections permit independent inspection.

| Bounded claim | Paper ID / location | Evidence | Tier / relationship | Context / confidence |
|---|---|---|---|---|
| Multi-level wavelet analysis, all-subband CNN processing and inverse-wavelet reconstruction are existing architecture choices. | arXiv:1805.07071v2, §§3.1–3.2 | Four DWT subbands are jointly input to CNN blocks between levels; IWT reconstructs high-resolution features. The transform pair permits exact reconstruction; learned nonlinear processing itself is not an identity guarantee. | A / supports | Image restoration; high |
| Full coefficients do not by themselves identify a new prior mechanism. | Same, §3.3 | The paper contrasts low-band-only pooling behavior with retaining all bands. | A / contextualizes | Interpretation of direct architectural overlap; high |
| Alternating learned-current and analytical-material updates already occur in nonlinear inverse-scattering unrolling. | arXiv:2209.03567v2, §II-B.1, Algorithm 1 | Four cascaded subnetworks update induced current, then reconstruct permittivity analytically. | A / supports | 2D TM scattering; high |
| Field consistency and material supervision in such an unrolled model are established. | Same, §II-B.2, Eqs.7–11 | Total loss includes induced-current, scattered-field, SSIM and material-MSE terms. The current target comes from MoM. | A / supports | Supervised physical-variable fitting; high |
| SOM-Net's displayed field loss is not the same operation as resolving a material state on every minibatch. | Same, Eq.2, Algorithm1 step8, Eq.9 | Predicted receiver field is G_S times the network-predicted induced current. | A / qualifies | Specific computational graph; high |
| No verified source here establishes equivalence to the proposed bounded-SPD learned metric and fixed-energy descent construction. | None retrieved | Both discovery queries failed; these two inspected mechanism sections cannot resolve the broader metric/fixed-point literature. | No substantive evidence / unresolved | **Do not turn missing evidence into an absence or novelty claim.** |

## Relation to proposed R5 components — analyst interpretation

- **Full-coefficient CNN prior:** cite MWCNN as direct architectural prior art. Wavelet completeness and multiscale nonlinear processing alone should not carry the claimed distinction.
- **Current/material unrolling and field-loss training:** cite SOM-Net. R5 instead declares no GT-current supervision and an actual material-state solve with retained-LU VJP. These differences do not establish novelty or superiority.
- **Learned bounded-SPD metric and fixed-energy descent:** assessment remains open. Neither the failed search nor this two-paper audit supports claiming the combination is unprecedented. Formal descent scope, prior-dependent energy changes and whole-iteration stability must be checked by the parent independently.
- **Final claim boundary:** the defensible comparison is a specific declared computational graph and training objective. Broad labels such as wavelet CNN, physics consistency, current/material alternation or unrolling already overlap the inspected literature.

This audit is non-exhaustive. It supplies two primary-source mechanism checks while explicitly failing Semantic Scholar discovery and batch-verification gates; it must not be described as a completed ScholarQA novelty review.
