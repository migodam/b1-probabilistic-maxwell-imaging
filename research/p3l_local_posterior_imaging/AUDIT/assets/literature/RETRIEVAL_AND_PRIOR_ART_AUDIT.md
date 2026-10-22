# Focused prior-art verification and retrieval outcomes

Search date 2026-10-22 UTC. Scope: conditional marginalization/material Bayes risk; Gaussian decoder moments/Rao-Blackwellization; derivative-informed ROM/DA posterior inference. No broad review, graph crawl, physical calls, training, or novel theorem derivation. Parent retains final novelty judgment.

## Executed retrieval and limits

Read scholarqa-research SKILL.md and retrieval/reporting guides completely. Executed one scholarqa-cli collect with exactly two formulations:

1. likelihood informed conditional marginalization Rao Blackwell posterior moments Bayesian inverse problems
2. derivative informed neural operator geometric MCMC delayed acceptance posterior potential error

All four underlying paper/snippet operations returned HTTP429 Too Many Requests; candidate_paper_ids=[]; raw errors preserved in collect.json and collect.stderr. No discovery retries or alternative broad search followed. The collect exit status itself was not captured separately; the surrounding shell completed after printing the preserved bundle, which must not be treated as successful retrieval.

Canonical-ID batch verification succeeded for 8/9 requested identities. The attempted DOI:10.1198/106186005X7695 did not resolve and is excluded from accepted citations; do not copy it. verify.json and verify.exit_status.txt preserve the failure (exit1). A targeted three-record batch requested abstracts for Tallis, Casella-Robert and Tierney-Kadane; all identities resolved but publisher-elided abstracts were null (foundational_records.json, exit0). A final two-record batch verified already-identified Zahm and Andrieu-Roberts (correction_records.json, exit0). No Semantic Scholar canonical record is marked resolved when it was not.

Already-identified primary arXiv pages and selected fulltext were inspected to fill the documented 429 gap. This verifies close-known sources, not search coverage or absence of closer work. Publisher access failed for Tallis DOI, Casella-Robert DOI/OUP, and Tierney-Kadane DOI (403). An attempted JMLR page was inaccessible; DINO publication status is not upgraded from its inspected 2024 arXiv v2. S2 years may be preprint years: Cui et al. record2015 versus journal2016; Spantini et al. record2016 versus journal2017. Preserve both rather than claiming a bibliographic conflict.

## Evidence-supported mature prior art

- [Cui et al.,2014](https://arxiv.org/abs/1403.4680) already identify likelihood-informed structure across posterior support. [Cui-Marzouk-Willcox,2016](https://arxiv.org/abs/1510.06053) combine parameter/state reduction and retain a prior complement. The latter is an approximate posterior factorization, not proof that a frozen complement reproduces the exact strong conditional marginal.
- [Spantini et al.,2017](https://arxiv.org/abs/1607.01881) explicitly optimize QoI posterior mean approximations using weighted Bayes risk in a linear-Gaussian setting. This rules out claiming task Bayes risk or mean-oriented reduction as new; it supplies no nonlinear P1 improvement guarantee.
- [Zahm et al.,2022](https://arxiv.org/abs/1807.03712) certify nonlinear posterior approximation using functional inequalities and integrated gradient information. One-anchor spectrum alone is an insufficient replacement for that contract.
- [Cao et al.,2024,v2](https://arxiv.org/html/2403.08220v2) already join derivative-informed operator learning, geometry, delayed acceptance, full-misfit correction, and training-cost accounting. Inspected sections2.8,4.3,5,AppendixC. Their misfit-error analysis concerns DA behavior, rather than enclosing bounds yielding identical decisions to a fixed Full-MH path. Examples are diffusion-reaction and hyperelasticity, not Maxwell. General architecture overlap is direct; measured TM24 benefit remains open.
- [Chen-Ghattas,2020 preprint](https://arxiv.org/abs/2002.10924) directly connect goal-oriented reduced potential/gradient errors, posterior KL error and transport sample errors. Generic goal-ROM-to-posterior error composition is therefore established.
- [Andrieu-Roberts,2009](https://arxiv.org/abs/0903.5480) establish pseudo-marginal stationary-target preservation and distinguish inexact Monte Carlo marginal approximations. P3-L must retain its own precise estimator/current-randomness/zero-weight proof obligations.

The structured ledger records exact source IDs, evidence levels, faithful source summaries, support/qualification relationships and limits. The preceding claims use TierA/B support; the three moment references below remain TierC.

## Moments and conditional integration: verified identity versus content gap

[Tallis,1961](https://doi.org/10.1111/j.2517-6161.1961.tb00408.x), [Casella-Robert,1996](https://doi.org/10.1093/biomet/83.1.81), and [Tierney-Kadane,1986](https://doi.org/10.1080/01621459.1986.10478240) all resolve to the package's intended titles/authors/years. Their fulltext/abstract was not recovered here. Therefore this run does not claim independent theorem-level verification of their truncated-normal recurrence, Rao-Blackwell guarantees, or Laplace remainder assumptions. No title-only entry is used as substantive evidence.

The package's own explicit P1 polynomial/truncated-normal derivation and math checks can be assessed on their mathematical content by Codex. This retrieval does not establish that the P1 formula is absent from prior literature, that generic analytic moments are novel, or that Gaussian-reference moments are exact under the true Maxwell conditional posterior. A defensible candidate contribution must concern the specific decoder/correction/physical cost mechanism and demonstrated risk, rather than the names 'Rao-Blackwell' or 'Laplace'.

## Additional inspected primary evidence, excluded from S2-verified citation list

The already-identified [Christen-Fox author preprint](https://physics.otago.ac.nz/data/fox/publications/ApproxMCMCPreprint.pdf) was accessible. Its title/authors, algorithm1(pp2), and theorem1(pp3) directly show surrogate screening followed by full-target correction, with conditions for invariance/irreducibility. It is an independent fulltext evidence source; the unrelated unresolved DOI attempt is not used to identify it. The paper does not license accepting a surrogate-promoted move without the second correction. No exact DOI or canonical S2 ID is asserted here.

## Contradictions and remaining gaps for final review

1. Strong/weak geometry, task-risk reduction, ROM posterior error control and surrogate geometric DA have close established predecessors. Combining their names is insufficient novelty evidence.
2. The inspected reduced-prior-complement method is not equivalent to exact conditional integration. Conversely P3-L's exact-integration algebra should not be presented as already empirically achieved because a LIS paper exists.
3. Linear-Gaussian risk optimality does not promise nonlinear material mean gains, individual-object MAE gains, or legal-decoder mean outputs.
4. DA changes the transition kernel while preserving its target under conditions. A same-proposal/same-uniform enclosing evaluator has a different coupling contract; this focused, rate-limited retrieval does not establish novelty of interval-MH itself.
5. None of these inspected papers validates P1 escape mass, boundary-Laplace quality, teacher mixing or actual Maxwell certificate tightness. Those are experimental/theoretical obligations of this project.
6. No closest-Maxwell-specific absence claim is justified: both discovery formulations were completely rate-limited. Claims such as 'no prior work combines these components' should remain UNVERIFIED, not promoted to novelty.

Accepted machine records: verify.json (8 resolved), correction_records.json (2 resolved), foundational_records.json (identity-only duplicates). CLAIM_CITATION_LEDGER.json is the bounded support ledger. Raw failed IDs/errors are retained solely as retrieval provenance, not references to cite.
