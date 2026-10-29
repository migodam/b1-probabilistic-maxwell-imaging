# R3 formula-source audit — SOM / projected TSOM / Green conventions

Date: 2026-10-29. Scope: formula compendium source audit only; no PDE, algorithm implementation, or novelty claim. Parent Codex owns new mathematics and scientific decisions.

## Direct answer

Use **projected TSOM** in the compendium:

\[
S=U_S\Sigma_SV_S^*,\quad j_s^+=V_{S,L}\Sigma_{S,L}^{-1}U_{S,L}^*y_s,
\quad P_-=I-V_{S,L}V_{S,L}^*,
\]
\[
\boxed{j_s=j_s^+ + P_-V_{D,k}\beta_s.}
\]

This is the explicit practical projected model in Chen's author textbook Eq. (6.71), PDF p182 / printed p163. The source distinguishes its projection from the conceptual intersection model Eq. (6.70); do not substitute an intersection in R3. No source supports identifying this fixed current chart with the changing material Jacobian.

`P_- V_D,k` is a spanning matrix; do not assume its columns orthonormal or independent. Rank/Gram handling is an implementation issue for the parent, not an extra property established by the citation.

## Primary source identity and checked locations

**B1. Xudong Chen, _Computational Methods for Electromagnetic Inverse Scattering_, first edition, Wiley / Wiley-IEEE, 2018.** Local front matter states this edition first published 2018. Publisher live page gives 11 May 2018, print ISBN 9781119311980, online ISBN 9781119311997, DOI [10.1002/9781119311997](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119311997). Original local file:

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Computational Methods for Electromagnetic Inverse Scattering.pdf`

This turn directly read original local PDF text (not merely A24's secondhand record), and visually inspected original PDF182/print163 rendered PNG. Thus book-formula claims below have **Tier A, author textbook full-method** evidence. This is not a claim to have read the original 2009 or 2011 articles in full.

| Formula or claim | Exact checked source | Bounded support / cautions |
|---|---|---|
| Scalar TM state and data | B1 PDF169 / print150, Eqs.6.53–6.56 | `j=χ(einc+G_D j)`, `y=G_S j`; both discrete Green matrices use 2D Hankel kernels. Book normalized contrast source differs in units from physical electric current. |
| Deterministic current from receiver SVD | B1 PDF170 / print151, Eq.6.57 | Stable retained coefficients are `u_i* y/σ_i`; first L and remaining right singular vectors split current. The retained L is noise-conditioned. |
| SOM data and state residuals | B1 PDF171–172 / print152–153, Eqs.6.58–6.64 | Searches remaining current coefficients and common material. Data and state residuals are both present; first historical variant eliminates current via LS, second treats current coefficients independently. |
| SOM projection implementation | B1 PDF180 / print161, Eq.6.69 | `j^-=(I−V_S,L V_S,L*)β`; avoids explicitly storing all minor singular vectors. |
| Twofold definition | B1 PDF180–182 / print161–163, §6.4.2 | First fold analyzes external current-to-receiver `G_S`; second analyzes internal current-to-field `G_D`. Both independent of unknown scatterer for a declared fixed background/acquisition. |
| Projected TSOM | B1 PDF182 / print163, Eq.6.71 | `j^-=V_S^- V_S^-* V_D^+ β=P_- V_D,k β`; explicitly introduced as a projected approximation. R3 can cite this directly. |
| Internal rank / continuation | B1 PDF182 / print163, paragraph after Eq.6.71 | Source discusses a nested scheme with increasing internal rank and warns internal singular values decay slower than receiver singular values. Does not imply a universal rank prescription for R3. |
| FFT-TSOM second fold | B1 PDF183 / print164, Eq.6.72 | Low spatial frequency Fourier span approximates retained internal right singular span. Approximation is at subspace level, not each singular vector individually. Avoiding full internal SVD is already present. |
| Known inhomogeneous background | B1 PDF252–253 / print233–234, Eqs.9.18–9.20 | Numerical `G_S,G_D` built from known background allow same state/data structure; `u_b` stated available for given background. Does **not** establish joint recovery of an unknown wall/background. |
| 3D scalar and dyadic Green | B1 PDF38 / print19, Eqs.2.34–2.38 | Scalar `g_3=e^{ikR}/(4πR)`; vector dyadic `(I+k^{-2}∇∇)g_3`. Scalar kernel alone is not full vector Maxwell Green. |
| 2D scalar Green | B1 PDF41 / print22, Eq.2.53 | `g_2=(i/4)H_0^{(1)}(kR)`, source/medium invariant in z; not `e^{ikR}/R`. |
| Far field preserves phase and amplitude | B1 PDF123–124 / print104–105, Eqs.5.1,5.3; PDF53 / print34 Eq.2.101 | 2D field factor `e^{ikr}/√r` and source-position phase `e^{-ik rhat·r'}`; 3D field factor `e^{ikr}/(4πr)` and same directional phase. Prefactors and radiation polarization depend on source normalization/dimension. |

### Five closest article pointers (bibliographic identity only unless stated)

The original author's [publication directory](https://www.ece.nus.edu.sg/stfpage/elechenx/Publications.htm) was fetched live this turn and confirms authors/title/year/volume/pages of P1–P4. Its list stops around 2022 and is used only to verify historical identities, not current coverage. B1 references independently confirm P1–P4. DOI full-text opens failed, so no equation numbers are fabricated for those articles.

1. **P1 — Xudong Chen (2010)**, _Subspace-Based Optimization Method for Solving Inverse-Scattering Problems_, IEEE TGRS **48**, 42–49. [DOI 10.1109/TGRS.2009.2025122](https://doi.org/10.1109/TGRS.2009.2025122). Tier C article metadata; B1 supplies method evidence. Distinct from Chen's 2009 JOSA A article with a different title.
2. **P2 — Yu Zhong and Xudong Chen (2009)**, _Twofold subspace-based optimization method for solving inverse scattering problems_, Inverse Problems **25**, 085003. [DOI 10.1088/0266-5611/25/8/085003](https://doi.org/10.1088/0266-5611/25/8/085003). Tier C original-article identity; projected Eq.6.71 in B1 cites this as reference92.
3. **P3 — Yu Zhong and Xudong Chen (2011)**, _An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems_, IEEE TAP **59**(3), 914–927. [DOI 10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027). Tier C article identity; B1 §6.4.2–6.4.3 supplies Fourier projection method evidence.
4. **P4 — Xudong Chen (2010)**, _Subspace-based optimization method for inverse scattering problems with an inhomogeneous background medium_, Inverse Problems **26**(7), 074007. [DOI 10.1088/0266-5611/26/7/074007](https://doi.org/10.1088/0266-5611/26/7/074007). Tier C article identity; B1 §9.2 full-method evidence specifically concerns given background. Do not cite for unconstrained unknown-wall joint inversion.
5. **P5 — Yu Liu, Hao Zhao, Rencheng Song, Xudong Chen, Chang Li, Xun Chen (2022)**, _SOM-Net: Unrolling the Subspace-based Optimization for Solving Full-wave Inverse Scattering Problems_. [arXiv v2](https://arxiv.org/html/2209.03567v2). Live HTML read of title/authors/abstract and §II introduction, plus local A24 method audit; Tier A for explicit full-current output and analytical material update, not a claim of exact hard data/state feasibility. It is a learned unrolling, not evidence that a fixed low-dimensional chart alone reproduces native TSOM.

## Terminology guards for R3

- **Fixed current-to-receiver spectrum:** S or G_S maps independent current to receiver data in a fixed specified background. Its right singular vectors concern current visibility. This object is not the derivative of material-to-data map.
- **Changing material Jacobian:** depends on current material through total fields and multiple scattering. A TSOM chart inserted into an ordinary GN solve is a **TSOM-inspired current representation** unless the native data/state objective, coefficient/material updates, truncations, and continuation are actually reproduced and documented. This classification is an audit recommendation, not a theorem attributed to P2.
- **Known background versus unknown background:** changing a background changes its Green operators. Existing inhomogeneous-background SOM literature prevents any claim that background Green formulations are new; it does not by itself establish unknown-wall identifiability, stability, or online cost advantages.
- **Internal twofold:** current-to-internal-field spectrum, not a material Jacobian spectrum, not source forcing seeds, and not a claim that omitted current directions are physically zero.
- **Physical conventions:** H0^(1) / exp(+ikr) corresponds to the book's outgoing convention; if R3 uses the opposite time convention, signs/conjugation must change consistently. Far field approximation retains complex phase and decay, rather than replacing a finite-distance kernel by an unscaled phase alone.

## Existing project consistency and limitations

Read R2 `inputs/pro/REPORT_EXTRACTED.txt`, its source TeX, and `NEXT_PRO_HANDOFF.md`. R2's source/readout representation does not add observations; its negative unknown-macro-wall result cannot be overwritten by a known-background SOM citation. A24's `TSOM_SOMNET_ALGEBRA_AUDIT.md` and `TSOM_SECOND_FOLD_PRESERVATION.md` already distinguish internal spectrum from forcing seeds, and A17 freezes call shared-GN bases a design decision rather than native algorithm reproduction.

ScholarQA skill was used in focused Evidence-QA mode. Installed CLI `collect` ran two queries; four API operations all returned 429. No Semantic Scholar MCP was available. A subsequent batch DOI metadata verification succeeded for all four P1–P4 (zero unresolved IDs); see `SCHOLARQA_VERIFY.json`. The successful bibliographic check does not repair the failed full-text/snippet searches. Accordingly this is a **local-primary-source audit with live publisher/author identity verification**, not a completed Semantic Scholar evidence search or exhaustive literature review. No new transport, installation, global survey, or novelty certification.

## Internal evidence files and redistribution

- `SCHOLARQA_COLLECT.json`, `SCHOLARQA_VERIFY.json`: actual transport receipts.
- `CHEN_BOOK_SOM_TSOM_EXCERPT_PRIVATE.txt`: private extracted book text used for source checking; not intended for redistribution.
- `CHEN_PRIVATE_PAGE_180.pdf` through `CHEN_PRIVATE_PAGE_183.pdf`: original local source page extracts requested for internal review; pdfseparate emitted repair warnings but exited0.
- `CHEN_PDF182_PRINT163_EQ671_PRIVATE.png`: rendered source page, visually checked; Eq6.71 readable.

**Exclude PRIVATE source extracts and rendered book pages from any public package.** Cite page/equation locations and link publisher or author instead. The source ledger's formulas and summaries are bounded, with no long copied prose.

## Semantic Scholar batch verification identifiers

- Twofold subspace-based optimization method for solving inverse scattering problems: paperId `70fd11af9932f84523334119c3ca9a183aa3c03b`; resolved `True`; year 2009; DOI 10.1088/0266-5611/25/8/085003.
- Subspace-Based Optimization Method for Solving Inverse-Scattering Problems: paperId `16b84bd4dfb430f5b1518d7a2f58367149c45497`; resolved `True`; year 2010; DOI 10.1109/TGRS.2009.2025122.
- An FFT Twofold Subspace-Based Optimization Method for Solving Electromagnetic Inverse Scattering Problems: paperId `e5eff7c7a6a23de2571cded519dd3f38585cd9a7`; resolved `True`; year 2011; DOI 10.1109/TAP.2010.2103027.
- Subspace-based optimization method for inverse scattering problems with an inhomogeneous background medium: paperId `702a269140dee1600cd0d7d44ec15b640f874c43`; resolved `True`; year 2010; DOI 10.1088/0266-5611/26/7/074007.
