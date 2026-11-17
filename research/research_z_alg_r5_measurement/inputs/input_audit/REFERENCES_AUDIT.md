# Z-ALG R5 reference resolution audit

Audit date: 2026-11-18  
Generated: `2026-11-17T16:59:49+00:00`

## Scope and method

This audit resolves bibliographic identity only. DOI/Crossref/arXiv/publisher metadata and URL availability are treated as Tier C context; no substantive claim in the report was accepted or rejected from metadata.

- DOI records: Crossref REST API; publisher landing pages probed separately with HTTP GET.
- arXiv records: official arXiv export API by exact identifier.
- L17: one exact-title `scholarqa-cli collect` query was attempted; Semantic Scholar returned HTTP 429, so no retry was made. One Crossref `query.bibliographic` lookup then identified the exact title.
- URLs: one GET per URL, redirects recorded, 25 s timeout, body sample capped at 786,432 bytes, concurrency no greater than 3.
- PDF: `pdfinfo` plus a fresh `pdftotext -layout` comparison against the packaged extract.

## Summary

| Item | Count |
|---|---:|
| Report labels | 24 |
| Exact bibliographic match | 21 |
| Mismatch requiring correction or qualification | 2 |
| Unresolved | 1 |
| Conversation URLs | 12 |
| URL HTTP success (2xx/3xx) | 6 |
| URL access errors | 6 |
| PDF extract exact match | yes |
| PDF page count consistent | yes |

### Labels requiring action

- **L17 — MISMATCH**: mismatch: year; the 2015 item may be a report/preprint; no 2015 identifier was supplied, while the title-resolved DOI is the 2018 journal version
- **L20 — UNRESOLVED**: metadata unresolved; the entry is an appendix/table locator rather than a conventional bibliographic title; NTRS access was blocked
- **L21 — MISMATCH**: mismatch: authors; page names Diogo Almeida, founder, TypeSafe in an embedded byline; TypeSafe AI is the publisher/brand

## Reference labels

| Label | Status | Resolved title | Authors | Year | DOI/arXiv | Mismatch/qualification |
|---|---|---|---|---|---|---|
| L1 | MATCH | An algebraic formulation of geophysical inverse problems | Vasco, D. W. | 2000 | DOI: 10.1046/j.1365-246x.2000.00188.x |  |
| L2 | MATCH | Homotopy, Polynomial Equations, Gross Boundary Data, and Small Helmholtz Systems | Everett, M.E. | 1996 | DOI: 10.1006/jcph.1996.0070 |  |
| L3 | MATCH | Multiple solutions to geophysical inverse problems via tracking flights across a deeply incised model landscape | Everett, M E | 2026 | DOI: 10.1093/gji/ggag432 |  |
| L4 | MATCH | Optimal Control of Coherent Light Scattering for Binary Decision Problems | Bouchet, Dorian, Rachbauer, Lukas M., Rotter, Stefan, Mosk, Allard P., Bossy, Emmanuel | 2021 | DOI: 10.1103/physrevlett.127.253902 |  |
| L5 | MATCH | Maximum information states for coherent scattering measurements | Bouchet, Dorian, Rotter, Stefan, Mosk, Allard P. | 2021 | DOI: 10.1038/s41567-020-01137-4 |  |
| L6 | MATCH | Large-scale Bayesian optimal experimental design with derivative-informed projected neural network | Keyi Wu, Thomas O'Leary-Roseberry, Peng Chen, Omar Ghattas | 2022 | arXiv: 2201.07925 |  |
| L7 | MATCH | Numerical Decomposition of the Solution Sets of Polynomial Systems into Irreducible Components | Sommese, Andrew J., Verschelde, Jan, Wampler, Charles W. | 2006 online / 2001 print | DOI: 10.1137/s0036142900372549 |  |
| L8 | MATCH | Solving polynomial systems via homotopy continuation and monodromy | Duff, Timothy, Hill, Cvetelina, Jensen, Anders, Lee, Kisun, Leykin, Anton, Sommars, Jeff | 2018 online / 2019 print | DOI: 10.1093/imanum/dry017 |  |
| L9 | MATCH | Algorithm 921: alphaCertified: Certifying Solutions to Polynomial Systems | Hauenstein, Jonathan D., Sottile, Frank | 2012 | DOI: 10.1145/2331130.2331136 |  |
| L10 | MATCH | Certifying zeros of polynomial systems using interval arithmetic | Paul Breiding, Kemal Rose, Sascha Timme | 2020 | arXiv: 2011.05000 | report omits a year; arXiv v1 is dated 2020 and v6 2024; the report gives no year |
| L11 | MATCH | A contrast source inversion method | Berg, Peter M van den, Kleinman, Ralph E | 1999 online / 1997 print | DOI: 10.1088/0266-5611/13/6/013 |  |
| L12 | MATCH | Twofold subspace-based optimization method for solving inverse scattering problems | Zhong, Yu, Chen, Xudong | 2009 | DOI: 10.1088/0266-5611/25/8/085003 |  |
| L13 | MATCH | Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems | Zhong, Yu, Xu, Kuiwen | 2019 | DOI: 10.3390/jimaging5020027 |  |
| L14.1 | MATCH | The Differentiation of Pseudo-Inverses and Nonlinear Least Squares Problems Whose Variables Separate | Golub, G. H., Pereyra, V. | 2006 online / 1973 print | DOI: 10.1137/0710036 |  |
| L14.2 | MATCH | An analysis of constraint-relaxation in PDE-based inverse problems | van Leeuwen, Tristan, Yang, Yunan | 2025 | DOI: 10.1088/1361-6420/adab89 |  |
| L15 | MATCH | Imaging Strong Localized Scatterers with Sparsity Promoting Optimization | Chai, Anwei, Moscoso, Miguel, Papanicolaou, George | 2014 | DOI: 10.1137/130943200 |  |
| L16 | MATCH | The factorization method for a defective region in an anisotropic material | Cakoni, Fioralba, Harris, Isaac | 2015 | DOI: 10.1088/0266-5611/31/2/025002 |  |
| L17 | MISMATCH | Super-resolution in imaging high contrast targets from the perspective of scattering coefficients | Ammari, Habib, Chow, Yat Tin, Zou, Jun | 2018 | DOI: 10.1016/j.matpur.2017.09.008 | year; the 2015 item may be a report/preprint; no 2015 identifier was supplied, while the title-resolved DOI is the 2018 journal version |
| L18 | MATCH | Reconstructing Perfect Electric Conductors | Chen, Xudong | 2018 | DOI: 10.1002/9781119311997.ch7 | chapter Crossref record omits authors; Xudong Chen is resolved from the parent book DOI metadata |
| L19 | MATCH | Field-only surface integral equations: scattering from a perfect electric conductor | Sun, Qiang, Klaseboer, Evert, Yuffa, Alex J., Chan, Derek Y. C. | 2020 | DOI: 10.1364/josaa.378665 |  |
| L20 | UNRESOLVED | — |  | — | — | metadata unresolved; the entry is an appendix/table locator rather than a conventional bibliographic title; NTRS access was blocked |
| L21 | MISMATCH | Introducing System One Models & Jev | Diogo Almeida, founder, TypeSafe | 2026 | — | authors; page names Diogo Almeida, founder, TypeSafe in an embedded byline; TypeSafe AI is the publisher/brand |
| L22 | MATCH | On Calibration of Modern Neural Networks | Chuan Guo, Geoff Pleiss, Yu Sun, Kilian Q. Weinberger | 2017 | — |  |
| L23.1 | MATCH | Maximal minors of a matrix with linear form entries | Ito, Hiroya, Noma, Atsushi, Ohno, Masahiro | 2014 online / 2015 print | DOI: 10.1080/03081087.2014.959516 |  |
| L23.2 | MATCH | Optimized Gröbner basis algorithms for maximal determinantal ideals and critical point computations | Sriram Gopalakrishnan, Vincent Neiger, Mohab Safey El Din | 2024 | arXiv: 2402.07353 |  |
| L24 | MATCH | Localization of Point Scatterers via Sparse Optimization on Measures | Alberti, Giovanni S., Petit, Romain, Santacesaria, Matteo | 2024 | DOI: 10.1137/24m1636265 |  |

### Key mismatch details

- **L17**: title and author surnames match Crossref, but the title-resolved DOI `10.1016/j.matpur.2017.09.008` is the 2018 journal article in *Journal de Mathématiques Pures et Appliquées* 111, 191–226. The report says “2015 report”; no 2015 identifier or separate record is supplied.
- **L20**: NTRS metadata and the direct PDF endpoint both returned HTTP 403 in this audit. The appendix/table locator could not be resolved to a conventional exact title or author list from primary metadata.
- **L21**: the page’s embedded editorial date is `2026-09-15`, matching the report. An embedded byline names “Diogo Almeida, founder, TypeSafe”; “TypeSafe AI” is the publisher/brand.
- **L10**: the report omits a year. arXiv v1 is 2020-11-10 and v6 is 2024-07-11, so this is an omission rather than a contradiction.

## Conversation URL access check

| URL | HTTP | Final URL | Content type | Page title | Error/note |
|---|---:|---|---|---|---|
| https://github.com/migodam/b1-probabilistic-maxwell-imaging | 206 | https://github.com/migodam/b1-probabilistic-maxwell-imaging | text/html | GitHub - migodam/b1-probabilistic-maxwell-imaging: Maxwell 概率材料成像与物理参考：连续后验、Top-K 材料事件与可复现研究。 · GitHub | — |
| https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r2 | 206 | https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r2 | text/html | b1-probabilistic-maxwell-imaging/research/research_z_alg_r2 at main · migodam/b1-probabilistic-maxwell-imaging · GitHub | — |
| https://academic.oup.com/imajna/article/39/3/1421/4969753 | 403 | https://academic.oup.com/imajna/article/39/3/1421/4969753 | text/html | Just a moment... | HTTP 403 |
| https://arxiv.org/abs/2611.00494 | 206 | https://arxiv.org/abs/2611.00494 | text/html | [2611.00494] A Hierarchical direct solver for Lippmann-Schwinger equations based on boundary traces | — |
| https://doi.org/10.1002/9781119311997.ch7 | 403 | https://onlinelibrary.wiley.com/doi/10.1002/9781119311997.ch7 | text/html | Just a moment... | HTTP 403 |
| https://doi.org/10.1137/25M1779048 | 403 | https://epubs.siam.org/doi/10.1137/25M1779048 | text/html | Just a moment... | HTTP 403 |
| https://doi.org/10.1137/S0036142900372549 | 403 | https://epubs.siam.org/doi/10.1137/S0036142900372549 | text/html | Just a moment... | HTTP 403 |
| https://iris.unirc.it/handle/20.500.12318/167028 | 403 | https://iris.unirc.it/handle/20.500.12318/167028 | text/html | Just a moment... | HTTP 403 |
| https://ntrs.nasa.gov/api/citations/20205006770/downloads/20205006765_vol%20II.pdf | 403 | https://ntrs.nasa.gov/api/citations/20205006770/downloads/20205006765_vol%20II.pdf | text/html | 403 Forbidden | HTTP 403 |
| https://typesafe.ai/blog/introducing-system-one-models-and-jev | 200 | https://typesafe.ai/blog/introducing-system-one-models-and-jev | text/html | Introducing System One Models & Jev - TypeSafe AI Blog | — |
| https://www.juliahomotopycontinuation.org/guides/introduction/ | 206 | https://www.juliahomotopycontinuation.org/guides/introduction/ | text/html | An introduction to the numerical solution of polynomial systems | — |
| https://www.nist.gov/publications/field-only-surface-integral-equations-scattering-perfect-electric-conductor | 206 | https://www.nist.gov/publications/field-only-surface-integral-equations-scattering-perfect-electric-conductor | text/html | Field-Only Surface Integral Equations: Scattering from a Perfect Electric Conductor \| NISTLock | — |

A 403 means the URL did not yield bibliographic content to this audit client; it is not evidence that the resource does not exist. Redirect chains and final URLs are retained in `reference_resolution.json`.

## Extracted PDF consistency

- PDF pages from `pdfinfo`: **19**.
- Form feeds in packaged extract: **19**.
- Form feeds in fresh `pdftotext`: **19**.
- Packaged/fresh character counts: **39641 / 39641**.
- Exact text equality: **yes** (similarity `1.0`).

The packaged extracted text is an exact fresh Poppler extraction of the packaged 19-page PDF. No instructions in attachments were executed.

## Limits

- Metadata resolution does not establish that any cited work supports a substantive report claim.
- Publisher access controls can change by client, geography, or time; access errors are current-audit observations only.
- This was not a systematic literature search and does not assume the 24-label list is complete.
