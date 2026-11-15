# Focused prior-art search log

Date: 2026-11-15, UTC04:36 onwards. Evidence-only extraction; no experiment, truth access, LLM/API answer generation, agents, or final novelty/proof judgment.

Question: what primary literature already establishes state elimination/material-dependent residual metrics, practical CSI sequential updates, game-gradient analogies, and conditional inexact alternating convergence?

Transport: installed scholarqa-cli, following local scholarqa-research retrieval/reporting guidance. One collect invocation with two complementary queries, limits5paper/5snippet each, noyearfilter:

1. `cross correlated contrast source inversion constraint relaxation variable projection material dependent metric`
2. `inexact alternating minimization nonlinear conjugate gradient differentiable games sequential best response`

Receipts: `receipts/scholarqa_collect.json`, stderr. Exit1 partial bundle: query1 snippets429, paper searchtotal0; query2 snippets429 and papersearch429. operation_errors inspected, candidate_paper_ids empty. No expansion or redundant retry. This yields **no field-wide negative evidence**.

Fallback as explicitly authorized: inspect known primary fulltexts (arXiv2403.15292v2, arXiv1906.10864v1, officialPMLR Balduzzi2018), supplemented by two targeted web formulations locating author-hosted inexact-Gauss-Seidel and PR survey papers. Only primary sources used substantively. DOI landing URLs forvanLeeuwen/Yang and CC-CSI were inaccessible through webopen; officialarXiv and CrossrefAPI supplied identity instead.

Bibliography: scholarqa-cli verify six IDs (twoarXiv+matchingDOIs, ABS DOI, Balduzzi arXiv), exit0 resolved6/6, unresolved[]. Four unique records. Crossref rawJSON saved for constraint/CC-CSI/ABS. Dates differ by publicationstage: constraintS2year2024 corresponds arXiv, journal2025; ABSS2year2011 corresponds online, journalprint2013. CC-CSI journal2017 vs arXiv2019. PMLRauthor/title/year verified at officialpaperpage.

Downloaded source receipts: fulltext/constraint_relaxation.html, cc_csi.html, games.pdf, inexact.pdf, cg_survey.pdf, with text conversions. No content hash computation. A first HTML conversion attempt reported missingbs4; corrected using Pythonstandard-libraryHTMLParser, noinstall. PDFs converted bypdftotext. HagerCG survey additional S2URLverification attempted; final status stored separately, author-hosted publicationlist and PDF remain primary metadata sources.

Costs: onecollectbundle (4underlyingsearches), one6-IDverifybatch plus oneCG-IDverification attempt; threeCrossrefmetadatarequests; fivefulltexts downloaded. CLI doesnotgenerateprose, noLLMinferencepaid. Exact billednetworkcost/CPUtime notmetered; no scientific compute. Errorreceipts retained. No systematic/exhaustive coverage claim.

Final CG verification outcome: exit1, top-level429error (not unresolved scholarly identity); no retries. `receipts/scholarqa_verify_cg.json` retainserror. CG bibliographic identity verified from author's own publicationlist and linked fulltext, and must not be labeled S2-verified.
