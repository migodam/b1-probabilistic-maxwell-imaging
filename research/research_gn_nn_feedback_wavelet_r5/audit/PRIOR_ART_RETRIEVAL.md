# Retrieval log

- Primary transport: `/Users/migodam/.local/bin/scholarqa-cli`.
- `collect` question and exactly two queries are preserved verbatim in `collect.json`; limits were three papers and three snippets per query. Exit status 1, all four operations HTTP429, no candidates.
- `verify ARXIV:1805.07071 ARXIV:2209.03567` requested paperId,title,authors,year,externalIds,abstract,url. Exit status1, top-level HTTP429; preserved `verified.json` unchanged.
- No API retry or graph expansion. Canonical S2 paper IDs unresolved.
- Direct inspection of task-supplied identified seeds used official arXiv abstract pages and linked v2 HTML full texts, through web open/click/find. No generic search. MWCNN fulltext inspected §§3.1–3.3; SOM-Net inspected §II-A and §II-B, including Eqs.2,7–11 and Algorithm1.
- Direct-record bibliographic checks: exact title/authors/year/arXiv IDs and comments/related DOI. This does not claim the unsuccessful CLI verification passed.
- No third paper or learned-SPD comparison source verified; bounded task stopped after concrete mechanism checks.
