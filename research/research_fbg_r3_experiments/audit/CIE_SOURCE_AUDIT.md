# CIE-I / Fourier-TSOM primary-source audit

Audit date: 2026-11-01. Scope: faithful algorithm structure for the present 2D real-material scalar-TM comparison, not novelty certification or reproduction of published 3D image errors.

## Acquired evidence

| Source | Current retrieval | Evidentiary use |
|---|---|---|
| Yu Zhong and **Kuiwen Xu**, *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*, Journal of Imaging 5(2), 27 (2019), DOI [10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027) | Complete publisher-deposited XML acquired from [Europe PMC](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8320915/fullTextXML). Stored `literature_cie2019.xml` and `literature_cie2019_body.txt`; full body read. | Tier A: formulation, twofold subspace, objective, update, continuation, stopping, dual beta, failures. XML front matter directly verifies title, authors and DOI. |
| Same article, [PMC landing page](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/) | Browser returned reCAPTCHA. | Failed route, not counted as body retrieval. |
| Same article, [MDPI](https://www.mdpi.com/2313-433X/5/2/27) | Browser returned HTTP 429. | Failed route; Europe PMC supplies the primary text. |
| Zhong, Lambert, Lesselier, Chen, 2016, *A New Integral Equation Method to Solve Highly Nonlinear Inverse Scattering Problems*, DOI [10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492) | Metadata verified against reference 29 of 2019 XML. Original full paper **not independently reacquired this round**. Supplied `memos/agent_cie.md` records prior full-text access. | The present algorithm follows the directly read 2019 description; inherited 2016 specifics are labeled as supplied-report evidence. |
| Zhong and Chen, 2011, *An FFT twofold subspace-based optimization method for solving electromagnetic inverse scattering problems*, DOI [10.1109/TAP.2010.2103027](https://doi.org/10.1109/TAP.2010.2103027) | Metadata verified against reference 8 of 2019 XML. Original full paper **SOURCE UNAVAILABLE in this round**. | Fourier-TSOM is reconstructed from 2019 Eqs. 9–15 and Secs. 2.2–2.3; no claim to independent line-by-line 2011 implementation verification. |

The supplied report abbreviates the 2019 second author as “J. Xu” in one bibliographic entry. The primary XML gives **Kuiwen Xu**; this audit uses that spelling.

## Claim-to-source ledger

| Algorithm claim | Primary location | Implementation consequence / scope |
|---|---|---|
| `j=m(E+Gj)` and `R=beta*m/(1+beta*m)` | Eqs. 1, 3–5 | Preserve the same scalar-TM physical state and positive pole-free branch. The source paper is vectorial 3D; this solver is an explicit 2D scalar adaptation. |
| DPCS from retained receiver singular vectors | Eq. 9, Sec. 2.3 step 1 | `jD=V_L diag(1/sigma_L) U_L* Y`; rank is frozen from public noise/stability information. No target support enters. |
| APCS obtained by projecting an IDFT current onto DPCS complement | Eqs. 10–11 | `Z=(I-V_L V_L*)F_mask`, complex Fourier coefficients independent for each source. Fourier columns are not generally orthonormal after projection. |
| Fourier restriction is a second regularization, not a material decoder | Sec. 2.2 | Full real material vector retained. The current mask controls internal propagation degrees of freedom. |
| Each 3D mask corner has edge length `MF` | Sec. 2.2 following Eq. 11 | 2D adaptation has four corners, each `mf×mf`, selecting `(2mf)^2` Fourier coefficients on an even grid. `mf` is not the number of modes. |
| Source-normalized data and current-state mismatches | Eqs. 12–15 | Fixed `1/||Y_s||²` and `1/||E_s||²` source weights. Implementation multiplies both squared terms by 1/2, with explicit common material prior added. |
| One Polak–Ribière current CG update, then LS contrast update | Sec. 2.3 step 3 | Sourcewise complex PR-CG and exact quadratic line minimization alternate with a shared real material block. PR+ and descent reset are disclosed numerical safeguards. |
| Smaller beta and larger mask, inheriting previous solution | Sec. 2.3 steps 4–6 | Native Fourier coefficients are zero-padded into a larger mask, and physical material is carried across stages. |
| Distinct beta for current and material in later stages | Sec. 2.4 | Recompute both modified contrasts from the **same physical m**. There is no claim of one common descending objective for the dual-beta method. |
| APCS stopping is source-average squared relative coefficient change | Eq. 16 and Sec. 3 | Default threshold `1e-3` on native Fourier coefficients, distinct from material convergence or physical reconstruction acceptance. |
| Published dual-beta may stop prematurely | Sec. 3.2, Example 2 | Keep premature stopping and original Maxwell residual in the record; no assumption that larger current beta always helps. |

## Search and coverage limits

The scholarqa-research evidence workflow was attempted with one exact-topic query. Both snippet and paper search returned Semantic Scholar HTTP 429. `literature_scholarqa.json` preserves both errors and the empty candidate set; it is **not** a successful scholarly search or bibliographic verification. The narrow baseline audit instead uses the directly acquired primary full text and its bibliography, as requested by the parent task. No broad prior-art or novelty conclusion follows.

The full supplied `report_sources/report_body.md`, `memos/agent_cie.md`, and `deliverables/CODEX_FBG_R3_HANDOFF_ZH.txt` were read. Their explicit warning is preserved: R2's 100-iteration current prototype is not the mature CIE-I/TSOM performance ceiling. The new module restores the mature structural ingredients; it does not reproduce the papers' acquisition, 3D discretization, complex material, or computational platform.
