# Preconditioned CSI (2000): targeted factual extraction

Status: **PARTIAL / formula-image access gap**. No originality verdict. No GT files read; no physics executed.

Identity: A. Abubakar and P. M. van den Berg, *Three-Dimensional Inverse Scattering Applied to Cross-Well Induction Sensors*, IEEE TGRS 38(4), 1669–1681 (July 2000), DOI [10.1109/36.851966](https://doi.org/10.1109/36.851966). Identity confirmed through [Crossref](https://api.crossref.org/works/10.1109/36.851966) and [TU Delft metadata](https://research.tudelft.nl/en/publications/three-dimensional-inverse-scattering-applied-to-cross-well-induct/). Primary full text: [TU Delft PDF](https://repository.tudelft.nl/file/File_e636219e-fe7f-4c5c-b870-558ee7dacd45).

| Item | Page / equation | Confirmed fact or explicit gap |
|---|---|---|
| Forward preconditioner | 1671, (21) | Extended-Born approximate inverse; a 3×3 matrix inverted separately at each position. Exact symbols unavailable in extracted text; not reconstructed. |
| Inverse CSI objective | 1673, (30) | Preconditioning enters the CSI cost functional itself, not merely a forward-solver acceleration. Exact norms/normalizers unavailable. |
| Material dependence | 1673, paragraph after (33) | (21) is recomputed per iteration because conductivity contrast changes; its underlying operator is computed once. |
| Second preconditioner | 1673, (33) | Computed once; described as approximate inverse of a positive-definite operator version. Exact formula unavailable. |
| Locality | 1671, after (21) | The 3×3 inversion is local in position, coupling vector components. This is not evidence for a retained off-diagonal spatial-cell block. |
| Appendix | 1679–1680, (54)–(62) | Explicit component formulas are located here; not visually verified. |
| Material gradient | 1674 onward | Whether the material update differentiates (21), freezes it, or uses a special algebraic identity is **UNRESOLVED**. No omission claim is warranted. |

R4 comparison uses only `THEORY_DERIVATIONS.md` §5: R4 explicitly mixes four distinct spatial cells through off-diagonal Green entries, uses one-exchange `R_one=(I+L)u`, and differentiates its material-dependent transport. The confirmed 2000 pointwise vector tensor differs from that spatial-cell coupling. However, this distinction alone establishes neither novelty nor superiority: material-dependent preconditioning inside a CSI penalty is already directly present in this paper. Exact residual-transform equivalence remains unassessed until equations are visually accessible.

Do not turn missing formula access into a negative prior-art finding. Local copies, if acquired later, belong under `local_only/`; no third-party full text is copied into this report.
