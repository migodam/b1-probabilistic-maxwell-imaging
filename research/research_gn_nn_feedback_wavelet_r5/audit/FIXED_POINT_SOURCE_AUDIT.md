# FP-T / inverse fixed-point source map

Read-only extraction, 2026-11-08. No experiment, theory replacement, architecture choice, or source edit was performed. All paths below are relative to `/Volumes/migodam's-external-brain/Research/Inv_SLAM`.

## Provenance verdict

**Concrete earlier local inverse fixed-point and transverse-corrector assets are accessible. The original independent FP-T report is not authenticated by these assets.**

`Gaussian/B1/b1-p3l-r12-bilateral-fixedpoint/research/p3l_r12_bilateral_fixedpoint/INPUT_AND_HISTORY.md:5` explicitly records that the referenced ChatGPT session supplied user questions but assistant placeholders, not a complete FP-T theory report; the R12 filename search did not find an independent FP-T-R1 file; its mathematics was derived from the task equations and existing code. `HIERARCHICAL_ALGORITHM.md:47` repeats that a named FP-T source was not found. These limitations must accompany any use of the following R12 work as prior local theory.

The later main-repository `research/research_z_alg_r4_basin/FORMULA_AND_CODE_REUSE.md:11,34-35` labels the R12 files as FP-T assets and `ACCESSIBLE_READ`. That corroborates availability of R12, **not identity with the missing original FP-T report**.

## Exact sources and useful line ranges

Let `B = Gaussian/B1/b1-p3l-r12-bilateral-fixedpoint/research/p3l_r12_bilateral_fixedpoint/`.

| File | Exact relevant locations | Evidence role |
|---|---|---|
| `B/FIXED_POINT_SPECTRUM_AUDIT.md` | 7-52; 54-77; 79-118; 120-180 | Original task-specified inverse map, full derivative, exact-root spectrum, nonzero residual correction, material metric, algorithmic unit modes and bounded numerical controls |
| `B/TRANSVERSE_CONTINUATION.md` | 7-40; 42-74; 76-99; 114-130 | Fixed weak-coordinate predictor/corrector, actual RML/native regularizers, version boundaries and qualification limits |
| `B/agent_b/engine.py` | 77-118; 225-380; 540-569 | Implemented canonical material lift, weak frame and C; refreshed-J constrained GN; explicitly supplied-map derivative compression |
| `B/INPUT_AND_HISTORY.md` | 3-11 | Prior provenance and physical/prior boundaries; missing original report acknowledged |
| `B/HIERARCHICAL_ALGORITHM.md` | 47 | Explicit missing FP-T-specific theory-source statement |
| `Gaussian/B1/b1-p3l-r12-bilateral-fixedpoint/research/maxwell_formula_atlas/SOURCE_FORMULA_INDEX.csv` | 97-99 | Stable formula IDs F-BF-006, F-BF-007, F-BF-008 map to the above report/code; says R12-added finite-model derivation, standard identity not novelty |
| `research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_basin/FORMULA_AND_CODE_REUSE.md` | 11; 34-35; 67 | Later acknowledgment of reading R12; adoption as local limiting theorem/optional standard corrector only |
| `research/delegated/z_alg_r2_publication/repository/research/research_gn_nn_preservation_r4/theory/A17_RECURRENT_FEEDBACK.md` | 138-146; 148-176; 178-216; 236-297 | Current R4's distinct current-profile feedback and limited fixed-point claims |

The B shorthand expands literally to the path above; it does not designate a new or inferred repository.

## Earlier inverse fixed-point equations (F-BF-006)

R12 assumes real parameters x; a differentiable realified, whitened residual r(x); J=Dr; fixed positive damping mu. Its map is

\[
 A=J^TJ+\mu I,\quad g=J^Tr,\quad T(x)=x-A^{-1}g.
\]

For h and dotJ=DJ[h], the report explicitly retains both curvature/inverse-derivative terms:

\[
 DT[h]=h-A^{-1}(J^TJh+\dot J^Tr)
 +A^{-1}(\dot J^TJ+J^T\dot J)A^{-1}J^Tr.
\]

Only at an **exact response root** r(x*)=0 does this reduce to

\[
 DT(x_*)=\mu(J^TJ+\mu I)^{-1},\qquad
 \lambda_i=\frac{\mu}{\sigma_i(J)^2+\mu}.
\]

This spectrum uses the same Euclidean parameter metric as mu I. Positive singular values yield eigenvalues below one; exact Jacobian null directions yield unit eigenvalues. Nearly unit eigenvalues may merely reflect weak sensitivity or strong damping. The report does not infer a global/continuous physical fiber from first-order nullity.

At a positive-residual stationary point J^Tr=0, r!=0,

\[
 DT=I-A^{-1}(J^TJ+S_r),\quad S_rh=(DJ[h])^Tr.
\]

The residual Hessian may change stability. The report expressly requires appropriate local regular-rank assumptions and exact-zero constraints for a continuous root manifold. Full column rank supports a locally isolated root; it does not forbid other isolated roots elsewhere.

## Damping versus RML/MAP (F-BF-007)

The damped data map has no prior gradient +x, so mu=1 does not make it MAP. The separate actual RML objective and map are

\[
 L_{\eta,\xi}=\tfrac12\|r-\eta\|^2+\tfrac12\|x-\xi\|^2,
\quad T_{RML}=x-(J^TJ+I)^{-1}\{J^T(r-\eta)+(x-\xi)\}.
\]

At its stationary point, DT_RML=-A^{-1}S_(r-eta). The report treats these as finite optimization candidates and explicitly rejects upgrading them to posterior samples/volume.

## Material metric, quotient and corrector (F-BF-008)

The P1 setup is Gaussianized66, with complex material measured on the full 576-cell grid by ||delta m||/24, including when the physical derivative uses an ROI256 reduction. With M=[Re D;Im D]/24 and retained SVD M=U Sigma V^T, define P=V_keep Sigma_keep^-1 and Q=U_keep. Then MP=Q and P^T M^T M P=I. Dropped decoder-zero directions are not asserted to be global physical equivalences.

Material orthonormalization changes the original damping to mu P^T P. Thus ordinary singular values of JP cannot simply be inserted into the Euclidean mu-I formula. The implemented effective derivative is Q^T M (DT) P, requiring DT to be explicitly supplied (`effective_spectrum`), never guessed from J.

At an exploration-frame reference W with W^T G W=I and G=M^TM, freeze C=W^TG throughout one corrector. Predict x_p, then maintain

\[
 C(x-x_p)=0,\quad Z=P(x)\,null(CP(x)),
\]

while refreshing the exact current J and material derivative at every GN step. The P1 subproblem is

\[
 \min_z \tfrac12\|r+JZz\|^2+\tfrac12\|x-\xi+Zz\|^2,
\]

with r replaced by r-eta for RML. Native-digit versions instead use their frozen non-Gaussian regularizer and explicit full-coordinate/box rules; these are not retroactively called P1. Backtracking accepts real objective descent and retains failures/costs.

If CT(x)=Cx then C DT=C, producing **algorithmic unit left modes**. Transverse contraction must be evaluated on the actual free domain ker C; these frozen coordinates do not establish Maxwell nullity or root-manifold dimension.

## What current R4 contains, and what it does not

R4 §5 defines a shared-weight learned composite T_theta on material and multisource currents. It says a declared invariant neighborhood/norm and sup||DT_theta||<1 would be needed for full local contraction; finite trained depth and decreasing updates are insufficient.

R4 §5.1 instead gives a deterministic current-profile objective

\[
 \Psi(j)=\|Sj-Y\|^2/(2a)+\rho\min_{m\in B}\|j-D_m(E+Gj)\|^2/(2b_0).
\]

Assumptions include positive row total-field energy (unique bounded row LS), positive fixed weights and compatible complex-real Hilbert convention. Material readout plus fixed-material exact line feedback plus rereadout decreases Psi. The NN proposal has no such guarantee.

R4 §5.2 assumes an interior physical solution, invertible L*, nonzero row field energy, and exact data fit. With N=(I-Pi_Q^R)L*, it gives

\[
 H_{prof}=S^*S/a+\rho N^*N/b_0,
\quad ker H_{prof}=ker S\cap ker N
 =\{T_{m_*}h:DF(m_*)h=0\}.
\]

For the **different fixed-step profile-gradient map**, DT=I-tau H_prof; with a declared positive nonzero spectral gap and constant rank, 0<tau<2/lambda_max yields linear transverse contraction and neutral kernel directions. The implemented direction-dependent exact line step does not automatically have this single Jacobian at zero gradient. R4 explicitly acknowledges this limitation. Its later noise hinge preserves qualifying fixed points for the physical subblock, not automatically for NN composition.

Therefore R4 contains a related local physical/profile-gradient argument and explicit cautions, but **does not reproduce the older inverse damped GN map/spectrum or authenticate original FP-T provenance**. It is inappropriate to substitute current-only damping, current energy decrease, or R4's profile-Hessian formula for the original inverse-map theorem. Whether/how to connect them is left to the parent scientific review.

## Historical evidence versus theorem scope

R12 reports two isolated two-real-parameter Maxwell roots with near-unit eigenvalue .9965852 despite full-rank J; a nonzero-residual control shows the omitted residual curvature changes DT. A separately declared three-real-DOF/two-real-data two-cell model has an explicit one-dimensional passive equal-response curve, with denominator/state invertibility bounds on a small interval. Those controls do not establish a fiber in the original Gaussianized66/P1 model, continuous PDE, 576-cell model, full material coverage or posterior calibration. The report labels these broader claims NOT_ESTABLISHED. These are extracted historical reports; no new run or independent numerical revalidation was performed.

## Search boundaries and unavailable assets

- Read applicable root user instructions, main repository AGENTS.md and formula-atlas AGENTS.md. Wrote only this new delegated output directory.
- Narrow searches covered main-repository relevant theory/handoff text, actual Gaussian/B1 source paths and atlas indices, and scoped Gaussian/A17 references. An initial research path inventory located existing atlas copies; authoritative source chosen here is the actual Gaussian/B1 R12 tree.
- Exact `/Users/migodam/Downloads/Research_FP_T_Report_ZH.pdf` is absent. A filename inventory of Downloads for FP-T/fixed-point/transverse also returned no hit. This does not prove absence from every possible external location.
- Original independent FP-T-R1 report, original full ChatGPT assistant derivation and a current download attachment path remain unavailable in the inspected scope. No whole-disk or unrelated-session scan was performed.
- Memory registry search yielded no relevant FP-T hit; this report relies on inspected project files, not memory assertions.
