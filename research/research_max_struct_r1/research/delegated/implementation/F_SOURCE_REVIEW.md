# Independent source review of F (before prescribed run)

Reviewed `src/tensor_families.py` and its adapter, without editing the F source.
This is a source review, not an additional physical experiment.

The finite-DDA physical self formulation is consistent: E is the internal
material field, j=(epsilon-I)E, and E=Einc+G j includes the fixed depolarization
and RR self terms. Thus V Eold=0 preserves old j and E under the same discrete
operator. H and every old-source external near/far dipole observation follow
from unchanged j. The passive-complex update uses q=E1 cross E2 and V=b q q^T;
the transpose dyad correctly implements the complex-symmetric, rather than
Hermitian, material null direction. The analytic real diagonal construction
has a direct continuum interpretation because the old incident electric
fields have zero z component while only epsilon_zz changes.

The finite-noise calculation uses the correct resolvent identity
`Rcur=I+X (I-GX)^-1 G` for `(I-XG)^-1`. Its ordering is valid:
`Lcur_t=(I-t V G Rcur)Lcur`, hence
`delta_j=t Rcur (I-t V G Rcur)^-1 V E0`. Frobenius bounds dominate operator
norms, so the denominator `1-T*||VG Rcur||_F` gives a conservative uniform
Neumann bound on that ray. This is ordinary floating arithmetic, not outward
rounded/interval certified. The weakest complex left singular vector is
conjugated before the symmetric dyad, correctly making q^T E a weak row.

Source mixtures use retained full factors with mixed RHS and charge four
additional triangular RHS. A unitary mixing preserves ideal sigma^2 I noise
covariance. The new-source SVD selects a normalized coherent combination of
the disclosed legal source library; it predicts separation and does not
acquire a new measurement. The ideal channel coefficient power metric is
explicitly distinguished from finite-antenna calibration.

Items communicated to the parent before running:

1. The refinement audit interpolates/extrapolates the coarse q field, then
   forms a new analytic-bump times q q^T extension. It is not direct linear
   interpolation of the coarse tensor V. Fine center coordinates ±.105 lie
   outside coarse center coordinates ±.09, and `fill_value=None` extrapolates.
   This is a legitimate *declared* fixed extension, but its semantics must be
   explicit. Fine-grid re-nulled tensors are different materials and cannot
   establish continuum convergence of one tensor family.
2. Raw `PASS` was initially unconditional. Legal/null and finite-noise
   bound/both-fit support conditions need separate recorded booleans or
   assertions; completion is not scientific support by itself.
3. Per-case timers initially stop before saving NPZ/public artifacts, while
   the main summary timer includes I/O. Report timer scopes explicitly.
4. Fine h=.03 cube has N=512 and 1536 complex electric unknowns per source;
   the multiple dense matrices are feasible on the current CPU machine but
   are not a new efficient Maxwell method. Report actual resources.
5. H-relative error is a finite-DDA cell-center diagnostic. No gap/interface
   H(curl) metric or compatible Hodge complex is supplied by this code.
6. The parent already identified the all-six-sources electric DOF metadata
   typo: complex electric states have 6*N*6 real scalar DOF, not 12*N*6.

No fatal issue was found in the finite-model null update or noise-ray algebra.
Final claims depend on recorded outcomes, lawful material margins, and the
scope limits above; the source review does not prejudge those results.
