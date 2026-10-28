# Local same-material scaling: conditional statement and validation limits

This addendum separates an exact decoder identity, a local existence result under stated rank conditions, and finite floating-point witnesses. It does not establish a universal Full66 gauge. No Maxwell operation is involved.

## Original P1 formulas and the strict stratum

Write the original decoder as

\[
 G(z,u)=\chi_h+\delta(u)p(z),\qquad
 p_i(z)=M_i S(f_i(z)/\tau),\qquad f(z)=\mu+A z.
\]

Here `S` is the unchanged clipped quintic, `M` the original ROI mask, and the two native property coordinates `u` lie in `(0,1)^2`. Gaussianized property coordinates are `xi=ndtri(u)`. On `(-1,1)`, `S` is smooth and strictly increasing with strictly positive derivative. The original real and imaginary inclusion increments are positive and bounded away from zero.

Let `T={i in ROI: -tau<f_i<tau}`, `H={i in ROI: f_i<-tau}`, and `P={i in ROI: f_i>tau}`. The following statements concern a strict stratum: no ROI field equals either seam. Seam examples require separate one-sided analysis.

## Exact inactive subspace, within a finite plateau margin

For any real `v` satisfying `A_T v=0`, and fixed properties, all transition fields are unchanged along `z(t)=z0+t v`. Plateau cells remain on their original plateau while their strict inequalities persist. Thus `G(z(t),u0)=G(z0,u0)` exactly in real arithmetic throughout that interval.

A conservative symmetric open interval is

\[
 |t|<\min_{i\in H\cup P:\,(Av)_i\ne0}
 \frac{|f_i|-\tau}{|(Av)_i|}.
\]

An empty minimum means no plateau restriction. At a seam this interval can be zero. The saved SVD basis only approximates `ker(A_T)` numerically; its residual and finite material differences are measured. A small singular value by itself is not an exact algebraic zero.

## Conditional local radial fiber: explicit construction / IFT

Assume:

1. There are no pure-phase cells and no seams: every ROI cell is in `T` or `H`.
2. Native properties are strictly interior.
3. `A_T` has full row rank. In particular this is possible when `0<|T|<=64`; the requested `|T|<64` case is included.

Choose any exact real right inverse `R` with `A_T R=I`. Let `p_T=S(f_T/tau)`, whose entries lie strictly between zero and one. For positive `alpha` sufficiently close to one, define

\[
 \delta_\alpha=\alpha\delta_0,\qquad
 f_T(\alpha)=\tau S^{-1}(p_T/\alpha),\qquad
 z_\alpha=z_0+R\bigl[f_T(\alpha)-f_T(1)\bigr].
\]

The two native properties are the unique affine coordinates giving the real and imaginary parts of `delta_alpha`. They remain interior for `alpha` close to one. The inverse quintic is well-defined and smooth because `p_T/alpha` stays in `(0,1)`. By the right-inverse identity, `A_T z_alpha+mu_T=f_T(alpha)` exactly. Therefore `delta_alpha p_T(z_alpha)=delta_0 p_T(z0)` on every transition cell. Host cells have positive initial plateau margins; finitely many continuous inequalities remain strict on some neighborhood of `alpha=1`. Outside the ROI the mask is zero. Hence `G(z_alpha,u_alpha)=G(z0,u0)` exactly on some open interval about one. Gaussianized property coordinates remain finite and smooth there.

This is also an immediate local implicit-function result: the derivative of the transition field constraints with respect to `z` is the surjective matrix `A_T`. The explicit right-inverse formula makes the local construction concrete. Adding sufficiently small `ker(A_T)` displacements gives additional local inactive freedom.

For nonempty `T` under these assumptions, `B=D_z p` has row rank `|T|`, `p` lies in `range(B)`, the real material derivative has rank `|T|+1`, and the local fiber dimension is `65-|T|`: `64-|T|` inactive spatial directions plus one coupled amplitude direction. The all-host case `T` empty is separate: the decoder is locally constant and its fiber dimension is 66.

## A sufficient interval formula, and what the numeric bounds mean

Property legality and the transition-fraction inversion require

\[
 \alpha\in(L,U),\quad
 L=\max\{0,a_r/\delta_r,a_i/\delta_i,\max_j p_j\},\quad
 U=\min\{b_r/\delta_r,b_i/\delta_i\}.
\]

These are necessary open bounds for this construction. They do not alone enforce the host plateaus.

For an interval `[1-eta,1+eta]` inside `(L,U)`, let `l_j` bound

\[
 |f'_{T,j}(\alpha)|=
 \frac{\tau p_j}{\alpha^2 S'(S^{-1}(p_j/\alpha))}.
\]

Then a sufficient host condition is

\[
 \eta\sum_j |(A_H R)_{ij}|l_j < -\tau-f_{H,i}(1)
 \quad\text{for every }i\in H.
\]

The derivative of the quintic is unimodal in its transition interval, so its minimum over each inverse-fraction interval occurs at an endpoint; this supplies the `l_j` bound. `E1_INDEPENDENT_VALIDATION.json` evaluates these expressions in float64 with a factor-two host-margin reserve. The mathematical implication above is symbolic; the reported numerical endpoints are diagnostics, **not** rigorous outward-rounded interval-arithmetic certificates.

## Why this does not become a universal rescaling gauge

A pure-phase cell has `p_i=1` and `B_i=0`. A nonzero radial property tangent `d delta=c delta` changes that material cell and cannot be compensated infinitesimally by its saturated spatial derivative. More generally, for nonzero fraction and interior properties,

\[
 \operatorname{rank}_{\mathbb R}DG
 =\operatorname{rank} B+1+\mathbf1\{p\notin\operatorname{range}B\}.
\]

Thus `rank(B)=64` and `p` outside its range give rank66. A rank66 point has no nontrivial local same-material fiber (select an invertible 66-output minor and apply the inverse-function theorem). This does not exclude remote disconnected preimages. Conversely, clipped priors contain open saturation strata with different ranks, so “generic no gauge” must be qualified by the stratum/conditions; it is not an almost-sure theorem for this entire prior. The observed 191/256 machine-scale full-column-rank `A_T` cases are finite empirical evidence, not a proof of a population statement.

Uniform transition controls can satisfy extra symmetry with `|T|>64`; their fiber does not establish the full-row-rank condition or the same gauge at nonuniform parents.

## Finite witnesses retained separately

The four sampled no-pure parents 092, 201, 233, and 251 have 19, 40, 10, and 14 transition cells. Their archived corrected curves include all four prescribed tangent steps `0.001,0.01,0.1,0.5`, complete Full66 `x`, complete 576-cell complex material, legal native properties, and residuals. Original arrays remain at `raw/decoder/NONLINEAR_EQUAL_MATERIAL_CORRECTIONS.npz`; the explicit 16-point subset is copied to `audit/E1_FOUR_RADIAL_CURVES.npz`.

Independent validation redecodes the saved coordinates with the original P1 prior, checks property legality, checks the conditional rank assumptions numerically, and compares NumPy spectra with a separate SciPy `gesvd` decomposition. These finite near-equalities support the stated cases. They are not global exact identities at floating-point coordinates, nor a symbolic certification of the numerical rank decisions.
