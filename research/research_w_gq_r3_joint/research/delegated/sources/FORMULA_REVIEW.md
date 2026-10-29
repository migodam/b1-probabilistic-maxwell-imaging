# R3 compendium algebra/source review

Reviewed: `FORMULA_COMPENDIUM_ZH.md`, read 2026-10-29. Scope F18–F31 plus source/notation/render checks. No PDE or numerical formula tests run. Parent Astra retains final scientific judgment. No main-file edits.

## Corrections needed

### 1. F23–F27: independent prior blocks must be an explicit assumption

F23 initially permits a general precision `Λ`, but then sets `C=AᵀB`, `H_w=AᵀA+Λ_W`, `H_x=BᵀB+Λ_X`. This block expression is only valid when `Λ=diag(Λ_W,Λ_X)` (and any step damping also has the declared compatible block structure). For a correlated wall/target prior, `C=AᵀB+Λ_WX`. Consequently F25's simple `S_X=BᵀQ_A B+Λ_X` and F26–27 as an error bound for that Schur object need the independent-block assumption.

Smallest correction: state before F23 that the present formulas assume block-diagonal fixed prior precision. General correlated priors retain the standard F24 block Schur identity with the true cross block; do not reuse F25 unchanged.

### 2. F20–F22: declare how the receiver whitening and realification enter

As written, `S_z=R L(z)^{-1}D^{1/2}` is the raw complex receiver map and its SVD/norms use raw complex Euclidean receiver data. F23 then announces real whitened/scaled coordinates. These are different metrics unless receiver noise is proportional to identity and scales are consistently absorbed. Therefore F22 as currently written is a valid **raw-data** operator bound but does not directly certify the whitened GN/Schur objects F23 onward.

Smallest correction: declare that `R` in F20 already denotes the fixed data weighting and consistently carry the realification convention, or define `S_z` as the real whitened map on realified forcing coordinates. `K_z` and the block projection must be changed consistently. When retaining a complex SVD, require a compatible complex-linear noise whitening (e.g. stated circular complex noise); arbitrary real/imaginary covariance cannot silently be represented by the same complex `R`.

Also explicitly include parameter scale normalization in `K_z` if F22 is used for the scaled Jacobian of F23. The positive fixed `D` construction itself is algebraically sound; it defines a load metric and need not be claimed a canonical exact FEM dual norm.

### 3. F12 source alignment: state residual normalization differs from the cited book

Chen Eq.6.64 (PDF172 / printed153) uses state normalization `||j_s^+||²`, not `||χ e_b,s||²`. The latter is material dependent during optimization and conflicts with “fixed normalization coefficients” in the introduction to F12. This is not a harmless transcription if it changes the objective being minimized.

Smallest correction: replace the state denominator with the declared nonzero fixed `||j_s^+||²` (and specify a floor if a source is zero), or label the present denominator as a different chosen objective and freeze its value outside the coefficient/material subproblem. Do not imply exact native SOM objective reproduction with the current expression.

### 4. Common active face and inverse assumptions should be made operational

F29 is correct for unconstrained systems. On a common active face it is correct **after projection to free coordinates**, with the same fixed active coordinate values and the reduced Hessian/gradient. The full `e` of F28 may have nonzero active components at an exact constrained optimum, so inserting full `e_w,e_x` into F29 is not a constrained identity. Existing prose is close but should explicitly say all blocks in F29–30 are reduced/free-coordinate blocks.

F24 and F29–30 require invertible `H_w,S_X`; F30 specifically needs `S_X≻0`, not merely its existence. `Λ_W≻0` guarantees wall invertibility and `Q_A≻0`, but does not alone guarantee target `S_X≻0` when `Λ_X` is singular and `B` lacks target rank. State the condition or use an explicitly declared damping/identified subspace; do not write a positive denominator without it.

## Typos / small clarifications

- Remove the stray literal `+` after “`H_x=B^TB+Λ_X`：” immediately before F24.
- F31 uses `μ` both for prior mean (a vector in F23/F28) and the strong-convexity scalar `H≽μI`. Rename the latter `m_H` or `α_H`.
- F31's gradient perturbation bound assumes the exact same current residual `r`, prior, and step damping in both QPs. If a surrogate also approximates the nonlinear forward residual, an additional residual error term is required. The intended chart appears to share exact current `r`; say so once.
- F20 prose “每频一次 receiver adjoint” could be read as one RHS total. Use “每频一组接收伴随（按实际接收通道计 RHS）” if constructing the full receiver map.
- If discussing lossy homogeneous background in F05–06, replace small-phase conditions by magnitude-aware complex-wavenumber statements or restrict this explanatory far-field paragraph to real `k_b>0`. The displayed formulas' outgoing signs match F01's declared time convention.

## Formula checks that passed under the stated assumptions

- F18 both resolvent identities have the correct sign/order for `L1=L0−ΔV`; field-response migration is not a material-coordinate transform.
- F19 follows from `||(I−G0ΔV)^{-1}||≤(1−η)^{-1}` and is correct for `η<1`.
- F20–21 factorization and chart error identity are correct for fixed positive `D` in a consistently declared metric.
- F22 follows from splitting `S_z=S_0+(S_z−S_0)` and `||S_0(I−P_r)||=σ_{r+1}(S_0)`; multiplication order is correct.
- F25 Woodbury identity and `0≺Q_A≼I` are correct for `Λ_W≻0`.
- F26 inverse-difference expansion and F27 target-Schur perturbation bound are correct with fixed independent-block prior; no eigenspace gap assumption is needed for these norm bounds.
- F28 normal residual and F29 sign are correct in the unconstrained or correctly reduced system.
- F31's strong-convexity/variational-inequality bound permits active-set changes and is algebraically correct with the same convex feasible set and exact minimizers; no requirement that the active set be unchanged.
- F09–F11 align with the book's projected Eq6.71. `B(B*B)^†B*` is the column-space orthogonal projector, including rank-deficient `B`. “SOM-inspired material-forcing chart” correctly avoids a native TSOM reproduction claim.

## Optional issue outside assigned F18–F31 range

F34 uses the minimized `c(v)` but states the nonlinear remainder along a generic `(h,v)`. Its displayed upper bound holds when `h` is a minimizer attaining `c(v)` (or with the actual `||Ah+Bv||` instead), and `t≥0` as written; otherwise use `|t|`. Also the scaled segment must remain feasible. This is a wording/quantifier correction worth making while editing.
