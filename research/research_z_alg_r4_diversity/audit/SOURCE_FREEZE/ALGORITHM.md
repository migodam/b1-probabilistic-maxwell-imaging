# R4 material diversity: four common-access solvers

`src.solver.solve(public, method, config=None, backend='cpu')` accepts a public
observation/operator dictionary or NPZ path. The methods are `native_trf`,
`reciprocal_trf`, `data_lbfgs`, and `globaldual`; CUDA is selected explicitly with
`backend='cuda'`. The solver reads observations and acquisition operators. It does
not import a truth file or the material generator.

## Common physical problem and access

Every forward uses the same full 576-cell disk-pulse scalar-TM VIE,
`j = diag(m)(E0 + G j)`, `F(m) = S j`, with 12 sources and 48 receivers.
The physical operators never become a coarse PDE. The material basis has
8×8, 12×12, and 24×24 coefficients; the last basis is the 576-dimensional
identity. Bounds are `[.02,20]`, with the single public initial material `.1`.
Coarse materials are prolonged by the frozen R3-DIGITS bilinear basis.
There are no truth-dependent starts, support masks, or conditional restarts.

The common score is

`L(m) = .5 ||F(m)-Y||² / ||Y||² + R(m)`.

`R` is the frozen R3-DIGITS regularizer:
`1e-4 mean(d²) + 1e-4 mean(sqrt(d²+.01²)-.01)`, where `d=Dm` comprises
horizontal/vertical neighbor differences divided by the original fixed scale
19.5. This normalization remains identical across methods even though the new
material box begins at .02. Both regularizer terms and their gradients are
included in every optimization phase and in final ranking.

## Methods

Native TRF uses scipy `least_squares(method='trf', tr_solver='exact',
x_scale='jac')`, with the analytic Full-J. Its regularization residuals square
to `2R`; thus the scipy half-squared residual objective equals `L`.

Reciprocal TRF uses the same solver and residuals with coefficient variable
`eta=1/c`, bounds `[1/20,1/.02]`, and `m=Bc`. Its exact residual Jacobian applies
`dm/deta = B diag(-c²)`. At the final identity basis this is the cellwise chain
rule `dm_i/deta_i=-m_i²`.

Data L-BFGS and GlobalDual share the same original-data L-BFGS-B optimization
at the 8 and 12 stages. L-BFGS-B optimizes the normalized box variable
`z=(c-.02)/19.98`; gradients apply the corresponding chain factor 19.98.
At the 24 stage both methods reset optimizer memory at the same two-phase
boundary. The first phase receives 45% and the second phase 55% of the final
stage's allocated wall/LU/RHS budget. Data L-BFGS uses `L` in both phases.
GlobalDual uses the source-field small-Gram objective plus `R` in phase one,
then returns to the original `L` in phase two. Each second phase starts from
the best original-score physical state visited in its own first phase.

For GlobalDual, put `Q=E0+Gj`, `Rdata=F(m)-Y`, and `H=QᴴQ`.
The proposal objective is

`alpha tr(Rdata H^-1 Rdataᴴ) + R(m)`,

where `alpha=.5 ||Q(.1)||_F²/(12 ||Y||_F²)` uses the shared initial physical
forward. The implementation reuses the mature R3 small-Gram/source-adjoint
formula. With `X=Rdata H^-1` and `C=XᴴX`, solve
`Aᴴ adj = SᴴX - Gᴴ(QC)` using the retained physical LU. Its gradient before
`alpha` is `2 Re sum_sources(conj(adj)*Q)`.
The existing thin-QR fallback is retained if the 12×12 Gram solve fails its
conditioning/residual criterion. No jitter changes the objective.

## Budgets, measurements, and candidate meaning

Defaults are 45 seconds, 2,000 full LU factorizations, and 60,000 full-system
complex RHS. The three stages have cumulative budget endpoints `.20/.45/1.00` for
wall time, LU, and full-system RHS. Unused allocation after early convergence
therefore remains available to later stages. At the final stage the first
L-BFGS phase has an endpoint 45% into the remaining stage budget; the second
phase ends at the common global endpoint and receives any unused first-phase
allocation. The same transfer rule applies to Data L-BFGS and GlobalDual.
Global limits remain active throughout. Guards run
before each counted physical block, following R3/R4 practice. An in-flight
block or opaque scipy optimizer operation can pass a deadline; measured budget
overshoots are reported, rather than hidden.

A full source forward factors one 576×576 complex matrix and solves 12 RHS.
An analytic Full-J reuses that LU and solves 48 receiver-transpose RHS.
An original-data or small-Gram source-adjoint gradient reuses it and solves
12 RHS. The local Gram factorization and its 48 right-hand sides are counted
separately as `Gram_factorizations` and `Gram_RHS`; they are not falsely added
to the count of 576-dimensional physical RHS. `costcounts` and `accounting`
report these quantities, operation wall/CPU times, and convergence or budget
reasons. Scipy TRF's internal dense SVD time is included in total wall/CPU time;
its opaque SVD call count is not invented.

Every completed physical trial is ranked with the original `L`. `m` is stored
once at the result root for the best trial and `top1` holds its scalar metadata.
The finite candidate list contains at most 16 phase-best material vectors.
The scalar `trace` is an evaluation history and never a count of distinct roots,
branches, or basin recovery. Candidate count is explicitly not branch count.

The default result is JSON-serializable. With
`retain_material_trajectory=True`, only `trajectory_arrays` additionally holds
NumPy arrays `m` and `seconds`; the runner should pop and save these as NPZ.
The JSON trace does not repeat material vectors or duplicate the trajectory.

## CPU and CUDA implementation provenance

`src/physics.py` mechanically reuses `Model`, operator construction, bilinear
basis, and `Regularizer` from the frozen R3-DIGITS source
`research/research_z_alg_r3_digits/src/{physics,multiscale}.py`. It adds the
mature R4 source-adjoint data gradient and a field-norm accessor. The original
R3/R4 files are unchanged. `src/solver.py` reuses the mature R4 analytic TRF
residual/Jacobian and matched source-adjoint/GlobalDual formulas.

`src/cuda_physics.py` retains the same operator values and uses torch CUDA
complex128 LU, source solves, receiver-transpose solves, source-adjoint solves,
and Gram/QR operations. Scipy optimizers and the common regularizer stay on
the CPU. CUDA synchronizes at solve entry, at the end of every counted physical
operation, and before reporting the final wall time. No float32 path, CUDA
installation, or new runtime service is introduced.

Implementation readiness, end-to-end execution, physical validity, and
scientific performance are separate statuses. Compilation alone does not
establish any of the latter three.
