# R4 finite-budget rational-valley accessibility candidate

The solver is a public-input-only algorithm candidate for the frozen finite 576-cell scalar-TM VIE. It imports the R3 `Model`, material interpolation, and `Regularizer` directly. Geometry, operators, measurements, and material bounds are unchanged. It does not read truth, digit labels, masks, support, historical endpoints, or offline diagnostic files. Multiple corrected candidates are not certified physical branches or probability samples.

## Interface and frozen common objective

```python
from research.research_z_alg_r4_basin.src.basin_solver import solve
result = solve(public, method, config, initial=None)
```

`public` is a dictionary or public NPZ path; the R3 public loader restricts readable keys. Operators may be merged into the dictionary, or supplied using the CLI `--operators`. The V1 method names are `data_lbfgs`, `globaldual`, `trf`, `basin_bridge`, and `ray_best_only`; V2 adds the three uniform methods described below. `initial=None` means uniform material 1. An explicit initial physical material is supported and must be disclosed when used.

All methods use bounds `[0.5,20]`, the fixed full physical grid, coefficient stages `8→12→24`, and exactly the same ranking objective

`0.5 * ||F(m)-Y||² / ||Y||² + smooth(m) + TV_epsilon(m)`.

Defaults are `smooth_weight=tv_weight=1e-4`, `tv_epsilon=.01`; the imported R3 regularizer defines the exact edge normalization. Stage wall fractions are `.20,.25,.55`. The global `wall_seconds` and `max_lu` caps apply to every physical, derivative, scan, and optimizer operation; an individual dense LAPACK operation is atomic and can cross a wall deadline. Overshoot is reported. `max_LU` and `maxLU` are accepted aliases for `max_lu`. The same frozen config must be supplied to each comparison arm.

The result contains `top1`, `m`, `candidates`, `trajectory`/`material_trajectory` (every full physical forward, including optimizer trials, with complete material vectors and objective), `raw_proposals`, `rays`, `events`, `costcounts`, and `accounting`. `all_unknowns_released` distinguishes actual final-stage execution. JSON conversion preserves arrays as lists.

## Methods

`data_lbfgs` runs bounded data-objective L-BFGS-B at all three stages. Its exact material gradient uses a source adjoint of the retained LU. `globaldual` shares its first two data-L-BFGS stages, then spends `.45` of its final stage on the frozen small-Gram objective and corrects its physically best visited proposal using the original objective. Its fixed dual scaling is computed from the paid uniform-m1 forward. It adds no measurement.

`trf` uses SciPy `least_squares(method='trf', tr_solver='exact')` and analytic Full-J at each stage. The least-squares residual concatenates normalized complex physical data, exact quadratic smoothness residuals, and exact smoothed-TV residuals. The TRF arm uses TRF at every stage; the other four arms share data-L-BFGS coarse stages. They share public m1 and priors, while their final-stage seeds need not coincide with TRF. Therefore its objective is the identical common objective, rather than a different quadratic surrogate for TV. SciPy's dense internal SVD/step linear algebra is included in correction wall/CPU time and identified as opaque-count work in events; physical LU/Full-J operations are separately counted.

`basin_bridge` shares the first two data stages. At the final full-material stage it first uses `.25` of the stage's remaining wall allocation for a local data correction, creating a public current anchor. It uses `.25` of the remaining time as a ray-scan deadline, retaining time for candidate correction. Default ray directions are:

1. Amplitude/shape amplification `p=m−1`, using only the public current material.
2. A regularized damped Gauss–Newton direction from the current analytic physical Full-J and measured residual, projected once into the material box to define a feasible straight direction.

Optional explicitly configured directions `shape=m−mean(m)`, uniform `offset`, data gradient, and small-Gram gradient are implemented; they must not be selected after inspecting holdout outcomes. A zero direction produces no ray. Every direction starts at the same anchor, with an explicitly paid fresh forward to retain its matching physical LU.

For each direction the code computes the **entire closed box-feasible real interval**, rather than a positive backtracking interval. It constructs

`K=M^{-1} diag(p) G`, `V=M^{-1} diag(p) Q`,

then a full complex Schur decomposition `K=Z T Zᴴ`. For every scan/refinement value the response is the exact finite-dimensional transfer

`F(t)=F0+t (S Z) (I−tT)^{-1} (Zᴴ V)`.

No rank truncation, Born approximation, Taylor step, or surrogate physical data is used. The formula requires the resolvent to exist; sampled points are finite floating-point evaluations and do not certify an unsampled interval free of singularities.

The default scan has 33 uniformly spaced values over the whole feasible interval. Every discrete valley (including a one-sided endpoint valley) is retained; interior valleys receive up to 10 bounded one-dimensional refinement iterations. All scan nodes, refined nodes, anchor materials, directions, intervals, raw scores, and minima are recorded. Finite sampling can miss a narrow valley; this is a proposal mechanism, not certified one-dimensional global optimization.

Current V2 material proposals are deduplicated before correction only at RMS `1e-7`; the archived V1 snapshot used `.25`. Distinct minima are enumerated in direction/grid order, **without raw-score sorting or raw-score top-k rejection**. Every retained distinct valley receives an equal share of the then-remaining correction time. Each correction runs the full original physical data objective with bounded L-BFGS-B; full physical evaluations determine all candidate scores. If a scan reaches its allotted deadline, its partial node record survives and already completed rays proceed to correction. A global wall/LU stop may leave explicitly uncorrected proposals. The final top1 is the minimum common objective over completed physical correction evaluations and the common local anchor. Candidate deduplication after correction also uses RMS `.25`.

`ray_best_only` is the mechanism ablation: it performs the identical anchor and exact scans but corrects only the raw-score-best proposed material. Discarded minima remain in the output with an explicit ablation selection status. The contrast with `basin_bridge` isolates the finite-budget decision to correct separated valleys before ranking their endpoints.

## Paid work and boundaries

A dense ray on 576 cells costs **576 full complex kernel RHS plus 12 illumination RHS**, in addition to its anchor forward. These are included in `actual_complex_RHS`, with separate `ray_kernel_RHS` and `ray_illumination_RHS`. A full 576-dimensional complex Schur factorization is counted and timed. Every scan/refinement transfer costs a 576-dimensional triangular system with 12 RHS; it is recorded separately as transfer linear algebra, not mislabeled as a free lookup or a full physical LU. All matrix products and factorization setup are inside timed operations.

The common forward, source adjoint, analytic Full-J receiver solves, small-Gram work, GN direction normal factorization, dense kernel, Schur, triangular transfer, and optimizer correction wall/CPU time are exposed. No hashes, derivative finite-difference checks, duplicate formula tests, smoke runs, or oracle initialization were added. Syntax compilation is the only pre-experiment executable check for each implementation version; actual performance and scientific acceptance require its authorized comparison execution.

Neither exact rational algebra nor a lower data objective establishes GT recovery, entry into every basin, a continuous-Maxwell theorem, or a new physical nonconvergence theorem. The experiment must separately assess imaging error and finite-budget accessibility, preserve failed/uncorrected proposals, and report the complete cost.

# V2 pivot: public uniform-amplitude multistart

V1 is preserved in its experiment source snapshot and original outputs. Its incumbent-based amplitude/GN rays did not reveal separated improving valleys on the first holdout7 comparison. V2 adds `uniform_bridge`, `uniform_best_only`, and `uniform_multistart`; it does not reinterpret the V1 records.

All V2 methods have the same global defaults: 60 seconds, 2,500 physical LU, and 60,000 full complex physical RHS. `max_rhs` is checked against actually paid `actual_complex_RHS` before each operation. A block already started can cross the cap; explicit RHS/LU/wall overshoot is returned. Kernel construction remains a real 588-RHS cost per target/run. Although the uniform kernel and Schur representation mathematically depend only on the common public operators and can be reused across targets, this implementation pays the cold cost in each independent solve and does not claim an unimplemented cross-target cache saving.

Each uniform method first runs the same public-m1 8→12→24 data-L-BFGS-B anchor, with a maximum 20% global wall allowance. Every stage still uses all 576 physical cells. A fixed empirical public-data gate permits branching only when the anchor relative data residual exceeds `1e-3` **and** the corrector has either projected gradient at or below its `gtol`, a relative-objective/line-search termination, or best-objective improvement no greater than `1e-6` over its last 10 full physical evaluations. Every component, optimizer message, threshold, and trigger decision is recorded. A wall deadline by itself is not the no-progress test. This empirical gate carries no basin certificate.

When the gate triggers, `uniform_bridge` constructs one exact full feasible rational ray from a fresh paid **uniform m1** state:

`m(t)=1+t·ones`, `t in [-0.5,19]`.

Its public direction is independent of the failed incumbent's digit shape or material support. The exact cold dense Schur response is scanned at 65 points and interior sampled valleys receive the same bounded refinement as V1. The method records every valley. Only floating-point-identical proposed materials, RMS below `1e-7`, are removed before correction; the `.25` RMS threshold applies only to corrected candidates. The same pre-correction duplicate threshold now applies to the legacy ray methods in V2; V1's earlier `.25` rule remains in its archived source.

Every retained uniform valley initializes a **new full 8→12→24 correction**, with bounded data L-BFGS-B and an exact retained-LU source-adjoint gradient. This is full nonlinear physical correction; it is not called a Full-J optimizer. Each seed receives an equal share of the then-remaining global time. An early-converged seed returns unused time automatically to subsequent seeds. All candidate scores use the original nonlinear data objective and the common `1e-4` smooth/TV priors. The final corrected-candidate deduplication remains RMS `.25`; all paid physical trajectories are retained.

`uniform_best_only` scans identically but corrects only its raw-score-best uniform valley, explicitly marking the discarded valleys. `uniform_multistart` uses the same anchor, gate, common budgets, and complete per-seed multiscale corrector, with the fixed ordinary seed list `[1,4,8,12,16]` and no rational-scan cost. Thus the ordinary multistart arm can spend the saved kernel/Schur time and RHS on correction. Both ablations are necessary for attributing any finite-budget benefit to data-conditioned uniform amplitude selection and correcting multiple valleys.

This candidate is **R10-style rational amplitude scanning plus physics/data-conditioned uniform multistart**. Ordinary multistart is direct prior art. Originality and any Maxwell-specific accessibility advantage remain to be established. The scan detects sampled/refined valleys; it does not certify global optimization, all valleys, GT inclusion, or new physical observations.

# V3 final freeze: resource-fair seeds and scaled L-BFGS variables

V2 source and actual records remain archived unchanged. V2's wall-only seed scheduling allowed a seed to consume the full RHS pool before other seeds reached the 24×24 parameter stage, so it did not provide a strong ordinary-multistart comparison. V3 fixes this implementation limitation before fresh evaluation. It does not discard V2's negative imaging evidence or upgrade a data-residual improvement to GT recovery.

All L-BFGS-B arms, including `data_lbfgs`, the data corrections, and `globaldual`, now optimize dimensionless coefficients `z=(c−0.5)/19.5` in `[0,1]`. They evaluate the original material `c=0.5+19.5z` and return the exact chain-rule gradient `19.5*g_c`. Bilinear interpolation preserves the box. Physical materials, nonlinear response, TV/smooth terms, data normalization, absolute `gtol`, and the permitted bounds are otherwise unchanged. The recorded `projected_gradient_inf` is in the optimizer's scaled coordinates; `material_projected_gradient_inf` is also retained. This affine variable change does not alter the objective or feasible set, but it changes optimization paths and makes the absolute stopping scale explicit. No tolerance sweep was performed. TRF retains its original material coordinates, analytic Full-J, and `x_scale='jac'`.

For each uniform arm the initial public-m1 anchor receives at most 20% of each global resource: wall, LU, and full complex RHS. Once the unchanged public-data gate triggers, each additional retained seed receives an equal share of the then-remaining **wall, LU, and RHS**, all three enforced simultaneously. Every seed's 8→12→24 stages receive cumulative resource fractions `.20,.45,1.00`; thus the stage increments are `.20,.25,.55`. Unused early-stage allowance remains available to later stages, and unused seed allowance remains available to subsequent seeds. Paid operations are never refunded. Atomic forward/adjoint/kernel blocks may cross a quota slightly; their actual costs and global overshoots remain exposed. Physical dimension stays 576 at every stage. The ordinary single-start baseline arms retain their stage wall schedule and common global LU/RHS caps; they receive the same scaled L-BFGS variables.

The fixed ordinary uniform-multistart set is still `[1,4,8,12,16]`, but seed 1 is **represented by the already paid common anchor** and is not rerun. Only seeds 4, 8, 12, and 16 are appended. Its raw proposal record labels the seed-1 representation explicitly.

The uniform rational scan still records every detected/refined valley. `uniform_bridge` now declares a maximum of **four** corrected ray candidates. It sorts distinct valleys by `t` and selects rounded uniformly spaced indices over that ordered list, including its first and last valleys. The selection uses no raw-score ranking. Valleys outside the declared finite candidate budget remain present and are marked `TRUNCATED_T_SPREAD_MAX_CANDIDATES`; floating-point duplicates alone use the `1e-7` pre-correction threshold. If there are at most four valleys, all are retained. `uniform_best_only` continues to correct exactly its raw-score-best valley as the explicit ablation. Post-correction candidate RMS deduplication remains `.25`.

Each attempted seed exposes its wall/LU/RHS allocation ends, start/end paid counts, stage quota ends, and correction status. A global stop preserves paid work and all actual material trajectories. V3 is the final declared code/configuration for fresh cases; no fresh-GT-based tuning is authorized by this implementation.
