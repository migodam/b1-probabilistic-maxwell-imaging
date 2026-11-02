# R2 actual baseline/API and blind scoring contract

Frozen before fresh generation: `configs/blind_cases_frozen.json`. Generation ran once through `experiments/generate_blind.py`; its public payload and exact cost are in `audit/GENERATION_RECEIPT.json`. This file records source-based recommendations and the parent's subsequently chosen matched-work schedule, not a new numerical outcome.

## Existing APIs and physical model

CPU model: `research.research_fbg_r3_experiments.src.physics.MaxwellModel(public)`, using the frozen full nonlinear disk-pulse VIE and retained LU. Public matrices are the existing `research_z_alg_r3_digits/public_data/operators.npz`, with G576×576, E0 576×12, S48×576 and unchanged centers/source/receiver coordinates. m is one shared real contrast per cell; final material dimension576, box[.02,20]. No support freeze, ROI restriction, glyph decoder or tensor adapter is used.

Actual mature current entry:

```python
from research.research_fbg_r3_experiments.src.cie_tsom_baseline import run_cie
out = run_cie(model, Y, np.full(576, .1), config)
```

It uses native projected Fourier coefficient coordinates, receiver DPCS/APCS, one exact quadratic current line step per alternation, sourcewise PR+/non-descent safeguards, shared-real legal material readout, and actual regularized material inner solves. Dual beta is a two-metric iteration, not monotone block descent on one common scalar objective. Stage memory resets and native zero-padding transport are part of the algorithm. It does not implement FFT complexity or generic correlated joint whitening.

Actual common physical corrector entry:

```python
from research.research_fbg_r3_experiments.src.full_maxwell_corrector import corrector
closed = corrector(model, Y, current_incumbent,
    {'method': 'data_lbfgs', 'max_nfev': 1000, 'box': (.02, 20.)})
```

For this entry, `max_nfev` means L-BFGS maxiter1000, with maxfun4000, ftol1e-12, gtol1e-8 and maxls30. It is not a cap of1000 function evaluations. Omit `wall_seconds` for no wall cutoff. Analytic data adjoints and the same original prior are used. Every evaluated physical incumbent is recorded; output `m` minimizes normalized original data plus prior over its evaluated candidates. It does not select using GT. Budget closure/re-evaluation and setup must remain paid.

The original helper imports the R4 physics module; R2 currently has that tracked dependency. Packaging must preserve those exact source paths or include the small dependency snapshots. Replacing `model.regularizer(m)` with its value/gradient-only equivalent may avoid building an unused dense Hessian, but only if the exact historical prior/gradient are retained; that is implementation efficiency, not scientific prior modification.

## Adequate-duration baseline

Historical45-second E3 runs frequently ended before all current stages; they cannot be treated as converged failures. A source-based adequate-duration recommendation was the strict native schedule (miniter50, coefficient_tol1e-8, maxiter500 per stage) plus the1000-iteration physical corrector. The parent selected a stronger **fixed-work** comparison to isolate update effects:

- Native sourcewise dual-beta CIE: mf3/beta1, mf6/beta_material.2/beta_current1, mf10/same dual beta; **500 steps each, coefficient_tol0**, all1500 steps run unless a real execution failure occurs. dpcs_rtol1e-3 and material_inner_maxiter60 are unchanged. No short wall cutoff.
- Global PR arm: same spaces, masks,1500 current steps, normalization, material blocks, initialization, prior, physical scoring and1000-iteration corrector. Only the explicitly declared global current iteration differs.
- Fixed-beta1 variable projection:250 declared material iterations plus the same1000-iteration corrector. This is a different optimization decomposition and work budget; record exact physical/reduced factorization work rather than calling it matched compute.
- Direct physical L-BFGS: the executed runner performs an auxiliary physical pass capped at1000 iterations from the public uniform m=.1, followed by the common physical corrector capped at1000 iterations again. Thus it is **two up-to1000 passes**, not1000 iterations in total. Actual early optimizer stops are retained. This stronger direct allowance is disclosed rather than retrospectively rerun. It is a control, not a replacement for mature CIE.

There is no parameter sweep or GT-triggered continuation. All MAXITER, line-search, material-inner failure and stopping statuses survive the report. Full scheduled work is still a finite budget, not a proof of convergence or a performance ceiling. Equal current steps are not equal FLOPs: sourcewise/global directions and line search differ. Report CPU/wall, current products/factorizations, original LU, each RHS type, preparation, candidate validation and failures separately.

## Frozen fresh cases and noise

Three independently realized shapes: asymmetric annulus, interrupted curved bar, unequal two-island control. Every material is `.1+12*alpha` with an8×8 per-cell area average; support is private `alpha>=.5`. The generator uses fresh128-bit entropy and UUID IDs, keeps the mapping private and writes public files sorted by opaque ID. Shapes are selected before inversion; no GT self-IoU filter or regeneration is allowed.

Each public observation has M=576 complex entries and iid circular complex Gaussian noise with `sigma=1e-3*||Yclean||_F/sqrt(M)`, real/imag variance sigma²/2. This is **not** the older fixed-norm noise direction. W=I/sigma; the per-observation99% radius is `sqrt(chi2.ppf(.99,2*M)/2)`. The exact realized noise can fall outside this radius with1% probability; no rejection or repair is permitted. This is a marginal per-case statement, not99% simultaneous coverage of all three objects.

The original objective is retained exactly:

`Phi = .5*||F(m)-Y||²/||Y||² + R_historical(m)`.

Let a=sigma²/||Y||². The same objective is `a*.5*||W(F-Y)||² + R_historical`. Dividing the entire objective by a gives `.5*||W(F-Y)||² + R_historical/a`; thus the prior multiplier in a **unit-weight whitened objective** is `||Y||²/sigma²`, not its reciprocal. Any implementation must choose one of these exactly equivalent forms consistently for values, gradients, selection and regularization.

## Independent scoring after seal only

The main/inversion worker must never open `sealed_truth` or generator-private labels/seeds/cleanY/materials. `data_public` contains only operators, coordinates/acquisition, Y, sigma/radius and known common priors. Before scoring, every frozen method×opaque-ID disposition and top1 material must be saved and sealed; the parent then explicitly signals that closure to this scoring worker. A failed run stays a sealed failed disposition and is not rerun using its score. No scoring code is imported by an inverse runner.

Primary per-case material/shape metrics retain historical semantics: **NRMSE=RMSE/12**, and fixed256-bin Otsu computed solely from the reconstruction, compared to private `alpha>=.5`. The joint gate is NRMSE≤.15, Otsu-IoU≥.7 and the frozen99% whitened data-radius gate. Also show the older noisy relative-residual≤.003 gate as a separately named comparability diagnostic. Record the material+data gate when IoU alone fails, and a GT self-Otsu diagnostic only after closure.

Additional frozen descriptive metrics: fixed known-midpoint IoU at m>6.1, foreground area error and physical support-centroid distance. Uniform prediction Otsu produces empty foreground, as in the actual historical scorer. IoU, empty support and undefined centroids must use explicit semantics; undefined centroid distance is null+reason. Candidate-pool `best@K` is oracle-only; top1 is the sealed online physical selection. Three objects cannot establish a population success rate, basin volume, universal superiority or a statistical significance claim.

## Initial paid work and closure update

The write-once generator paid3 original LU and36 forward complex RHS; wall .150064s, including operator loading, forward/noise and public/private output inside the generation function. Imports and the receipt's own serialization are outside that timer. Historical operator assembly was not repeated and its old cost is unknown rather than zero. Reading R1 arrays/reporting costs are separate uninstrumented overhead.

All12 predictions were independently sealed at06:09:07.838542UTC before the scorer's first private truth access at06:09:12.846938UTC. The actual gates give native/global2/3 joint passes (ring and two-island control), known profile1/3 (ring), direct0/3. Curved bar fails every arm. These are three specific objects; no extra coverage or new global-PR imaging superiority was established. Detailed scores and plots are in `results/BLIND_SCORES.*` and `figures/BLIND_RECONSTRUCTIONS.*`; the independent final review and cost ledger supersede the initial accounting paragraph for full R2 totals.
