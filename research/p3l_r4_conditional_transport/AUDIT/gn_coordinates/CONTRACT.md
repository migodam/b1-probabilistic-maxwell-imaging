# Full66 conditional GN coordinates — implementation contract

This is an invertible observation-derived coordinate family for an otherwise learned full66 conditional density. It does not declare the GN Gaussian to be the posterior, and no superiority has been established. No campaign fit, optimizer update, physical state, Jacobian, or new oracle teacher was launched by this worker.

## Map and density

`src/p3r4/conditional_coordinates.py` computes the Gaussianized Jacobian from the existing whitened native Jacobian:

- Jx = Jq diag(1,...,1,φ(anchor_eta1),φ(anchor_eta2)).
- U is positive-diagonal upper triangular, with UᵀU = I + JxᵀJx.
- U is obtained from QR of stacked[Jx;I]. This is algebraically chol(I+JxᵀJx)ᵀ, without explicitly forming the normal matrix, adding a numerical H floor, or discarding weak directions.
- GN: r = U(x−context_x), and x = context_x + solve_triangular(U,r).
- Random: W = U Rfixedᵀ; r = W(x−context_x), and x = context_x + Rfixed solve_triangular(U,r).

Rfixed is one real orthogonal matrix, generated independently of TRAIN and GT using default seed20261025113, shared across parents and training seeds and saved in the checkpoint. It preserves the covariance spectrum while rotating its raw-coordinate alignment. Its actual floating-point log absolute determinant is retained alongside sum(log diag U). No full physics inverse is formed.

The density is with respect to the original Gaussianized dx66:

`log q_x = log q_r + log|det W|`.

TRAIN-only r standardization and any supplied fixed orthogonal frame are composed afterwards; their determinants are also included. Validation raw NLL subtracts log|det W| and is commensurable with the `none` family. CFM energy-score selection undoes both TRAIN normalization/frame and W before measuring raw latent error. TRAIN labels may fit the declared global normalizer and diagnose representation error; they never define W, its rotation, or an extra H floor.

All66 coordinates remain stochastic and modeled. Raw512+context_x66 conditioner inputs and network counts are unchanged; cached J can enter the coordinate map without adding NN features. Optional J summaries remain a separate explicit feature control. GN/random require saved context_J even when no J-summary feature is selected. The online evaluator explicitly allowlists that field for these checkpoints and still never reads ambiguous x or GT fields.

## CLI and CFM base

`src/p3r4/train_anchor_density.py` adds:

```text
--coordinate-preconditioner none|gn|random
--alignment-seed 20261025113
--preconditioner-base coordinate|raw
```

Default `none` retains the previous translation rule and raw covariance semantics. Old checkpoints without new metadata load as `none`. GN/random use context_x as the center and require residual centering. They default to `--preconditioner-base coordinate`: CFM starts from a fullrank Gaussian in r coordinates, with supplied Cholesky or anchor_scale*I. Thus an isotropic r base induces raw covariance proportional to H^-1, or its Rfixed conjugate for the random control. This is a computational source, not a GN-posterior claim. Supplied base_mean remains a y-derived raw-coordinate mean, mapped consistently into r coordinates.

`--preconditioner-base raw` is an explicit alternate control: supplied raw covariance is transformed through W. That option retains the same raw base law and therefore does not test the intended GN-versus-random source alignment. Cholesky units are explicitly declared in each checkpoint contract. Covariance factors transform with QR of the factor transpose, avoiding A Aᵀ and associated loss of conditioning. CFM log density remains an opt-in numerical continuous-ODE density; these affine coordinate Jacobians do not make the finite RK4 sampler an exact-MH proposal.

Examples, for later parent-budgeted runs only:

```text
python -m p3r4.train_anchor_density ... --kind spline --coordinate-preconditioner gn
python -m p3r4.train_anchor_density ... --kind spline --coordinate-preconditioner random --alignment-seed 20261025113
python -m p3r4.train_anchor_density ... --kind cfm --coordinate-preconditioner gn --preconditioner-base coordinate
```

Use PYTHONPATH=src or `python -m src.p3r4.train_anchor_density` from repository root. Same TRAIN/validation parent sets, model capacity, three training seeds and declared budgets are needed for meaningful comparisons. No actual comparison has been run here.

## Numerical evidence and limit

`NUMERICAL_COORDINATES.json` records:

- Gram-factor error3.11e-15; full66 float64 inverse errors<=3.9e-16 on moderate cached Jacobians.
- Genuine float32 triangular/map inverse errors<=1.2e-7 in the moderate test.
- Full66 autodiff determinant agrees with the analytic determinant.
- Random covariance is the intended orthogonal conjugate, with matching full66 eigenvalues.
- Gaussian/diagonal/spline/triangular and coordinate/raw-base CFM density references agree to<=2e-5 for float32 networks; all sample coordinates are retained.
- **Retained stress limit:** an extreme rank1 mixed-scale J at1e8 remains fullrank after the identity prior, but float32 inverse error reaches1.016, while float64 error is2.6e-9. This is not hidden by an H floor, truncation, or a universal numerical PASS claim.

Production coordinate construction and inverse solves use CPU float64. NN coordinates remain float32 on CPU/MPS. Every fit now writes actual TRAIN/validation `COORDINATE_AUDIT.json` representation replay: normalize r, cast to float32, invert normalization/frame/W and report maximum/median/p95 parent error and latent RMS. These are diagnostics, not automatically relaxed acceptance thresholds. Parent must inspect actual cached-data precision before interpreting gains. The moderate synthetic update0 fixtures have maximum raw representation errors about1.6e-8.

`PIPELINE_AUDIT.json` records a three-parent, four-variant synthetic freeze→online→separate-PID score replay, K128, required Jacobian context, zero failures, no forbidden-label reads and no physical diagnostics claimed. `NONE_DEFAULT_REGRESSION.stdout` reruns the earlier centered/absolute, J-feature and all-five-family algebra audit. These are implementation tests, not new-parent scientific evidence.

## Why a learned coordinate family is being tested

Historical failures remain in the unchanged R2 reports: full66 GN/Gaussian strong conditionals were unreliable; independent Student proposals had true median ΔΦ≈40203 versus a local prediction≈10.98; learned10×prior56 correction failed with large/all rejections. Decoder saturation seams, nonlinear finite Maxwell response and coupled high-dimensional conditional mass invalidate a naive global use of local curvature. A normalized covariance or a cheaper evaluator does not establish posterior quality or mixing.

Here U is only an invertible chart. A full66 network learns the nonlinear residual conditional law in that chart; it is free to learn non-Gaussian shapes and widths rather than treating isotropic r as the final posterior. This is a candidate inductive bias, with a random-alignment control and proper raw NLL. Calibration, new-parent material scores, paid residual diagnostics, finite-precision audits and final scientific judgment remain parent-owned.
