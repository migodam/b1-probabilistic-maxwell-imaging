# Illumination, observability and host error

These synthetic tests use full nonlinear disk-pulse 2D scalar TM current solves for the physical observations. The material is real. All feed banks use the radiation Gram matrix \(Q_{ab}=J_0(k|s_a-s_b|)\). Three experiments examine different questions and geometries.

## One illuminated cell

Six active cells are embedded in the nonzero six-cell H12 host. The prescribed value is 15.25. For each active cell, 256 scrambled Sobol points vary the five nuisance values over [0.5,30]. Two observations at the same nuisance values compare the prescribed material to 17.25. Each arm uses one unit-radiated-power beam. Fixed-target drift is the difference from the reading with all six materials at 15.25. The paired signal changes only the target.

| Full eight-feed ring, one receiver | Maximum fixed-target drift | Median paired signal |
|---|---:|---:|
| Matched beam | 1.22e−16 | 0.00502 |
| Minimum-power host field focused at the target | 0.0972 | 0.00463 |
| Maximum local target sensitivity | 0.142 | 0.00731 |
| Maximum sensitivity under complex nuisance-Jacobian nulling | 8.72e−17 | 0.00502 |
| Random complex feed | 0.0796 | 0.00206 |

At the nominal state every nuisance receiver factor is nonzero. The independently constructed Jacobian-null drive matches the analytic drive to 1.51e−15 in relative Euclidean norm after a global phase alignment. Its nuisance Jacobian sets the total field to zero at every off-target cell. This equality identifies an existing local design that attains the global identity here. Unconstrained local sensitivity gives a larger median paired signal while admitting substantial nuisance drift.

With eight equal-power emitters on 120° and 60° arcs, the matched median signal falls to 0.000512 and 0.000192. The 60° optimal radiation-only drives require Euclidean norms from 104 to 230. Imposing the additional drive cap of 4 and reoptimizing each matched beam lowers the median to 0.000120. The capped region includes radiation powers below one when the drive constraint binds first.

A deliberate bounded search varied two nuisance cells while holding the other three at 15.25. Six of 60 cell/pair combinations yielded a first-receiver reciprocal factor below 1e−12. For example, set active materials (2,5) to (29.019046807605, 6.352452407794) in one-based indexing and vary cell 1 from 0.5 to 30. At the original receiver the largest response difference over eight primitive feeds is 1.00e−16, versus 0.0501 at the opposite receiver. This is a constructed numerical darkness witness. The archive includes failed searches and the primitive complex responses.

## Uniform patch under rank-deficient illumination

A fresh jittered 3×4 active grid is surrounded by five known scattering cells. Ten physical feeds give rank 10 in the 12-dimensional active space. One receiver records each exposure. For each nested, connected patch of size 1–7, 128 Sobol states vary the common patch material in [2,28] and the independent exterior materials in [0.5,30].

Write \(B\) for the map from feeds to active host fields, \(H\) for the reduced Green matrix, \(R\) for the receiver row and \(w\) for the unknown reciprocal response row. For a patch mode \(p\), one prescribed-value test needs \((I-qH)p\in{\rm ran}(B)\). The constant-mode two-shot readout needs both \(p,Hp\in{\rm ran}(B)\). It gives \(q=wp/(Rp+wHp)\) when the denominator is nonzero. We selected the mode using only nominal material 12 by maximizing nominal squared reciprocal visibility per predicted variance of the two individually power-normalized shots.

| Patch size | One-value mode dimension | Constant two-shot mode dimension | Largest clean complex error, optimized | Median real error at 1e−4 RMS per shot |
|---:|---:|---:|---:|---:|
| 1 | 0 | 0 | — | — |
| 2 | 0 | 0 | — | — |
| 3 | 1 | 0 | — | — |
| 4 | 2 | 0 | — | — |
| 5 | 3 | 1 | 4.88e−12 | 1.111 |
| 6 | 4 | 2 | 8.76e−13 | 0.0802 |
| 7 | 5 | 3 | 6.87e−13 | 0.0374 |

The maximum full-solve fixed-value drift for the reachable one-shot tests is 5.36e−16. For sizes one and two, selecting the mode with smallest absolute off-range field and projecting it to the reachable space gives maximum drifts of 0.00367 and 0.00252. For the two-shot projection of a constant spatial template at sizes 5, 6 and 7, the median clean complex errors are 19.67, 22.04 and 17.42. Here “constant-mode” refers to independence of the template from the material being tested. Its complex spatial coefficients can vary within the patch. Higher-degree material-dependent modes are studied in the separate pencil experiment.

On the six-cell patch, nominal optimization lowers the median real error from 0.820 for one random reachable mode to 0.0802 with 1e−4 circular complex RMS noise per shot. The optimized 95th-percentile error is 1.001. Global phase was canonicalized before applying common per-shot noise. For five cells, both designs span the same one-dimensional space and give identical noisy estimates up to roundoff. With zero-mean within-patch perturbations of peak 0.1, 0.5 and 1.0, the six-cell optimized noiseless real quotient misses the arithmetic mean by medians 0.0263, 0.137 and 0.279. The exact heterogeneous quotient is a complex weighted average with weights \((R+wH)_j p_j\). At peak 1 its real part lies outside the physical patch range in 2/128 six-cell cases.

For a uniformly 12 nominal object, force the six-cell patch to 13 and fit all six exterior materials to preserve the two readings. The projected-uniform control admits a pair with full nonlinear data distance 1.72e−13. For the reachable optimized mode, the two-shot identity gives the lower bound 0.00212322974520030 on the Euclidean complex data distance for every such exterior. A bounded fit attains 0.00212322974520031. With the optimized design, an exterior cell can change from 12 to 16 at fixed patch material while a change of the other exterior cells restores the data within 3.31e−15. The complete paired materials and complex readings are in the JSON. On 24 separate clean states, one nominal-start trust-region fit of the same two optimized readings gives median patch error 2.75e−9 and worst error 6.81. Six fits reached the 300-evaluation limit.

## Three-shot pencil under known-host calibration error

All 24 nonactive contrasts in the 36-cell square experiment were multiplied by 1.01. Both layouts use the archived 16 materials, three physical beams, original noise vectors with norm bounded by 2e−4, nominal readout, and the same error radius. A separate full solve supplies the noiseless empty-active actual-host baseline. Templates and transfer rows remain nominal in both perturbed arms.

| Subtraction and material | Coverage, primary | Coverage, confirmation | Empty sets, combined | Worst hull midpoint error among nonempty sets |
|---|---:|---:|---:|---:|
| Nominal reference and material | 16/16 | 16/16 | 0/32 | 0.128 |
| Nominal reference, changed host | 9/16 | 5/16 | 9/32 | 0.345 |
| Exact actual-host reference, changed host | 2/16 | 6/16 | 17/32 | 0.426 |

Model error before noise exceeds the allowed 2e−4 norm on all 32 cases for both perturbed arms. Its median norm falls from 6.42e−4 to 4.51e−4 upon subtracting the actual empty-host reference. Its projection onto the identifying complex direction still exceeds the radius at truth after adding noise in 24 cases, compared with 18 using nominal subtraction. Relative changes of the reduced Green matrix, active illumination and receiver row are 0.283%, 0.570%, 0.659% in primary and 0.248%, 0.540%, 0.621% in confirmation. This selected perturbation shows that empty-host measurement alone does not compensate for changed illumination and material-to-receiver transfer.

### One additional calibration exposure per fixture

Under the explicit model that every nonactive contrast shares one unknown multiplicative scale, a single empty-active measurement estimates this scale without target data. The external calibration drive maximizes the magnitude of the nominal host-only derivative with respect to scale, under radiation energy ≤1 and feed norm ≤4. Its energy is 1 and the respective feed norms are 1.336 and 1.428. We used one fixed random noise phase per fixture at magnitude 2e−4, fit the real scale by bounded complex least squares on [0.97,1.03], and resynthesized all three target shots and the host subtraction with the fitted host. Each target measurement set has total radiated energy 2. The calibration exposure was reused for all 16 objects on its fixture.

| Layout | Fitted scale (true 1.01) | Covered | Empty | Worst midpoint error |
|---|---:|---:|---:|---:|
| Primary | 1.009606406 | 16/16 | 0 | 0.0824 |
| Confirmation | 1.011348834 | 16/16 | 0 | 0.1716 |

The combined median midpoint error is 0.0161. The target model-error norm peaks at 1.10e−4. After target noise is added, five confirmation objects exceed the assumed 2e−4 data-vector radius, although their projected errors at truth all fit the radius. The observed 32/32 coverage therefore applies to these synthetic cases and two calibration draws. In the diagnostic using the true scale for redesigned shots and subtraction, coverage is also 32/32 and the worst midpoint error is 0.1291. The one-parameter calibration result leaves spatially varying host errors open.

### Propagating calibration uncertainty uniformly over the unknown material

The one calibration reading also bounds its own estimation error when the uniform-scaling model and both noise radii hold. A 512-bin subdivision of [0.97,1.03] bounds the change of the complex empty-host trace by a resolvent derivative envelope on each bin. For each fitted fixture, the resulting scale uncertainty and the passive scattering inequality bound the contribution of the host error to the focused polynomial for every admissible real nuisance material. We evaluated a conservative polynomial outer set using the square root of twice the sum of squared systematic and additive error envelopes. The only per-object input to the set calculation is the archived three-shot observation.

| Layout | Retained scale interval | Radius about fitted scale | Neumann product | Median hull width, noise only | Median hull width, uniform outer |
|---|---|---:|---:|---:|---:|
| Primary | [1.00890625, 1.01031250] | 0.0007061 | 0.08891 | 0.1399 | 1.0638 |
| Confirmation | [1.00996094, 1.01277344] | 0.0014246 | 0.11102 | 0.2146 | 2.8096 |

Both calibration intervals contain the simulated scale 1.01. The uniform outer contains the true patch value in all 32 cases; none spans the full material box of width 14. Its combined median width is 1.665, with maximum 7.853. The largest measured projected systematic residuals reach 1.20% and 2.97% of their envelopes on the two fixtures. The largest sampled passive scattering norms are 0.99956 and 0.99980. Computing the two calibration lookup tables required 1,024 linear solves of order 24. The postprocessing used zero further full loaded solves and took 0.12 seconds. The uniform bound assumes the stipulated single host-scale family, passive real active loads and bounded noise. Kernel and interval endpoints were evaluated in floating point.

The paper figure shows coverage and empty sets before and after calibration, with paired per-target interval widths for both fixtures. [PDF](figures/host_calibration.pdf) · [PNG](figures/host_calibration.png) · [Plotting script](src/plot_host_calibration.py)

## Reproduction and raw data

From the R7 study directory, with NumPy, SciPy, Matplotlib and SymPy available, run:

```sh
python src/geometry_controls.py
python src/geometry_patch.py
python src/geometry_mismatch.py
python src/geometry_calibrated.py
python src/host_error_bound.py
python src/plot_host_calibration.py
```

The first two accept `--physics /path/to/research_z_alg_r2/src/physics.py` for sparse checkouts. The mismatch study reads `results/pencil_raw.json` and `results/pencil_observations.npz`. JSON outputs carry per-arm summaries, selected pairs and observed timings. Corresponding `*_raw.npz` files carry individual complex responses, physical feed vectors or archived noise. In the mismatch NPZ, the three clean data rows correspond in order to nominal, nominal-host subtraction on the changed object, and actual-host subtraction. The control and patch runs took about two seconds each including plots in this environment. The host-error and one-shot calibration runs took approximately 0.10 and 0.13 seconds after imports, respectively. Figures are available as both PDF and PNG under `figures/geometry_controls` and `figures/geometry_patch`.
