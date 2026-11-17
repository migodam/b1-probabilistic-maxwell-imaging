# Coordinated illumination of a uniformly loaded TM patch

A coordinated three-shot design identified the uniform contrast of a four-cell patch in two independent geometries, while independently optimized local anchor illuminations left ambiguity between anchors. For the three-cell L patch, the coordinated and local methods produced the same four-dimensional feed space and the same relaxed profile. A noisy L-patch example demonstrates that the coordinated method can still admit widely separated, physically feasible contrasts.

## Physical model and matched comparisons

Both layouts contain 36 circular TM pulse cells of radius $0.045\lambda$, including 24 fixed heterogeneous host cells and 12 actual variable cells in a 3-by-4 subarray. Each patch has one common real contrast $q$, and the remaining active cells vary independently. The primary host has pitch $0.31\lambda$. The confirmation host has pitch $0.34\lambda$, a rotation and translation, a different host contrast pattern, angularly and radially perturbed transmitters, and a different receiver position. Both geometries use ten coherent external line sources, one receiver, and $k_0=2\pi/\lambda$. Patch and nuisance contrasts lie in $[1,15]$. There are 16 independently drawn materials per patch and geometry, for 64 material states. Geometry, target regions, source construction, and settings were set before evaluating the materials. The initial material for both designs fixes nuisance values to 8. Measurements come from independent full 36-cell VIE solves followed by subtraction of the exact known-host response.

The coordinated design builds the least-degree polynomial internal field supported on the patch, subject to physical incident range constraints. The calibration objective at $q=8$ selects a template by nominal reciprocal gain per mean radiated power. Its incident preimage coefficients span the acquired feeds. The local design selects a separate maximum calibrated gain per unit power at each of equally spaced contrasts, with the patch supported at that anchor and the same nuisance calibration. The local shots span those anchor feeds. Random acquisition uses an unselected complex orthonormal span of the matched dimension. Completion measures all ten whitened feed modes. All shot sets have total radiated energy two under the far-field Gram. Both Gram matrices are positive definite, with smallest eigenvalues 0.447 and 0.348. The maximum physical Euclidean feed norm among all shots is 1.116.

Each matched row receives detector noise at total complex $\ell_2$ radius $\epsilon=0,2\times10^{-5},2\times10^{-4}$. The same random functional on the underlying ten-dimensional whitened feed space is projected to each shot basis, normalized to unit $\ell_2$ norm and scaled by the same sampled fraction of $\epsilon$. When two acquisitions differ only by a unitary basis change, their data and noise transform by that same unitary. The code adds an absolute $8\times10^{-13}$ floating-point margin in profile and interval comparisons. A bounded real-material TRF baseline uses eight starts, a summed patch derivative and the same data for each acquisition. Full VIE validation uses the same pulse Green operator and a separate full-domain factorization.

## Source geometry and observability

| Layout | Patch cells | Degree-zero chain smallest singular value | Rank of $[R_T;H_{NT}]$ | Shortest template degree | Acquired shots | Largest angle to local span |
|---|---:|---:|---:|---:|---:|---:|
| Primary | 3 | 0.30695 | 3 | 2 | 4 | Below $10^{-12}$° |
| Primary | 4 | 0.22711 | 4 | 1 | 3 | 29.80° |
| Confirmation | 3 | 0.28118 | 3 | 2 | 4 | Below $10^{-12}$° |
| Confirmation | 4 | 0.26594 | 4 | 1 | 3 | 34.20° |

The ten physical incident directions have rank ten on the active cells in each layout. The regional stacked observability matrix has full column rank in all four cases, with smallest singular values 0.0620, 0.0394, 0.0702 and 0.0454 in table order. The zero-degree chains have full column rank and supply no constant template in these regions. The two-cell diagnostic had no polynomial template in either layout. The five-cell diagnostic first admitted a degree-zero template with two shots. These chain diagnostics apply to the specified class of focused identities.

We also formed the full incident coefficient matrix of all independent shortest template chains. The center evaluation has rank one for the L chain and rank two for the square chains, matching the generic focused-field kernel dimensions. The coefficient matrix has rank four for the L patch and rank six for the square. The smallest singular values below quantify the independence of those columns in the two physical systems.

| Layout and patch | Incident coefficient rank / number of columns | Smallest singular value | Column-normalized smallest singular value | Smallest singular value of the template evaluation at $q=8$ |
|---|---:|---:|---:|---:|
| Primary L | 4 / 4 | 0.10902 | 0.18118 | 0.16722 |
| Primary square | 6 / 6 | 0.14814 | 0.18763 | 0.29866 |
| Confirmation L | 4 / 4 | 0.09803 | 0.15892 | 0.28679 |
| Confirmation square | 6 / 6 | 0.13426 | 0.23994 | 0.40403 |

These rank observations support the four- and three-shot lower bounds for a fixed acquisition span that nulls arbitrary complex off-patch currents at every generic contrast. The bounded real-material inverse problem permits methods outside that class.

For the L patches, local and coordinated feed projector distances were $1.89\times10^{-15}$ and $4.60\times10^{-15}$. The maximum difference after unitary transport of their measured data was $7.64\times10^{-17}$ and $1.87\times10^{-16}$, including the equivariant noise. For the square patches, projector distances were 0.497 and 0.562. Whitened polynomial synthesis relative defects were below $1.19\times10^{-15}$.

## Inversion and ambiguity

| Layout | Patch | Median interval width at $\epsilon=2\times10^{-4}$ | Interval midpoint error | Pencil TRF error | Local-anchor TRF error | Random TRF error | Ten-shot TRF error |
|---|---|---:|---:|---:|---:|---:|---:|
| Primary | L | 0.1682 | 0.00611 | 0.00600 | 0.00600 | 2.371 | 0.00567 |
| Primary | Square | 0.1394 | 0.01457 | 0.01118 | 0.10109 | 2.012 | 0.00849 |
| Confirmation | L | 0.2530 | 0.01568 | 0.02307 | 0.02307 | 1.233 | 0.00838 |
| Confirmation | Square | 0.2154 | 0.02251 | 0.01829 | 0.10041 | 2.505 | 0.00919 |

All point errors are median absolute errors in the patch contrast. At $\epsilon=2\times10^{-5}$, the interval medians were 0.01693, 0.01407, 0.02359 and 0.02092 in table order. All 192 reported coordinated sets contained the true contrast. The maximum clean midpoint error was $6.04\times10^{-14}$. The 64 clean central polynomials yielded zero additional real roots in the box within the stated numerical tolerance. The minimum observed clean template visibility at the true contrast was $1.78\times10^{-3}$ for the coefficient normalization in the code.

The profile minimizes data mismatch over arbitrary complex off-patch reciprocal currents. On the 113-value candidate grid, the off-anchor statistic excludes the anchor values and candidates within 0.25 of the true contrast. At $\epsilon=2\times10^{-4}$, the coordinated square design had median rejection fraction 1, with minimum 0.991 by layout. The independently localized square design and the random design had zero-dimensional nuisance left null spaces throughout the off-anchor grid and rejected zero candidates through this relaxed test. The L designs produced identical relaxed profiles. The full ten-shot acquisition had median rejection fraction 1 in every layout and patch. Relaxed rejection measures necessary conditions on the physical model. The finite-dimensional, real, box-constrained model is tested separately by TRF and full VIE checks.

At the numerical radius $8\times10^{-13}$, independent full VIE re-evaluation validated distinct legal square target values for the clean data from local anchors in 16/32 sampled cases. Primary square case 3 admitted witnesses at $q=12.10151$ and $q=13.15373$, with residual norms $4.44\times10^{-14}$ and $1.30\times10^{-15}$. Confirmation square case 3 admitted witnesses at $q=11.69361$ and $q=13.73812$, with residuals $3.74\times10^{-13}$ and $5.43\times10^{-14}$. These are numerical feasibility examples. The corresponding target truths were 12.71329 and 12.95317. The witness starts and full contrast vectors are preserved.

One coordinated noisy L case in the primary geometry had two components, $[1.44452,1.85122]$ and $[6.57532,6.73438]$, at $\epsilon=2\times10^{-4}$. The true contrast was 1.63588. Pointwise constrained searches independently validated three interior values per component. Two illustrative full-VIE feasible contrasts were $q=1.64787$ and $q=6.65485$, with residuals $4.45\times10^{-5}$ and $2.35\times10^{-5}$ against the same data. This verifies physical ambiguity in the bounded real model for this noisy case. An interval hull spanning the two components is conservative and its midpoint can lie in the excluded gap.

The maximum full/reduced library response discrepancy was $1.02\times10^{-16}$. Median relative first-order Born response disagreement about the host was 1.448. In 21,696 direct grid evaluations, the original scalar inequality and component membership agreed. The focused scalar residual divided by its feed norm agreed with the rank-one complex nuisance profile to $7.03\times10^{-17}$. The patch TRF Jacobian matched centered differences with at most $1.11\times10^{-9}$ relative error across a material point per layout and region. Eight analytic interval regressions passed, including an isolated irrational clean root and positive-noise tangent contact. Re-running the corrected interval solver on the 192 saved rows left all observed endpoints unchanged.

At $\epsilon=2\times10^{-4}$, median single-target interval evaluation time was 0.30–0.36 ms CPU, including polynomial construction and component extraction. Same-data TRF took 73–78 ms CPU and 591–654 median reduced factorizations over eight starts. Setup for each host factorizes one 36-cell system and uses 22 complex right-hand sides. All 768 measurement rows were processed in 59 CPU seconds. The raw results retain every TRF endpoint, including starts at the evaluation limit.

## Reproduction and evidence

Run from the repository root in a Python environment with NumPy, SciPy, Matplotlib and SymPy. The forward operator is shared with R2.

```bash
PYTHONDONTWRITEBYTECODE=1 uv run --with numpy --with scipy --with matplotlib --with sympy \
  python research/research_z_alg_r7/src/pencil_experiment.py
PYTHONDONTWRITEBYTECODE=1 uv run --with numpy --with scipy --with matplotlib --with sympy \
  python research/research_z_alg_r7/src/pencil_followup.py
PYTHONDONTWRITEBYTECODE=1 uv run --with numpy --with scipy --with matplotlib --with sympy \
  python research/research_z_alg_r7/src/plot_pencil.py
```

The [comparison figure](../figures/pencil_comparison.pdf) shows the physical source geometry, the relaxed profile for confirmation cases 0, the median grid rejection and median errors, and bounded-real witness endpoints for the first noisy square case in both layouts. The [disconnected example](../figures/pencil_disconnected.pdf) plots the original focused residual inequality and the independently verified contrast witnesses. Individual [CSV rows](pencil_individual.csv), [raw fit records](pencil_raw.json), [complex observations and feeds](pencil_observations.npz), [summary](pencil_summary.json) and [full VIE follow-up](pencil_feasible_followup.json) are available.
