# ε9 sphere: material/geometry/convention audit and one authorized refinement

This review uses qualification outputs and analytic Mie only; no inversion truth or material scores were opened. εr=9 means χ=8 and n=sqrt(9)=3. The inherited analytic routine accepts εr, not χ or refractive index. Both reference and DDA use k=2π, radius .145λ, incident k=z/E=x, exp(−iωt), outgoing exp(ikr)/r amplitude. The vector reference is i S2/k times e_theta on phi=0, identical to the previously inherited ε3 qualification convention. No exact convention/normalization wiring error was found by inspection; that does not prove the coarse DDA accurate.

|Grid|h/λ|Cells/internal λ|Far complex-vector error|Unapplied best gain amplitude|Gain phase|Shape residual after unapplied gain|
|---|---:|---:|---:|---:|---:|---:|
|8³|.05|6.667|61.583%|.52009|−5.334°|.38004|
|10³|.04|8.333|55.526%|.57840|−5.428°|.35409|
|16³|.025|13.333|43.057%|.68379|−5.363°|.28181|

The gains above are diagnostic complex least-squares projections of the meridian vector curves. They are **never applied** to data, materials or physical outputs. A substantial angular-shape error remains even under that diagnostic scalar projection. Backscattering also disagrees strongly in the saved complex curves; the result cannot be repaired by a claimed scalar normalization adjustment. Measured errors decrease under refinement but remain far outside an accurate continuous-physics qualification. Three coarse levels do not establish an asymptotic convergence order or exclude further physics/discretization defects.

For example, the analytic backward/forward complex-field amplitude ratio is approximately .0176, whereas grid16 gives approximately .358. The reference's strong backward cancellation is not reproduced by this coarse discretization. This is an observed angular-mode mismatch, not proof of its unique cause; it cannot be removed by applying one global gain.

The same fixed24³ midpoint MG rule is used at every level: χ_eff=3fχ/[3+(1−f)χ]. Geometric volumes are approximately .0127646, .0127730 and .0127688λ³, versus exact .0127701λ³. Small volume error does not imply accurate high-contrast dynamic scattering. This MG sphere control also differs from linear-χ box-volume fractions used for cube/hollow data; it is not an error bound for their fields.

The unique added16³ solve eliminated strictly zero-contrast background cells using the existing TensorDDA arbitrary-occupied-points API. This leaves1136 active cells, vector dimension3408. A complex dense matrix is185.8MB; a conservative eight-matrix estimate was1.487GB. Measured peakRSS was1.250GB, CPU3.2135s, wall3.544s, with exactly1 LU/1 source RHS. A full16³ single matrix alone would be2.416GB. No new Green/self-term/factorization backend was introduced. The grid16 cross section is .132947λ² versus analytic .192433λ², a30.913% error. Linear residual7.6e−15 certifies only the finite linear solve.

The sphere has external size parameter ka=.911062 and internal nka=2.733186. At the declared frequency the exact electric/magnetic dipole magnitudes are .45577/.44142, with dipoles accounting for99.887% of the analytic total scattering weight. The root-authorized five-frequency **analytic-only** diagnostic is:

|k/k0|Internal size parameter|Magnetic dipole magnitude|Electric dipole magnitude|Exact cross section λ²|
|---|---:|---:|---:|---:|
|.95|2.59653|.25707|.38898|.115152|
|.975|2.66486|.33457|.42173|.145728|
|1|2.73319|.44142|.45577|.192433|
|1.025|2.80152|.58686|.49096|.266332|
|1.05|2.86984|.76765|.52713|.375875|

This establishes strong frequency sensitivity and an increasing magnetic-dipole response within the inspected interval. It does **not** show a peak in that interval, prove that the measured DDA error is solely a resonance shift, or authorize changing frequency/material/radius to improve agreement. No further sweep or fit is performed.

The high-contrast two-cube inversion may therefore support only a finite-DDA algorithm comparison at its declared grid/model mismatch. It does not establish accurate continuous ε9 imaging. All original8/10 failure evidence remains under `RESULTS/vector_extension/validation`; the added result and full complex samples are under `validation_refinement16`.
