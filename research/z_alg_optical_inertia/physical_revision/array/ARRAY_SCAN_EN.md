# Full-ring illumination and reception on the native grid

At a tested peak contrast of 0.75, adding illuminations improved the number of strict exclusions from three to five of six materials at 1% noise, and from two to five at 2% noise. At 2% noise, the fifth case required 24 illuminations with 20 receivers, or 16 illuminations with 32 receivers. The outcomes at 24 and 32 illuminations were identical. Every tested cap of 1 or above remained open. [Case results](array_scan.csv) · [Array figure](array_scaling.pdf)

| Receivers | Noise | 10 illuminations | 16 illuminations | 24 illuminations | 32 illuminations |
|---:|---:|---:|---:|---:|---:|
| 20 | 1% | 3 | 5 | 5 | 5 |
| 20 | 2% | 2 | 4 | 5 | 5 |
| 32 | 1% | 3 | 5 | 5 | 5 |
| 32 | 2% | 2 | 5 | 5 | 5 |

Each entry counts the same six materials and uses the full 9,216-cell model. Material 06 remains open across the scan. With 16 illuminations at 2% noise, the 32-receiver ring additionally resolved material 02. The six generating contrasts, with peaks 3.346–3.632, provide feasible upper witnesses for all error balls. [Forward checks](forward_audit.json) · [Observation manifest](observations/manifest.json)

## Angular bandwidth and the remaining source spectrum

At 400 MHz, the 2 m square is 2.667 wavelengths wide and its enclosing circle gives \(ka=11.848\), so \(2ka=23.695\). The 32-column incident matrix \(B\) and 32-row receiver matrix \(S\) retain 23 and 24 directions at 10% relative singular amplitude, and 29 directions each at 1%. Both operators have full numerical rank. The last relative singular amplitudes are 0.00247 and 0.00266. The drop near order 24 agrees with the angular bandwidth scale. The effective count depends on the accuracy threshold. [Spatial bandwidth](https://doi.org/10.1109/TAP.1987.1144024) · [Degrees of freedom of scattered fields](https://doi.org/10.1109/8.29386)

For scalar weights \(\mu>0\) and \(t=\lambda/\mu\), write
\[
 A=K/\mu=I-u\operatorname{Re}_H G+t\operatorname{Im}_H G+S^*S/\mu,\quad
 L=S^*Y/\mu+(u+it)B/2,\quad
 H_Y=Y^*Y-\mu L^*A^{-1}L .
\]
The joint lower distance is the square root of the sum of the positive source eigenvalues after discarding the first \(\nu=n_-(A)\) in decreasing order. When \(\nu\ge s\) the sum is empty. When \(\nu<s\), the retained positive terms can still be too small to exceed the radius.

| Negative state directions | \(u=0.75\) | \(u=1\) | \(u=1.25\) | \(u=1.5\) | \(u=2\) |
|---|---:|---:|---:|---:|---:|
| Real sector, without reception | 19 | 21 | 25 | 36 | 45 |
| Stable optical weight, 20 receivers | 4 | 16 | 20 | 22 | \(\ge39\) |
| Stable optical weight, 32 receivers | 4 | 13 | 16 | 20 | \(\ge34\) |

These are native interval counts or one-sided lower counts. They are independent of the material that generated the measurements. The receiver term is positive semidefinite with rank at most \(p\), and can remove at most \(p\) negative directions for fixed weights.

At 0.75 and 2% noise, the stable scalar optical weight's median distance/radius rose from 0.189 to 1.613 between 10 and 24 illuminations with 20 receivers, then to 1.645 at 32. With 32 receivers, the medians were 0.345, 1.766 and 1.798. This is an observed leveling-off at one electrical size.

At \(u=1\), the leading screened scalar profile for \(s=32\) had two retained positive eigenvalues with 20 receivers and seven with 32 receivers. The corresponding native ratios \(d/r\) were 0.13703 and 0.30106, both for material 04 at 2% noise. With 32 receivers the leading screened profiles at \(u=1.25\) and \(u=1.5\) replayed at 0.11331 and 0.01494. A tested 20-receiver profile had \(\nu=20<32\) at \(u=1.25\), with a strictly zero positive tail. At \(u=2\), the stable optical weight exhausts the source count on both rings. These cases separate insufficient illuminations from insufficient positive source spectrum. The numerical screening results concern 825 fixed scalar proposals and a second comparison with six inherited designs. A low score for any such profile leaves material feasibility open. [Results and scope](summary.json) · [Raw native source spectra](source_spectra.csv)

## Acquisition and figure notes

The 3.75 m rings start at \(0^\circ\). Plane-wave amplitudes scale as \(\sqrt{10/s}\), holding summed incident fluence fixed for equal dwell, background impedance and reference aperture. Receiver and data coordinates use \(\sqrt{20/p}\) quadrature weights. Noise is paired in each material, shares a complex Gaussian direction at coincident angular pairs, and is normalized to the desired relative Frobenius magnitude. The radius adds \(10^{-6}\|Y\|_F\) for the forward solve. The original \(10\times20\) data and operator arrays are reused unchanged. [Frozen design](FREEZE.md)

**Array figure.** (a) Relative singular amplitudes of the incident and receiver matrices. (b) Native negative counts, with upward arrows for lower bounds. (c) Counts of objects whose peaks are shown to exceed 0.75. (d) Numerically screened scalar maxima over objects and both noise levels at 32 illuminations. Filled points are full-grid interval replays of selected profiles. [Array figure](array_scaling.pdf)

**Source figure.** Both panels use 32 receivers and 2% noise. Empty circles mark the discarded leading eigenvalues, and filled circles mark retained positive eigenvalues. Panel (a) varies illumination count using the same stable optical weight at \(u=0.75\). Panel (b) compares that weight at \(u=0.75\) with the leading screened scalar profile at \(u=1\), both with 32 illuminations. [Source figure](source_spectra.pdf)
