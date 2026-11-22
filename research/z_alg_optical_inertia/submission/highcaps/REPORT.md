# Quadrant multipliers at higher full-aperture caps

On the native 96 × 96 full-aperture model, two continuous families of genuinely local quadrant multipliers have at least ten negative directions for every cap \(1\le u\le1.5\). All six materials share this 20-receiver, ten-illumination operator. These parameter boxes therefore exhaust the retained spectrum for all six materials at any noise condition. The proof does not read a measurement or a reference material.

Let \(D_b\) indicate the four disjoint 48 × 48 quadrants. Write \(D_\mu=\mu\sum_b w_bD_b\) and \(D_\lambda=\mu\sum_b(t+\delta_b)D_b\), where \(\mu>0\) and the four \(w_b\) and four \(\delta_b\) can vary independently within the specified boxes.

| Independent sector weights \(w_b\) | Optical centre \(t\) | Independent jitter \(\lvert\delta_b\rvert\) | Vertices | Verified upper bound for largest eigenvalue of the projected \(K/\mu\) |
|---|---|---|---:|---:|
| \([0.875,1.125]\) | \([-4,16]\) | \(\le0.125\) | 1,024 | \(-0.039393812878198\) |
| \([0.5,1.5]\) | \(0\) | \(\le0.06\) | 512 | \(-0.03777103969586489\) |

An eigensolve first proposes ten currents in the numerical receiver kernel for each of nine scalar centres. The interval verifier corrects the chosen basis \(V_0\) through
\[
V=V_0-S^*(SS^*)^{-1}SV_0.
\]
An interval solve with the saved binary64 receiver matrix validates the correction, and \(SV=0\) holds algebraically. The rounded Gram lower bounds exceed \(0.99999999999998\); the interval norm of the computed receive defect is below \(1.36\times10^{-48}\). Applying the full saved 192 × 192 Fourier multiplier to each projected column uses separable 160-bit Arb transforms, native padding and cropping. Define \(A_b=V^*D_bV\), \(H_b=\operatorname{Re}_H(V^*D_bGV)\) and \(R_b=\operatorname{Im}_H(V^*D_bGV)\). On this subspace the receive term vanishes and
\[
V^*(K/\mu)V=\sum_b\left[w_b(A_b-uH_b)+(t+\delta_b) R_b\right].
\]
The Hermitian matrix is affine in each scalar in isolation. Its values over the rectangle are convex combinations of vertex values, including the two cap endpoints. Independent interval congruences establish the displayed negative upper bounds at every vertex. The full-space inertia therefore satisfies \(\nu\ge10\).

The numerical box screen did not find a shared proof subspace in the wider \(\lvert w_b-1\rvert\le0.25\), \(\lvert\delta_b\rvert\le0.25\), \(t\in[-4,16]\) box. Its best proposed worst vertex was \(+0.0672814\). This is an unresolved box, since a different subspace may work at that vertex or over smaller boxes. An additional native numerical check at \(u=1\), \(\mu=10^{-3}\) examined six more skewed sector patterns at each of three constant optical ratios, \(0,0.4,1\). All 18 tests returned at least 12 negative eigenvalues of \(K\), with the largest returned eigenvalue of \(K/\mu\) at most \(-0.0145648\). This is a numerical screen with discrete parameters and does not extend the interval result. More general spatial partitions and parameter ranges remain open.

**Reproduction.** From the project root, run `screen.py`, `verify.py`, then `skew.py` in this folder using Python with NumPy, SciPy, python-flint and threadpoolctl. The numerical screen takes 8.64 s wall / 8.56 s CPU and saves the two chosen proposal arrays. Interval verification takes 6.50 s wall / 6.44 s CPU. The optional skew check takes 15.32 s wall / 14.99 s CPU. These are measured on a single-threaded Apple M4. [Numerical screen](screen.json) · [Interval certificate](certificates.json) · [Skewed weights](skew.json) · [Existing project and native geometry](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/)
