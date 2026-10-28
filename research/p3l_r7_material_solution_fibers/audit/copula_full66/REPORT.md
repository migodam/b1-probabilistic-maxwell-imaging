# Controlled discrete nonuniform Full66 fibre coupling experiment

**CONTROLLED_DISCRETE_FIBRE:** four deterministic matched points per parent, each weighted 1/4, versus the 4×4 product of their field/property marginals, each weighted 1/16. These are not original-prior conditional draws or posterior samples. All 64 crossed coordinates are legal under the unchanged original Full66 decoder. New Maxwell RHS / states / LU / NN training = 0.

Rows use spatial coordinates from curve point i; columns use Gaussianized property coordinates from curve point j. The diagonal retains the numerically same-material coupling. Off-diagonal crossings are fixed before decoding. The finite decoder prediction is chi_ij = host + (alpha_j/alpha_i)(chi_reference-host), with the finite equality residual measured separately. Thus the elementary finite-sum mean factor is E[alpha] E[1/alpha]. No new general theorem is claimed.

| Parent | Ealpha E1/alpha − 1 | Mean bias, material RMS | Product variance, RMS² | Bias prediction error, RMS |
|---|---:|---:|---:|---:|
| R7_E1_092 | 6.33503518e-05 | 4.41351173e-06 | 6.15046105e-07 | 1.62e-17 |
| R7_E1_201 | 5.27417115e-09 | 1.0211934e-09 | 3.9545018e-10 | 3.43e-17 |
| R7_E1_233 | 0.00265788182 | 5.68146157e-05 | 2.44253946e-06 | 1.46e-17 |
| R7_E1_251 | 2.70532361e-06 | 1.32549928e-07 | 1.29888962e-08 | 1.47e-17 |

The matched-curve material variance remains at numerical roundoff. The independent product produces nonzero material variance and the predicted mean bias in these finite constructions. This demonstrates what the specified destruction of coupling does in this controlled ensemble; it does not establish original-prior conditional weights, a posterior effect, a causal diagnosis of R6 results, or any scattering consequence.

Full x, native q, 576-cell complex materials, field/fraction arrays, means, per-cell complex variance, Re/Im covariance, weights, and prescribed pairs are in `CONTROLLED_DISCRETE_FIBRE_ARRAYS.npz`. `FIXED_ENDPOINT_CANDIDATES.json` contains the predeclared first/last crossings and one matched control per parent for possible later separately billed physics. Their Maxwell status is NOT_RUN.

The beta−1 calculation uses a nonnegative pairwise finite-sum expression to avoid cancellation for narrow alpha ranges. All raw alpha values, direct Ealpha E1/alpha, pointwise ratio residuals and covariance checks remain in `RESULTS.json`. The source curves and original E0/E1 evidence are read-only.
