# First competition: numerical scope and possible later amendment

No first-competition model/training code was changed during the running competition.

1. Fixed Gauss-Legendre rule exactly integrates polynomial h' and constant moments, but actual labels delta*1{s>t} are discontinuous in t. Thus exact labels at finite8/16 thresholds integrate to a staircase approximation h_K(s), not original h(s). A statistical network minimizing weighted threshold loss converges to the conditional mean of this approximate material, with a persistent discretization term even at infinite training data. Three historical validation parents yielded max component label-integration errors0.279(K8),0.136(K16), with weights summing1. This is separate from finite-regression error and does not prove observed trained risk failure is caused by discretization.

2. Phase initial zero logits are uniform positive masses, yielding constant half-of-maximum legal increments. It differs substantially from the exact anchor initialization of direct residual heads. Update0 selection preserves negative training evidence but does not make initialization/optimization identical. Same5842 local CNN architecture also uses unscaled residual outputs whereas historical R2 learned residuals used a train-fitted residual_scale; frozen R2 inference retains its original scale exactly. Family comparisons must describe these actual training contracts.

3. direct_moment and local_cnn currently have the same local residual material architecture and inputs. They are aliases as an output-role control, not two scientifically independent architectures. A future shared-vs-unshared density encoder experiment is not implemented by those names.

A later, separately frozen remedy can remove deterministic quadrature bias without changing original prior/decoder or increasingK. Partition standardized field threshold domain by bins[a_j,b_j], with positive exact weights w_j=h(b_j)-h(a_j). Define each conditional bin-moment as h'-weighted average of the event moment within the bin. Its single-joint-sample label is

L_j = delta * clamp(h(s)-h(a_j), 0, w_j) / w_j.

Then sum_j w_j L_j = delta*h(s) exactly for every training parent/pixel, including plateaus. Bins ordered in t give nonnegative decreasing labels, so a monotone head remains suitable. This is a bin-integrated event moment, not point-threshold probability; preserve the distinction in naming/calibration. Equal h-mass bins use w_j=1/K and label delta*clip(K*h(s)-j,0,1),j=0,...,K-1; all decoder evaluations must call original source saturation. No independent E(delta)*phase factorization is introduced. Fit real/imag components jointly as before. This proposal is mathematics for a future amendment; it was not substituted into the first frozen competition.

Alternative randomized thresholds T~density h'(t) give unbiased stochastic integration of delta*1{s>T} when integrated over T, but a fixed finite inference rule still requires explicit quadrature-error accounting. The integrated-bin design above is the more direct deterministic identity.
