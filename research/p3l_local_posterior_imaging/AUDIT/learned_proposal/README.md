# Supplied-array learned full proposal

`src/p3l/learned_proposal.py` has no campaign file reads, Book, physical backend or trainer. `LearnedFullProposal` accepts a fixed normalized weak/probit network, one supplied scalar/image context, a complete orthogonal frame `Vw,Vs`, anchor `q0=(z0,u0)`, paid whitened Jacobian `J` and residual `r`.

All proposal methods use **x=(z,eta)** and density measure **dz d_eta**, with `u=ndtr(eta)`. The network receives raw `(a,eta)`, where `a=Vw.T@(z-z0)`. `b=Vs.T@(z-z0)`. No additional probit or property Jacobian is applied to the network's eta density. The full frame has absolute determinant one. No clipping or changes to priors/network are performed.

The strong reference has precision `H=I+Jb.T@Jb` and mean `-H^-1@(Vs.T@z0+Jb.T@(r+Jw@a+Ju@(u-u0)))`. It retains the complete correlated covariance `H^-1`. `learned_log_prob(x)` evaluates the normalized product; `sample(count,rng=...)` and `log_prob(x)` use the same defensive mixture, with default weights **0.9 learned product + 0.1 N(d+p)**. The default campaign dimensions are d=64, k=8, p=2; generic dimensions enable analytic controls. All methods return NumPy arrays/scalars; torch network calls run without gradients. Put the supplied torch network in evaluation mode first.

`conditional_field_moments(a,u,matrix,offset=...)` and `conditional_material_moments(...)` provide pixel marginal moments conditional on fixed outer a,u under this Gaussian strong reference, using the existing moments module. The material increment is supplied by the caller at that u. These functions retain b correlation through matrix propagation; they do not supply spatial covariance or posterior-certified moments. Outer integration must preserve a/u dependence and use total variance.

`to_chart(x)`/`from_chart(a,b,eta)` provide the complete orthogonal chart. `GaussianWeakDensity` is an independent NumPy Gaussian control with the network's single-context shape API. NumPy proposal randomness and network randomness use separate streams: callers seed both their NumPy generator and torch (or Gaussian control RNG) when reproducibility is required.

The parent wrapper remains responsible for paid full likelihood and exact MH ratio. In these coordinates the true unnormalized target is `exp(-Phi(x))*N(d+p)(x)`; proposal-reference moments and normalization controls do not establish posterior accuracy.

Run `analytic_controls.py` with the existing campaign NumPy/scipy/torch environment. It writes `RESULTS.json`, comparing exact Gaussian and defensive-mixture densities, inverse chart, synthetic sampling moments, correlated strong field/material moments, torch interface and a synthetic 66D prior identity. It reads no campaign data and runs no physics/training.
