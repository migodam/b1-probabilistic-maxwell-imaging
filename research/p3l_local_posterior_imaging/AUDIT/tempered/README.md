# Generic adaptive tempered SMC

`src/p3l/tempered.py` imports NumPy and SciPy only. It never reads campaign data, instantiates a solver/Book, imports an old ledger, starts training or changes the prior. The parent supplies all densities, coordinates, wrappers, cache ownership, budgets and stopping decisions.

```python
from p3l.tempered import tempered_smc
result = tempered_smc(sampler_prior, logprior, loglikelihood,
    rng=rng, N=128, mutation_steps=2,
    mutation="mala", step_size=.1,
    gaussian_dimension=66,
    gradloglikelihood=gradient_likelihood,
    precision=deterministic_precision,
    callback=save_event)
```

`sampler_prior(rng,N)` returns finite `(N,d)` states from the normalized declared prior law. `logprior(x)` and `loglikelihood(x)` return scalar log densities, finite or `-inf`; `NaN`/`+inf` raise. A sampled point outside the supplied prior support is an explicit error. No implicit transformation occurs: exact probit coordinates require the parent to provide the correctly transformed prior, target and gradients. The reported evidence is the likelihood integral against the sampled normalized prior, with likelihood constants exactly as supplied. The module does not assert that omitted likelihood constants are free evidence constants across different models.

Optional `gaussian_dimension=d` declares a standard Gaussian reference in all coordinates, including exact probit Gaussian coordinates. `gaussian_dimension=d-2` declares Gaussian field coordinates plus two explicit native independent uniform properties in `[0,1]`. This declaration permits a known prior gradient `[-z,0,0]`; otherwise MALA requires an explicit `gradlogprior(x)`. Prior sampling remains the parent's responsibility.

The tempering target is `prior(x) * likelihood(x)^beta`. For normalized incoming weights w, the next beta is the largest step satisfying conditional ESS `N*(sum(w*G))²/sum(w*G²) >= .8*N`, `G=exp(delta_beta*loglikelihood)`. This is **conditional ESS**, distinct from ordinary final weight ESS. Bisection never invents extra samples or lowers the threshold. Log offsets are removed for stable CESS; beta zero avoids `0*(-inf)`. If the largest safe increment is smaller than `min_delta_beta` (default `1e-8`), the run returns `TEMPERING_STALLED_BELOW_BETA_FLOOR` at its actual beta. A safe final step to one may be smaller than that floor. A maximum of 256 stages is enforced by default; exhausting it returns `MAXIMUM_TEMPERING_STAGES_REACHED`. Neither failure is a posterior result at beta one.

Each stage records the normalizing log-evidence increment `logsumexp(logw + delta_beta*loglikelihood)`, adds it to cumulative evidence, and normalizes new logweights. Ordinary weight ESS triggers systematic resampling when below `.5*N` by default; the threshold is configurable and recorded through call configuration/history. A single systematic offset produces parent indices. Root prior ancestor IDs are carried through every stage and mutation; history records both immediate parent indices and root IDs. Resampling resets weights to uniform, while mutation preserves the assigned weights.

Mutations are fixed within a run:

- `random_walk`: a symmetric Gaussian increment, with symmetric reflection only for explicitly declared native uniform properties.
- `pcn`: Gaussian-reference proposal plus a symmetric native property reflection when present. Its full MH ratio includes the Gaussian-reference proposal-density term; likelihood-only acceptance is not assumed for an undeclared arbitrary prior.
- `mala`: drift `(h/2)*H^-1*(gradlogprior + beta*gradloglikelihood)` and covariance `h*H^-1`. `precision(x,beta)` must be deterministic, finite, symmetric and SPD; identity is used if absent. Both endpoint normalized proposal log determinants and quadratic forms enter MH. No divergence drift or ideal Riemannian-diffusion claim is made. Native-box MALA uses full Gaussian proposals and rejects points outside prior support; it never reflects/clips a Langevin proposal. A Gaussian-reference declaration does not replace the supplied logprior in acceptance ratios.

Zero-likelihood particles at positive beta have zero target mass and cannot require an undefined log-likelihood gradient. Those particles receive a symmetric random-walk move with exact zero-target handling; this is not a state-dependent mixture on the positive target support. Proposed zero likelihood from a positive-target state is rejected before reverse gradients/precision. When all weighted likelihoods vanish, tempering fails explicitly rather than fabricating a normalized population.

Only one particle's current likelihood gradient/metric is retained during its consecutive mutation steps. Arrays of state coordinates, scalar log values, normalized weights and ancestry are retained globally; no per-particle physical Point/LU cache or giant metric array is owned here. Parent wrappers may implement and account for their own physical caches.

Every likelihood invocation emits `LIKELIHOOD_EVALUATION_ATTEMPT` and `LIKELIHOOD_EVALUATED` events, and every gradient/precision invocation has an attempt event. Counters separately record prior evaluations, likelihoods, likelihood gradients, supplied prior gradients and precision calls. Callbacks also receive initialization, mutation outcomes, stage completion and final status. Initialization/stage completion includes numeric particle/log-value/weight arrays for parent persistence; arrays may contain `-inf` for zero support and should be written to array artifacts. Events expose cumulative snapshots rather than additive cost estimates. These counters measure callable invocations, not RHS; all actual paid costs remain with the supplied wrappers. Callback/budget exceptions propagate without an automatic replay or free fallback. The last completed stage can already have been persisted by the callback.

The returned dictionary contains terminal beta, status, weighted particles, normalized weights, cumulative evidence, genealogy/history and counters. A completed beta-one run is labeled `COMPLETED_FINITE_N_UNQUALIFIED`, with `posterior_qualification=False`. There is no inference that N particles are N independent posterior samples. Adaptive schedules/resampling introduce finite-N error and can affect evidence bias; **unbiasedness of this adaptive evidence/log-evidence is not claimed**. Large ESS does not exclude missed modes, shared bias or poor mutation mixing. Independent repeated populations and parent-owned material-function comparisons are required for any physical reference decision.

## Independent controls

`check_tempered.py` runs eight independent 1,024-particle Gaussian linear-observation populations with variable-metric MALA, comparing weighted moments and likelihood evidence to closed forms. It separately checks prior invariance under random walk, pCN and variable-metric MALA, native uniform properties, normalizer increment accumulation, normalized weights, ancestry, systematic resampling, callback/invocation counts, beta-floor failures, zero likelihoods and stage-cap failures. These are small probability controls with zero campaign or physical calls. Their fixed tolerances are not scientific qualification gates.

Reproduce from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 ../../.venv_nn/bin/python research/p3l_local_posterior_imaging/AUDIT/tempered/check_tempered.py
```

`results.json` records the independent analytic evidence and limitations. Passing does not establish Maxwell posterior exploration, material-risk improvement, solver cost savings or equality of two physical mechanisms' budgets.
