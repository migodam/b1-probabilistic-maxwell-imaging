# Noise-locked contrast sources and auxiliary material paths in high-contrast Maxwell inversion

**Research manuscript draft.** This document reports mechanism evidence and limitations, not a claim of a publication-ready universally superior solver. The Chinese report and separate classical/vector benchmark reports are the authoritative experimental ledger.

## Abstract

Spatial transport of constitutive residuals can change the trajectories of contrast-source inversion, but a cheaper conditional material minimizer does not by itself remove nonconvexity of the physical inverse problem. We examine this distinction using A-material, its Möbius-based A+ readout, ordinary sourcewise and global CIE/FFT-TSOM, and direct physical optimization. In a finite two-dimensional Maxwell model, noise amplification in the hard deterministic-current subspace explains failures shared by all material readouts. A public-noise discrepancy rule increases recovery from two to four objects on a frozen six-object confirmation set for both ordinary CIE and A+. A+ provides no additional successful object and has higher average computational cost. A full-current relaxed initialization, motivated by established CSI/WRI formulations, recovers one previously failed development shell; a single weak-penalty stage followed by physical correction preserves that recovery at lower cost. Its advantage does not transfer to the confirmation set. These observations identify distinct representation, noise, and optimization failures and limit claims based on isolated high-contrast successes. All candidate selection uses the same physical objective and no test truth.

## 1. Introduction

High contrast is not a single experimental axis. Electrical size, internal wavelength, topology, illumination, noise, and numerical discretization all influence the inverse problem. Consequently, a successful difficult curve and an unsuccessful ring cannot establish a general ordering of algorithms. We ask which part of CIE is changed by a transported material residual, whether more accurate material readout addresses the limiting failure, and whether changing the auxiliary optimization path provides a different capability.

Our contributions are an auditable failure decomposition, a noise-consistent repair shared with strong CIE controls, and a negative transfer result for complementary auxiliary paths. We do not claim a new discrepancy principle, variable-projection theorem, or global convexification algorithm.

## 2. Related work

The [2016 NIE formulation](https://doi.org/10.1109/TAP.2016.2535492) and [2019 vector CIE study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8320915/) combine contraction-form residuals with FFT-TSOM subspaces and staged updates. Dual choices of the current and material contraction parameters are already present in that literature. Their published comparisons to LSIE are not interchangeable with comparisons to a modern fully converged physical optimizer under another acquisition contract.

Relaxing the wave equation and eliminating the state are established strategies; see the [penalty-method analysis](https://doi.org/10.1088/0266-5611/32/1/015007). A [recent constraint-relaxation analysis](https://doi.org/10.1088/1361-6420/adab89) also gives special convex finite-dimensional constructions. The matching source/receiver-span assumptions of that construction do not hold for our full internal Maxwell state space. A source-resolved claim ledger is provided in `delegated/nonconvex_prior_art/`.

## 3. Methods

For real contrast m and fixed permeability, let L(m)=I−diag(m)G, J(m)=L(m)⁻¹diag(m)E and F(m)=SJ(m). Final candidates minimize the common physical objective Φ=||F(m)−Y||²/(2||Y||²)+R(m), with the same box and prior. A-material transports the fixed-current constitutive residual through local Green couplings. A+ uses the conditional Möbius coordinate to reduce the cost of this material readout. Global A/A+ is therefore compared to ordinary Global-Dual current updates, while sourcewise A+ requires a sourcewise Dual control.

For S=UΣV*, a hard rank-r current readout has expected squared noise error Ns σ² Σ(i≤r)σᵢ⁻² under independent circular complex noise with E|n|²=σ². We select the smallest contiguous prefix for which the discarded measured energy lies within one declared 99% full-data noise ball. All prefixes use the same event; the complete data remain in the subsequent objective. This changes the hard current constraints, not the number of material unknowns or the physical observations.

The auxiliary-path comparison minimizes ||SJ−Y||²/(2yd)+λ||L(m)J−diag(m)E||²/(2ed)+R(m). The current is eliminated by an exact positive-definite solve. A three-stage schedule is compared with a fixed penalty, equally repeated fixed-penalty restarts, data-target continuation, and direct physical optimization. A later development-only simplification hands the first weak-penalty material directly to the physical corrector.

## 4. Theory and scope

A strictly increasing scalar transformation preserves local minima. A local diffeomorphic change of material coordinates preserves Hessian inertia at stationary points. Thus a Möbius coordinate or a monotone transformation of the physical loss cannot alone remove its erroneous local minima. Auxiliary residual transport can nevertheless change intermediate objectives and finite-step trajectories.

Eliminating the relaxed current yields a residual metric proportional to (yd I+(ed/λ)AA*)⁻¹, A=SL(m)⁻¹. This metric depends on material; state-wise convexity is not material-wise convexity. The small-λ limit can also suppress information. Detailed assumptions, derivatives and boundary caveats are in `THEORY_DERIVATIONS.md`.

## 5. Finite-model experiments

The first phase comprises 47 development inversions and 66 confirmation inversions with 576 free real material parameters. Six independent confirmation objects cover ellipse, multiple disks, annulus, comb, star and perforated geometry, with different declared contrasts, noise and apertures. Predictions and public-objective portfolio choices were sealed before confirmation truth was read. Six objects, not 66 runs, are the statistical units. This is a small same-discretization confirmation, not experimental-data or continuum evidence.

The fixed-rank single CIE, Dual, Global-Dual and A+ each recover 2/6. Their noise-adapted versions each recover 4/6. Mean NRMSE is .18666, .18882, .18766 and .18552 respectively; corresponding mean inclusive CPU times are 38.84, 40.81, 42.98 and 46.37 seconds. A+ adds no recovered object. Its mean improvement over noise-adapted Global-Dual is dominated by a still-unsuccessful annulus, while its two ordinary successful controls are slightly worse in continuous error.

The relaxed-current continuation recovers the public 40ac shell with NRMSE .05807. A weak-stage handoff gives .05835 at 19.63 CPU seconds including initialization, versus 28.10 seconds for the full schedule. Neither yields an additional confirmation recovery: both recover 2/6. Paid two-branch portfolios remain at 4/6. Restart-only and negative-curvature probes do not reproduce the development recovery.

## 6. Classical and vector benchmark extension

The companion extension distinguishes original Austria acquisition/grid contracts, scaled full-vector tests, and the original 3λ vector configuration. It explicitly corrects two inherited implementation confounds: A must not also change the current metric when testing material-only transport, and overlapping Fourier coefficients must be transported when the basis grows. The published CIE stopping quantity and error both contain square roots; an audit transcription error was corrected before any new vector inversion. Completed results and execution limits are recorded separately in `CLASSIC_BENCHMARKS.md` and the vector reports; no unexecuted example is claimed here.

## 7. Discussion

The results distinguish three phenomena. The historical curve demonstrates a useful A-material trajectory under its development contract. The noisy confirmation cases demonstrate a repair to the common deterministic-current constraint that ordinary CIE also benefits from. The shell demonstrates a complementary path from state relaxation. None establishes dominance across all three phenomena.

A smaller residual can coexist with a worse material image. Confirmation truth supplies offline witnesses that some failed endpoints have much larger physical objective than an existing feasible material; this rules out explaining those endpoints solely by absence of a better feasible point. It does not establish uniqueness, certify a local minimum on the full constraint cone, or prove that an optimizer will reach the witness.

## 8. Limitations and conclusion

The confirmation sample is small, the first phase uses the same discretization for generation and inversion, and the proposed repairs have established methodological antecedents. Material-only A+ cannot be directly generalized from a real scalar arc to an unknown complex reciprocal tensor. CPU, full-wave solves and auxiliary dense factorizations must all be charged; fewer physical LU calls alone do not establish lower total cost.

The evidence currently supports a reproducible mechanism and failure study. It does not yet support an algorithm paper claiming broad superiority over mature CIE. The next substantive question is how to preserve observable material directions while relaxing erroneous current constraints, rather than how to make one conditional material minimizer increasingly exact.
