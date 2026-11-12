# Physical conditions and evidence boundaries

This record describes the frozen physical inputs and previously completed qualification controls. It contains no reconstruction-quality review from the new C2 campaigns. No new 3D truth arrays, measured heldout values or fresh reconstruction results were opened. Final scientific classification remains with the parent after the required prediction seals.

`PHYSICS_CONDITIONS.csv` contains 14 campaign acquisition-condition rows: six fresh scalar conditions, four vector conditions, two measured conditions and two classical conditions. Its additional 12 rows are pre-campaign analytic Mie controls, explicitly excluded from inversion-object counts. Each row records acquisition, real material/current dimensions, normalization, ball/noise convention, Green/self term, prior, independence and limitations.

## Physical model and acquisition

All inverse materials are real scalar susceptibility chi=epsilon_r−1. The scalar-TM model uses one complex current component per cell. The full-vector 3D model has three coupled x/y/z current components and one shared real isotropic chi per cell; 216 material variables must not be reported as 648 independent material variables. Its J is unscaled susceptibility current chiQ; the electric-current convention in the original vector backend is −ikJ. Propagation and receiver maps retain their original volume factors.

Radiation is represented by outgoing free-space integral Green functions. **No PML is used.** FFT zero-padding evaluates nonperiodic convolution; it is not an absorbing layer. Scalar disk-pulse diagonal entries are `i*pi*k*a*H1^(1)(ka)/2−1`; vector CM+RR self-density is `−1/3+i*k³h³/(6*pi)`. The scalar fresh-cohort disks use a = .45h. Measured/classical models use the distinct equal-area quadrature a = h/sqrt(pi), whose disks are integration surrogates rather than disjoint physical objects.

| Family | Actual acquisition | Unknown material | Independent-object accounting |
|---|---|---:|---|
| Fresh scalar C2 | Cylindrical **line sources**, E0=(i/4)H0^(1)(kr); 6 or 12 sources at 1.4 lambda0; 24 or 48 receivers at 1.1 lambda0; .96 lambda0 square | 576 real cells, 24² | 48 distinct objects; six acquisition/noise conditions × 8. Mechanism/local-reference subsets reuse these objects. |
| New vector 3D | **60 external electric dipoles**, 20 on each of three offset circles of radius .6 lambda0; three electric components at 60 receivers; .24 lambda0 cube | 216 real cells, 6³; 648 complex current components/source | Four conditions but **three independent nonsphere objects**. The 8³→6³ frame repeats the first geometry. |
| Fresnel2001 single/twin | Actual 2 GHz horn measurements; empty-target outgoing-multipole continuation; 27 training sources with 39 actual receiver rows/source | 576 real cells, 24²/.16 m square | Two measured specimens. **1053 training + 711 heldout** actual channels; no pixel material GT. |
| Classical Austria | **20 plane waves**, 40 receivers at 3.75 m; 400 MHz ; 2 m square | 1024 real cells, 32² | Two conditions of **one object**: 150²→32² and matched 32². The matched control is post-exposure diagnosis. |

The fresh scalar sources and vector dipoles are not plane waves. Plane-wave acquisition belongs to the Austria contract and to the separate analytic Mie reference. The fresh cohort contains 36 matched 24² observations and 12 generated 36²→24² observations. Crossgrid placement is not balanced across acquisition/noise: all eight 12-source/48-receiver/1% cases are crossgrid, four of eight 6-source/48-receiver/2% cases are crossgrid, and the other conditions are matched. Different condition groups use different objects; group differences cannot isolate a causal aperture, noise or discretization effect.

## Scaling, noise and feasibility

Scalar candidate selection uses the original `Phi=.5||Yhat−Y||F²/||Y||F²+R`. Vector selection retains the different inherited `Phi=mean_s ||Yhat_s−Y_s||²/||Y_s||²+R`, without a factor 1/2. Scalar prior parameters remain smooth=tv=1e−4, epsilon=.01; vector prior remains 2e−4 mean-neighbour squared difference plus 1e−6 mean chi². The lifted projection metric is `||J/6||F²+||Q||F²`; its unwhitened global data ball is distinct from these selection objectives.

Fresh scalar noise is iid circular complex Gaussian at 1% or 2% clean-data RMS; new vector noise is 1%. Their raw noise-only 99% radii use the declared sigma and complex chi-square event. Such radii do not include discretization error. The vector frame's recorded 8.006% clean crossgrid difference already exceeds its 1% noise. This does not certify coarse true-material feasibility or justify a posthoc radius change.

Fresnel uses the fixed **15% training-norm model-discrepancy ball, not CN99**. Actual noise covariance/sigma and receive-horn response remain unresolved. The stored sigma field is only a discrepancy-per-channel scale. Previously checked affine numerical-rank floors are 1.865% and 3.252% of training norm; 15% balls are affine-reachable and zero balls fail. Affine reachability is not feasibility for a real material. Empty-field calibration validation medians are 2.515% and 2.408%. The 711 closed heldout channels comprise 441 on withheld sources and 270 withheld receiver channels on training sources. They may be evaluated only after authorization, without absent-channel zero filling, new gain fitting or candidate selection.

Austria uses historical 10% fixed-norm complex noise and the unchanged 12% observed-norm discrepancy ball, not an iid-CN99 event. The previously recorded coarse target residual is 22.526% for 150²→32² and 9.951% for matched 32². Thus the original mismatch condition and matched feasible-target diagnostic must remain separate. Lower residual or membership of either ball does not establish material recovery.

## Qualification and algorithm names

The original dyadic backend has independent analytic Mie controls, not a general continuum certificate. At fixed radius .09 lambda0, epsilon 9, staircase h=.02/.015 far-vector errors were 9.015%/7.064%; corresponding volume-averaged susceptibility errors were 19.76%/16.93%. Coarser negative controls remain in the CSV. Near-machine-precision linear solves establish discrete consistency; they do not remove boundary-material or continuum error. These figures must not be converted into a seven-percent accuracy promise for the new nonsphere inverse 6³ experiments.

`global_dual` denotes **global-PR Dual-CIE**: the auxiliary current update combines source columns in its PR step. It is distinct from local-reference `physical_globaldual`, the existing **physical GlobalDual total-field source-Gram objective at beta=0**, with initial physical reference held fixed. The latter retains original-Phi candidate selection and the common original corrector. Beta=.2 feedback is another objective. Local-reference jobs reuse 16 fresh objects on a separate host/budget and do not enlarge object counts or supply CPU times that can be pooled with the main host.

Measured CIE/Dual retain each source's receiver SVD/DPCS and native projected Fourier APCS. Vector controls reuse the existing native 3D implementation. Classical controls retain DPCS 15, mf 5/8/10, material beta 8/2/.5 and Dual/global current beta 8; original 150→32 allows 100 iterations/stage and matched diagnosis 300. Cold mature classical runs preserve paperloop stage stopping. G19 uses its declared physical-warm loop/global stopping with the same equal-area physics; full 32² masks are 16. Neither scalar common-port G19 transplanted into measured data nor scalar G19 relabeled as vector is claimed.

## Existing reference provenance

The cached original full text identifies Yu Zhong and Kuiwen Xu, *Contraction Integral Equation for Three-Dimensional Electromagnetic Inverse Scattering Problems*, Journal of Imaging 5(2),27, **8 February 2019**, [DOI 10.3390/jimaging5020027](https://doi.org/10.3390/jimaging5020027). The 2019 paper's3D grid/domain, continuation and two-beta experiments precede this work; current reduced grids/domain/acquisition are adaptations, not a literal reproduction. Date/source were checked in the already available R7 `delegated/public_literature/CIE.txt` and R8 `PAPER_CONTRACTS.json`; no new literature retrieval was performed.

The cached NIE 2016 contract is [DOI 10.1109/TAP.2016.2535492](https://doi.org/10.1109/TAP.2016.2535492):150² forward/64² inverse, rank 15, zero initial contrast/APIC and no additional prior. R13's 32² inverse, common prior/corrector, finite real box and equal-area quadrature remain explicit adaptations. Actual 2001 measured package provenance is the official supplementary record [DOI 10.1088/0266-5611/17/6/301](https://doi.org/10.1088/0266-5611/17/6/301), with original local paths retained in the portable-input metadata. No claim of author-certified code or measured pixel truth follows from reuse.

Sources for the table are the public fresh48manifest, frozen physics manifest, portable input metadata, original acquisition code and pre-campaign qualification receipts. New C2 material/heldout scoring and final scientific conclusions are intentionally absent.
