# Data and scientific scope

Development uses already-published R2 observation obs_b20d5a55677d4eca.
Those results are not new blind samples. Inverse code reads only the R2 public allowlist
(G,E0,S,Y,geometry,noise,material box and prior parameters), never truth/shape/mask.
Independent post-seal scoring can read the already-disclosed R2 truth and must
remain separate from online candidate selection. The separate six-object
confirmation was generated after freezing the one-exchange algorithm.

Physics: exp(−i omega t), lambda0=1, k=2pi, fixed mu/background1, real scalar
contrast m in [.02,20]. j=diag(m)(Einc+Gj), Y=Sj. Existing24×24 disk-pulse Hankel
VIE,12external point sources,48external receivers, complex128 CPU1. This is the
existing finite scalar-TM model, not a 3D vector or continuous-PDE certificate.

Original selection: Phi=.5||F(m)−Y||F²/||Y||F²+R(m). R is the same historical
smooth/TV: weights1e-4 each,width19.5,epsilon.01. State/Gram/port surrogate losses
are recorded separately. Source metrics are explicit. No material or support GT
is used to choose a candidate, direction, cutoff, step or stopping threshold.
The historical normalized objective and noise whitening are distinguished:
Wnoise=I/sigma, while Wphi=I/||Y||. Equivalently the objective is sigma²/||Y||²
times the whitened data term plus R. In a unit-weight whitened convention the
identical prior is ||Y||²/sigma² R, shared by all methods. Noise qualification
uses Wnoise; an auxiliary normalized loss is never substituted for it.

Noise: independent circular CN entries with real/imagvariance sigma²/2. Retain
original per-case sigma and99%white radius25.1654. Scoring retains NRMSE/12,
Otsu256-bin IoU, existing thresholds. These summarize each development object;
number of methods is not number of independent materials.

Candidate masks, geometric blocks and singular cutoffs are declared before their
runs. They are development choices, not probability guarantees or GT knowledge.
Dense and small factorizations, complex RHS, Schur decompositions, Gram/SVD work,
rejected proposals and startup failures are accounted separately. No SHA checks.

## Fresh confirmation boundary

The separate six-object confirmation was frozen and generated once.
See CONFIRMATION_PROTOCOL.md and configs/confirmation_frozen.json; the draft is
retained as preparation history, not the executed protocol.
New public observations use data_public/confirmation; truth remains only in
sealed_truth/confirmation and is excluded from all inversion/remote packages.
Current development evidence above remains GT-disclosed development.
The solver package has no truth, mapping, generator or scoring files. It contains
only explicit source paths and public observation/operator keys. All18 frozen
dispositions must close before the separate scorer opens truth. Two objects per
ring/curved-bar/two-islands family are six material clusters, not18 independent
samples. A three-endpoint same-Phi retained set was declared before generation
and is charged all three arms. Its result is not attributed to onehop alone.

The32x32 replay is separately marked DISCLOSED_DEVELOPMENT_GRID_DIAGNOSTIC.
It retains the same192-raster geometry, sources and receivers but changes pulse
grid, block physical width and noise realization. Prior numerical parameters
are retained, not asserted to have an invariant continuum scale. It tests
representation sensitivity and is not a fresh confirmation object or rigorous
continuous-PDE convergence bound. No3D vector experiment is implied by this2D work.
