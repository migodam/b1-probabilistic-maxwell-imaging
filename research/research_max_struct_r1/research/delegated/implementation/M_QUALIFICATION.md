# M representation mechanism: NOT_RUN

The credible implemented backend is an outgoing finite vector CM+RR DDA. For
the radius 0.12 lambda, epsilon=3, mu=1 sphere, refining h=.04 to .03 lambda
(123 to 257 dipoles) reduces absolute complex electric far-amplitude error
against independent analytic Mie from 10.468% to 3.655%, and scattering
cross-section error from 21.808% to 5.989%. These are useful bounded physical
controls; they do not qualify internal H(curl) or Maxwell weak residual norms.
The respective staircase volume biases are +8.76% and -4.13%.

The adapter provides cell-center E/H/J and external electric/magnetic kernels.
Its analytic H kernel is the curl of the reused dipole electric radiation
kernel. The symmetric cell self curl is set to zero. A globally reconstructed
point-dipole field has singularities at source locations; off-source analytic
derivatives do not determine a consistent material-interior weak curl norm.
Consequently H snapshots cannot simply be interpreted as samples of a
compatible finite-element H(curl) space.

The following required contract is missing:

1. Compatible edge/face spaces and discrete gradient/curl/divergence, with
   separate primal and dual degrees of freedom, boundary traces, and any
   harmonic modes retained.
2. Electric and magnetic mass matrices representing the required physical
   weighted norms, including material interfaces and small gaps.
3. A weak Maxwell residual and its dual norm on the stated H(curl) test space;
   H(curl) alone does not control a strong second-curl L2 residual.
4. An external boundary/receiver map consistent with those edge/face traces,
   incident sources, radiation condition, and the derivative/residual norm.
   The existing external dipole receiver map is a valid finite-DDA map, but
   does not complete that compatible-space contract.

For planned gaps 0.06/0.12 lambda, h=.03 gives only 2/4 nominal cell samples
across a gap. Center-inclusion geometry and point-field collocation at that
resolution do not demonstrate trustworthy interface/gap H(curl) behavior.
The sphere's far-field refinement result cannot validate those missing norms.

The M primary experiment therefore remains **NOT_RUN**. No L2 singular-value
tail is substituted for the requested metric, no claim that H is easier to
represent is made, and no imaging mechanism extension follows. Replacing the
backend with compatible edge/face FEM or a qualified finite-volume vector VIE
would be a new implementation and qualification task, beyond this bounded run.

The F tensor-fiber experiment can still use this backend with its exact
finite-DDA contract: epsilon is physical tensor input; Eint is distinguished
from Eexc; j=(epsilon-I)Eint and Eint=Einc+G j. A tensor update annihilating all
local Eint columns preserves the discrete j, Eint, and external dipole field
exactly. Such a result is an algebraic ambiguity of the finite model. It is
not, by itself, a continuous material-interface Maxwell certificate or proof
of a continuum indistinguishable tensor family. Parent scientific review
must maintain that distinction.
