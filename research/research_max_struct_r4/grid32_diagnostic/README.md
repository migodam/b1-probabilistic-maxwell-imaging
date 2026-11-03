# 32×32 necessary-grid diagnostic — prepared only

Classification: `DISCLOSED_DEVELOPMENT_GRID_DIAGNOSTIC`.
No32 observation, forward solve, inversion, new result or remote work has been
performed during preparation. The existing24 geometry replay is confirmed by
`audit/GEOMETRY_REPLAY_PREPARATION.json`: exact alpha equality, maximum error0,
correct original family order and exact original-generator RNG consumption.
No frozen confirmation SOURCE_FILES or runner are changed.

`generate.py` is the only GT-reading component. It reads the already disclosed
R2 mapping/seed and the existing curved-bar alpha for identity. It replays
family shuffle `two_islands→curved_bar→ring`; after each geometry the original
noise consumes two48×12 normal arrays. The copied192 geometry implementation is
checked against the original function including RNG state, and its24 area
average must equal the original target exactly. If identity fails, generation
stops rather than guessing or searching seeds. Current replay passed.

The SAME192 binary support is area averaged in6×6 cells into32×32 material,
`.1+12*alpha`. Original domain endpoints are±.48 (old centers±.46, h.04,
radius.018). New h=.03, centers±.465, radius.0135=.45h. Original12 sources,
48 receivers and k are reused. `build_operators` is the existing disk-pulse
Hankel scalar-TM implementation; no alternative backend or mesh model is added.
The new model has1024 independent real material cells. The generator accounts
for operator assembly and one retained-LU forward/all12 source RHS explicitly.

New observation noise is CN with independent real/imag variance sigma²/2,
sigma=.001||Yclean||/sqrt(576), fixed declared development seed202611030032.
Both32 methods use exactly this one observation. Truth/alpha/mask/support/cleanY/
noise realization remain in `development_truth/`; inverse code reads only
`data_public/grid32_curved_bar.npz`, whose exact public keys are checked.
This is a disclosed same-object development diagnostic, never a fresh blind test.

`solver.py` is an independent32 adapter using only the existing MaxwellModel,
`run_onehop_cie`, R2 `global_cie`, and R4 `refine`; it never imports the generator
or reads original/new truth. API is `solve(public,method,config,backend='cpu')`.
Only onehop/global are configured. Both retain mf3/6/10 for500 updates each,
beta1 then current1/material.2, no coefficient early stop, original material
inner settings, same common host.1, box[.02,20], and original prior weights,
width19.5/epsilon.01. Candidate selection uses original normalized physical Phi
plus prior, common host and true stage endpoints, followed by one shared
corrector1000 (same original scaled-Phi R4 implementation). Retain host/seed/
refined incumbent. There is no mf12 full-space competitor or parameter search.
Fourier masks have the same declared mode capacity; actual APCS ranks and stage
history are recorded. No scientific wallcap is introduced.

Physical limits of interpretation:

- Two finite disk-pulse discretizations are compared; this is not a continuous
  PDE convergence study or an electromagnetic material truth certificate.
- Area averaging, pulse radius and G/E/S quadrature change with grid. The fixed
  two-by-two-cell onehop blocks shrink in physical width from.08 to.06.
- Original numerical prior weights AND width19.5 are kept, so no claim is made
  that the regularizer defines a matched continuum prior across grids.
- New fixed-seed noise differs from historical24 noise; this necessary diagnostic
  measures sensitivity of the whole stated setup and cannot isolate a pure
  discretization effect from every coupled change.
- Auxiliary iterations and mode masks are matched; dense1024-cell physical work
  is more expensive than576-cell work. Report actual setup/LU/source RHS,
  auxiliary scalar/small-block operations, corrector cost and failed attempts.
- If the24 advantage disappears, report sensitivity or a possible discretization
  artifact. Do not enlarge masks, tune beta/prior/noise, retry seeds or expand
  budgets to rescue the advantage. Independent parent scoring/claim judgment
  follows prediction sealing; no scorer is supplied here.

Later executable commands, from this study root and existing research Python:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" grid32_diagnostic/generate.py --generate
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" grid32_diagnostic/solver.py --method onehop
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" grid32_diagnostic/solver.py --method global
```

All generation/result directories are write-once. These commands are prepared,
not executed. CPU1 per solve; no dependencies installed, CUDA contexts opened,
remote resources reserved, or persistent services started. Local standalone
caller receipts include data-read/output time; imports are outside their timer
and explicitly excluded rather than implicitly charged. An unchanged existing
remote worker can later call the module API if the parent chooses to package it.
