# Crossgrid8: implementation-ready plan, pending root review

Status: **DRAFT / NOT EXECUTED**. No new observation, fine material or inverse result has been generated. No main48 private metadata/GT has been read. This independent supplementary8-case cohort does not replace the frozen48 primary cohort, and no scientific method/phase/hyperparameter is selected from it.

## Exact APIs and physical scope

Existing geometry generator is `research/research_max_struct_r10/src/twod_materials.py::draw_material(rng,family,regime,high_chi,side=48)`. Its grid-size parameter already exists. It samples fixed family parameters once, then renders the material at cell centers; no forward-based rejection/reroll. Reuse the existing `research/research_z_alg_r4_diversity/src/physics.py::Model(geometry)` and `.forward(m)['y']`: internally `build_operators` uses the original disk-pulse Hankel off-diagonal kernel and analytic disk self term, and each full forward uses one2304-complex LU with all12 incident RHS. This is the same physical implementation behind the existing MaxwellModel, with no new backend.

Coarse G/E0/S/centers/source/receiver geometry is copied verbatim from `research/research_max_struct_r10/data_public/twod/operators.npz` (the original48 generator uses this source). Fine source/receiver coordinates and k are copied from that public package; only the material grid centers and radius change. Both grids cover0.96×0.96, wavelength1 (k=2π), sources12 and receivers48 at their inherited positions.

Coarse spacing h=.04 and pulse radius a=.45h=.018; fine h=.02 and a=.009. Coarse disk area is 0.001017876020, fine disk area 0.000254469005; four fine disks exactly have the same total area as one coarse disk. Disk fill fraction is π(.45)²=0.636172512352, not a filled square or an overlapping equal-square-area disk.

The scored target is explicitly `m_fine.reshape(24,2,24,2).mean(axis=(1,3)).ravel()`. Equal fine-disk areas make this the same-total-polarization-volume first-order material projection. Pulse support positions and self/interactions still differ. This is **fine disk-pulse data versus coarse disk-pulse model mismatch**, not identical continuous geometry, an exact coarse truth for scattering, a strict discretization-error bound, or continuous-PDE consistency evidence.3D remains a separate experiment/model.

## Frozen candidate contract

Eight new case identities are `r11_crossgrid8_0000`..`0007`, ordered annulus high/mixed, slotted_shell high/mixed, bridge high/mixed, multi_island high/mixed. SeedSequence202611071118 spawns8 independent case children. Each case child spawns separate geometry and noise substreams. Thus high/mixed shapes are independently drawn, not one geometry relabeled; all8 are independent of earlier seeds/cohorts. Fixed high12, mixed high12/low2.5, host.1, box.02..20. No reroll or choice by residual/truth.

Observed Y is generated only on48×48; sigma=.01||Y_fine_clean||F/sqrt(12×48), additive CN noise sigma/sqrt2*(Nreal+iNimag). Public sigma remains the raw measurement-noise sigma, not an inflated model-mismatch allowance. CN99 failure is allowed and is not a reason to redraw/tune. Each method receives only copied24-grid operators/geometry, fine-generated noisyY, rawsigma/noise-radius, legal24-grid/box metadata. Fine materials/labels/parameters/seed are generator-private.

Four fixed inherited methods: CIE single, Dual CIE, Global Dual, original real-G19 material readout (`prepare_full/full_material_readout`, scalar current). No BothBlocks/absorptive algorithm is added. Solver config copies frozen48 prior/init/FFT/APCS/source weights/beta schedules,3000aux+max3000 corrector,480whole-childCPU and60/120/240/480 paid prefixes. Case-local method order rotates by case index; all32 outcomes retained. No new policy/threshold, no GT online.

## Artifacts and executable entry

- `GENERATOR_CONFIG.json`: explicit8-case generation/physics/independence/gating contract; statusDRAFT blocks execution.
- `SOLVER_CONFIG.json`: exact48 solver settings with only cohort identity, four-method list and workers4 changed; statusDRAFT blocks execution.
- `generate_crossgrid8.py`: root-gated executable, refuses existing public/private/receipt/TASKS paths and preserves failed partial generation. It generates data plus portable forward-slash32-task `TASKS.json`; no inverse or coarse-projected-truth forward is in this script.
- Public future output: `research/research_max_struct_r11/data_public/crossgrid8/{operators.npz,case.npz,manifest.json}`.
- Private future output: `sealed_truth/crossgrid8/{case.npz,FINE_OPERATOR_GEOMETRY.npz,generator_metadata.json,FROZEN_GENERATOR_CONFIG.json}`; each case stores projected `m`, `m_fine` and fine `Y_clean`.
- Future inverse output: `RESULTS/crossgrid8/`; future queue receipts must be under `delegated/twod/crossgrid8/queue` to satisfy unchanged existing queue's ownership guard. Protocol/generation receipt stays in this directory.

After root review changes generation status toFROZEN_FOR_GENERATION and solver toFROZEN_FOR_INVERSION, use the existing numerical interpreter (do not resolve its symlink):

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r11/delegated/diagnostics/crossgrid8_protocol/generate_crossgrid8.py --execute
```

On the existing i9 Windows tree, the same file/config uses its existing interpreter/vendor environment; script `--cpu` and `--vendor-path` mirror inherited generation options. No environment installation is requested. Inversion is a separate later action, only after primary48 queues naturally finish and root releases the same four i9 P cores (logical0/2/8/10), using unchanged `experiments/twod_queue.py --manifest <this>/TASKS.json --config <this>/SOLVER_CONFIG.json --receipt-dir <R11>/delegated/twod/crossgrid8/queue --cpus 0 2 8 10 --execute`. Do not occupy GPU,3D cores, or current48 slots. No queue is started by this plan or generator.

## Costs, seals and later scoring

Generation planned work is8 physical LU at order2304,96 complex source RHS,0adjoint, plus one fine Green setup. Dense single complex2304² array is84,934,656bytes; resident G/G*/identity/forward factor and setup temporaries make a several-hundred-MB generation process, so run only one generator. Setup/physics/export/failed attempts require an external complete child CPU receipt; per-object/self timers are nested. No seconds estimate is asserted before measurement.

The32 inverse allowances total at most15,360nominal CPU seconds before charged atomic/output overhead; this is a cap sum, not predicted actual spend/wall. Four-lane ideal cap-only division is3,840seconds but is not a throughput promise; early stops and shared-host overhead remain actual measured values. Original24-grid physics counts and small19-LU auxiliary counts remain separate. Do not merge this cohort with48 or3D as rawCPU superiority evidence.

Only after all32 terminal public outcomes and all available eligible predictions are sealed may the separate offline scorer read these8 private materials. Score against the2×2 projected target, with inherited material NRMSE/strict support>1/box gate, and disclose fractional projected boundaries. Coarse and fine phase masks are not identical sets. Keep this8-object table separate; no32-run or noise pseudo-replicate sample count.

After the32seal, and **only then**, execute exactly one576-order forward at projected truth for each8 case:8LU/96sourceRHS,0adjoint. Compare Y_coarse(projected truth) with saved fineY_clean to report relative model-error norm and model-error/sigma. Also show raw noise size and candidate residuals, allowing CN99 failure. That operation is offline model-mismatch characterization, never an inverse input, candidate choice, noise inflation, retuning or extra inverse solve. Its concrete sealer/scorer is the next implementation stage after root review; none is silently called by generation.

Current cost: no new physical work and no scientific data generated. This draft used source/public-header inspection and ordinary code writing only; planning CPU was not externally measured and remains unknown.
