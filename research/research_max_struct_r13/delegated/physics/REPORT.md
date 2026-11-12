# R13 physical adapter and new input delivery

Completed new callable public-only vector/measured adapters, four actual new3D observations, two portable actual measured Fresnel training packages, bounded physical qualification/diagnosis and public constructor checks. This worker did not run the remote inversion comparison, score inversion material recovery, change legacy results or publish. Parent owns final scientific acceptance and actual comparison execution.

## Exact physics and scientific scope

Reused R1 `src/vector_backend.py::TensorDDA`, R8 `src/vector_cie3d.py::VectorCIE3D`, and R8 `src/vector_pr3d.py::VectorPR3D`. G is the complete dyadic CM+RR DDA operator, including self-density `-1/3+i k³h³/(6pi)`, with k=2pi and lambda0=Z0=1. All three vector components couple through G. The real isotropic material m has N coordinates and repeats each coordinate on its three consecutive component rows. No known support is given to inversion.

States solve `(I-D_m G)J=D_m E`, with `Q=E+GJ`, `Yhat=SJ`. Here J is unscaled susceptibility/polarization current j=chi E, matching existing field/current projectors; dimensional electric current in R1 convention is `-ik j`. One LU serves all60 source RHS and its analytic adjoint. The exact physical gradient sums every source and all three components for each cell.

Vector physical objective is `mean_s ||Yhat_s-Y_s||²/||Y_s||² + 2e-4 mean_edges(m_i-m_j)² + 1e-6 mean_cells m²`, matching R8 exactly. The physical value/gradient has zero numerical difference from original R8 on the checked state. A new directional-derivative check differs by1.043e-8 relative. This checks adapter wiring; it does not certify an inverse solution.

The original R12 highcontrast6³ cases had h=.1 and max internal cell phase1.885 at chi8. New bounded cases use width=.24, h=.04 on6³ and h=.03 on8³; phase decreases to.754/.565. This is a more resolved finite dyadic-DDA study, not a universal continuum qualification. Qualifications below retain meaningful negative evidence.

## Newly generated actual observations

Frozen before observation generation: `new3d/FROZEN_PROTOCOL.json`. Three independent analytic rectilinear composites were fixed without searching outcomes: a slotted high-chi frame, connected unequal inclusions with chi8/3/5, and a mixed-chi crossed bracket with a cavity. The fourth observation uses the SAME first geometry on8³; it is a crossgrid pair, not an independent geometry. Piecewise box volumes are analytically integrated by the existing R12 disjoint-brick volume-average routine.

Each public observation has180x60 complex channels and new iid circular-complex1% noise with declared sigma and fixed seed. No clean prediction, target coefficients, geometry parts or labels enter public observation/runner files. Inverse chi and generation chi, clean observations and noise are under separate `new3d/sealed_truth/`, generator-only until all inversion predictions and choices are sealed. The online adapters never reference this directory.

Generation paid4 full Maxwell LU/240 forward source RHS. All new physical state residuals are near1e-15. The independent-grid frame changes clean data by **8.0063%** relative, compared with1% observation noise. Therefore a strict noise-only radius is not a true-material feasibility certificate for crossgrid inversion. Preserve this mismatch and any strict-noise failure; any larger diagnostic discrepancy radius requires explicit labeling and parent decision. No radius was changed here.

Generation receipt: `new3d/GENERATION_RECEIPT.json`, CPU.527766s, wall.547176s excluding process import before timer, with inverse/fine setup metadata included. Portable operator-enriched files are in `new3d/runner_public/`, generated only from public observation files and acquisition operators.

## Independent analytic Mie diagnosis

All controls use the exact same TensorDDA Green/self normalization and independent analytic Riccati-Bessel Mie amplitudes. Linear solve residuals are near1e-15; those residuals do not establish continuum accuracy. Radius is fixed at.09lambda, eps9 for highcontrast refinement. Both negative coarse and later finer controls are retained.

| Boundary rule | h | Far vector relative Mie error |
|---|---:|---:|
| Staircase sphere | .04 | .33731 |
| Staircase sphere | .03 | .18672 |
| Volume-averaged chi | .04 | .36627 |
| Volume-averaged chi | .03 | .27430 |
| Volume-averaged chi | .02 | .19761 |
| Volume-averaged chi | .015 | .16934 |
| Staircase sphere | .02 | .09015 |
| Staircase sphere | .015 | .07064 |

Coarse staircase sphere volume differs from analytic sphere by19.46% and8.76%, explaining a substantial geometry contribution. Volume-averaged chi uses48x48 Gauss x/y integration and exact z intersections; at h=.02/.015 its integrated volume agrees within about1e-5 relative, but highcontrast far error remains19.8%/16.9%. On those same finer grids the homogeneous staircase error is9.0%/7.1%. Thus discretization and boundary material representation remain material biases even after volume matching. The observed convergence and same-backend boundary comparison do not suggest a simple normalization bug; parent owns final interpretation. A diagnostic fitted complex gain was recorded only to expose amplitude/shape error and was never used to alter a prediction.

At lower eps3, coarse volume-average errors are11.16%/9.35%. This qualification cannot support a chi8 continuum accuracy claim for the6³ inverse cube. It does support real dyadic3D implementation and bounded discretization evidence with explicit accuracy limits. New nonsphere inverse runs remain meaningful as finite-DDA comparisons; lower residual must not be presented as continuum/material recovery.

Receipts: `QUALIFICATION.json`, `VOLUME_MIE_QUALIFICATION.json`, `MIE_REFINEMENT_VOLUME_ONLY_FIRST_ATTEMPT.json`, `MIE_REFINEMENT_DIAGNOSIS.json`. Finer volume-only first attempt paid3.667062CPU s; adding staircase controls repeated those two volume solves and paid4.749531CPU s total for the expanded diagnostic. Both costs are retained, including duplicate work.

## Actual Fresnel calibrated inputs

Read actual existing official2001 local packages in sibling R8 worktree: `data_local/fresnel/dielTM_dec8f.exp` and `twodielTM_8f.exp`, using existing R8 `fresnel_adapter.prepare`. No simulated replacement, gain fit to scattered data or invented antenna calibration was used. Frequency2GHz,24²/0.16m square,27 training sources with39 actual receiver rows each =1053 training channels,711heldout actual channels. Order2 outgoing source multipoles are fit only to empty-target measurements; same original heldout-source/channel split is preserved.

Empty-field calibration validation medians are.0251476 single/.0240807 twin. Receive-horn response and actual noise covariance remain unknown; source-field continuation is not a certified receive-antenna model. Source-dependent receiver maps and matched adjoints are preserved. The portable model reproduces original source-specific predictions, objective and gradients exactly (all differences0).

The fixed .15||Ytrain|| discrepancy balls are actually affine-reachable. Numerical-rank affine floors are.0981254 absolute for single (1.865% relative) and.236506 for twin (3.252% relative); zero balls correctly fail. Projected residuals match fixed delta.789134/1.090861. These are receiver/propagation graph reachability statements, not real-material feasibility or a noise certificate. See `EXTERNAL_INPUT_RECEIPT.json`.

Portable training packages: `fresnel_public/*.npz` with kind='fresnel', S_ports shape27x39x576 and Y shape39x27. Heldout measured values are kept in explicitly offline-only files. `sigma` is a required-storage discrepancy scale only, with noise_sigma_known=False and explicit sigma_role; it must never be treated as known CN sigma or a 99% statistical event. No pixel material truth exists. Permitted measured outcome is predictive fit on fixed heldout data after predictions are sealed.

Existing R8 source-dependent `fresnel_solver.run` supplies mature sourcewise CIE/Dual and global-PR controls, with each source's own SVD DPCS/native projected Fourier APCS. The new seam preserves that implementation and adds caller-owned physical stage audits/budget checks through a local proxy. Scalar G19 has no demonstrated source-dependent-port adapter here and is explicitly rejected; do not silently run the scalar common-port method under a measured G19 label. Direct/AP/DR and actual mature measured CIE/Dual are callable.

## Classic Austria reuse

Loaded actual R8 `RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz`, no GT read:20 sources,40 receivers,32² inverse,400MHz. Historical forward grid is150², real chi=epsilon-1=3,10% fixed-Frobenius complex noise. Reused equal-area disk-pulse VIE and prior. R12's actual offline coarse-target physical residual was.225256 against a.12 diagnostic ball; this existing negative fact is retained rather than calling that target noise-feasible. `prepare_austria` exposes actual public observation/model for parent use. No new Austria run was requested from this worker.

## Implementation/reproduction boundaries

Read `INTEGRATION.md` for exact APIs and source-specific ports, objective ledger normalization, metadata and bounded local commands. `INTEGRATION_CHECK.json` records a two-cycle native-vector Dual interface exercise only; no scientific comparison is claimed from it. Local package construction/verification, imports and qualification costs are distinct from parent remote inversion budgets.

Early attempts: default Python lacked NumPy; first uv generation and qualification imports failed on an old R12 generator's sibling import. The entrypoint now supplies that legacy import path. These failures produced no forward solves or datasets. No driver, CUDA, permanent environment service or legacy file was changed. The one-shot generation/refinement artifacts, including failures and duplicate diagnostic cost, remain reviewable.

## Second handoff: actual classic Austria inputs and runner audit

Prepared `classic_public/austria_150to32.npz` and `austria_matched32.npz` from the actual existing public observation archives. Original public files/source receipts are copied alongside; metadata records exact provenance, physical acquisition, all-zero air-host initial material, fixed[0,20] box, fixed.12||Y|| ball, inherited prior and source-stage settings. Existing truth archives were copied opaquely to `classic_sealed_truth/`, with no material arrays opened or new truth/clean observation generated. This handoff paid0newforward/0LU/0RHS/0validation; acquisition assembly/file-copy costs are recorded separately.

`classic_auxiliary_config(public,method,model,Y,out=None)` preserves DPCS15, restricted mf5/8/10, materialbeta8/2/.5, CIEcurrentbeta8/2/.5 versusDual/globalcurrentbeta8, original100/stage and matched300/stage. The originalNIE2016 physical source contract and adaptation differences are explicit in metadata and `QUALIFICATION_DIFFERENCES.md`. Actual primitive operators are equal-area disk quadrature with radius h/sqrt(pi); generic nonoverlap-disk constructors cannot accept that radius and must not shrink it. The full-current worker is implementing the explicit equal-area private-clone seam.

Runner review found and addressed standalone Fresnel import reexports in vector_adapter; physics_inputs.py must also be bundled. Vector objective decomposition/prior attributes and source-specific projector ports are currently wired correctly. Parent must branch kindclassic to the classic factory instead of generic fixed-host/noise-rank/mf12 configuration. Mature paperloop supports perstage stop thresholds and zero initial APIC, while original G19loop supports physical-current warm initialization and global stop settings: preserve/report each declared contract and obtain the parent's scientific decision before claiming identical stopping. No scalarG19 method is silently relabeled as vector or measuredsource-specific G19.

## Final adapter audit amendment before source freeze

Reviewer identified that the first measured seam ignored audit_every and would lose auxiliary counters on a budget exception. Corrected by a per-call private copy of the unchanged R8 Fresnel function: wrap its original material readout for each50th paid candidate, retain stage-end validations, capture Counter/SVD setup counts, and always publish model.r13_auxiliary_cost plus optional cost_callback from finally. No beta, readout, PR coefficient/memory or physical-model change. Wrong/missing source-port shapes reject before executing the loop.

One new bounded actual measured public-only50-update wiring check intentionally interrupted at periodic audit50 and retained all counters/callbacks: CPU.403842s, inherited receiverSVD27, completedmaterialupdates50, inheritedPRcounter49 (last increment occurs after returned readout), physicalaudits2. The pending last paid update is explicitly counted separately rather than pretending the old PR counter includes it. No heldout or pixeltruth read. Receipt `FRESNEL_PERIODIC_AUDIT_CHECK.json` and integration contract amendment are retained before parent freezes the final shared source.

The finite `configs/physics_v1.json` manifest is READY with60 fixed tasks, one4CPUqueue, noGPU, shared proposal.6/steps3000/audit50/AP128. Vector/Fresnel caps180CPU; classic240CPU. Truth/heldout specs contain only path/key metadata; all public paths exist, no truth arrays were read for the manifest, and this worker did not execute the campaign. Parent must freeze final reviewed sources after these audit corrections before launch.

## C2 shared CPU-cadence addition before any physics campaign run

Added optional audit_due callable to vector and measured auxiliary APIs. When present it replaces periodic iteration50 scheduling with the parent's common absolute-processCPU2-second proposal clock, checked at completed material-update boundaries. Native initial/stage/final audit extras remain; all completed work and audits are paid. Corrector objective evaluations keep existing paid auditing. Both vector receipts and Fresnel finally receipts state the effective schedule; interruption cost callbacks are preserved. The60-task manifest/ready record now fixes physical_audit_CPU_interval=2.0 and records scope. This change received syntax/basic API review only, without new numerical execution or changes to an existing development runtime. Parent will freeze C2 before launching any physics_v1 task.
