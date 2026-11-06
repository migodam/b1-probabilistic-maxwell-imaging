# R9 Phase2 strict real-vector execution

Frozen contract: `../../configs/vector_extension_phase2.json`. The implementation used at launch is copied under `implementation_frozen/`; this is provenance, not a new physics backend. Numerical work uses inherited R1 TensorDDA and R8 vector conventions. R8 and R9 PhaseI are preserved.

Entry from repository root, using the existing environment:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r9/experiments/vector_extension.py all
```

This is the original launch command, not a command to rerun into existing directories. `all` is deliberately write-once and originally schedules24 tasks. A root-authorized pre-GT amendment adds five resolved-hollow controls; the actual execution uses an explicit scoring HOLD to stop that already-running original process before its scoring reads. `complete_amendment.py` then runs the unique grid16 qualification and the five-task queue before `finalize29.py` seals/scores all29. See `configs/vector_extension_all29.json` and `AMENDMENT_COORDINATOR.json` for exact task identities and sequence.

For self-contained public-only inversion in a fresh output tree, the existing worker entry is `worker --case CASE --arm ARM`: it reads only the bundled R9 public Y/noise/acquisition files and does not need any R8 untracked array. The resolved-hollow worker is `delegated/3d_extension/resolved_hollow_control.py worker --arm ARM`. Both preserve existing outputs by refusing overwrite. The historical-copy portion of original `all` was generation-only provenance; it is not required by these public workers. `finalize29.py` performs no physics and requires all29 terminal public outputs before any private scoring read.

The user subsequently explicitly requested a remote P-core handoff. Actual main cohorts are **24 completed Mac tasks** and **five complete XINAN i9 resolved-hollow tasks**, each with the same frozen physics/budgets. The first two local resolved-hollow attempts were already running when steering arrived; they were stopped only for this explicit resource move and archived under `local_aborted_resolved_hollow`, with sampled CPU and available LU/RHS lower bounds in `LOCAL_REMOTE_MIGRATION_ABORT.json`. They are neither algorithm failures nor part of the remote paired CPU table. The completed Mac Mie16 control is kept and not repeated.

The portable remote entry is `windows_resolved.py queue`, using only the small public/code bundle, logical P cores8/10, one pre-import BLAS thread, foreground connected SSH and an inherited Windows kill-on-close JobObject. It reuses `StrictVectorInverse` unchanged. There is no GPU, truth, `resource` or `threadpoolctl` dependency in that portable wrapper. Whole-process CPU is Windows psutil user+system/GetProcessTimes; peak working set is explicitly measured, not replaced by an unknown zero. See `REMOTE_TRANSPORT.json`. All29 outputs are imported and sealed before scoring; CPU comparisons remain within each object/platform.

The legacy four-cube Y is copied unchanged into R9 public files. Fine-grid clean fields, reference materials, generator seeds and target geometry contracts remain under `sealed_truth/vector_extension`. Public observations contain exactly one array, Y. No worker loads private truth, clean fields, target support or a private noise realization.

## Strict split and real-isotropic A+

All current steps use ordinary CIE. A/A+ replace only the material residual. In a local H block, the scalar-isotropic self term is s I3 and Hii=0. For one real material coordinate x, let b=beta-s and z=x/(1+b x). With J,Q fixed, u=beta J+z[-beta(Q+bJ)]. The full block residual is affine in this single complex scalar z, including all three components and sources: R=C+z D. Therefore its state loss is q|z|²+2 Re(h z)+constant, with q>=0. Its real-coordinate derivative has a quadratic numerator. Legal stationary roots, box endpoints and the current x are finite candidates. The original quadratic spatial/material prior is evaluated for each candidate; the whole unchanged objective is checked before accepting the sweep. Eight colors have disjoint H blocks and no same-color prior edges, permitting equivalent independent conditional proposals within each color.

This applies to one independent **real isotropic** scalar per voxel. It is not a tensor-material or complex-permittivity A+ implementation. The prior is not silently folded into the state-only quadratic roots; endpoints/current and unchanged-prior acceptance remain explicit.

The source-wise PR and global PR variants share the same finite data, material box, prior, Fourier spaces and physical corrector. Native projected-IDFT coordinates are transported by frequency/component identity across MF1→2; there is no APCS whitening. The paper Eq16 diagnostic includes its square root. The fixed 240-update budget remains primary and is not reported as paper convergence. Offline χ relative-L2 NRMSE is not paper Eq17.

## Scope and electrical size

All inversions have 8³ independent real unknowns and fine10³ generation. The .6λ legacy four-cube side is .16λ (2.13 inverse cells), whereas the original paper side is .8λ. At ε5 these are respectively .358 and 1.789 internal wavelengths. The ε9 target cubes have side .12λ = .36 internal wavelengths, gap .10λ = two inverse cells, and inverse h=.05λ = .15 internal wavelengths (6.67 cells/internal wavelength); the fine grid gives 8.33 cells/internal wavelength. Higher contrast alone does not make this smaller target harder than the paper target.

The .6λ hollow cube has outer side .30λ, inner side .18λ, shell thickness .06λ = .8 inverse cell. Each rod has .07λ transverse width = .93 inverse cell. These thin boundaries are explicitly underresolved. Cube/hollow cells use linear contrast times exact geometric box-volume fraction inherited from R8. The independent sphere validation instead uses MG polarizability mixing from fixed 24³ midpoint fractions. Its Mie error is a backend/geometry control, not a bound on the box-object fields.

The new cases use raw iid circular complex 3% expected Frobenius noise and publish scalar sigma with E|n|²=sigma². Qualified DPCS uses the smallest supported contiguous SVD prefix whose discarded Y energy is at most sigma² chi2(.99,2*nr*ns)/2. The legacy exact-energy 10% observation uses the deterministic bound (.1/.9)||Y||F instead. Both concern noise only; neither absorbs the independently measured fine/coarse model discrepancy. No per-rank significance or GT selection is used.

## Evidence and accounting

- `phase2_all.log`: preparation plus every child start/end/error.
- `phase2_queue/RECEIPT.json`: finite parent queue and per-job wall/exit status.
- `../../RESULTS/vector_extension/validation/RECEIPT.json`: fixed ε9 sphere 8→10 Mie result, including failures and setup/solve costs.
- `../../RESULTS/vector_extension/generation/`: independently assembled fine/coarse mismatch and generation costs; no repeated legacy forward cost.
- Per-arm directories: inherited-physics provenance, whole-worker CPU/wall, setup/online nested costs, LU/all-source/adjoint RHS counts, solver stop, paid original-Phi incumbent trajectory, terminal fields and complete stage joint m/J/Q/gamma.
- Accepted-step80 and objective-evaluation80 checkpoints copy already paid incumbents; terminal finalization adds no LU/RHS.
- `PUBLIC_PREDICTIONS_SEAL.json` is written before offline private scoring for all29. The archived `SCORE_HOLD_ROOT_AMENDMENT_29.json` is explicitly **not** a prediction seal; it is the pre-GT flow-control record. Evidence stamps use size/mtime, no SHA.
- `PHASE2_COMPLETE.json` adds parent CPU to all child CPU once. Parent includes imports, qualification, generation, queue, scoring and figures; child costs include their imports/setup/online/export. Nested setup/online costs must not be added a second time. Direct300 and joint240+corrector300 are different budgets; compare measured paid-cost curves as well as terminal values.

The final cost record explicitly retains lower-bound/unknown overhead where a process was stopped before its normal receipt: it does not turn those costs into zero or claim a falsely exact total. Original24-parent validation/generation receipts and the local control-generation receipt are nested in their observed parent CPU and are not double-counted. Cross-hardware sums are consumption bookkeeping only.

No continuous-physics accuracy, original-paper physical reproduction, optimizer convergence or material-recovery claim follows merely from execution completion or a small linear-system residual.
