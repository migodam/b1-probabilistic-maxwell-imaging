# R13 physical inputs integration

Online modules: `src/vector_adapter.py`, `src/physics_inputs.py`. They have no target generation, GT file read or sealed-truth path. Dependencies are the existing R1/R8/R12 modules, NumPy/SciPy; no new physics library, GPU runtime or service.

## Vector observations

Portable common-runner input files: `delegated/physics/new3d/runner_public/*.npz`. Each contains acquisition G/E0/S, 3D grid metadata, Y, sigma and m_init. All have components=3,216 real material unknowns,648 complex current rows,60 sources,180 electric receiver channels, bounds[0,8].

```python
from research.research_max_struct_r13.src.vector_adapter import VectorPublicModel,run_vector_auxiliary
model=VectorPublicModel(public)
D=model.data_objective(state['Y']-Y,Y) # use this in the record ledger
R=model.reg.evaluate(m)[0]
# Otherwise generic model.state/data_value_gradient and 3-component projector work.
result=run_vector_auxiliary(model,Y,m0,method,check,audit_callback,audit_every=50)
```

Vector `reg.smooth/tv/epsilon` are None, since this prior has different parameters. `reg.name` and `model.metadata` describe the actual inherited prior. The vector objective has NO factor1/2 and is source-normalized; generic scalar Frobenius records must not overwrite it. Auxiliary supports `cie`, `dual_cie`, `global_dual`; the last is a declared global-PR/Dual-beta adaptation. Scalar G19/A+ are explicitly rejected as vector substitutes. Caller owns common corrector, global budgets, pre/post audit eligibility, and retention of partial-run paid incumbents. No auxiliary corrector runs inside this adapter.

## Fresnel measured observations

Portable training input files: `delegated/physics/fresnel_public/{dielTM_dec8f,twodielTM_8f}.npz`. `kind='fresnel'` selects `SourceSpecificPublicModel(public)`. Use `model.projector_S`, a list of27 matrices, in the A projector; common `model.S` is only a compatibility storage field. The authoritative source maps are S_ports with shape(27,39,576).

```python
from research.research_max_struct_r13.src.physics_inputs import SourceSpecificPublicModel,run_fresnel_auxiliary
model=SourceSpecificPublicModel(public)
projector=DenseFieldCurrentProjectors(model.G,model.E,model.projector_S,Y,float(public['delta']),
    components=1,lower=0,upper=20,scale=6,grouped_ports=True)
result=run_fresnel_auxiliary(model,Y,m0,method,check,audit_callback)
```

Mature measured CIE/Dual/Global use original R8 source-specific SVD DPCS/native projected Fourier APCS and matched receive adjoints. They have source-normalized auxiliary metrics but original physical candidate score is the shared `.5||residual||²/||Y||² + prior`. Supported strings are `cie`, `dual_cie`, `global_dual`; no measured G19 adapter is claimed. AP/DR/direct are fully compatible with source-specific model and projector ports.

Actual noise covariance/sigma is unknown. The existing .15||Ytrain|| ball is a declared model discrepancy diagnostic. Required `sigma` in these portable files is ONLY `delta/sqrt(Y.size)` to satisfy scalar storage, with `sigma_role` and `noise_sigma_known=False`; never use it as iid-CN noise, a 99% event, or certified DPCS noise rank. `delta` overrides any Gaussian-derived radius. Receipts and figures must preserve this distinction.

711 heldout measured channels are separate `*_HELDOUT_OFFLINE_ONLY.npz`; never upload/read these in online workers. Those files also contain full empty-calibrated E0 and original S for postseal physical heldout prediction. Heldout masks correspond to all actual source/receiver channels; do not fill absent channels with zeros. Metadata records original local raw-file provenance without copying ASCII files.

## Reproduction commands

Use an existing Python environment with NumPy/SciPy (local default `python` lacks NumPy). The bounded local runs used `uv run --with numpy --with scipy --with threadpoolctl python <entrypoint>`. Entrypoints are `generate_inputs.py` (new directory required), `package_public.py`, `package_fresnel.py`, `qualify_adapter.py`, `qualify_volume_sphere.py`, `diagnose_mie_refinement.py`, `check_external_inputs.py`, `check_integration.py`. Output-generating entrypoints preserve existing data by refusing overwrite. No remote job or push was performed by this worker.

## Classic Austria inputs (second handoff)

Use `classic_public/austria_150to32.npz` and `classic_public/austria_matched32.npz`. Both retain actual existing public Y byte-for-byte in accompanying ORIGINAL_PUBLIC_OBSERVATIONS archives, and use m_init=exactly0 on all1024 unknowns, box[0,20], the inherited prior, fixed .12||Y|| discrepancy radius, and kind='classic'. Separate `classic_sealed_truth/` contains opaque copies of existing truth archives, never loaded for this preparation. Do not upload these archives to an online worker.

`physics_inputs.classic_auxiliary_config(public,method,model,Y,out=None)` returns the actual inherited Austria contract. Do not call R6/R8 generic arm_config for this domain: it explicitly enforces host.1/box[.02,20]/mf3,6,10 and will reject or change the classic contract. Do not infer a new DPCS rank from noise: inherited classic rank is15. Restricted stages are mf5/8/10, material beta8/2/.5; CIE current beta follows material while Dual/global current beta remains8. Original150→32 allowance is100/stage; matched32 diagnostic is300/stage, per actual R12 contracts. Caller may declare an additional whole-process budget but must retain these perstage ceilings/stops rather than silently shortening a mature comparator.

The returned `physics_family='equal_area_disk'` requires the compatible R8 paper_cie_fft loop/EqualAreaGreen in the private inherited-loop clone. Generic DiskPulseToeplitzGreen rejects the equal-area quadrature radius h/sqrt(pi) because it interprets pulse disks as nonoverlapping. This is a constructor-assumption mismatch; reducing radius would change the physical operator and is not allowed. `paper_physics.EqualAreaGreen` and `PaperOneHopAbsorption` implement the exact inherited classic quadrature. Full native Fourier masks on32² are16; restricted stages retain5/8/10. Parent and the full-current worker own that explicit private-clone seam.

Classic `sigma` is actual per-complex-channel RMS of the historical fixed-norm noise, with explicit `sigma_role`; it is not iid CN standard deviation. For150→32 the scalar norm is recovered from publicY and the known original generator seed/direction and fixed-norm equation, without computing/loading cleanY or GT. For matched32 it is copied from the existing generation receipt's recorded noise-to-observed ratio. The fixed ball and fixed DPCS rank do not depend on any Gaussian 99% event.

## Measured periodic audit and interruption receipt amendment

`run_fresnel_auxiliary(...,audit_every=50,*,cost_callback=None)` now honors periodic paid audits after every50 completed material updates and preserves all original initial/stage-end audits. A private `FunctionType` copy wraps the existing `_baseline_readout`, SVD and Counter globals; it does not change the material update, coefficients, beta schedule or PR memory. No shared module globals are mutated. Actual source-specific port dimensions are checked before the legacy loop.

Every exit, including caller BudgetEnd (a BaseException), leaves `model.r13_auxiliary_cost` with captured inherited operation counters, material-update count, SVD setup count, paid physical-audit counts, physical-model counter increments and elapsedCPU/wall. Optional `cost_callback(receipt)` receives the same finally receipt; caller should retain it after a budget exception. Completed calls additionally return `result['r13_audit_receipt']`. Legacy PR_CG_steps increments after readout; an interrupt during periodic audit can leave one paid material update recorded only by `material_updates_completed`, explicitly stated in receipt.

`FRESNEL_PERIODIC_AUDIT_CHECK.json` records a single actual measured public-only check:50 unchanged Dual material updates, periodic audit50, intentional immediate interruption, finally counters retained and callback delivered. CPU.403842s,2 physical audits (inherited initial plus periodic), no heldout/GT reads or full solver comparison.

## C2 proposal audit cadence (supersedes iteration cadence when supplied)

Both `run_vector_auxiliary(...,*,audit_due=None,...)` and `run_fresnel_auxiliary(...,*,audit_due=None,cost_callback=None)` accept the parent's common absolute-process-CPU clock. When supplied, `audit_due()` replaces the periodic iteration-modulo condition; it is evaluated at each completed material-update boundary. Native initial/stage/final audits remain, and the returned/finally receipt declares the effective periodic schedule. Corrector objective evaluations remain audited whenever paid.

`physics_v1` fixes `physical_audit_CPU_interval=2.0` seconds for proposals across all physical families. The caller owns the clock and budget; expensive in-progress material updates may cross a due time and are charged in full before the next legal audit boundary. Fallback audit_every50 remains only for callers that do not supply audit_due. No existing development runtime was modified or rerun. This API/manifest addition was syntax checked only; no additional numerical experiment was performed.
