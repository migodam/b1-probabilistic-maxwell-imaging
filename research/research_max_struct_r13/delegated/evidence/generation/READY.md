# Fresh48 generation ready

Status: **48 public input files ready**, generation only; no inversion or reconstruction quality scoring. Parent must freeze production methods after development before dispatch.

Public manifest: `research/research_max_struct_r13/data_public/fresh48/MANIFEST.json`.

Each entry contains an `input_file` path relative to the R13 root, directly compatible with `src/run_experiment.py`. Example: `data_public/fresh48/r13_fresh_0000.npz`. Case IDs contain only opaque indices. Public NPZ contains G/E0/S, actual cell/source/receiver coordinates, physical radius/k, noisy Y, sigma/noise-only delta, full576 uniform m_init=.1, box[.02,20] and grid/acquisition information. It contains no family, material seed, truth, support mask, clean response or realized noise array.

The sealed generator manifest, geometry records, material averages, clean responses, noise and model-mismatch diagnostics are exclusively under `sealed_truth/fresh48/`. Do not stage that directory to online solvers. No future solver NRMSE has been read or computed.

## Frozen design

Eight material families × six distinct geometry seeds: multi-island rotated ellipses, continuous random fields, random-level patchwork, perforated slabs with closed holes, open/notched porous bodies, thin-bridge networks, eccentric multilayers, and mixed-scale disconnected inclusions. Family/variant mapping is privately permuted across opaque IDs. Material seeds are131300+opaque_index; independent noise seeds are231300+opaque_index. No object is repeated with another noise draw.

Six acquisition/noise groups contain eight objects each:12 sources/full48receivers at1% or2%;6 sources/full48receivers at1% or2%;12 sources/half24receivers at1% or2%. Six sources select alternating members of the actual twelve-source circle. Half aperture selects the upper half of the actual48receiver circle. CN noise has real/imag variance sigma²/2 with sigma=relative clean-observation RMS; draws are not renormalized to a fixed Frobenius length.

All inverse models have576 independent real contrasts on the inherited0.96λ24² disk-pulse grid, k=2π, host/initial material.1 and allowed box[.02,20]. Geometry creates high contrasts in12–18 and low components in2–5; no GT support is an inverse constraint.

Twelve predeclared objects use36² forward observations with24² inversion. One216² sampling of each continuous geometry is averaged to24²/36² square-cell contrasts. The sealed coarse truth is explicitly the cell average, while the physical kernel retains the inherited disk radius. A separate coarse forward at the averaged material measures cross-grid mismatch. It is not added to delta:delta remains a **noise-only** radius and is not a certificate that coarse truth fits cross-grid observations.

## Actual generation cost and verification

`GENERATION_RECEIPT.json`:60 full LU /600 complex RHS; two operator setups; measured generator interval4.598228 CPU seconds /4.933500 wall seconds. The interval includes package imports, generation and per-case serialization up to writing the final receipt; final receipt write/process exit are outside this interval. Memory peak was not measured.

`CASE_GENERATION_COSTS.json` preserves every case's counters and interval; sealed diagnostics preserve physical solver residual and coarse/fine mismatch. `PUBLIC_INPUT_VALIDATION.json` checks all48 actual public files for finite arrays, dimensions, uniform initialization, full unknown count and absence of truth fields. Total public NPZ size34,127,079 bytes. No extra physics was called for validation.

Two failed launches preceded generation and performed zero physical solves:one syntax typo, corrected;one default Python without NumPy. Their exact logs remain `generator_stdout.log` and `generator_stdout_attempt2.log`; their process CPU times were not measured and are unknown. The successful run used the existing `Gaussian/.venv_nn/bin/python`, one BLAS thread, and installed nothing.

Existing Austria has a separate R13 public-input adapter in `src/physics_inputs.py`; this fresh-generator task did not regenerate or rename that historical benchmark as new data or literature reproduction.
