# Final data and anchor contract

The final generation is `data/GENERATION_RECEIPT.json`. The first generation
used every even test index for stress, which inadvertently selected only two
families. Before inspecting any reconstruction truth, the design was changed to
`test index % 8 < 4`, balanced across all four families and ID/OOD. Both measured
generation costs remain in `GENERATION_RECEIPT_INITIAL_STRESS_DESIGN.json` and
`GENERATION_RECEIPT_BALANCED_STRESS_DESIGN.json`; their combined costs are 848
full VIE LU and 10,176 complex source RHS, with 9.292 s generation wall. No GT
scoring was used to choose the amendment. DEV records and inputs are unchanged.

## Parent distribution and independent generation

| Split | Parents | Digits | Primitive composites | Islands | Smooth |
|---|---:|---:|---:|---:|---:|
| Train | 256 | 64 | 64 | 64 | 64 |
| Val | 64 | 16 | 16 | 16 | 16 |
| Test ID | 32 | 8 | 8 | 8 | 8 |
| Test OOD | 32 | 8 | 8 | 8 | 8 |
| DEV timing only | 8 | 2 | 2 | 2 | 2 |

Each parent has independent SeedSequence `[20261107, split_code, parent_index]`.
All variants of a parent retain its split. There is no parent-derived train/test
overlap. This is procedural independence rather than a claim that random
geometries can never be visually similar. The seed, digit/font, shape generator,
background, alpha and amplitudes are evaluator-only metadata.

Each parent gets one low amplitude, cyclic .1/.2/.4, with an independently drawn
constant background uniform on [.08,.15]. Counts per split differ by at most one
between the three amplitudes. Material is `m=background+amplitude*alpha` on all
576 cells; no background is supplied to inference or fixed by the optimizer.
Alpha lies in [0,1]. Smooth alpha is privately normalized by its positive peak
after subtracting its 40th percentile. Rendered alpha uses pixel-area averaging
from 192² geometry. Neither normalization uses a reconstruction or GT score.

Train/val/ID digit fonts are Arial Bold and Georgia Bold, with fresh scale,
anisotropic scale, rotation and position. OOD digits exclusively use Courier New
Bold and Times New Roman Bold. Digits cycle 1–9. Primitive ID scenes use ellipse
and rectangle unions; OOD scenes use random polygons. ID islands use separate
ellipses; OOD islands use arcs. ID smooth scenes use Gaussian random fields;
OOD smooth scenes use random Fourier fields. FAEMIS inspired geometry choices
only; its linear scattering and reconstructed channels are never imported.

Exactly 32 test parents were predesignated for A=2.5 stress: 16 ID and 16 OOD,
four from each family within each distribution. Stress shares the parent alpha,
background and normalized noise direction with its low-amplitude counterpart.
Every low/stress material has paired clean and noise=1e-3 observations. Noise is
complex Gaussian direction rescaled to exactly `1e-3*||Y_clean||`, independently
drawn per parent, with the same direction used for its paired variants.

The final dataset has 848 observations: train512, val128, test192, DEV16.
The final generation performs 424 full nonlinear R4 solves, 424 LU / 5088 source
RHS. Truth `m,j,alpha,amplitude` uses the full nonlinear VIE state, not Born or
FAEMIS. Public observation NPZs contain only `Y,noise`. Operators are separate;
`V48` is one shared complex receiver-right basis of S, computed by one economy
SVD. Family metadata is absent from `manifests/online.json` and anchor APIs.

## Anchors and fixed feature accounting

Reuse the existing R4 `Model`, `Regularizer`, `material_basis` and
`resample_material`; original sources remain unchanged. The base anchor is
bounded native SciPy TRF, full original normalized data objective without
regularization, stages 8→12→24 with .20/.25/.55 shares of the 45 s / 500 LU /
30,000 complex RHS cap, material box [.02,20], public uniform start .1 and
ftol=xtol=gtol=1e-12. The physical solve remains 576-dimensional at every stage.
The separately requested optional weak anchor uses smooth=TV=1e-5 under the
same R4 regularizer; it is saved with `_weakreg1e-5` suffix.

`m0,j0,Q0` belong to the actual last-stage accepted TRF point; rejected residual
trials cannot replace that endpoint. The best valid physically evaluated point
is preserved separately. Accepted endpoint state arrays are retained from their
actual forward operation, so final features require no uncounted forward solve.
Any incomplete fine-stage release is explicitly reported.

The shared reference at public m=.1 is `cache/reference_jacobian.npz`, with its
own group setup receipt (one LU, 12 source + 48 receiver-transpose RHS). Packing
is source-major, receiver-major, real/imaginary interleaved. The online feature
is exactly `gb=J_ref.T pack(Y-F(m0))/||Y||²`; cached J_ref is fixed. Its matvec
wall is recorded independently in `feature_wall_seconds`, and its feature
counts are zero LU/RHS plus one fixed-J transpose matvec. Base solver costs and
reference setup must not be charged again as a new method's 10 s incremental
budget. Input I/O and child process wall are reported separately.

`run_anchors.py` uses one fresh subprocess per observation, BLAS threads=1,
at most two workers locally or four on Windows. It consumes only online.json,
public NPZ and operators, never evaluation.json/truth. Every observation has a
JSON receipt, and the finite batch has aggregated costs and stop histograms.
Failed processes retain stderr and unknown-cost status. No shared service or
GPU is started. Windows `--core-mask` binds only this process's descendants.

Initial DEV timing: eight unrelated clean parents completed in 12.474 s total
with two fresh local workers; per-anchor solve 0.591–4.483 s, all final stage24,
data residual 5.394e-8–6.370e-7. These residuals are optimization diagnostics,
not GT material or scientific acceptance. No DEV GT scoring was performed.

## Predesignated append-only training expansion (implemented, not executed)

From the study directory, the later parent decision can invoke:

```sh
python experiments/prepare_data.py --output data --expand-train 1024
python experiments/prepare_data.py --output data --expand-train 2048
```

Each command increases the training parent count to the requested target using
the same independent `[20261107,1,index]` streams for previously unused indices.
The original first 256 train parents, all other split parents, their public and
truth files, and all existing manifest records remain unchanged. New records
are appended at the end of both manifests, retaining the old record order;
existing files are never opened for rewriting. New scenes remain ID training
geometries with the same cyclic amplitude, background, and paired noise rules.
The command does not read existing material arrays or select by GT/model score.

Each expansion writes its own
`GENERATION_EXPAND_TRAIN_<before>_TO_<target>_RECEIPT.json`, including all actual
generation LU/RHS and wall costs, completed parent/observation counts, and
failure status. The original generation receipt is retained as the initial
generation record. Partial failures preserve written artifacts and report
whether the manifest commit completed; existing collision artifacts are never
silently overwritten. Expansion has not been executed as part of implementation.

The read-only cost utility is:

```sh
python experiments/prepare_data.py --output data --summarize-anchors anchors --summary-output audit/ANCHOR_COST_SNAPSHOT.json
```

It reads online.json and scalar anchor receipts only, grouping base and optional
weak anchors by split, stop/status and backend. It reports solver counts,
feature counts, separate wall sums, incomplete base observations and shared
reference setup. Per-task wall sums are computational work totals; overlapping
batch wall records are retained without pretending their sum is elapsed time.
This command performs no physics and opens no truth/NPZ payloads.
