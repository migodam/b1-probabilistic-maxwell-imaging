# Reproduce frozen Track B inference

Restore the read-only R4 source/assets at the paths in `SELECTION_FROZEN.json:source_assets`; keep this folder under the R5 repository path. Runtime: existing `Gaussian/.venv_nn/bin/python` with NumPy, SciPy, PyTorch, matplotlib. No additional installations/cloud are needed. All paths and data regimes are in ASSET_INVENTORY.json and SELECTION_FROZEN.json. No SHA checks were used.

```python
import sys, torch
sys.path.insert(0, "/absolute/path/to/research/p3l_r5_parallel_ab/stage_b_material_risk")
from inference import Action
torch.set_num_threads(1)
action = Action("distilled_w128_n32768", 2026102601)
prediction = action.predict(y, online, batch_size=64)
```

`y` is real source-major shape [N,512]. `online` contains only matching original R4 anchor arrays: observation[N,32,8], residual[N,32,8], q0[N,66], chi0[N,576], J[N,512,66]. The action never acquires physics or reads labels. Valid finalists are deterministic_matched and distilled_w128_n32768, seeds2026102601/02/03. Original_K32_frozen additionally requires parent identifiers. Default deterministic quadrature is robustQ64; Q16 is only a measured numerical fast path for qualified checkpoint inputs. Replaying oldVAL is a reproducibility check, not a fresh blind evaluation.

The raw timing benchmark, original/MPS predictions, source checkpoints, all fit losses, training split memberships, teacher caches, failures and metadata remain in this folder. Teacher-generation cost and privileged pretrainedGaussian dependency are charged/disclosed. Cached teacher32768 was reused for the capacitypair: it is an offline generation cost, not a runtime provider requirement. It does not remove sunk originalteacher training/physics costs. Compact deployment needs one PROVIDER.pt and one HEAD.pt perseed, not teacher cache or originalGaussian encoder.
