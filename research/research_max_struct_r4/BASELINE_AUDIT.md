# Baseline audit (delivery draft; no scientific judgment)

Evidence inspected: actual experiment entry points/configs, immediate development run `RECEIPT.json` files, shared port setup receipts, existing cost audit, and public-only confirmation adapter/package code. No confirmation GT, private NPZ, new solver execution, or third-party full-text scan was used.

| Baseline/control | Actual scope | Cost/comparability boundary |
|---|---|---|
| Historical R2 global | Saved mature `research_max_struct_r2/src/cie_global.py` auxiliary schedule and historical result | Historical 689 LU / 16452 RHS and 9.22095 s are reference history, not new R4 work. Historical unscaled corrector differs. |
| `global_common_corrector` | Reuses saved global auxiliary candidate; new common scaled physical corrector | Charge new receipt 1013 physical LU / 24288 model RHS. Do not call it a freshly timed complete auxiliary+corrector run. Historical auxiliary 3 LU / 36 RHS is separately identified, not added again. |
| `native_extended` | Existing native current solver; 8000 updates and separate common corrector | Auxiliary `model` delta of 5 LU / 60 RHS is already inside top physical model counters. Not an equal-compute comparison with 1500-update arms. |
| Generic block 2x2 | Mature CIE extension, fixed blocks; 1500 updates | Small 4x4 factorizations/RHS and physical 576-cell LU/RHS are different operations. Embedded one-direction wiring cost is included once. |
| Cached block 2x2 | Arithmetic efficiency replay, same configured scientific schedule | Entire replay is chargeable. Inverse multiplication may change floating-point trajectory; no exact-trajectory or free-cache claim. |
| Self-only | Same generic material readout, scalar diagonal application | Scalar division counts are separate; no small LU. Not a specialized exact-rational scalar readout implementation. |
| Onehop / phase pi | Fixed one-exchange operator and its fixed phase control | Both runs charged; phase control reuses existing derivative evidence without another test. Counts include scalar divisions and block matrix applications, not small LU. |
| T-matrix | DCTMC-style real-diagonal box variation, explicitly not faithful paper reproduction | Top counter LU includes 30 nonlocal locality factorizations: split 1112 into 1082 physical plus 30 locality. Split 60744 RHS into 43464 retained physical/model plus 17280 locality RHS. |
| Illuminated localization | Separate known-method current/localization control | Its top counters cover its own completion/readout/corrector work; auxiliary readout counts are nested metadata. |
| Port lambda0 / lambda1 | Known WRI/state-profile control / finite-eta port augmentation candidate | One shared reference charged once (1 LU/60 RHS); independent 4x4 envelope check charged once (5 LU/106 RHS, dimension16, not 576). Arm auxiliary/corrector counters are subsets of top cost. Small48 profile Cholesky/RHS are separate from physical LU. |
| Cooperative / box paths | Full geometry plus finite rays; box reuses paid full geometry | `curvature` free-face diagnostic is an independent paid run. Full geometry is charged in cooperative receipt only; box's own forward/ray work remains chargeable. Schur/eigen/triangular counts are separate operations. |
| Current invariants | Saved-current algebra/diagnostic only | Receipt explicitly records zero new physical LU/RHS; elapsed work is not zero. |

All arms use original physical Phi for candidate ranking within their stated protocol, but that does not make their budgets, starts, arithmetic, or objective histories equal. This file does not determine novelty, convergence, recovery, superiority or confirmation acceptance. R3 base `1bd07e49` is ordinary provenance; hashes were not verified. Monetary costs remain unknown (`null`).
