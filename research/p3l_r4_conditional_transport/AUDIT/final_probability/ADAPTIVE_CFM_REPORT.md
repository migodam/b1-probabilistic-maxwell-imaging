# Adaptive CFM numerical-density diagnostic: bounded partial result

Status: PARTIAL_BUDGET_STOP. No rerun or expansion. The completed matched-raw DEV0 tolerance pair **fails** the frozen0.05nat agreement criterion:103.88341198 vs104.01940092, absolute difference0.13598895nat. All other model/object agreement conclusions are unestablished.

| Model/object | Loose NLL | Tight NLL | Loose raw-cycle maxabs | Tight raw-cycle maxabs | Status |
|---|---:|---:|---:|---:|---|
| Matched raw CFM / DEV0 |103.883412 |104.019401 |0.0475475 |0.00208430 |Agreement criterion failed |
| Matched raw CFM / DEV1 |83.402244 |Unavailable |0.0381988 |Unavailable |Tight solve budget-interrupted |
| Anchor CFM_none / DEV0,1 |Unavailable |Unavailable |Unavailable |Unavailable |Not attempted within cap |
| Geometry CFM / DEV0,1 |Unavailable |Unavailable |Unavailable |Unavailable |Not attempted within cap |

Fixed same-model Float64 RK4 NFE32 at true DEV0 input gives NLL104.36160165,0.34220072nat above the tight adaptive estimate. Original public Float32 vs same-rounded-input Float64 fixed32 logq differs only4.01e-6nat; standardized-input rounding adds3.66e-7nat. Those precision terms do not explain the much larger fixed/adaptive and tolerance differences. The adaptive estimate itself is not established as converged.

Full trace is deterministic66-coordinate autograd, CPU Float64, same stored network coefficients and source condition/base values promoted. The completed loose/tight DEV0 backward solves used16070/42146 velocity evaluations and1,060,620/2,781,636 diagonal VJPs. Across completed and interrupted ODE work:156173 velocity calls,5,614,422 trace VJPs. Forward cycle tests used velocity only and the actual frozen raw-coordinate inverse. See ADAPTIVE_CFM_DENSITY_PROTOCOL.md one directory above for equations, logdet signs and provenance.

The single analysis-only Book closed at895.090029s wall /885.654420s CPU, peakRSS264814592B, zero GPU/RHS/optimizer updates/heldout GT reads. Three solve+cycle pairs completed; one tight solve was interrupted by the fixed global budget; eight declared solves were not reached. Final raw records are preserved in RECORDS.jsonl/RESULT.json; CONVERGENCE.csv distinguishes budget interruption from numerical solver failure and unattempted cases. All3 checkpoint and original source byte/mtime identities are unchanged.

Only the two already exposed DEV GT q arrays were read, inside this isolated diagnostic; prepared arrays are explicitly marked numerical-audit-only. No locked or confirmation GT was read. Existing fresh NFE32 sampling/model/normalizer choices and watcher remain unchanged. Earlier fixed16/32/64 failure evidence is retained; its small unmatched raw CFM is a different checkpoint. This test establishes neither density accuracy nor poor distribution quality, exact finite-map/MH density, locked NLL, or posterior qualification. No stochastic trace substitution, tolerance relaxation, new object, new fit or new physical call occurred.
