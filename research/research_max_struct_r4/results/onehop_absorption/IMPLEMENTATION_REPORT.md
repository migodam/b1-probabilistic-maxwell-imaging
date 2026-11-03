# One-hop residual control: implemented and single development run sealed

No GT was read by this implementation or runner. Scoring and final scientific judgment remain with the parent Codex thread.

The candidate uses `u + diag(a) Goff u`, with `u=beta*(J-mQ)/(1+(beta-Gii)*m)` and `a=m/(1+(beta-Gii)*m)`. `Goff` is the fixed offdiagonal Green matrix inside disjoint 2x2 spatial blocks. There is no four-cell inverse or small LU. The exact material derivative includes both the derivative of `u` and the derivative of `a`.

Whole-box sufficient invertibility bounds on [.02,20] are 0.16809971545079827 at beta 0.2 and 0.03784328908973763 at beta 1. The single material directional check passed with relative error 1.8949077143400813e-8; fixed J/Q, no physical solve.

The frozen run retains the scalar CIE seed, global PR, source metrics, Fourier/DPCS, 3x500 stages, material inner maxiter 60, material bounds/prior and the R4 corrector maxiter 1000. All 1500 PR steps ran. The corrector hit its 1000-iteration limit; convergence is not established.

| Variant | Auxiliary wall (s) | Total wall (s) | Auxiliary original Phi | Final original Phi | Final relative residual |
|---|---:|---:|---:|---:|---:|
| self_absorption | 26.06145 | 37.28599 | 0.01080189472 | 0.005832503897 | 0.1078206516 |
| block_absorption | 43.24644 | 59.42373 | 0.0007066885667 | 6.506278835e-06 | 0.0009670144345 |
| block_absorption_cached | 35.42837 | 46.37474 | 0.0007045806703 | 6.49499025e-06 | 0.0009695368273 |
| onehop_absorption | 31.14277 | 42.16572 | 0.0007039845411 | 6.492062459e-06 | 0.0009678872549 |

Against cached full-block absorption, the observed auxiliary time ratio is 1.1376 and whole-run ratio is 1.0998. These are single-run CPU1 observations with fixed iteration budgets; arithmetic and optimizer trajectories differ.

Counted one-hop transformed-residual and derivative arithmetic (inherited scalar-seed arithmetic is outside these new counters): 1,064,386,944 scalar complex divisions; 248,840,064 four-cell matrix RHS applications. No small LU is charged. Total physical work: 1,032 LU and 24,684 complex RHS.

All three stage NPZ files retain material/current/Q/coefficients. `prediction.npz` retains selected physical material, predicted data, current and Q. Source/config snapshots and receipts are preserved in this directory.

No imaging/originality/generalization claim is made. Lower original Phi or residual alone does not establish GT recovery. No tuning, additional cases, GPU/remote work, commit or push occurred.

A prelaunch executable path check exited127 before Python started because of a literal path-quoting typo; it consumed zero physics and is disclosed in `IMPLEMENTATION_SUMMARY.json`. The scientific runner was launched once and completed.
