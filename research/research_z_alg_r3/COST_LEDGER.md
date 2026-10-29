# R3 cost ledger

All RHS below are actual complex columns for 12-dimensional full systems. Reduced and algebraic least-squares RHS are reported separately, never converted to full RHS equivalents. Every original DEV and holdout receipt is charged once; frozen copies are not new runs.

| Stage | Full RHS | Reduced RHS | Algebraic LS RHS | Full factors | Reduced factors | Inner QR | solver/validation wall s |
|---|---:|---:|---:|---:|---:|---:|---:|
| dev_solver | 249318 | 9096 | 45009 | 44086 | 2542 | 15003 | 21.746799 |
| holdout_solver | 532989 | 14220 | 108126 | 94318 | 4260 | 36042 | 27.617566 |
| confirm_v7_solver | 324207 | 0 | 136167 | 58109 | 0 | 45389 | 13.644571 |
| all_solver | 1106514 | 23316 | 289302 | 196513 | 6802 | 96434 | 63.008937 |
| generation | 30 | 0 | 0 | 10 | 0 | 0 | 0.008526 |
| independent_validation | 13551 | 0 | 0 | 4517 | 0 | 0 | 0.334443 |
| scoring | 0 | 0 | 0 | 0 | 0 | 0 | 2.380141 |

holdout_v6 end-to-end subprocess time: 42.805927 s including interpreter/imports and result IO. This contains the numerical wall above and is not added to it.

confirm_v7 end-to-end subprocess time: 18.631361 s including interpreter/imports and result IO. This contains the numerical wall above and is not added to it.

Generation operation time is separately retained in JSON. Shared geometry setup is charged there. Solver time includes cached operator import and all numerical work; Python startup/import and file IO are excluded. The failed initial generation consumed zero RHS; its elapsed time and legacy DEV scoring times were not retained.

The JSON provides per-run method/revision, termination status, all operation counts/times, fit invocations, candidate counts and rejection counts. No R2 experiment is rerun or charged as new work.

Grand paid physical full RHS including all DEV/V6/V7 solver work, target generation and independent validation: 1120095. Separately, reduced RHS: 23316; algebraic complex least-squares RHS: 289302. The stage subtotals and all-solver total are overlapping views and must not be added together.
