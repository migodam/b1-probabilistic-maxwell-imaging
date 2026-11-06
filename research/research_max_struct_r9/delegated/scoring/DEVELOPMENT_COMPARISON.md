# R9 fixed four-case development comparison

Each cell: amplitudeNRMSE / OtsuIoU / inclusive CPU seconds. PASS is the unchanged material success rule. Fixed rows prefer fresh R9 same-machine controls; [old] denotes reused R8 Mac0000 only. Old Windows0027 CPU is retained in full scores, never used as a local control here.

| Method | b20d | 40ac | 0000 | 0027 |
|---|---:|---:|---:|---:|
| fixed .001 Dual | 0.3896 / 0.423 / 46.7 | 0.6335 / 0.269 / 52.2 | 0.0817 / 0.939 / 45.7 PASS [old] | 1.0139 / 0.234 / 66.3 |
| fixed .001 GlobalDual | 0.4880 / 0.303 / 27.9 | 0.6190 / 0.283 / 55.0 | 0.0816 / 0.936 / 47.8 PASS [old] | 1.0128 / 0.234 / 65.4 |
| fixed .001 A+ global | 0.0406 / 0.952 / 44.0 PASS | 0.6260 / 0.277 / 66.2 | 0.0811 / 0.936 / 59.3 PASS [old] | 1.0151 / 0.227 / 70.8 |
| fixed .001 A+ sourcewise | 0.0404 / 0.940 / 46.7 PASS | 0.6259 / 0.271 / 67.7 | 0.0813 / 0.939 / 54.3 PASS [old] | 1.0188 / 0.236 / 61.5 |
| noise99 Dual | 0.3896 / 0.418 / 49.0 | 0.6335 / 0.269 / 53.3 | 0.0845 / 0.922 / 28.3 PASS | 0.1191 / 0.856 / 46.3 PASS |
| noise99 GlobalDual | 0.3895 / 0.423 / 40.7 | 0.6190 / 0.283 / 54.5 | 0.0812 / 0.979 / 30.4 PASS | 0.1197 / 0.850 / 43.8 PASS |
| noise99 A+ global | 0.0427 / 0.940 / 45.2 PASS | 0.6260 / 0.277 / 68.0 | 0.0809 / 0.958 / 35.5 PASS | 0.1178 / 0.856 / 49.0 PASS |
| noise99 A+ sourcewise | 0.0429 / 0.940 / 46.2 PASS | 0.6259 / 0.271 / 70.4 | 0.0835 / 0.959 / 35.2 PASS | 0.1184 / 0.856 / 51.9 PASS |
| Data homotopy | 0.7321 / 0.000 / 11.1 | 0.4602 / 0.181 / 34.9 | 0.0788 / 0.898 / 30.2 PASS | 0.4989 / 0.062 / 10.3 |
| Profile fixed | 0.3813 / 0.080 / 12.1 | 0.4577 / 0.219 / 24.6 | 0.0818 / 0.938 / 15.0 PASS | 0.5583 / 0.160 / 16.3 |
| Profile continuation | 0.3771 / 0.289 / 27.6 | 0.0581 / 0.962 / 28.1 PASS | 0.0807 / 0.917 / 16.5 PASS | 0.5144 / 0.155 / 13.7 |
| Direct physics | 0.5345 / 0.228 / 37.0 | 0.4510 / 0.212 / 27.6 | 0.0805 / 0.938 / 22.6 PASS [old] | 0.4609 / 0.324 / 25.7 |

1. 40ac: known ordinary-LSIE profile continuation recovers this exposed development case: NRMSE .05807, IoU .96183, residual .0009556, original Phi8.05112e-6, CPU28.10s. Fixed profile, direct physics, and all fixed/noise CIE/A+ controls fail. This is one local known-control result, not a new algorithm or universal basin theorem.
2. 0027: noise99 selects rank11 in all four arms; all pass material and original noise metrics. Fresh local fixed controls all fail and reach the3000-corrector-iteration limit. This is not an A+-specific gain. Its three landscape arms and direct physics also fail.
3. b20d: both A+ couplings pass under fixed/noise rank, while ordinary Dual/GlobalDual, three landscape arms and direct fail. Normal control0000 passes across displayed methods.

New43 task inclusive CPU1748.976415s; task wall sum1788.842416s is concurrent workload, not queue elapsed time. Full-wave LU72717; complexRHS1743240. Original28 import failures are preserved separately: zero physical calls, known archived CPU .098573s is a lower bound; baseline failure CPU was unmeasured. Inherited R8 costs are historical and not charged again.

Successful offline scoring CPU1.849007s/wall2.056370s, zero new physical calls. Initial KeyError scoring CPU/wall unknown, separately retained in FAILED_SCORING_ATTEMPT_01.json. Saved-array original Phi agrees for all59 final rows within2.776e-17.

All59 rows include Phi, stop, actual steps, selected stage, costs and platform in RESULTS/development_scoring/SCORES.json and CSV.10727 already-paid incumbent snapshots are graded offline;287 use an explicit UPPER_SHIFT_including_output_tail. No GT iterate selection or interpolation.

Actual saved receiver ranks (fixed/noise99): b20d17/15;40ac17/17;000017/12;002717/11. These are read from saved stage metadata, not recomputed SVDs.
