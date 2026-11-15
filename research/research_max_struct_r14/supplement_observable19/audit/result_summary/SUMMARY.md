# Observable19 saved-output mechanical summary

Source: `research/research_max_struct_r14/supplement_observable19/analysis/observable19_v1`. Receipt `COMPLETE`, 0 issues;144 tasks, all144 external OS CPU measured, all24 objects retained per method. No new scorer/physics or GT-array access in this summary.

Exact domain: exposed24 scalar-TM observations,576 independent real material variables, fixed19 auxiliary denominator basis; G/S/O0 × PGP/GP.20 matched-discrete cases and4 cross-grid diagnostics. Raw CN99 noise-only feasibility is not a model-mismatch certificate. No blind/independent-cohort/3D or scientific superiority claim.

## Primary / raw / corrected / secondary

|Arm|View|Recovered/24|NRMSE mean / median|IoU mean / median|Noise feasible/24|
|---|---|---:|---:|---:|---:|
|G-PGP|primary|10|0.198433 / 0.100683|0.637959 / 0.709692|18|
|G-PGP|raw|5|0.207667 / 0.140351|0.638159 / 0.624488|3|
|G-PGP|corrected|10|0.198433 / 0.100683|0.637959 / 0.709692|18|
|G-PGP|secondary|11|0.197173 / 0.094628|0.633487 / 0.688983|18|
|S-PGP|primary|12|0.158919 / 0.097155|0.695155 / 0.798655|20|
|S-PGP|raw|9|0.173227 / 0.134876|0.675187 / 0.699283|1|
|S-PGP|corrected|12|0.158919 / 0.097155|0.695155 / 0.798655|20|
|S-PGP|secondary|13|0.158715 / 0.104512|0.688710 / 0.815263|20|
|O0-PGP|primary|12|0.159954 / 0.099093|0.696943 / 0.815718|20|
|O0-PGP|raw|10|0.172511 / 0.133170|0.671266 / 0.690659|1|
|O0-PGP|corrected|12|0.159954 / 0.099093|0.696943 / 0.815718|20|
|O0-PGP|secondary|13|0.158079 / 0.099093|0.690174 / 0.827406|20|
|G-GP|primary|10|0.197667 / 0.098848|0.639345 / 0.703843|18|
|G-GP|raw|4|0.206821 / 0.148877|0.634629 / 0.625440|3|
|G-GP|corrected|10|0.197667 / 0.098848|0.639345 / 0.703843|18|
|G-GP|secondary|11|0.195943 / 0.094487|0.631083 / 0.693220|18|
|S-GP|primary|12|0.161074 / 0.102036|0.691570 / 0.792785|20|
|S-GP|raw|7|0.176243 / 0.134435|0.685456 / 0.680438|1|
|S-GP|corrected|12|0.161074 / 0.102036|0.691570 / 0.792785|20|
|S-GP|secondary|13|0.162163 / 0.102036|0.679823 / 0.779347|20|
|O0-GP|primary|12|0.162209 / 0.101813|0.690316 / 0.799011|20|
|O0-GP|raw|6|0.178752 / 0.135142|0.673359 / 0.666891|1|
|O0-GP|corrected|12|0.162209 / 0.101813|0.690316 / 0.799011|20|
|O0-GP|secondary|13|0.162471 / 0.099310|0.679648 / 0.780800|20|

Recovery is the existing descriptive NRMSE<=.15 and IoU>=.70 rule, never a selector. All views have24 accuracy-available cases. Primary and corrected coincide in these saved rows. Secondary is a distinct original-R criterion with explicitly retained primary deployment fallback.

## Primary membership, relative to fresh G-PGP

- G-PGP success IDs: r14_new_0000, r14_new_0005, r14_new_0006, r14_new_0009, r14_new_0011, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0020.
- S-PGP success IDs: r14_new_0000, r14_new_0003, r14_new_0005, r14_new_0006, r14_new_0011, r14_new_0012, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0021, r14_new_0022.
  Gain: r14_new_0003, r14_new_0012, r14_new_0021, r14_new_0022; loss: r14_new_0009, r14_new_0020.
- O0-PGP success IDs: r14_new_0000, r14_new_0003, r14_new_0005, r14_new_0006, r14_new_0011, r14_new_0012, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0021, r14_new_0022.
  Gain: r14_new_0003, r14_new_0012, r14_new_0021, r14_new_0022; loss: r14_new_0009, r14_new_0020.
- G-GP success IDs: r14_new_0000, r14_new_0005, r14_new_0006, r14_new_0009, r14_new_0011, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0020.
  Gain: ; loss: none.
- S-GP success IDs: r14_new_0000, r14_new_0003, r14_new_0005, r14_new_0006, r14_new_0011, r14_new_0012, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0020, r14_new_0022.
  Gain: r14_new_0003, r14_new_0012, r14_new_0022; loss: r14_new_0009.
- O0-GP success IDs: r14_new_0000, r14_new_0003, r14_new_0005, r14_new_0006, r14_new_0011, r14_new_0012, r14_new_0014, r14_new_0015, r14_new_0016, r14_new_0017, r14_new_0020, r14_new_0022.
  Gain: r14_new_0003, r14_new_0012, r14_new_0022; loss: r14_new_0009.

## Actual paid cost

|Arm|Whole-child OS CPU s|Physical LU (576-order)|Actual complex physical RHS|Actual factory setup CPU s|
|---|---:|---:|---:|---:|
|G-PGP|2838.046875|16003|313860|12.234375|
|S-PGP|2882.562500|16703|323748|7.296875|
|O0-PGP|2889.609375|16578|322512|7.593750|
|G-GP|2866.187500|16490|317958|7.078125|
|S-GP|2853.328125|15845|310854|7.218750|
|O0-GP|2869.656250|16768|324600|7.578125|

Total whole-child OS CPU=17199.390625s. Offline remote scoring CPU=87.390625s, wall=137.503906s, separate from solver cost. Queue wall=1614.144759s and controller CPU=4.250000s; wall is not summed with child times.
Whole-child GetProcessTimes includes imports/final serialization. Setup CPU is nested in whole-child CPU, not added. Actual factory and nominal R13 counters are separate in COMPACT. O0 receiver-adjoint solve RHS are already in physical aggregate; dense Green/basis actions and small-LU work are auxiliary, not free physical RHS. Core type fields are null; no P/E performance attribution follows from logical IDs alone.

## Normalized saved-candidate views

|Arm|Common native CPU120 recovery / mean NRMSE|Common componentwise physical recovery / mean NRMSE|
|---|---:|---:|
|G-PGP|10/24 / 0.198156|10/24 / 0.197906|
|S-PGP|11/24 / 0.161249|11/24 / 0.161942|
|O0-PGP|11/24 / 0.162230|11/24 / 0.162312|
|G-GP|10/24 / 0.197328|10/24 / 0.197208|
|S-GP|11/24 / 0.164530|11/24 / 0.164551|
|O0-GP|11/24 / 0.165609|11/24 / 0.165799|

COMPACT retains every15/30/60/90/120 paid-prefix and common-native-CPU summary, per-arm success IDs, and normalized gain/loss sets. Common physics uses componentwise minima of actual final same-case physical LU/RHS across all six; selections are min eligible completed public Phi under both limits. Native candidate process timestamps are not exact external OS candidate timestamps. These are already-paid pool views, not reruns or refunds; full costs and auxiliary actions remain charged.

## r14_new_0021 (primary)

|Arm|NRMSE|IoU|Phi|Recovered|
|---|---:|---:|---:|---|
|G-PGP|0.49928930|0.20164609|0.035384846|False|
|S-PGP|0.05105661|0.93253968|5.3332699e-05|True|
|O0-PGP|0.05078431|0.93253968|5.3311627e-05|True|
|G-GP|0.48801800|0.22083333|0.035079892|False|
|S-GP|0.47222178|0.21338912|0.03605411|False|
|O0-GP|0.48702924|0.27935223|0.018444517|False|

This is r14_new_0021, not r13_fresh_0021. All four view values are retained in COMPACT; no cross-ID comparison is implied.

## Selector accounting and historical baseline difference

Selector audit: 720/720 consistent; ISSUES empty. Raw initial fallbacks=0; corrected pre-corrector fallbacks=0. Secondary→primary fallback counts: G-PGP=6, S-PGP=4, O0-PGP=4, G-GP=6, S-GP=4, O0-GP=4.

The frozen historical R14 G19 primary is11/24; fresh G-PGP is10/24. The only membership difference is r14_new_0019: old NRMSE=.0512005628, IoU=.7875, Phi=.0001834785882; fresh NRMSE=.0597677298, IoU=.6896551724, Phi=.0001830363462. Both NRMSE pass; fresh IoU is below.70. A lower saved Phi does not imply the same threshold membership.

This is not an algorithm-change comparison: both record PR/restricted/physical-current warm start, inherited audit schedule, nominal original G19 setup and absolute120-with-unused-proposal-reallocated contract. Actual old/fresh external CPU120.609375/120.53125s; internal final CPU120.375/120.390625s; selected candidate CPU119.984375/119.953125s. Logical cores16/2 (core type unknown), physical LU432/881 and RHS4962/10344, candidate IDs433/881, setup CPU.4375/.21875s differ. Fresh factory has private provenance/snapshot/spectral observers; identical nominal math/budget labels do not freeze the number of completed updates under a CPU cap. No causal explanation of this change is inferred here. Exact cost/current/setup records are in COMPACT.

Historical CIE/Dual/A+-Arc(pair01)/11/G19/DR references remain separately labeled and are excluded from the fresh six-arm strict normalizations. They are not new equal-stop runs.

Mechanical summary CPU=0.112428s, wall=0.118950s (reading existing scalar JSON and writing these summaries only; initial interactive inspection not instrumented). Final scientific interpretation remains with parent.
