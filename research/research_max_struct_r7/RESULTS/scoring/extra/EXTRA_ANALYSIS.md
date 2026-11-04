# Offline extra-run recovery diagnostics

These repeated/disclosed objects add zero independent objects to the R7 new30 denominator. All GT scoring follows complete primary sealing.

|Cohort|Arm/readout|Runs|Existing objects|Material successes|Median NRMSE|CPU s|LU|
|---|---|---:|---:|---:|---:|---:|---:|
|DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20|original_a|4|4|3|0.0559|284.21|6967|
|DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20|mobius|4|4|3|0.0551|347.44|7021|
|DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20|mobius_arc_candidates|4|4|3|0.0552|197.78|6709|
|DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20|mobius_guarded|4|4|3|0.0557|268.02|7014|
|DISCLOSED_FOUR_OBJECT_APLUS_DEVELOPMENT20|projected|4|4|2|0.2275|323.03|7511|
|DISCLOSED_OLD_DIGIT9_REGRESSION|a_plus|9|9|8|0.0832|345.31|15352|
|REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC|a_material|6|6|5|0.0641|435.58|10865|
|REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC|a_plus|6|6|5|0.0632|338.20|10497|
|REPEATED_NEW_OBJECT_INITIALIZATION_DIAGNOSTIC|cie_dual|6|6|4|0.0754|286.50|12015|

Init18 repeats six existing new30 objects at m=1 against their fixed primary m=.1 endpoints for three arms. Both predictions remain visible; no GT-selected multi-start replacement enters the primary table.

Development20 is five readout variants on four disclosed historical objects. Recovery is scored only after the Aplus publicPhi/cost-based choice was frozen; post-hoc GT does not authorize algorithm reselection.

The9 Aplus digit cases are disclosed historical regression. Old R6 dense endpoints are copied as labeled references without rerun. FFT action roundoff, readout changes, and resulting long-trajectory path changes prevent interpreting the new-vs-old difference as a pure readout causal effect. These old digits are not new confirmation.

Scoring reads complete576-variable endpoint materials and paid Y_prediction only;0new physical forward/adjoint calls.
