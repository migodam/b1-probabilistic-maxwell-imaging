# Full-state catalog for theory2 / Pro

Only typed JSONmetadata andNPYheaders were inspected. No m/J/Q/A payload,GT,forward or optimizer was read/run. CSV/JSON preserve every path andtyped field.

|Local source group|Files|Kinds|Material dimensions|Source counts|Known stage rows|
|---|---:|---|---|---|---:|
|RESULTS/dualmetric_pilot12_v1|435|SCALAR_TM_NATIVE_FOURIER_APCS|576|12,6|435|
|RESULTS/shared_readout4_v1|160|SCALAR_TM_NATIVE_FOURIER_APCS|576|12,6|160|
|RESULTS/vector4_v1|168|UNKNOWN_OR_PARTIAL,VECTOR_CELL_MAJOR_XYZ|216,None|60|144|
|audit/coupled_mode_trajectory_v1|65|SCALAR_TM_NATIVE_FOURIER_APCS|64|6|65|
|audit/finite_mode_trajectory_v1|65|SCALAR_TM_NATIVE_FOURIER_APCS|64|6|65|
|audit/interventions|665|SCALAR_TM_NATIVE_FOURIER_APCS|576|12,6|665|
|audit/local_closed_loop_v1|1|SCALAR_TM_NATIVE_FOURIER_APCS|64|6|1|
|audit/wiring_v1|3|SCALAR_TM_NATIVE_FOURIER_APCS|576|6|3|

Scalar state: m(realN),A(complexnativecoefficients×sources),J/Q(complexN×sources),oldprecurrentgradient/direction,DPCS/source metrics,actual4cell/Green caches andstage loopmetadata. A count is ambient,not independent current rank. Vector state is cell-majorXYZ withone real material percell,vectorrangebasis/prstate,andrequires thecompatibleexistingmodel.

A full auxiliary iterator checkpoint isnot a whole-process restart: physicalLU cache,timers,controller/selector state andruntimeprior objects areoutside thisstructural certificate. PRmemory explicitNone may be a legalstageentry/reset;never silently substitute zero memory. Use SAME fullstate path for paired interventions;material-only endpoints areinsufficient. Cache presence isstructural,not a recomputed physics/identity check. Missing stage/pair metadata remainsunknown;path-derived pair labels areexplicitlyqualified.
