# R9-S exploration E0/E1/E2/E3 evidence

All eight R7 locked objects are now exposed development objects. New Maxwell work is complete; fresh confirmation has not been run.

| Arm | r | available | qualified | operational RHS | audit RHS | operational wall(s) | initial-state-only |
|---|---:|---:|---:|---:|---:|---:|---:|
|A_NATIVE_FULLJ|0.01|8|8|696|32|1.6028|0|
|B_FIXED_WEAK|0.01|8|8|368|32|1.0542|7|
|C_DYNAMIC_WEAK|0.01|8|8|368|32|1.0239|7|
|D_FOURIER|0.01|8|6|704|256|1.4368|0|
|D_WAVELET|0.01|8|2|704|256|1.5645|0|
|E_SOM|0.01|0|0|0|0|0.0000|0|
|F_FROZENJ|0.01|8|8|192|264|0.7232|0|
|R_RANDOM|0.01|8|4|704|256|1.4039|0|
|S2_GLOBAL|0.01|8|7|656|224|1.7239|1|
|U_UNCORRECTED|0.01|8|8|64|256|0.1639|0|
|A_NATIVE_FULLJ|0.02|8|8|696|32|1.5555|0|
|B_FIXED_WEAK|0.02|8|7|360|0|1.0965|7|
|C_DYNAMIC_WEAK|0.02|8|7|368|32|0.9940|7|
|D_FOURIER|0.02|8|0|704|256|1.5511|0|
|D_WAVELET|0.02|7|0|656|224|1.3643|0|
|E_SOM|0.02|0|0|0|0|0.0000|0|
|F_FROZENJ|0.02|8|8|192|256|0.6182|0|
|R_RANDOM|0.02|8|1|704|256|1.2250|0|
|S2_GLOBAL|0.02|8|4|656|224|1.3075|1|
|U_UNCORRECTED|0.02|8|6|64|256|0.1938|0|

B/C paired endpoint comparisons:
- r=0.01: exact same x in 7/8; B/C operational RHS=368/368; summed operational wall=1.0542/1.0239 s; B relative wall saving=-2.959%; B/C qualified=8/8.
- r=0.02: exact same x in 7/8; B/C operational RHS=360/368; summed operational wall=1.0965/0.9940 s; B relative wall saving=-10.314%; B/C qualified=7/7.

Corrector stops (all retained):
- A_NATIVE_FULLJ: NONPOSITIVE_PREDICTED_REDUCTION x 2.
- B_FIXED_WEAK: NONPOSITIVE_PREDICTED_REDUCTION x 15.
- C_DYNAMIC_WEAK: NONPOSITIVE_PREDICTED_REDUCTION x 14.
- S2_GLOBAL: NO_FEASIBLE_RETRACTED_STEP x 2.

Total actual new physics: 13976 complex RHS, 387 states and LU factorizations, 0 solver failures. Reserved and callback counts match: True. Controlled wall=36.440s; CPU=32.477s; physics wall=5.909s; peak RSS=173293568 bytes. All method, rejected-trial, audit-restoration and path work is charged.

Native A reproduced R7 endpoint Phi to max absolute difference 0. This is a reproducibility check, not independent evidence.

Limits and interpretation boundaries:
- eight exposed historical DEV objects; no fresh confirmation.
- A uses historical radial retraction; B/C/S2 use fixed weak-coordinate decoder Newton retraction.
- B vs C same initial weak direction; dynamic refresh only current correction complement.
- a failed step or radius projection is retained, no extra optimizer rescue.
- SOM failure scoped to chosen dark-current lift, not all SOM null directions.
- Fourier/Wavelet mapped basis projections change target directions; projection errors disclosed.
- all current endpoint J diagnostic costs charged after method-selection clock.
- five-node actual correction path samples do not prove continuous fibers.
- no posterior or exact/continuous material-fiber claim.
