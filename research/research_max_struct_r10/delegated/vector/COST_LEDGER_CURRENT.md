# Vector R10 live cost ledger

Status: FINAL_KNOWN_LOWER_BOUND; sealed 18/18; GT scoring NOT_RUN.

Known CPU lower bounds: Mac 4095.307625s; XINAN 5962.796875s. Active cumulative counters are replaced, never added to final same-process CPU.

Whole wait4/GetProcessTimes children are additive. Setup/arm/forward/adjoint/device-event views are nested. CUDA events include host/pause gaps and are not active kernel time. Different vector/scalar dimensions are not reduced to one speed metric.

Historical GPU waiting and engineering failures remain charged/unknown. Publication and old R9 costs excluded. Raw source ledgers are unchanged.

Running: 

Unmeasured helper/failed-receipt CPU remains null, not zero.
