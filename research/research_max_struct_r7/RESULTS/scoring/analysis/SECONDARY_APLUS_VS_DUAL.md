# Secondary requested comparison: Aplus versus mature sourcewiseDual

SECONDARY_USER_REQUESTED. This append does not alter the original three comparison ledgers. It compares complete algorithms on the same30 objects. No new threshold, noninferiority margin, optimizer run or physical call is introduced.

|Endpoint metric|Mean Aplus-minus-Dual|95% paired object bootstrap CI|
|---|---:|---|
|material_nrmse|-0.0377787|[-0.0929187, -0.00308697]|
|material_image_success|0.0333333|[-0.0666667, 0.133333]|
|cpu_seconds|5.73959|[1.91586, 9.54081]|
|LU|-41.7333|[-207.806, 101.335]|
|actual_complex_RHS|-1134|[-4440.85, 1665.22]|

|Paired common maximum axis|Aplus material successes|Dual material successes|Mean NRMSE delta|95% CI|
|---|---:|---:|---:|---|
|CPU|22/30|20/30|-0.0368365|[-0.0914344, -0.00269315]|
|LU|22/30|20/30|-0.0385691|[-0.0930752, -0.00376244]|
|complex_RHS|22/30|20/30|-0.0385691|[-0.0932036, -0.00380701]|

Each pair/object cost window is [max(first observed cost), min(last observed cost)] separately for actual process-prefix CPU, physical LU and complex RHS. Eleven fixed window fractions use the last saved paid original-Phi incumbent at or below budget. No stopped-method extrapolation. Auxiliary work is counted in CPU; equal LU is not equal CPU. Float32/every-ten-record plus sparse auxiliary-stage resolution limitations remain.

Endpoint continuous means and paired medians answer different questions. Success-rate intervals include finite-object uncertainty; no formal noninferiority or population-wide superiority claim follows.
