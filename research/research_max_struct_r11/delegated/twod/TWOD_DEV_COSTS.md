# Exposed DEV actual CPU costs

Known CPU lower bound: 1165.234148 seconds. Status: COMPLETED_KNOWN_CPU_LOWER_BOUND.

|Task|Actual child CPU seconds|LU|Complex RHS|Stop|
|---|---:|---:|---:|---|
|r10_twod_0001__full_replay__480CPU|284.109375|3114|74688|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|r10_twod_0004__lagged__480CPU|113.796875|500|11940|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|r10_twod_0001__lagged__480CPU|266.328125|3127|75000|STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT|
|r10_twod_0022__lagged__480CPU|213.109375|2063|49452|CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH|
|r10_conditional_0000__lagged__480CPU|175.21875|1500|35952|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|
|r10_conditional_0003__lagged__480CPU|105.234375|358|8532|CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL|

No GT scoring. No budget padding. Child actual CPU includes setup/failure/finalization. Unknown control/setup overheads remain explicit in JSON. Nested internal timers are not charged twice.
