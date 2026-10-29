# Actual fixed-shape probability workload replay

Status: **NUMERICAL_FIXED_SHAPE_TASK_REPLAY**. Same exposed DEV shapeR9V_FRESH_000, exactly836 saved revision2 event_split_original_GL level1 nodes and weights. No refit, new rank selection, material GT access, or new quadrature reference.

Rank64 passed836/836 residual/potential gates, zero fallback. Maximum computed white bound=1.8820536e-08; all sampled bounds cover actual response error.

Original Uniform2 density=1/.15. Exact/reduced log evidence=-233.7226700023664/-233.7226700023636. Relative evidence difference=2.7853275e-12, maximum amplitude mean difference=8.437695e-15, event probability difference=9.5368158e-14; all pass original tolerances2e-4/2e-5/5e-4. Events: Re(delta)>1.3 and Im(delta)>.105.

Actual cold kernel plus quadrature: exact Schur0.287601074s, reduced plus full residual/gate0.117269774s; saving59.2248%. Both shape operators were rebuilt cold. Historical profile/Hessian/event-grid construction is explicitly reused common history. Its original transfer69784RHS256/270LU includes construction, but no separately measured pure-grid time exists; none is invented or subtracted.

Actual added cost=13376 cross-order RHS:6688order256+6688order64;0LU256,0RHS576. Job wall0.865565s, CPU0.842065s. Responses, weights and all paid attempts remain in results/reduced_response_probability_workload/.

Independent audit: `audit/independent_qa/PROBABILITY_WORKLOAD_AUDIT.json`, PASS,0newphysics. This is a DEV fixed-shape numerical task replay, not independent material evidence, a new conditional reference, a full-support tail/rounding certificate, or Full66 posterior qualification.
