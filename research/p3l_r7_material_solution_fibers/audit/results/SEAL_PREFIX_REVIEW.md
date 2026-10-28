# Independent R7 sealed-result audit

PASS. Zero new physics; no NN/Maxwell imports and no offline truth loaded. Only fresh locked8 appear in outcome tables; development and preserved revisions remain billed in campaign costs. This audit supplies factual checks and descriptive comparisons, not novelty judgment.

## Accounting

{
  "totals": {
    "RHS": 36752,
    "states": 1285,
    "LU": 1285,
    "forward_RHS": 10280,
    "derivative_RHS": 26472,
    "failures": 0
  },
  "sequence": 6567,
  "event_counts": {
    "JOB_BEGIN": 16,
    "PHYSICS_BEGIN": 2176,
    "STATE_RESERVED": 1285,
    "PHYSICS_END": 2176,
    "HISTORICAL_GEOMETRY_REUSE": 7,
    "DERIVATIVE_RESERVED": 891,
    "JOB_END": 16
  },
  "CPU_seconds": 65.251603,
  "controlled_wall_seconds": 69.81946049863473,
  "physics_wall_seconds": 14.141339980298653
}

Historical geometry reuse is recorded separately in METRICS.json; cold-equivalent derivative RHS is not a new paid charge. All journal charges, including failures and preserved development revisions, remain in the total. CPU and wall clocks are separately replayed.

## Main outcomes

Success requires distance, signed progress, noise ball, DeltaPhi <=2 and legal decoder; the prior heuristic is an additional column. Missing/infeasible endpoints count as failure for all8 success comparisons.

| arm | radius | executed | success_no_prior | success_with_prior | RHS_total |
| --- | --- | --- | --- | --- | --- |
| A_FROZEN | 0.0025 | 8 | 8 | 8 | 192 |
| A_FROZEN | 0.005 | 8 | 8 | 8 | 192 |
| A_FROZEN | 0.01 | 8 | 8 | 8 | 192 |
| A_FROZEN | 0.02 | 8 | 8 | 8 | 192 |
| B_FULLJ | 0.0025 | 8 | 8 | 8 | 656 |
| B_FULLJ | 0.005 | 8 | 8 | 8 | 656 |
| B_FULLJ | 0.01 | 8 | 8 | 8 | 696 |
| B_FULLJ | 0.02 | 8 | 8 | 8 | 696 |
| C_JETS | 0.0025 | 8 | 8 | 8 | 656 |
| C_JETS | 0.005 | 8 | 8 | 8 | 688 |
| C_JETS | 0.01 | 8 | 8 | 8 | 688 |
| C_JETS | 0.02 | 6 | 6 | 6 | 520 |
| D_RANDOM | 0.0025 | 8 | 8 | 8 | 704 |
| D_RANDOM | 0.005 | 8 | 8 | 8 | 704 |
| D_RANDOM | 0.01 | 8 | 4 | 4 | 704 |
| D_RANDOM | 0.02 | 8 | 1 | 1 | 704 |

## Paired descriptive bootstrap

Parent resampling, seed 2026102811, 20,000 replicates, percentile95. No large-population inference. Endpoint-only RHS excludes shared setup. The second cost comparison charges the full shared 24 RHS jet setup to C at .02 on every parent, including unavailable endpoints; this is conservative for a four-radius amortized use. Both costs remain in the campaign ledger. Phi comparisons condition on both endpoints being available and must not be interpreted as all8 superiority.

```json
{
  "success_no_prior_C_minus_B_all8": {
    "n": 8,
    "mean": -0.25,
    "percentile95": [
      -0.625,
      0.0
    ],
    "unit": "parent",
    "replicates": 20000
  },
  "success_with_prior_C_minus_B_all8": {
    "n": 8,
    "mean": -0.25,
    "percentile95": [
      -0.625,
      0.0
    ],
    "unit": "parent",
    "replicates": 20000
  },
  "C_setup_RHS_per_parent": [
    24,
    24,
    24,
    24,
    24,
    24,
    24,
    24
  ],
  "paid_endpoint_plus_full_shared_setup_RHS_C_minus_B_all8": {
    "n": 8,
    "mean": 2.0,
    "percentile95": [
      -31.0,
      24.0
    ],
    "unit": "parent",
    "replicates": 20000
  },
  "paid_endpoint_RHS_C_minus_B_all8": {
    "n": 8,
    "mean": -22.0,
    "percentile95": [
      -55.0,
      0.0
    ],
    "unit": "parent",
    "replicates": 20000
  },
  "both_executed_parents": [
    "R7_LOCKED_000",
    "R7_LOCKED_001",
    "R7_LOCKED_002",
    "R7_LOCKED_005",
    "R7_LOCKED_006",
    "R7_LOCKED_007"
  ],
  "C_unavailable_explicit_failures": [
    "R7_LOCKED_003",
    "R7_LOCKED_004"
  ],
  "Phi_C_minus_B_CONDITIONAL_BOTH_EXECUTED": {
    "n": 6,
    "mean": -0.011539306680210379,
    "percentile95": [
      -0.06021921132051489,
      0.02132922320575877
    ],
    "unit": "parent",
    "replicates": 20000
  }
}
```

## All16 statuses per parent

| parent | A_FROZEN | B_FULLJ | C_JETS | D_RANDOM |
| --- | --- | --- | --- | --- |
| R7_LOCKED_000 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_001 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_002 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_003 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: INFEASIBLE_LATENT_RADIUS | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_004 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: INFEASIBLE_LATENT_RADIUS | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_005 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_006 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |
| R7_LOCKED_007 | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED | 0.0025: EXECUTED, 0.005: EXECUTED, 0.01: EXECUTED, 0.02: EXECUTED |

## Paths and source groups

Two predetermined paths per parent are tested at five nodes only. Path failure does not show disconnectedness; node success does not certify between-node continuity. Compensation identity is independently recomputed using real arrays and the same lower-order path. Geometry and source spectra use the original disjoint source groups, not a new acquisition layout.

| parent | path | status | qualified |
| --- | --- | --- | --- |
| R7_LOCKED_000 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_000 | C_RAW_CUBIC | EXECUTED | True |
| R7_LOCKED_001 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_001 | C_RAW_CUBIC | EXECUTED | True |
| R7_LOCKED_002 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_002 | C_RAW_CUBIC | EXECUTED | True |
| R7_LOCKED_003 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_003 | C_RAW_CUBIC | UNAVAILABLE | False |
| R7_LOCKED_004 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_004 | C_RAW_CUBIC | UNAVAILABLE | False |
| R7_LOCKED_005 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_005 | C_RAW_CUBIC | EXECUTED | True |
| R7_LOCKED_006 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_006 | C_RAW_CUBIC | EXECUTED | True |
| R7_LOCKED_007 | B_FULLJ_CHORD | EXECUTED | True |
| R7_LOCKED_007 | C_RAW_CUBIC | EXECUTED | True |

## Anchors and regularity

Locked centers are measurement-only Gaussianized-coordinate optimizer points. Optimizer status, nonzero gradient, fresh-state attempts, and costs are retained in METRICS.json; these points are not material MAPs. All six available C raw paths cross decoder seams by their final node (1 to 6 cells). Their finite jet errors do not establish exact Taylor validity across seams.

## Trial selection

{
  "paid_main_trials": 365,
  "infeasible_trials": 0,
  "not_selected_trials": 239,
  "stop_events": 13
}

Not selected includes feasible paid trials rejected by objective selection; infeasible and not-selected overlap. Counts here concern main arms; ledger failures count physical action errors. Detailed check records, conditional means, costs and source identities are in METRICS.json and CSV companions.
