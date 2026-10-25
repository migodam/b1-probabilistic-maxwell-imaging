# A0 separate scorer implementation

Implemented `tools/a0_score.py`; not run, no physical calls spent by this worker. No fitting changes.

The main entrypoint first requires a terminal durable online optimization report and E0 PASS. Only then does it read the independently stored evaluation labels and outer-aperture holdout. Every saved online record survives into `results/a0/scoring/FULL_TABLE.csv` and `RESULTS.json`, including failed/unrun records and selected nonconverged fits. Labels never enter optimization, candidate selection, nuisance profiling, or stopping.

Each selected candidate is evaluated with the qualified-family FEM settings h=.17, n=256, field grid=2 at the two frozen outer Tx and 16 Rx; the saved shared complex gain is applied unchanged. All actual solves charge the same A0 BudgetBook. Saved base predictions are used to recompute training energy. Both training and independent holdout energy are compared only with the existing 64-real-component conservative one-sided noise envelope. This is not a fitted-degrees-of-freedom test or posterior-calibration gate.

Offline outputs include target-center Euclidean/squared error, area-weighted public-24x6-grid epsilon-contrast and loss L2/RMSE errors, support IoU and symmetric difference, exact grid boundary precision/recall/F1 and length disagreement, component/hole-count agreement, and declared-geometry feasibility. No new accuracy cutoff is imposed.

After primary scoring, fixed parent IDs 0,3,6,9, all saved noise repeats, all methods receive serial independent h=.17/n=384/grid=3 outer-channel prediction checks. Selection is fixed by ID, never by error or label. Prediction differences and finer holdout energies are descriptive; no new numerical gate is invented. Field state is cleared between calls. A budget exception stops subsequent physical work and preserves remaining statuses as NOT_RUN; no cap reset/refund occurs.

Per-candidate primary and independent NPZs contain parameters/support, predictions, whitened residuals/differences, Tx/Rx and solver metadata. Incremental JSON/CSV include costs, errors, nonconvergence and audit paths. `results/a0/scoring/gate.json` remains INCONCLUSIVE pending parent scientific judgment, even when code completes. This scorer provides no formal posterior or unrestricted-topology claim.

Execution is the parent's responsibility through the existing budgeted A0 job wrapper; import alone performs no physical work.
