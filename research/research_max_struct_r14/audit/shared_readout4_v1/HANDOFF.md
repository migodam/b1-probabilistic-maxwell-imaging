# Common material solver control handoff

16/16 attempt02 tasks completed with3000nativeupdates each; allcheckpointtags match shared_LBFGS60, no seed, identicalchi-coordinate L-BFGS-B60/ftol1e-12/gtol1e-8/maxls25, old-m nonincrease safeguard, zeroArc readouts. Private stepcode is the identical corecodeobject with only the two material globals rebound; originalcoreglobals remain untouched. Native/current/stage/physicalcorrector code and allfrozen/independent algorithms were not changed.

Attempt01 failed at psutil import before model/physics setup,0LU/RHS; all16 failurelogs/seal and.514267OSchildCPU retained. Attempt02 reuses existingABI3psutil7.2.2 as private package-only PYTHONPATH under the requestednumericalPython, with genuine ProcessRSS; no installation/provider/newbackend. Runtime reusecost is PSUTIL_RUNTIME_RECEIPT.json.

Actual16predictions are sealed at predictions/shared_readout4_v1/attempt_02/PREDICTION_SEAL.json; queue receipt at queues/shared_readout4_v1/attempt_02/RECEIPT.json. NoGTaccess/scoring yet; require parentauthorization after reviewing this seal. Allraw/corrected/primary/secondaryavailability/history/pool/PRstates/costs follow unchangedparentrunner; missingsecondary remains explicit, and raw retains the parent'sphase-specific semantics.

One inheritedcounter reconciliation: core.step adds2 whitening/source-product counters after everyM1readout, but directonehopvaluegradient usesOmega without JLs/QLs. Rawcounters remain unchanged; TASK_COSTS_AND_COUNTER_RECONCILIATION subtracts2perM1readout for those two counts only. No numericalstates, CPU, LU/RHS orselection depend on those counters. No reruns or hiddenoverwrites.

Four exposedcases only; this commonmaterialinnerprogram is a mechanismcontrol, not a newmainalgorithm/full48result or blindconfirmation. Compare actualCPU within this local16 group; noMac-versus-Windows speedclaim. No scorer was added or run for this control.
