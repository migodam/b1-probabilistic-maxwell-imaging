# R15 public integration check

Status: COMPLETE; attempts 7, failures 0. Only public r14_new_0000; no GT, main experiment, scoring or corrector.

Native G19 uses original frozen R13/R8 program with diagnostic max_iterations3; new stateful uses the same config and three step calls. Callback reads existing R8 frame arrays without replacing formulas. Physical warm state and original readout are preserved. Cost/trajectory equivalence is not assumed.

Per-attempt receipts preserve every failure and full cumulative physical/auxiliary/projector counters. Restored historical counters are not fresh-work totals; attempt process CPU includes all continuation/replay and serialization work.

Recorded attempt CPU 3.763848s; whole process CPU 4.105439s; physical LU 13, actual complex RHS 156.

## Native vs stateful G19
Step 1: m rel=0.0, J rel=0.0, Q rel=0.0, A rel=0.0, old_grad rel=0.0, old_direction rel=0.0; PR rel=0.0
Step 2: m rel=0.0, J rel=0.0, Q rel=0.0, A rel=0.0, old_grad rel=0.0, old_direction rel=0.0; PR rel=0.0
Step 3: m rel=0.0, J rel=0.0, Q rel=0.0, A rel=0.0, old_grad rel=0.0, old_direction rel=0.0; PR rel=0.0

## twofold_hop
17 steps, cycle events=1; all finite=True, all material legal=True.
Checkpoint next-step differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=0.0, old_direction abs=0.0, memory_X abs=0.0, memory_Q abs=0.0.

## twofold_dr
17 steps, cycle events=1; all finite=True, all material legal=True.
Checkpoint next-step differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=0.0, old_direction abs=0.0, memory_X abs=0.0, memory_Q abs=0.0.
Pending-before-cycle restore differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=None, old_direction abs=None, memory_X abs=0.0, memory_Q abs=0.0.

## twofold_ap
17 steps, cycle events=1; all finite=True, all material legal=True.
Checkpoint next-step differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=0.0, old_direction abs=0.0, memory_X abs=0.0, memory_Q abs=0.0.

## twofold_full_dr
17 steps, cycle events=1; all finite=True, all material legal=True.
Checkpoint next-step differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=0.0, old_direction abs=0.0, memory_X abs=0.0, memory_Q abs=0.0.

## twofold_reset
17 steps, cycle events=1; all finite=True, all material legal=True.
Checkpoint next-step differences: m abs=0.0, J abs=0.0, Q abs=0.0, A abs=0.0, old_grad abs=0.0, old_direction abs=0.0, memory_X abs=0.0, memory_Q abs=0.0.

This checks software integration only. Parent owns interpretation and numerical scope. No parameter was retuned to force equality.
