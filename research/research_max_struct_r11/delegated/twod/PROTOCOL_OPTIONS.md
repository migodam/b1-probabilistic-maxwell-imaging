# R11 2D budget contract: root decision received

Root authorized DEV wiring at 480 whole-child CPU seconds and checkpoints 60/120/240/480. This is **equal allowed CPU budget**, not exact equal consumed work. Natural algorithm convergence retains its actual spent CPU and unused allowance. No idle filling, wall cap, automatic retry, or scientific schedule alteration.

- Each method has a separate, single-thread child. `time.process_time()` is absolute since child creation; imports, public data, noise99 DPCS, inherited setup, candidate prepare, all auxiliary stages, corrector, failures, observation and finalization are charged. External Windows GetProcessTimes is the authoritative completed child total, including tail after the last internal timer.
- Only already evaluated original-Phi states completed at/before a checkpoint are eligible. The best of those states is selected without GT. Atomic auxiliary/material/physical calls can overrun; their full work is charged, their late state is excluded from earlier checkpoints. No extra forward is invoked to export a prediction.
- Preserve the mature inherited Fourier schedule and physical corrector ceilings. Budget termination is cooperative at the original outer-loop `budget_check`/callback, and before/after each paid physical evaluation; it is not an iteration-count surrogate for cost. A candidate prepare is atomic and charged.
- Keep original Phi, continuous material scores after sealing, actual-cost curves and separate original physical LU/RHS views. Auxiliary small LU/FFT costs must remain distinct operation types. CPU and LU cannot be declared equivalent.
- Failed runs and setup-only exhaustion stay in the ledger. NO_PAID_STATE_WITHIN_CAP is a missing cap prediction, not zero cost or success.
- 48 independent 8-family × high/mixed × 3-geometry cases remain NOT_GENERATED. Contrast8/12/18 occurs within each family. No fresh GT or confirmation run until root freezes the candidate/protocol. Cross-grid testing needs a separate compatible runner: inherited R8 explicitly requires24x24/576 cells.

## Alternative not selected

Equal iteration/3000+3000 schedules are useful only as matched algorithm allowance, not equal actual CPU. Reusing late states for earlier checkpoints creates lookahead. Forcing an already-converged solver to consume480 CPU seconds changes the algorithm and is prohibited.
